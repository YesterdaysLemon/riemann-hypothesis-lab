"""Finite rejection certificate for one arithmetic Nyman trial subspace.

The canonical aggregate space consists of the natural dilates
``rho_1, ..., rho_256`` together with eight functions supported on the new
block ``257 <= j <= 512``.  If ``V`` is the corresponding 256-by-8
coefficient matrix, the eight functions are the columns of

``sum_(j=257)^512 V[j, q] rho_j``.

This module certifies a strict lower bound for the squared distance ``F_V``
from the target to that 264-dimensional aggregate space.  The decisive
machine check is fixed-order interval ``LDL^T`` positivity of

``[[G_V, -b_V], [-b_V^T, 1-L]]``

at ``L=(9/10)U_256``.  Consequently ``F_V > (9/10)U_256``.  Since the
certified source endpoint gives ``d_256^2 <= U_256``, the gain captured by
this trial subspace satisfies ``Gamma_V < d_256^2/10``.

This rejects one finite ansatz.  It neither rejects every possible trial
subspace nor proves or disproves the Riemann Hypothesis.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

import flint
from flint import arb

from . import nyman as core
from . import nyman_beta2 as beta2
from .artifacts import content_sha256
from .balls import arb_record, arb_to_dyadic
from .nyman_summary import SUMMARY_SCHEMA


AUDIT_SCHEMA = "rh-lab/nyman-trial-subspace-audit/v1"
LOWER_CERTIFICATE_SCHEMA = "rh-lab/nyman-trial-subspace-lower/v1"
FROZEN_AUDIT_ID = "nyman-trial-subspace-v1"

FROZEN_OLD_N = 256
FROZEN_NEW_N = 512
FROZEN_GENERATION_BITS = 768
FROZEN_REPLAY_BITS = 1536

FROZEN_U256 = beta2.FROZEN_U256
REJECTION_FACTOR = Fraction(9, 10)
GAIN_CEILING_FACTOR = 1 - REJECTION_FACTOR
FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND = REJECTION_FACTOR * FROZEN_U256

FROZEN_AUDIT_PAYLOAD_SHA256 = (
    "467d6d819700a87f656e917bcb3e63008d7243fa5f412f75b6eec4ade5a67956"
)
FROZEN_BASIS_COEFFICIENTS_SHA256 = (
    "ec9ea0b1a3cb7439f2e30c552c38e2e5a4620bf5a5ec1285d68901a2a29a7f99"
)
FROZEN_AGGREGATE_SYSTEM_SHA256 = (
    "dd9a571997cc10a32fe332ad2127d58913c021054f058aee90253fdb85422fec"
)
FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256 = (
    "a58186360f701a26fede8872e463e28d57f4e2fd6f06a001137495216b29eea0"
)

TRIAL_COLUMN_IDS = (
    "mobius_log_taper",
    "mobius",
    "mobius_log_taper_squared",
    "squarefree_log_taper",
    "squarefree",
    "squareful",
    "divisible_by_2",
    "divisible_by_3",
)

TRIAL_COLUMN_FORMULAS = (
    "mu(j)*log(512/j)",
    "mu(j)",
    "mu(j)*log(512/j)^2",
    "abs(mu(j))*log(512/j)",
    "abs(mu(j))",
    "1-abs(mu(j))",
    "1_{2|j}",
    "1_{3|j}",
)

AUDIT_OUTCOME = "CANONICAL_EIGHT_COLUMN_TRIAL_SUBSPACE_REJECTED_AT_N512"
TERMINAL_STATEMENT = (
    "F_V > (9/10)*U_256, hence Gamma_V < d_256^2/10; "
    "finite canonical trial subspace rejected"
)
AUDIT_LIMITATION = (
    "This CERTIFIED_FINITE artifact rejects one explicit eight-column trial "
    "subspace at the single 256-to-512 extension. Its verifier performs a "
    "same-backend Arb replay at exactly 1536 bits, not a clean-room replay. "
    "It does not bound all trial subspaces, does not prove an all-scale "
    "recurrence or an infinite-tail floor, and does not prove or disprove "
    "the Riemann Hypothesis; the global status remains UNRESOLVED."
)


class NymanTrialError(ValueError):
    """Raised when the finite trial-subspace audit cannot be generated."""


class NymanTrialVerificationError(NymanTrialError):
    """Raised when a supplied trial-subspace artifact fails reproduction."""


@dataclass(frozen=True)
class TrialAggregateSystem:
    """Interval Gram system for the old span plus the trial columns."""

    gram: tuple[tuple[Any, ...], ...]
    target: tuple[Any, ...]
    old_dimension: int
    trial_dimension: int
    precision_bits: int
    basis_coefficients_sha256: str

    @property
    def dimension(self) -> int:
        return self.old_dimension + self.trial_dimension


def _backend_record() -> dict[str, str]:
    return {
        "python_flint": str(flint.__version__),
        "flint": str(flint.__FLINT_VERSION__),
    }


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _require_payload_hash(
    artifact: Mapping[str, Any],
    error_type: type[NymanTrialError],
    *,
    label: str,
) -> None:
    try:
        supplied = artifact.get("payload_sha256")
        body = {
            key: value
            for key, value in artifact.items()
            if key != "payload_sha256"
        }
        expected = content_sha256(body)
    except (OverflowError, RecursionError, TypeError, ValueError) as exc:
        raise error_type(f"{label} is not canonical JSON") from exc
    if supplied != expected:
        raise error_type(f"{label} payload hash mismatch")


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _load_object(
    source: Mapping[str, Any] | Path,
    error_type: type[NymanTrialError],
    *,
    label: str,
) -> dict[str, Any]:
    # The beta=2 loader rejects duplicate keys, NaN/Infinity, floats, mutable
    # mapping aliases, non-string keys, and all non-JSON values.
    try:
        return beta2._load_object(source, error_type, label=label)
    except error_type:
        raise
    except (OverflowError, RecursionError, TypeError, ValueError) as exc:
        raise error_type(f"{label} is not a finite canonical JSON tree") from exc


def _canonical_integer(
    value: Any,
    name: str,
    *,
    minimum: int | None = None,
    error_type: type[NymanTrialError] = NymanTrialError,
) -> int:
    return beta2._canonical_integer(
        value,
        name,
        minimum=minimum,
        error_type=error_type,
    )


def _source_fraction(
    value: Any,
    name: str,
    *,
    error_type: type[NymanTrialError] = NymanTrialError,
) -> Fraction:
    return beta2._source_fraction(value, name, error_type=error_type)


def _arb_from_fraction(value: Fraction) -> Any:
    return arb(value.numerator) / value.denominator


def _mobius_sieve(limit: int) -> tuple[int, ...]:
    """Return exact Mobius values from zero through ``limit``."""

    if isinstance(limit, bool) or not isinstance(limit, int):
        raise TypeError("limit must be an integer")
    if limit < 1:
        raise ValueError("limit must be positive")
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes: list[int] = []
    composite = [False] * (limit + 1)
    for value in range(2, limit + 1):
        if not composite[value]:
            primes.append(value)
            mu[value] = -1
        for prime in primes:
            product = value * prime
            if product > limit:
                break
            composite[product] = True
            if value % prime == 0:
                mu[product] = 0
                break
            mu[product] = -mu[value]
    return tuple(mu)


def _build_trial_basis(
    old_n: int = FROZEN_OLD_N,
    new_n: int = FROZEN_NEW_N,
) -> tuple[tuple[Any, ...], ...]:
    """Build the ordered new-block coefficient matrix with Arb logs.

    Zero Mobius values and the endpoint ``j=new_n`` are handled before any
    logarithm multiplication or squaring.  This both expresses the exact
    zero coefficients and prevents future backend changes from turning a
    mathematically dead ``0*log`` expression into a non-finite ball.
    """

    if isinstance(old_n, bool) or not isinstance(old_n, int):
        raise TypeError("old_n must be an integer")
    if isinstance(new_n, bool) or not isinstance(new_n, int):
        raise TypeError("new_n must be an integer")
    if old_n < 1 or new_n <= old_n:
        raise ValueError("trial block endpoints are invalid")
    mu = _mobius_sieve(new_n)
    rows: list[tuple[Any, ...]] = []
    for j in range(old_n + 1, new_n + 1):
        mu_j = mu[j]
        abs_mu_j = abs(mu_j)
        if mu_j == 0 or j == new_n:
            mobius_log = arb(0)
            mobius_log_squared = arb(0)
            squarefree_log = arb(0)
        else:
            log_taper = (arb(new_n) / j).log()
            mobius_log = mu_j * log_taper
            mobius_log_squared = mu_j * log_taper * log_taper
            squarefree_log = abs_mu_j * log_taper
        rows.append(
            (
                mobius_log,
                arb(mu_j),
                mobius_log_squared,
                squarefree_log,
                arb(abs_mu_j),
                arb(1 - abs_mu_j),
                arb(1 if j % 2 == 0 else 0),
                arb(1 if j % 3 == 0 else 0),
            )
        )
    return tuple(rows)


def _is_exact_zero(value: Any) -> bool:
    return bool(value.is_exact() and value == 0)


def _basis_manifest(
    basis: Sequence[Sequence[Any]],
    *,
    old_n: int,
    new_n: int,
) -> dict[str, Any]:
    rows = tuple(tuple(row) for row in basis)
    expected_rows = new_n - old_n
    if len(rows) != expected_rows:
        raise NymanTrialError("trial basis row count changed")
    if any(len(row) != len(TRIAL_COLUMN_IDS) for row in rows):
        raise NymanTrialError("trial basis column count changed")
    for row in rows:
        for value in row:
            if not hasattr(value, "is_finite") or not value.is_finite():
                raise NymanTrialError("trial basis contains a non-finite coefficient")

    mu = _mobius_sieve(new_n)
    coefficient_payload = {
        "row_indices": [str(j) for j in range(old_n + 1, new_n + 1)],
        "ordered_column_ids": list(TRIAL_COLUMN_IDS),
        "coefficients": [
            [arb_to_dyadic(value) for value in row]
            for row in rows
        ],
    }
    endpoint_row = rows[-1]
    endpoint_log_columns_zero = all(
        _is_exact_zero(endpoint_row[index]) for index in (0, 2, 3)
    )
    zero_mobius_log_columns_zero = all(
        all(_is_exact_zero(rows[j - old_n - 1][index]) for index in (0, 2, 3))
        for j in range(old_n + 1, new_n + 1)
        if mu[j] == 0
    )
    if not endpoint_log_columns_zero or not zero_mobius_log_columns_zero:
        raise NymanTrialError("trial basis lost its exact-zero log policy")
    return {
        "new_block": {
            "first_j": str(old_n + 1),
            "last_j": str(new_n),
            "row_count": str(expected_rows),
        },
        "ordered_columns": [
            {
                "index": str(index),
                "column_id": column_id,
                "formula": formula,
            }
            for index, (column_id, formula) in enumerate(
                zip(TRIAL_COLUMN_IDS, TRIAL_COLUMN_FORMULAS, strict=True)
            )
        ],
        "coefficient_encoding": "Arb exact-dyadic enclosure per coefficient",
        "coefficient_matrix_content_sha256": content_sha256(coefficient_payload),
        "mobius_values_content_sha256": content_sha256(
            {
                "first_j": str(old_n + 1),
                "last_j": str(new_n),
                "values": [
                    str(mu[j]) for j in range(old_n + 1, new_n + 1)
                ],
            }
        ),
        "checks": {
            "column_order_frozen": True,
            "endpoint_j_equals_512_log_columns_are_exact_zero": (
                endpoint_log_columns_zero and new_n == FROZEN_NEW_N
            ),
            "zero_mobius_rows_skip_log_multiplication": (
                zero_mobius_log_columns_zero
            ),
            "binary_float_used": False,
        },
    }


def _validate_basis_shape(
    basis: Sequence[Sequence[Any]],
    *,
    new_block_size: int,
) -> tuple[tuple[Any, ...], ...]:
    rows = tuple(tuple(row) for row in basis)
    if len(rows) != new_block_size:
        raise ValueError("basis row count must equal the new-block dimension")
    if not rows:
        raise ValueError("basis must have at least one row")
    trial_dimension = len(rows[0])
    if trial_dimension < 1 or any(len(row) != trial_dimension for row in rows):
        raise ValueError("basis rows must have one common positive dimension")
    for row in rows:
        for value in row:
            if not hasattr(value, "is_finite") or not value.is_finite():
                raise TypeError("basis entries must be finite Arb enclosures")
    return rows


def _build_aggregate_system(
    system: core.NaturalSystem,
    basis: Sequence[Sequence[Any]],
    *,
    old_n: int = FROZEN_OLD_N,
    basis_coefficients_sha256: str | None = None,
) -> TrialAggregateSystem:
    """Form ``W^T G W`` and ``W^T b`` without an interval inverse."""

    size = core._validate_system(system)
    if isinstance(old_n, bool) or not isinstance(old_n, int):
        raise TypeError("old_n must be an integer")
    if old_n < 1 or old_n >= size:
        raise ValueError("old_n must split the supplied system")
    rows = _validate_basis_shape(basis, new_block_size=size - old_n)
    trial_dimension = len(rows[0])
    if basis_coefficients_sha256 is None:
        basis_coefficients_sha256 = content_sha256(
            [[arb_to_dyadic(value) for value in row] for row in rows]
        )

    # First compute G times every trial column.  This makes the construction
    # O(size * new_block * trial_dimension), rather than evaluating each
    # trial Gram entry as a fresh double sum.
    transformed = [
        [arb(0) for _ in range(trial_dimension)] for _ in range(size)
    ]
    for row in range(size):
        for trial_column in range(trial_dimension):
            total = arb(0)
            for block_index, coefficient_row in enumerate(rows):
                coefficient = coefficient_row[trial_column]
                if _is_exact_zero(coefficient):
                    continue
                total += system.gram[row][old_n + block_index] * coefficient
            transformed[row][trial_column] = total

    aggregate_dimension = old_n + trial_dimension
    gram = [
        [arb(0) for _ in range(aggregate_dimension)]
        for _ in range(aggregate_dimension)
    ]
    for row in range(old_n):
        for column in range(row, old_n):
            entry = system.gram[row][column]
            gram[row][column] = entry
            gram[column][row] = entry
        for trial_column in range(trial_dimension):
            aggregate_column = old_n + trial_column
            entry = transformed[row][trial_column]
            gram[row][aggregate_column] = entry
            gram[aggregate_column][row] = entry

    for left in range(trial_dimension):
        for right in range(left, trial_dimension):
            entry = arb(0)
            for block_index, coefficient_row in enumerate(rows):
                coefficient = coefficient_row[left]
                if _is_exact_zero(coefficient):
                    continue
                entry += coefficient * transformed[old_n + block_index][right]
            aggregate_left = old_n + left
            aggregate_right = old_n + right
            gram[aggregate_left][aggregate_right] = entry
            gram[aggregate_right][aggregate_left] = entry

    target = list(system.target[:old_n])
    for trial_column in range(trial_dimension):
        entry = arb(0)
        for block_index, coefficient_row in enumerate(rows):
            coefficient = coefficient_row[trial_column]
            if _is_exact_zero(coefficient):
                continue
            entry += system.target[old_n + block_index] * coefficient
        target.append(entry)

    precision_bits = system.precision_bits
    if precision_bits is None:
        raise ValueError("aggregate source system must record its precision")
    return TrialAggregateSystem(
        gram=tuple(tuple(row) for row in gram),
        target=tuple(target),
        old_dimension=old_n,
        trial_dimension=trial_dimension,
        precision_bits=precision_bits,
        basis_coefficients_sha256=basis_coefficients_sha256,
    )


def _aggregate_system_content_sha256(aggregate: TrialAggregateSystem) -> str:
    dimension = aggregate.dimension
    if (
        len(aggregate.gram) != dimension
        or any(len(row) != dimension for row in aggregate.gram)
        or len(aggregate.target) != dimension
    ):
        raise ValueError("aggregate system dimensions are inconsistent")
    return content_sha256(
        {
            "schema": "rh-lab/nyman-trial-aggregate-system/v1",
            "precision_bits": str(aggregate.precision_bits),
            "old_dimension": str(aggregate.old_dimension),
            "trial_dimension": str(aggregate.trial_dimension),
            "basis_coefficients_sha256": aggregate.basis_coefficients_sha256,
            "target": [arb_to_dyadic(value) for value in aggregate.target],
            "gram_upper_triangle": [
                {
                    "row": str(row),
                    "column": str(column),
                    "value": arb_to_dyadic(aggregate.gram[row][column]),
                }
                for row in range(dimension)
                for column in range(row, dimension)
            ],
        }
    )


def _build_augmented_matrix(
    aggregate: TrialAggregateSystem,
    lower_bound: Fraction,
) -> tuple[tuple[Any, ...], ...]:
    if not isinstance(lower_bound, Fraction):
        raise TypeError("lower_bound must be an exact Fraction")
    if lower_bound < 0 or lower_bound >= 1:
        raise ValueError("lower_bound must lie in [0,1)")
    dimension = aggregate.dimension
    matrix = [
        [arb(0) for _ in range(dimension + 1)]
        for _ in range(dimension + 1)
    ]
    for row in range(dimension):
        for column in range(row, dimension):
            entry = aggregate.gram[row][column]
            matrix[row][column] = entry
            matrix[column][row] = entry
        off_diagonal = -aggregate.target[row]
        matrix[row][dimension] = off_diagonal
        matrix[dimension][row] = off_diagonal
    matrix[dimension][dimension] = 1 - _arb_from_fraction(lower_bound)
    return tuple(tuple(row) for row in matrix)


def _lower_certificate_labels(
    old_dimension: int,
    trial_dimension: int,
    lower_bound: Fraction,
) -> tuple[str, str]:
    """Return truthful canonical or parameterized certificate labels."""

    is_canonical = (
        old_dimension == 256
        and trial_dimension == 8
        and lower_bound == Fraction(9, 10) * beta2.FROZEN_U256
    )
    if is_canonical:
        return (
            (
                "rho_1,...,rho_256; eight frozen trial columns in manifest "
                "order; target Schur coordinate"
            ),
            "F_V > (9/10)*U_256",
        )
    return (
        (
            f"{old_dimension} old aggregate columns; {trial_dimension} trial "
            "columns in supplied order; target Schur coordinate"
        ),
        (
            "restricted squared distance > "
            f"{lower_bound.numerator}/{lower_bound.denominator}"
        ),
    )


def _certify_aggregate_lower_bound(
    aggregate: TrialAggregateSystem,
    lower_bound: Fraction,
    *,
    aggregate_system_sha256: str | None = None,
) -> dict[str, Any]:
    if aggregate_system_sha256 is None:
        aggregate_system_sha256 = _aggregate_system_content_sha256(aggregate)
    augmented = _build_augmented_matrix(aggregate, lower_bound)
    ldlt = core.fixed_order_interval_ldlt(augmented)
    certified = ldlt.get("classification") == "POSITIVE_DEFINITE"
    matrix_order, strict_statement = _lower_certificate_labels(
        aggregate.old_dimension,
        aggregate.trial_dimension,
        lower_bound,
    )
    body = {
        "schema": LOWER_CERTIFICATE_SCHEMA,
        "classification": "CERTIFIED_FINITE" if certified else "INCONCLUSIVE",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(aggregate.precision_bits),
        "normalization": core.NYMAN_NORMALIZATION,
        "aggregate_system_content_sha256": aggregate_system_sha256,
        "aggregate_dimension": str(aggregate.dimension),
        "old_dimension": str(aggregate.old_dimension),
        "trial_dimension": str(aggregate.trial_dimension),
        "claimed_lower_bound": _fraction_record(lower_bound),
        "matrix_order": matrix_order,
        "strict_statement": strict_statement,
        "decision": (
            "LOWER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
        ),
        "augmented_ldlt": ldlt,
        "checks": {
            "fixed_order_interval_ldlt_positive": certified,
            "aggregate_gram_positive_as_leading_principal_submatrix": certified,
            "trial_columns_independent_modulo_old_span": certified,
            "approximate_solve_used_as_evidence": False,
        },
        "limitation": (
            "A strict finite lower bound for one aggregate subspace is not an "
            "all-scale Nyman result and does not resolve RH."
        ),
    }
    return _with_payload_hash(body)


def _record_interval(
    record: Any,
    name: str,
    error_type: type[NymanTrialError],
) -> tuple[Fraction, Fraction]:
    if not isinstance(record, Mapping) or set(record) != {
        "dyadic",
        "display",
        "is_exact",
    }:
        raise error_type(f"{name} Arb record changed")
    if not isinstance(record["display"], str) or type(record["is_exact"]) is not bool:
        raise error_type(f"{name} Arb metadata changed")
    dyadic = record["dyadic"]
    if not isinstance(dyadic, Mapping) or set(dyadic) != {
        "mid_mantissa",
        "mid_exponent",
        "radius_mantissa",
        "radius_exponent",
    }:
        raise error_type(f"{name} dyadic enclosure changed")
    midpoint_mantissa = _canonical_integer(
        dyadic["mid_mantissa"],
        f"{name} midpoint mantissa",
        error_type=error_type,
    )
    midpoint_exponent = _canonical_integer(
        dyadic["mid_exponent"],
        f"{name} midpoint exponent",
        error_type=error_type,
    )
    radius_mantissa = _canonical_integer(
        dyadic["radius_mantissa"],
        f"{name} radius mantissa",
        minimum=0,
        error_type=error_type,
    )
    radius_exponent = _canonical_integer(
        dyadic["radius_exponent"],
        f"{name} radius exponent",
        error_type=error_type,
    )

    def dyadic_fraction(mantissa: int, exponent: int) -> Fraction:
        if exponent >= 0:
            return Fraction(mantissa << exponent, 1)
        return Fraction(mantissa, 1 << -exponent)

    midpoint = dyadic_fraction(midpoint_mantissa, midpoint_exponent)
    radius = dyadic_fraction(radius_mantissa, radius_exponent)
    if record["is_exact"] != (radius_mantissa == 0):
        raise error_type(f"{name} exactness flag changed")
    return midpoint - radius, midpoint + radius


def _validate_lower_certificate(
    record: Mapping[str, Any],
    *,
    precision_bits: int,
    aggregate_dimension: int,
    aggregate_system_sha256: str,
    error_type: type[NymanTrialError] = NymanTrialError,
) -> None:
    _require_payload_hash(record, error_type, label="trial lower certificate")
    expected_fields = {
        "schema",
        "classification",
        "hypothesis_status",
        "precision_bits",
        "normalization",
        "aggregate_system_content_sha256",
        "aggregate_dimension",
        "old_dimension",
        "trial_dimension",
        "claimed_lower_bound",
        "matrix_order",
        "strict_statement",
        "decision",
        "augmented_ldlt",
        "checks",
        "limitation",
        "payload_sha256",
    }
    if set(record) != expected_fields:
        raise error_type("trial lower certificate fields changed")
    matrix_order, strict_statement = _lower_certificate_labels(
        FROZEN_OLD_N,
        len(TRIAL_COLUMN_IDS),
        FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND,
    )
    expected_values = {
        "schema": LOWER_CERTIFICATE_SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(precision_bits),
        "normalization": core.NYMAN_NORMALIZATION,
        "aggregate_system_content_sha256": aggregate_system_sha256,
        "aggregate_dimension": str(aggregate_dimension),
        "old_dimension": str(FROZEN_OLD_N),
        "trial_dimension": str(len(TRIAL_COLUMN_IDS)),
        "matrix_order": matrix_order,
        "strict_statement": strict_statement,
        "decision": "LOWER_BOUND_CERTIFIED",
        "limitation": (
            "A strict finite lower bound for one aggregate subspace is not an "
            "all-scale Nyman result and does not resolve RH."
        ),
    }
    for key, expected in expected_values.items():
        if record.get(key) != expected:
            raise error_type(f"trial lower certificate {key} changed")
    if _source_fraction(
        record.get("claimed_lower_bound"),
        "trial lower bound",
        error_type=error_type,
    ) != FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND:
        raise error_type("trial lower certificate endpoint changed")
    expected_checks = {
        "fixed_order_interval_ldlt_positive": True,
        "aggregate_gram_positive_as_leading_principal_submatrix": True,
        "trial_columns_independent_modulo_old_span": True,
        "approximate_solve_used_as_evidence": False,
    }
    if record.get("checks") != expected_checks:
        raise error_type("trial lower certificate checks changed")

    ldlt = record.get("augmented_ldlt")
    order_size = aggregate_dimension + 1
    if not isinstance(ldlt, Mapping) or set(ldlt) != {
        "classification",
        "fixed_order",
        "failed_pivot_index",
        "pivots",
    }:
        raise error_type("trial lower certificate LDL transcript changed")
    if (
        ldlt.get("classification") != "POSITIVE_DEFINITE"
        or ldlt.get("failed_pivot_index") is not None
        or ldlt.get("fixed_order") != [str(index) for index in range(order_size)]
    ):
        raise error_type("trial lower certificate LDL decision changed")
    pivots = ldlt.get("pivots")
    if not isinstance(pivots, list) or len(pivots) != order_size:
        raise error_type("trial lower certificate pivot count changed")
    for index, pivot in enumerate(pivots):
        if (
            not isinstance(pivot, Mapping)
            or set(pivot) != {"index", "value"}
            or pivot.get("index") != str(index)
        ):
            raise error_type("trial lower certificate pivot record changed")
        lower, _ = _record_interval(
            pivot.get("value"),
            f"trial lower pivot {index}",
            error_type,
        )
        if lower <= 0:
            raise error_type("trial lower certificate stores a nonpositive pivot")


def _exact_threshold_record() -> dict[str, Any]:
    lower_bound = REJECTION_FACTOR * FROZEN_U256
    if lower_bound != FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND:
        raise NymanTrialError("restricted-distance threshold changed")
    return {
        "source_upper_bound_u256": _fraction_record(FROZEN_U256),
        "rejection_factor": _fraction_record(REJECTION_FACTOR),
        "gain_ceiling_factor": _fraction_record(GAIN_CEILING_FACTOR),
        "restricted_distance_lower_bound": _fraction_record(lower_bound),
        "checks": {
            "lower_bound_equals_9_over_10_times_u256": True,
            "source_certifies_d256_squared_at_most_u256": True,
            "all_comparisons_exact_rational": True,
        },
    }


def _theorem_bridge() -> dict[str, Any]:
    return {
        "definitions": {
            "aggregate_space": (
                "span(rho_1,...,rho_256) plus the eight frozen new-block "
                "trial functions"
            ),
            "restricted_squared_distance": (
                "F_V=dist(chi,aggregate_space)^2"
            ),
            "captured_gain": "Gamma_V=d_256^2-F_V",
        },
        "certificate_implication": [
            (
                "fixed-order interval LDL positivity of "
                "[[G_V,-b_V],[-b_V^T,1-(9/10)U_256]] gives "
                "F_V>(9/10)U_256"
            ),
            "the frozen source endpoint gives d_256^2<=U_256",
            (
                "Gamma_V=d_256^2-F_V"
                "<d_256^2-(9/10)U_256<=d_256^2/10"
            ),
        ],
        "rank_implication": (
            "positive definiteness of the augmented matrix makes its "
            "aggregate Gram leading principal submatrix positive definite; "
            "therefore the eight trial columns are independent modulo the "
            "old span"
        ),
        "finite_rejection": (
            "this canonical eight-column space captures strictly less than "
            "one tenth of the N=256 residual energy"
        ),
        "all_scale_trial_subspace_lemma_proved": False,
        "all_trial_subspaces_rejected": False,
        "resolves_rh": False,
    }


def _validate_source_summary(
    summary: Mapping[str, Any],
    checkpoint_dir: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        return beta2._validate_source_summary(summary, checkpoint_dir)
    except beta2.NymanBeta2Error as exc:
        raise NymanTrialError("source Nyman summary verification failed") from exc


def _source_record(
    verification: Mapping[str, Any],
    source_cell: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "summary": {
            "schema": SUMMARY_SCHEMA,
            "payload_sha256": beta2.FROZEN_SUMMARY_PAYLOAD_SHA256,
            "plan_sha256": beta2.FROZEN_PLAN_SHA256,
            "index_payload_sha256": beta2.FROZEN_INDEX_PAYLOAD_SHA256,
            "verification": dict(verification),
        },
        "n256": {
            "cell_id": "n-0256",
            "cell_payload_sha256": beta2.FROZEN_N256_CELL_PAYLOAD_SHA256,
            "cell_contract_sha256": beta2.FROZEN_N256_CELL_CONTRACT_SHA256,
            "candidate_sha256": beta2.FROZEN_N256_CANDIDATE_SHA256,
            "lower_certificate_payload_sha256": (
                beta2.FROZEN_N256_LOWER_CERTIFICATE_SHA256
            ),
            "upper_certificate_payload_sha256": (
                beta2.FROZEN_N256_UPPER_CERTIFICATE_SHA256
            ),
            "numeric_prefix_sha256": beta2.FROZEN_N256_PREFIX_KERNEL_SHA256,
            "lower_bound": _fraction_record(beta2.FROZEN_L256),
            "upper_bound": _fraction_record(FROZEN_U256),
            "terminal_statement": source_cell.get("terminal_statement"),
        },
        "source_verification_scope": (
            "structural regeneration of frozen v1; the frozen certified "
            "N=256 upper endpoint is the theorem-bridge input"
        ),
    }


def _check_optional_freeze(
    actual: str,
    frozen: str,
    label: str,
    error_type: type[NymanTrialError],
) -> None:
    if frozen and actual != frozen:
        raise error_type(f"frozen {label} hash changed")


def _is_sha256(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _build_trial_evidence(precision_bits: int) -> dict[str, Any]:
    if (
        isinstance(precision_bits, bool)
        or not isinstance(precision_bits, int)
        or precision_bits < 96
    ):
        raise ValueError("precision_bits must be an integer of at least 96")
    try:
        with beta2._clean_precision(precision_bits):
            system = core.build_natural_system(
                tuple(range(1, FROZEN_NEW_N + 1))
            )
            kernel_sha256 = core.natural_system_content_sha256(system)
            prefix_256_sha256 = beta2._prefix_kernel_sha256(
                system, FROZEN_OLD_N
            )
            prefix_512_sha256 = beta2._prefix_kernel_sha256(
                system, FROZEN_NEW_N
            )
            if precision_bits == FROZEN_GENERATION_BITS:
                if kernel_sha256 != beta2.FROZEN_GENERATION_KERNEL_SHA256:
                    raise NymanTrialError(
                        "generation max-512 kernel changed"
                    )
                if (
                    prefix_256_sha256
                    != beta2.FROZEN_N256_PREFIX_KERNEL_SHA256
                ):
                    raise NymanTrialError(
                        "generation N=256 numeric prefix changed"
                    )
                if (
                    prefix_512_sha256
                    != beta2.FROZEN_GENERATION_PREFIX_512_SHA256
                ):
                    raise NymanTrialError(
                        "generation N=512 numeric prefix changed"
                    )

            basis = _build_trial_basis()
            basis_manifest = _basis_manifest(
                basis,
                old_n=FROZEN_OLD_N,
                new_n=FROZEN_NEW_N,
            )
            basis_sha256 = basis_manifest[
                "coefficient_matrix_content_sha256"
            ]
            aggregate = _build_aggregate_system(
                system,
                basis,
                old_n=FROZEN_OLD_N,
                basis_coefficients_sha256=basis_sha256,
            )
            aggregate_sha256 = _aggregate_system_content_sha256(aggregate)
            lower_certificate = _certify_aggregate_lower_bound(
                aggregate,
                FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND,
                aggregate_system_sha256=aggregate_sha256,
            )
    except NymanTrialError:
        raise
    except Exception as exc:  # pragma: no cover - backend-specific failures
        raise NymanTrialError(
            "canonical trial-subspace interval calculation failed"
        ) from exc

    _validate_lower_certificate(
        lower_certificate,
        precision_bits=precision_bits,
        aggregate_dimension=FROZEN_OLD_N + len(TRIAL_COLUMN_IDS),
        aggregate_system_sha256=aggregate_sha256,
    )
    if precision_bits == FROZEN_GENERATION_BITS:
        _check_optional_freeze(
            basis_sha256,
            FROZEN_BASIS_COEFFICIENTS_SHA256,
            "basis coefficient matrix",
            NymanTrialError,
        )
        _check_optional_freeze(
            aggregate_sha256,
            FROZEN_AGGREGATE_SYSTEM_SHA256,
            "aggregate system",
            NymanTrialError,
        )
        _check_optional_freeze(
            lower_certificate["payload_sha256"],
            FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256,
            "lower certificate payload",
            NymanTrialError,
        )
    return {
        "precision_bits": str(precision_bits),
        "backend": _backend_record(),
        "max_512_kernel": {
            "core_system_content_sha256": kernel_sha256,
            "prefix_256_numeric_sha256": prefix_256_sha256,
            "prefix_512_numeric_sha256": prefix_512_sha256,
            "prefix_nesting_check": True,
        },
        "trial_basis": basis_manifest,
        "aggregate_system": {
            "content_sha256": aggregate_sha256,
            "dimension": str(aggregate.dimension),
            "old_dimension": str(aggregate.old_dimension),
            "trial_dimension": str(aggregate.trial_dimension),
            "ordering": (
                "rho_1,...,rho_256 followed by the eight manifest columns"
            ),
        },
        "restricted_distance_lower_certificate": lower_certificate,
    }


def _derive_audit(
    summary: Mapping[str, Any],
    checkpoint_dir: Path,
) -> dict[str, Any]:
    verification, source_cell = _validate_source_summary(summary, checkpoint_dir)
    evidence = _build_trial_evidence(FROZEN_GENERATION_BITS)
    body = {
        "schema": AUDIT_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "policy": {
            "old_n": str(FROZEN_OLD_N),
            "new_n": str(FROZEN_NEW_N),
            "generation_precision_bits": str(FROZEN_GENERATION_BITS),
            "replay_precision_bits": str(FROZEN_REPLAY_BITS),
            "trial_column_order": list(TRIAL_COLUMN_IDS),
            "failed_generation_behavior": (
                "fail closed; no adaptive precision, basis, or matrix ordering"
            ),
        },
        "source": _source_record(verification, source_cell),
        "generation": evidence,
        "exact_thresholds": _exact_threshold_record(),
        "theorem_bridge": _theorem_bridge(),
        "audit_outcome": AUDIT_OUTCOME,
        "terminal_statement": TERMINAL_STATEMENT,
        "limitation": AUDIT_LIMITATION,
    }
    return _with_payload_hash(body)


def generate_nyman_trial_subspace_audit(
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Generate the canonical finite trial-subspace rejection artifact."""

    source_summary = _load_object(
        summary,
        NymanTrialError,
        label="source Nyman summary",
    )
    artifact = _derive_audit(source_summary, checkpoint_dir)
    _check_optional_freeze(
        artifact["payload_sha256"],
        FROZEN_AUDIT_PAYLOAD_SHA256,
        "audit payload",
        NymanTrialError,
    )
    return artifact


def _validate_source_binding(source: Any) -> None:
    if not isinstance(source, Mapping) or set(source) != {
        "summary",
        "n256",
        "source_verification_scope",
    }:
        raise NymanTrialVerificationError("audit source binding changed")
    summary = source.get("summary")
    n256 = source.get("n256")
    if not isinstance(summary, Mapping) or not isinstance(n256, Mapping):
        raise NymanTrialVerificationError("audit source records changed")
    if set(summary) != {
        "schema",
        "payload_sha256",
        "plan_sha256",
        "index_payload_sha256",
        "verification",
    }:
        raise NymanTrialVerificationError("audit summary record fields changed")
    if set(n256) != {
        "cell_id",
        "cell_payload_sha256",
        "cell_contract_sha256",
        "candidate_sha256",
        "lower_certificate_payload_sha256",
        "upper_certificate_payload_sha256",
        "numeric_prefix_sha256",
        "lower_bound",
        "upper_bound",
        "terminal_statement",
    }:
        raise NymanTrialVerificationError("audit N=256 record fields changed")
    if (
        summary.get("schema") != SUMMARY_SCHEMA
        or summary.get("payload_sha256")
        != beta2.FROZEN_SUMMARY_PAYLOAD_SHA256
        or summary.get("plan_sha256") != beta2.FROZEN_PLAN_SHA256
        or summary.get("index_payload_sha256")
        != beta2.FROZEN_INDEX_PAYLOAD_SHA256
    ):
        raise NymanTrialVerificationError("audit summary binding changed")
    if (
        n256.get("cell_id") != "n-0256"
        or n256.get("cell_payload_sha256")
        != beta2.FROZEN_N256_CELL_PAYLOAD_SHA256
        or n256.get("cell_contract_sha256")
        != beta2.FROZEN_N256_CELL_CONTRACT_SHA256
        or n256.get("candidate_sha256")
        != beta2.FROZEN_N256_CANDIDATE_SHA256
        or n256.get("lower_certificate_payload_sha256")
        != beta2.FROZEN_N256_LOWER_CERTIFICATE_SHA256
        or n256.get("upper_certificate_payload_sha256")
        != beta2.FROZEN_N256_UPPER_CERTIFICATE_SHA256
        or n256.get("numeric_prefix_sha256")
        != beta2.FROZEN_N256_PREFIX_KERNEL_SHA256
        or n256.get("lower_bound") != _fraction_record(beta2.FROZEN_L256)
        or n256.get("upper_bound") != _fraction_record(FROZEN_U256)
    ):
        raise NymanTrialVerificationError("audit N=256 binding changed")
    if source.get("source_verification_scope") != (
        "structural regeneration of frozen v1; the frozen certified "
        "N=256 upper endpoint is the theorem-bridge input"
    ):
        raise NymanTrialVerificationError("audit source scope changed")


def _preflight_audit(supplied: Mapping[str, Any]) -> None:
    """Reject structural and hash mutations before rebuilding a 512 kernel."""

    _require_payload_hash(
        supplied,
        NymanTrialVerificationError,
        label="Nyman trial-subspace audit",
    )
    _check_optional_freeze(
        str(supplied.get("payload_sha256")),
        FROZEN_AUDIT_PAYLOAD_SHA256,
        "audit payload",
        NymanTrialVerificationError,
    )
    expected_fields = {
        "schema",
        "audit_id",
        "classification",
        "hypothesis_status",
        "policy",
        "source",
        "generation",
        "exact_thresholds",
        "theorem_bridge",
        "audit_outcome",
        "terminal_statement",
        "limitation",
        "payload_sha256",
    }
    if set(supplied) != expected_fields:
        raise NymanTrialVerificationError("audit fields changed")
    expected_preamble = {
        "schema": AUDIT_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": AUDIT_OUTCOME,
        "terminal_statement": TERMINAL_STATEMENT,
        "limitation": AUDIT_LIMITATION,
    }
    for key, expected in expected_preamble.items():
        if supplied.get(key) != expected:
            raise NymanTrialVerificationError(f"audit {key} changed")
    expected_policy = {
        "old_n": str(FROZEN_OLD_N),
        "new_n": str(FROZEN_NEW_N),
        "generation_precision_bits": str(FROZEN_GENERATION_BITS),
        "replay_precision_bits": str(FROZEN_REPLAY_BITS),
        "trial_column_order": list(TRIAL_COLUMN_IDS),
        "failed_generation_behavior": (
            "fail closed; no adaptive precision, basis, or matrix ordering"
        ),
    }
    if supplied.get("policy") != expected_policy:
        raise NymanTrialVerificationError("audit policy changed")
    _validate_source_binding(supplied.get("source"))
    if supplied.get("exact_thresholds") != _exact_threshold_record():
        raise NymanTrialVerificationError("audit exact thresholds changed")
    if supplied.get("theorem_bridge") != _theorem_bridge():
        raise NymanTrialVerificationError("audit theorem bridge changed")

    generation = supplied.get("generation")
    if not isinstance(generation, Mapping) or set(generation) != {
        "precision_bits",
        "backend",
        "max_512_kernel",
        "trial_basis",
        "aggregate_system",
        "restricted_distance_lower_certificate",
    }:
        raise NymanTrialVerificationError("audit generation evidence changed")
    if generation.get("precision_bits") != str(FROZEN_GENERATION_BITS):
        raise NymanTrialVerificationError("audit generation precision changed")
    if generation.get("backend") != _backend_record():
        raise NymanTrialVerificationError("audit backend contract changed")
    kernel = generation.get("max_512_kernel")
    if not isinstance(kernel, Mapping) or kernel != {
        "core_system_content_sha256": beta2.FROZEN_GENERATION_KERNEL_SHA256,
        "prefix_256_numeric_sha256": beta2.FROZEN_N256_PREFIX_KERNEL_SHA256,
        "prefix_512_numeric_sha256": beta2.FROZEN_GENERATION_PREFIX_512_SHA256,
        "prefix_nesting_check": True,
    }:
        raise NymanTrialVerificationError("audit max-512 kernel binding changed")
    basis = generation.get("trial_basis")
    aggregate = generation.get("aggregate_system")
    if not isinstance(basis, Mapping) or not isinstance(aggregate, Mapping):
        raise NymanTrialVerificationError("audit trial construction changed")
    basis_sha256 = basis.get("coefficient_matrix_content_sha256")
    aggregate_sha256 = aggregate.get("content_sha256")
    if not _is_sha256(basis_sha256):
        raise NymanTrialVerificationError("audit basis hash changed")
    if not _is_sha256(aggregate_sha256):
        raise NymanTrialVerificationError("audit aggregate hash changed")
    expected_basis_fields = {
        "new_block",
        "ordered_columns",
        "coefficient_encoding",
        "coefficient_matrix_content_sha256",
        "mobius_values_content_sha256",
        "checks",
    }
    if set(basis) != expected_basis_fields:
        raise NymanTrialVerificationError("audit basis manifest fields changed")
    expected_columns = [
        {
            "index": str(index),
            "column_id": column_id,
            "formula": formula,
        }
        for index, (column_id, formula) in enumerate(
            zip(TRIAL_COLUMN_IDS, TRIAL_COLUMN_FORMULAS, strict=True)
        )
    ]
    mu = _mobius_sieve(FROZEN_NEW_N)
    expected_mobius_sha256 = content_sha256(
        {
            "first_j": str(FROZEN_OLD_N + 1),
            "last_j": str(FROZEN_NEW_N),
            "values": [
                str(mu[j])
                for j in range(FROZEN_OLD_N + 1, FROZEN_NEW_N + 1)
            ],
        }
    )
    if (
        basis.get("new_block")
        != {
            "first_j": str(FROZEN_OLD_N + 1),
            "last_j": str(FROZEN_NEW_N),
            "row_count": str(FROZEN_NEW_N - FROZEN_OLD_N),
        }
        or basis.get("ordered_columns") != expected_columns
        or basis.get("coefficient_encoding")
        != "Arb exact-dyadic enclosure per coefficient"
        or basis.get("mobius_values_content_sha256")
        != expected_mobius_sha256
        or basis.get("checks")
        != {
            "column_order_frozen": True,
            "endpoint_j_equals_512_log_columns_are_exact_zero": True,
            "zero_mobius_rows_skip_log_multiplication": True,
            "binary_float_used": False,
        }
    ):
        raise NymanTrialVerificationError("audit basis manifest changed")
    _check_optional_freeze(
        basis_sha256,
        FROZEN_BASIS_COEFFICIENTS_SHA256,
        "basis coefficient matrix",
        NymanTrialVerificationError,
    )
    _check_optional_freeze(
        aggregate_sha256,
        FROZEN_AGGREGATE_SYSTEM_SHA256,
        "aggregate system",
        NymanTrialVerificationError,
    )
    expected_aggregate = {
        "content_sha256": aggregate_sha256,
        "dimension": str(FROZEN_OLD_N + len(TRIAL_COLUMN_IDS)),
        "old_dimension": str(FROZEN_OLD_N),
        "trial_dimension": str(len(TRIAL_COLUMN_IDS)),
        "ordering": (
            "rho_1,...,rho_256 followed by the eight manifest columns"
        ),
    }
    if aggregate != expected_aggregate:
        raise NymanTrialVerificationError("audit aggregate record changed")
    certificate = generation.get("restricted_distance_lower_certificate")
    if not isinstance(certificate, Mapping):
        raise NymanTrialVerificationError("audit lower certificate is missing")
    _validate_lower_certificate(
        certificate,
        precision_bits=FROZEN_GENERATION_BITS,
        aggregate_dimension=FROZEN_OLD_N + len(TRIAL_COLUMN_IDS),
        aggregate_system_sha256=aggregate_sha256,
        error_type=NymanTrialVerificationError,
    )
    _check_optional_freeze(
        str(certificate.get("payload_sha256")),
        FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256,
        "lower certificate payload",
        NymanTrialVerificationError,
    )


def verify_nyman_trial_subspace_audit(
    artifact: Mapping[str, Any] | Path,
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
    *,
    replay_precision_bits: int = FROZEN_REPLAY_BITS,
) -> dict[str, Any]:
    """Regenerate at 768 bits and independently replay the LDL proof at 1536."""

    if (
        isinstance(replay_precision_bits, bool)
        or not isinstance(replay_precision_bits, int)
        or replay_precision_bits != FROZEN_REPLAY_BITS
    ):
        raise NymanTrialVerificationError(
            f"replay precision must equal the frozen {FROZEN_REPLAY_BITS} bits"
        )
    supplied = _load_object(
        artifact,
        NymanTrialVerificationError,
        label="Nyman trial-subspace audit",
    )
    _preflight_audit(supplied)
    source_summary = _load_object(
        summary,
        NymanTrialVerificationError,
        label="source Nyman summary",
    )
    try:
        expected = _derive_audit(source_summary, checkpoint_dir)
    except NymanTrialError as exc:
        raise NymanTrialVerificationError(str(exc)) from exc
    if supplied != expected:
        raise NymanTrialVerificationError(
            "Nyman trial-subspace audit does not canonically regenerate"
        )

    try:
        replay = _build_trial_evidence(replay_precision_bits)
    except NymanTrialError as exc:
        raise NymanTrialVerificationError(
            "higher-precision trial-subspace replay failed"
        ) from exc
    replay_aggregate = replay["aggregate_system"]
    replay_certificate = replay["restricted_distance_lower_certificate"]
    _validate_lower_certificate(
        replay_certificate,
        precision_bits=replay_precision_bits,
        aggregate_dimension=FROZEN_OLD_N + len(TRIAL_COLUMN_IDS),
        aggregate_system_sha256=replay_aggregate["content_sha256"],
        error_type=NymanTrialVerificationError,
    )
    return {
        "classification": (
            "REPRODUCED_CERTIFIED_FINITE_NYMAN_TRIAL_SUBSPACE_REJECTION"
        ),
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": AUDIT_OUTCOME,
        "artifact_payload_sha256": supplied["payload_sha256"],
        "generation_precision_bits": str(FROZEN_GENERATION_BITS),
        "replay_precision_bits": str(replay_precision_bits),
        "generation_kernel_sha256": supplied["generation"]["max_512_kernel"][
            "core_system_content_sha256"
        ],
        "replay_kernel_sha256": replay["max_512_kernel"][
            "core_system_content_sha256"
        ],
        "generation_basis_coefficients_sha256": supplied["generation"][
            "trial_basis"
        ]["coefficient_matrix_content_sha256"],
        "replay_basis_coefficients_sha256": replay["trial_basis"][
            "coefficient_matrix_content_sha256"
        ],
        "generation_aggregate_system_sha256": supplied["generation"][
            "aggregate_system"
        ]["content_sha256"],
        "replay_aggregate_system_sha256": replay_aggregate["content_sha256"],
        "generation_fixed_order_ldlt_positive": True,
        "replay_fixed_order_ldlt_positive": True,
        "restricted_distance_lower_bound_replayed": True,
        "exact_theorem_bridge_rechecked": True,
        "same_backend_replay_only": True,
    }


__all__ = [
    "AUDIT_LIMITATION",
    "AUDIT_SCHEMA",
    "FROZEN_AGGREGATE_SYSTEM_SHA256",
    "FROZEN_AUDIT_ID",
    "FROZEN_AUDIT_PAYLOAD_SHA256",
    "FROZEN_BASIS_COEFFICIENTS_SHA256",
    "FROZEN_GENERATION_BITS",
    "FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256",
    "FROZEN_NEW_N",
    "FROZEN_OLD_N",
    "FROZEN_REPLAY_BITS",
    "FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND",
    "FROZEN_U256",
    "GAIN_CEILING_FACTOR",
    "LOWER_CERTIFICATE_SCHEMA",
    "NymanTrialError",
    "NymanTrialVerificationError",
    "REJECTION_FACTOR",
    "TRIAL_COLUMN_FORMULAS",
    "TRIAL_COLUMN_IDS",
    "TrialAggregateSystem",
    "generate_nyman_trial_subspace_audit",
    "verify_nyman_trial_subspace_audit",
]
