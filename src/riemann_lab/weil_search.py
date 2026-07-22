"""Deterministic, checkpointed exploratory search over finite Weil matrices.

This module deliberately sits outside the claim-bearing certificate path in
``riemann_lab.weil``.  It may locate useful finite positive examples or an
integer negative-witness candidate, but every artifact it emits remains
``EXPLORATORY`` and the global Riemann Hypothesis status remains
``UNRESOLVED``.

The only terminal positive gate is fixed-order interval ``LDL^T``.  Approximate
eigenvectors and Rump spectra are search diagnostics, never proof inputs.  A
negative result is accepted only for an exact primitive integer vector whose
full ``P-R-S`` Rayleigh quotient has a strictly negative Arb upper bound at two
precisions.  Even then it is quarantined for independent mathematical and
backend audit rather than promoted to a disproof.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
import json
from math import gcd
import os
from pathlib import Path
import tempfile
from typing import Any, Callable, Iterable, Mapping, Sequence

import flint
from flint import arb, arb_mat, ctx

from .artifacts import content_sha256
from .balls import arb_from_dyadic, arb_record, arb_to_dyadic
from .weil import (
    _prime_transcript,
    _serialize_entries,
    _symmetric_components,
    certify_positive_ldlt,
    certify_rump_spectrum,
    evaluate_integer_rayleigh as _evaluate_integer_rayleigh,
)


PLAN_SCHEMA = "rh-lab/weil-search-plan/v1"
RUN_SCHEMA = "rh-lab/weil-search-run/v1"
ATTEMPT_SCHEMA = "rh-lab/weil-search-attempt/v1"
CELL_SCHEMA = "rh-lab/weil-search-cell/v1"
INDEX_SCHEMA = "rh-lab/weil-search-index/v1"
ENGINE_ID = "finite-weil-grid-search-v1"

EXPLORATORY_LIMITATION = (
    "This is a bounded search over finitely many cutoffs and Fourier modes. "
    "Positive cells do not prove RH, an inconclusive cell is not evidence of "
    "a negative direction, and any interval-negative candidate remains "
    "quarantined until independent sign, normalization, source-equation, and "
    "backend audits reproduce it."
)


class WeilSearchError(ValueError):
    """Base class for invalid search inputs or artifacts."""


class WeilSearchPlanError(WeilSearchError):
    """Raised when a search plan is malformed or incompatible."""


class WeilSearchCheckpointError(WeilSearchError):
    """Raised when checkpoint state is missing, corrupt, or incompatible."""


class WeilSearchVerificationError(WeilSearchError):
    """Raised when an exploratory search artifact fails replay."""


def _parse_integer(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise WeilSearchPlanError(f"{name} must be an integer")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            parsed = int(value)
        except ValueError as exc:
            raise WeilSearchPlanError(f"{name} must be an integer") from exc
        if value != str(parsed):
            raise WeilSearchPlanError(f"{name} must use canonical decimal syntax")
        return parsed
    raise WeilSearchPlanError(f"{name} must be an integer")


def _validated_integer(value: Any, name: str, minimum: int) -> int:
    parsed = _parse_integer(value, name)
    if parsed < minimum:
        raise WeilSearchPlanError(f"{name} must be at least {minimum}")
    return parsed


@dataclass(frozen=True)
class WeilSearchCell:
    """One exact rational cutoff and symmetric Fourier degree."""

    cutoff_numerator: int
    cutoff_denominator: int
    degree: int

    def __post_init__(self) -> None:
        numerator = _validated_integer(
            self.cutoff_numerator, "cutoff numerator", 1
        )
        denominator = _validated_integer(
            self.cutoff_denominator, "cutoff denominator", 1
        )
        degree = _validated_integer(self.degree, "degree", 0)
        if numerator <= denominator:
            raise WeilSearchPlanError("cutoff must be strictly greater than one")
        common = gcd(numerator, denominator)
        object.__setattr__(self, "cutoff_numerator", numerator // common)
        object.__setattr__(self, "cutoff_denominator", denominator // common)
        object.__setattr__(self, "degree", degree)

    @property
    def cell_id(self) -> str:
        return (
            f"c-{self.cutoff_numerator}-over-{self.cutoff_denominator}"
            f"-n-{self.degree}"
        )

    @property
    def sort_key(self) -> tuple[Fraction, int, int, int]:
        return (
            Fraction(self.cutoff_numerator, self.cutoff_denominator),
            self.degree,
            self.cutoff_numerator,
            self.cutoff_denominator,
        )

    def to_record(self) -> dict[str, Any]:
        record = {
            "cutoff_c": {
                "numerator": str(self.cutoff_numerator),
                "denominator": str(self.cutoff_denominator),
            },
            "degree": str(self.degree),
            "dimension": str(2 * self.degree + 1),
        }
        return {
            "cell_id": self.cell_id,
            "input_sha256": content_sha256(record),
            **record,
        }


@dataclass(frozen=True)
class WeilSearchPolicy:
    """Frozen deterministic precision and rationalization schedules."""

    attempt_bits: tuple[int, ...] = (96, 192, 384)
    confirmation_bits: int = 768
    eigenpair_count: int = 3
    scale_bits: tuple[int, ...] = (8, 12, 16, 24, 32, 48)

    def __post_init__(self) -> None:
        attempt_bits = tuple(
            _validated_integer(value, "attempt precision", 96)
            for value in self.attempt_bits
        )
        if not attempt_bits:
            raise WeilSearchPlanError("attempt precision schedule cannot be empty")
        if any(left >= right for left, right in zip(attempt_bits, attempt_bits[1:])):
            raise WeilSearchPlanError(
                "attempt precision schedule must be strictly increasing"
            )
        confirmation_bits = _validated_integer(
            self.confirmation_bits, "confirmation precision", 96
        )
        if confirmation_bits <= attempt_bits[-1]:
            raise WeilSearchPlanError(
                "confirmation precision must exceed every search precision"
            )
        eigenpair_count = _validated_integer(
            self.eigenpair_count, "eigenpair count", 1
        )
        scale_bits = tuple(
            _validated_integer(value, "rational scale bits", 1)
            for value in self.scale_bits
        )
        if not scale_bits:
            raise WeilSearchPlanError("rational scale schedule cannot be empty")
        if any(left >= right for left, right in zip(scale_bits, scale_bits[1:])):
            raise WeilSearchPlanError(
                "rational scale schedule must be strictly increasing"
            )
        object.__setattr__(self, "attempt_bits", attempt_bits)
        object.__setattr__(self, "confirmation_bits", confirmation_bits)
        object.__setattr__(self, "eigenpair_count", eigenpair_count)
        object.__setattr__(self, "scale_bits", scale_bits)

    def to_record(self) -> dict[str, Any]:
        return {
            "attempt_bits": [str(value) for value in self.attempt_bits],
            "confirmation_bits": str(self.confirmation_bits),
            "positive_terminal_gate": "fixed-order-interval-ldlt",
            "candidate_confirmation": (
                "same primitive integer witness is interval-negative at a "
                "strictly higher precision"
            ),
            "witness_policy": {
                "eigenpair_count": str(self.eigenpair_count),
                "eigensolver_role": "approximate-untrusted-hint-only",
                "parts": ["real", "imaginary"],
                "projections": ["even", "odd", "raw"],
                "scale_bits": [str(value) for value in self.scale_bits],
                "rounding": "exact-rational-round-ties-to-even",
                "normalization": "primitive-gcd-and-first-nonzero-positive",
            },
        }


def _backend_record() -> dict[str, str]:
    return {
        "python_flint": str(flint.__version__),
        "flint": str(flint.__FLINT_VERSION__),
    }


def _policy_from_plan(raw: Mapping[str, Any]) -> WeilSearchPolicy:
    precision = raw.get("precision_policy", {})
    witness = raw.get("witness_policy", {})
    if not isinstance(precision, Mapping) or not isinstance(witness, Mapping):
        raise WeilSearchPlanError("precision and witness policies must be objects")
    return WeilSearchPolicy(
        attempt_bits=tuple(precision.get("attempt_bits", (96, 192, 384))),
        confirmation_bits=precision.get("confirmation_bits", 768),
        eigenpair_count=witness.get("eigenpair_count", 3),
        scale_bits=tuple(witness.get("scale_bits", (8, 12, 16, 24, 32, 48))),
    )


def canonicalize_plan(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and canonicalize a v1 search plan, including its hash."""

    if not isinstance(raw, Mapping):
        raise WeilSearchPlanError("search plan must be an object")
    if raw.get("schema") != PLAN_SCHEMA:
        raise WeilSearchPlanError("unexpected search plan schema")
    if raw.get("classification", "EXPLORATORY") != "EXPLORATORY":
        raise WeilSearchPlanError("search plans must remain EXPLORATORY")

    raw_cutoffs = raw.get("cutoffs")
    raw_degrees = raw.get("degrees")
    if not isinstance(raw_cutoffs, list) or not raw_cutoffs:
        raise WeilSearchPlanError("cutoffs must be a nonempty list")
    if not isinstance(raw_degrees, list) or not raw_degrees:
        raise WeilSearchPlanError("degrees must be a nonempty list")

    cutoffs: list[tuple[int, int]] = []
    seen_cutoffs: set[tuple[int, int]] = set()
    for index, record in enumerate(raw_cutoffs):
        if not isinstance(record, Mapping):
            raise WeilSearchPlanError(f"cutoff {index} must be an object")
        cell = WeilSearchCell(
            record.get("numerator"), record.get("denominator"), 0
        )
        cutoff = (cell.cutoff_numerator, cell.cutoff_denominator)
        if cutoff in seen_cutoffs:
            raise WeilSearchPlanError("duplicate reduced cutoff")
        seen_cutoffs.add(cutoff)
        cutoffs.append(cutoff)
    cutoffs.sort(key=lambda value: (Fraction(*value), value[0], value[1]))

    degrees: list[int] = []
    seen_degrees: set[int] = set()
    for index, value in enumerate(raw_degrees):
        degree = _validated_integer(value, f"degree {index}", 0)
        if degree in seen_degrees:
            raise WeilSearchPlanError("duplicate degree")
        seen_degrees.add(degree)
        degrees.append(degree)
    degrees.sort()

    policy = _policy_from_plan(raw)
    backend = raw.get("backend_contract", _backend_record())
    if not isinstance(backend, Mapping):
        raise WeilSearchPlanError("backend contract must be an object")
    backend_record = {
        "python_flint": str(backend.get("python_flint", "")),
        "flint": str(backend.get("flint", "")),
    }
    if not all(backend_record.values()):
        raise WeilSearchPlanError("backend contract versions are required")

    supplied_engine = raw.get("engine", {})
    if supplied_engine and not isinstance(supplied_engine, Mapping):
        raise WeilSearchPlanError("engine must be an object")
    if supplied_engine and supplied_engine.get("algorithm_id", ENGINE_ID) != ENGINE_ID:
        raise WeilSearchPlanError("unsupported search algorithm")

    policy_record = policy.to_record()
    body = {
        "schema": PLAN_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "cutoffs": [
            {"numerator": str(numerator), "denominator": str(denominator)}
            for numerator, denominator in cutoffs
        ],
        "degrees": [str(value) for value in degrees],
        "precision_policy": {
            "attempt_bits": policy_record["attempt_bits"],
            "confirmation_bits": policy_record["confirmation_bits"],
            "positive_terminal_gate": policy_record["positive_terminal_gate"],
            "candidate_confirmation": policy_record["candidate_confirmation"],
        },
        "witness_policy": policy_record["witness_policy"],
        "engine": {
            "algorithm_id": ENGINE_ID,
            "matrix": "A=P-R-S",
            "matrix_builder": "complete exact-cutoff finite Weil components",
            "source_oracles_required": True,
            "cell_order": "cutoff ascending, then degree ascending",
            "parallelism": "serial; FLINT ctx.prec is shared mutable state",
        },
        "backend_contract": backend_record,
        "limitation": EXPLORATORY_LIMITATION,
    }
    plan_sha256 = content_sha256(body)
    supplied_hash = raw.get("plan_sha256")
    if supplied_hash is not None and supplied_hash != plan_sha256:
        raise WeilSearchPlanError("search plan hash mismatch")
    return {**body, "plan_sha256": plan_sha256}


def load_search_plan(path: Path) -> dict[str, Any]:
    """Load and canonicalize a search plan from JSON."""

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise WeilSearchPlanError(f"cannot read search plan: {path}") from exc
    return canonicalize_plan(raw)


def _cells_from_plan(plan: Mapping[str, Any]) -> list[WeilSearchCell]:
    cells = [
        WeilSearchCell(
            _parse_integer(cutoff["numerator"], "cutoff numerator"),
            _parse_integer(cutoff["denominator"], "cutoff denominator"),
            _parse_integer(degree, "degree"),
        )
        for cutoff in plan["cutoffs"]
        for degree in plan["degrees"]
    ]
    return sorted(cells, key=lambda cell: cell.sort_key)


def _policy_from_canonical_plan(plan: Mapping[str, Any]) -> WeilSearchPolicy:
    return WeilSearchPolicy(
        attempt_bits=tuple(
            _parse_integer(value, "attempt precision")
            for value in plan["precision_policy"]["attempt_bits"]
        ),
        confirmation_bits=_parse_integer(
            plan["precision_policy"]["confirmation_bits"],
            "confirmation precision",
        ),
        eigenpair_count=_parse_integer(
            plan["witness_policy"]["eigenpair_count"], "eigenpair count"
        ),
        scale_bits=tuple(
            _parse_integer(value, "rational scale bits")
            for value in plan["witness_policy"]["scale_bits"]
        ),
    )


@contextmanager
def _working_precision(precision_bits: int) -> Iterable[None]:
    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        yield
    finally:
        ctx.prec = previous_precision


def _dyadic_fraction(mantissa: Any, exponent: Any) -> Fraction:
    mantissa_int = int(mantissa)
    exponent_int = int(exponent)
    if exponent_int >= 0:
        return Fraction(mantissa_int << exponent_int, 1)
    return Fraction(mantissa_int, 1 << (-exponent_int))


def _midpoint_fraction(value: Any) -> Fraction:
    midpoint = value.mid()
    return _dyadic_fraction(*midpoint.man_exp())


def _fraction_to_dyadic(value: Fraction) -> tuple[str, str]:
    denominator = value.denominator
    if denominator & (denominator - 1):
        raise ValueError("value is not dyadic")
    mantissa = value.numerator
    exponent = -(denominator.bit_length() - 1)
    if mantissa:
        while mantissa % 2 == 0:
            mantissa //= 2
            exponent += 1
    return str(mantissa), str(exponent)


def _record_bounds(record: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    encoded = record["dyadic"]
    midpoint = _dyadic_fraction(
        encoded["mid_mantissa"], encoded["mid_exponent"]
    )
    radius = _dyadic_fraction(
        encoded["radius_mantissa"], encoded["radius_exponent"]
    )
    return midpoint - radius, midpoint + radius


def _stable_arb_record(ball: Any, precision_bits: int) -> dict[str, Any]:
    """Outward-round an enclosure onto a cache-stable dyadic grid."""

    if ball.is_exact():
        return arb_record(ball)
    midpoint = _dyadic_fraction(*ball.mid().man_exp())
    radius = _dyadic_fraction(*ball.rad().man_exp())
    retained_bits = max(precision_bits - 32, 16)
    quantum = Fraction(1, 1 << retained_bits)
    lower_index = (midpoint - radius) // quantum - 1
    upper_index = -((-(midpoint + radius)) // quantum) + 1
    lower = lower_index * quantum
    upper = upper_index * quantum
    stable_midpoint = (lower + upper) / 2
    stable_radius = (upper - lower) / 2
    midpoint_pair = _fraction_to_dyadic(stable_midpoint)
    radius_pair = _fraction_to_dyadic(stable_radius)
    encoded = {
        "mid_mantissa": midpoint_pair[0],
        "mid_exponent": midpoint_pair[1],
        "radius_mantissa": radius_pair[0],
        "radius_exponent": radius_pair[1],
    }
    reconstructed = arb_from_dyadic(encoded)
    return {
        "dyadic": encoded,
        "display": reconstructed.str(40),
        "is_exact": stable_radius == 0,
    }


def _stabilize_ball_records(value: Any, precision_bits: int) -> Any:
    if isinstance(value, dict):
        if set(value) == {"dyadic", "display", "is_exact"}:
            return _stable_arb_record(
                arb_from_dyadic(value["dyadic"]), precision_bits
            )
        return {
            key: _stabilize_ball_records(item, precision_bits)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_stabilize_ball_records(item, precision_bits) for item in value]
    return value


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _round_ties_even(value: Fraction) -> int:
    sign = -1 if value < 0 else 1
    numerator = abs(value.numerator)
    denominator = value.denominator
    quotient, remainder = divmod(numerator, denominator)
    twice_remainder = 2 * remainder
    if twice_remainder > denominator or (
        twice_remainder == denominator and quotient % 2 == 1
    ):
        quotient += 1
    return sign * quotient


def _primitive_witness(values: Sequence[int]) -> tuple[int, ...]:
    if not values or not any(values):
        raise ValueError("integer witness must be nonzero")
    common = reduce(gcd, (abs(value) for value in values if value), 0)
    primitive = tuple(value // common for value in values)
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = tuple(-value for value in primitive)
    return primitive


def _project_direction(
    direction: Sequence[Fraction], projection: str
) -> tuple[Fraction, ...]:
    reverse = tuple(reversed(direction))
    if projection == "even":
        return tuple(left + right for left, right in zip(direction, reverse))
    if projection == "odd":
        return tuple(left - right for left, right in zip(direction, reverse))
    if projection == "raw":
        return tuple(direction)
    raise ValueError("unknown direction projection")


def extract_integer_witnesses(
    matrix: Sequence[Sequence[Any]],
    policy: WeilSearchPolicy = WeilSearchPolicy(),
) -> list[dict[str, Any]]:
    """Extract deterministic primitive integer candidates from low eigenmodes.

    The approximate eigensolver is used only to propose vectors.  Every returned
    vector must subsequently be checked with :func:`evaluate_integer_rayleigh`.
    No approximate eigenvalue or eigenvector is certificate evidence.
    """

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    eigenvalues, right_vectors = arb_mat(matrix).eig(
        right=True, algorithm="approx"
    )
    order = sorted(
        range(len(eigenvalues)),
        key=lambda index: (_midpoint_fraction(eigenvalues[index].real), index),
    )[: min(policy.eigenpair_count, len(eigenvalues))]

    candidates: list[dict[str, Any]] = []
    seen: set[tuple[int, ...]] = set()
    for rank, eigen_index in enumerate(order):
        eigenvalue_midpoint = _midpoint_fraction(eigenvalues[eigen_index].real)
        for part_name, attribute_name in (
            ("real", "real"),
            ("imaginary", "imag"),
        ):
            direction = tuple(
                _midpoint_fraction(
                    getattr(right_vectors[row, eigen_index], attribute_name)
                )
                for row in range(size)
            )
            if not any(direction):
                continue
            for projection in ("even", "odd", "raw"):
                projected = _project_direction(direction, projection)
                if not any(projected):
                    continue
                pivot = max(
                    range(size),
                    key=lambda index: (abs(projected[index]), -index),
                )
                pivot_value = projected[pivot]
                normalized = tuple(value / pivot_value for value in projected)
                for scale_bits in policy.scale_bits:
                    scale = 1 << scale_bits
                    rounded = tuple(
                        _round_ties_even(value * scale) for value in normalized
                    )
                    witness = _primitive_witness(rounded)
                    if witness in seen:
                        continue
                    seen.add(witness)
                    candidates.append(
                        {
                            "eigen_rank": rank,
                            "eigen_original_index": eigen_index,
                            "eigenvalue_real_midpoint": _fraction_record(
                                eigenvalue_midpoint
                            ),
                            "source_part": part_name,
                            "projection": projection,
                            "scale_bits": scale_bits,
                            "pivot_index": pivot,
                            "witness": list(witness),
                        }
                    )
    return candidates


def evaluate_integer_rayleigh(
    matrix: Sequence[Sequence[Any]], witness: Sequence[int]
) -> dict[str, Any]:
    """Certify the sign of an integer-vector Rayleigh quotient with Arb balls."""

    if any(isinstance(value, bool) or not isinstance(value, int) for value in witness):
        raise TypeError("witness entries must be integers")
    primitive = _primitive_witness(witness)
    evaluation = _evaluate_integer_rayleigh(matrix, primitive)
    classification = (
        "NEGATIVE"
        if evaluation["sign"] == "NEGATIVE"
        else "POSITIVE_FOR_WITNESS"
        if evaluation["sign"] == "POSITIVE"
        else "INCONCLUSIVE"
    )
    return {
        **evaluation,
        "classification": classification,
        "primitive_gcd": "1",
    }


def _stable_evaluation(
    evaluation: Mapping[str, Any], precision_bits: int
) -> dict[str, Any]:
    stable = _stabilize_ball_records(dict(evaluation), precision_bits)
    quotient = arb_from_dyadic(stable["rayleigh_quotient"]["dyadic"])
    negative = quotient < 0
    positive = quotient > 0
    stable["sign"] = (
        "NEGATIVE" if negative else "POSITIVE" if positive else "INCONCLUSIVE"
    )
    stable["classification"] = (
        "NEGATIVE"
        if negative
        else "POSITIVE_FOR_WITNESS"
        if positive
        else "INCONCLUSIVE"
    )
    stable["upper_bound_is_negative"] = negative
    stable["lower_bound_is_positive"] = positive
    return stable


def _matrix_evidence(components: Mapping[str, Any]) -> dict[str, Any]:
    indices = components["indices"]
    pair_count = len(indices) * (len(indices) + 1) // 2
    evidence = {
        "normalization": "A=P-R-S",
        "indices": [str(value) for value in indices],
        "dimension": str(len(indices)),
        "prime_power_transcript": _prime_transcript(dict(components)),
        "matrix_entries_upper_triangle": _serialize_entries(dict(components)),
        "checks": {
            "exact_prime_power_cutoff": True,
            "symmetric_enclosures_by_construction": True,
            "pole_oracle_pairs": str(components["pole_oracle_pairs"]),
            "archimedean_oracle_pairs": str(
                components["archimedean_oracle_pairs"]
            ),
            "all_pole_oracle_pairs_checked": (
                components["pole_oracle_pairs"] == pair_count
            ),
            "all_archimedean_oracle_pairs_checked": (
                components["archimedean_oracle_pairs"] == pair_count
            ),
        },
    }
    return _stabilize_ball_records(evidence, ctx.prec)


def _stable_ldlt(
    certificate: Mapping[str, Any], precision_bits: int
) -> dict[str, Any]:
    stable = _stabilize_ball_records(dict(certificate), precision_bits)
    if stable.get("classification") == "POSITIVE_DEFINITE":
        for pivot in stable["pivots"]:
            lower, _ = _record_bounds(pivot["value"])
            if not (lower > 0):
                stable["classification"] = "INCONCLUSIVE"
                stable["failed_pivot_index"] = pivot["index"]
                stable["failure_code"] = (
                    "STABLE_SERIALIZATION_DID_NOT_PRESERVE_POSITIVE_PIVOT"
                )
                break
    return stable


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _require_payload_hash(
    artifact: Mapping[str, Any], error_type: type[WeilSearchError]
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type("payload hash mismatch")


def _rump_diagnostic(matrix: Sequence[Sequence[Any]]) -> dict[str, Any]:
    try:
        return certify_rump_spectrum(matrix)
    except (ArithmeticError, ValueError):
        return {
            "classification": "INCONCLUSIVE",
            "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
        }


def _candidate_artifact_record(candidate: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "eigen_rank": str(candidate["eigen_rank"]),
        "eigen_original_index": str(candidate["eigen_original_index"]),
        "source_part": candidate["source_part"],
        "projection": candidate["projection"],
        "scale_bits": str(candidate["scale_bits"]),
        "pivot_index": str(candidate["pivot_index"]),
        "witness": [str(value) for value in candidate["witness"]],
    }


def _build_components(cell: WeilSearchCell) -> dict[str, Any]:
    return _symmetric_components(
        cell.cutoff_numerator, cell.cutoff_denominator, cell.degree
    )


def _search_attempt(
    cell: WeilSearchCell,
    policy: WeilSearchPolicy,
    precision_bits: int,
    has_next_precision: bool,
) -> dict[str, Any]:
    base = {
        "schema": ATTEMPT_SCHEMA,
        "classification": "EXPLORATORY",
        "kind": "SEARCH",
        "cell": cell.to_record(),
        "precision_bits": str(precision_bits),
        "limitation": EXPLORATORY_LIMITATION,
    }
    with _working_precision(precision_bits):
        try:
            components = _build_components(cell)
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE",
                    },
                    "decision": (
                        "ESCALATE_PRECISION"
                        if has_next_precision
                        else "INCONCLUSIVE_MAX_PRECISION"
                    ),
                }
            )

        matrix = components["q"]
        evidence = _matrix_evidence(components)
        ldlt = _stable_ldlt(certify_positive_ldlt(matrix), precision_bits)
        rump = _stabilize_ball_records(
            _rump_diagnostic(matrix), precision_bits
        )
        if ldlt["classification"] == "POSITIVE_DEFINITE":
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {"classification": "COMPLETE"},
                    "matrix_evidence": evidence,
                    "positive_definiteness": ldlt,
                    "secondary_rump_spectrum": rump,
                    "candidate_generation": {
                        "classification": "SKIPPED_LDLT_POSITIVE",
                        "role": "untrusted-search-hint-only",
                        "tested": [],
                    },
                    "decision": "FINITE_POSITIVE_CERTIFIED",
                }
            )

        try:
            candidates = extract_integer_witnesses(matrix, policy)
        except (ArithmeticError, ValueError):
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {"classification": "COMPLETE"},
                    "matrix_evidence": evidence,
                    "positive_definiteness": ldlt,
                    "secondary_rump_spectrum": rump,
                    "candidate_generation": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": "APPROXIMATE_EIGENVECTOR_EXTRACTION_FAILED",
                        "role": "untrusted-search-hint-only",
                        "tested": [],
                    },
                    "decision": (
                        "ESCALATE_PRECISION"
                        if has_next_precision
                        else "INCONCLUSIVE_MAX_PRECISION"
                    ),
                }
            )

        tested: list[dict[str, Any]] = []
        negative_index: int | None = None
        for candidate in candidates:
            evaluation = _stable_evaluation(
                evaluate_integer_rayleigh(matrix, candidate["witness"]),
                precision_bits,
            )
            tested.append(
                {
                    "candidate": _candidate_artifact_record(candidate),
                    "evaluation": evaluation,
                }
            )
            if evaluation["classification"] == "NEGATIVE":
                negative_index = len(tested) - 1
                break

        candidate_generation: dict[str, Any] = {
            "classification": "COMPLETED",
            "role": "untrusted-search-hint-only",
            "eigensolver": "flint-arb_mat-eig-approx",
            "tested": tested,
        }
        if negative_index is not None:
            candidate_generation["negative_candidate_index"] = str(negative_index)
            decision = "NEGATIVE_WITNESS_DISCOVERED"
        else:
            candidate_generation["negative_candidate_index"] = None
            decision = (
                "ESCALATE_PRECISION"
                if has_next_precision
                else "INCONCLUSIVE_MAX_PRECISION"
            )
        return _with_payload_hash(
            {
                **base,
                "matrix_construction": {"classification": "COMPLETE"},
                "matrix_evidence": evidence,
                "positive_definiteness": ldlt,
                "secondary_rump_spectrum": rump,
                "candidate_generation": candidate_generation,
                "decision": decision,
            }
        )


def _entry_map(evidence: Mapping[str, Any]) -> dict[tuple[str, str], Any]:
    return {
        (entry["m"], entry["n"]): entry
        for entry in evidence["matrix_entries_upper_triangle"]
    }


def _precision_containment(
    original_evidence: Mapping[str, Any], replay_evidence: Mapping[str, Any]
) -> tuple[bool, dict[str, str] | None]:
    original = _entry_map(original_evidence)
    replay = _entry_map(replay_evidence)
    if original.keys() != replay.keys():
        return False, {"component": "index_schedule", "m": "", "n": ""}
    for key in sorted(original, key=lambda item: (int(item[0]), int(item[1]))):
        for component in ("pole", "archimedean", "prime_power", "q"):
            original_ball = arb_from_dyadic(original[key][component]["dyadic"])
            replay_ball = arb_from_dyadic(replay[key][component]["dyadic"])
            if not original_ball.contains(replay_ball):
                return False, {
                    "component": component,
                    "m": key[0],
                    "n": key[1],
                }
    return True, None


def _confirmation_attempt(
    cell: WeilSearchCell,
    precision_bits: int,
    witness: Sequence[int],
    discovery_attempt: Mapping[str, Any],
) -> dict[str, Any]:
    base = {
        "schema": ATTEMPT_SCHEMA,
        "classification": "EXPLORATORY",
        "kind": "CONFIRMATION",
        "cell": cell.to_record(),
        "precision_bits": str(precision_bits),
        "limitation": EXPLORATORY_LIMITATION,
    }
    with _working_precision(precision_bits):
        try:
            components = _build_components(cell)
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE",
                    },
                    "witness": [str(value) for value in _primitive_witness(witness)],
                    "decision": "NEGATIVE_CONFIRMATION_FAILED",
                }
            )
        evidence = _matrix_evidence(components)
        evaluation = _stable_evaluation(
            evaluate_integer_rayleigh(components["q"], witness),
            precision_bits,
        )
        contained, failed_entry = _precision_containment(
            discovery_attempt["matrix_evidence"], evidence
        )
        decision = (
            "NEGATIVE_WITNESS_CONFIRMED"
            if evaluation["classification"] == "NEGATIVE" and contained
            else "NEGATIVE_CONFIRMATION_FAILED"
        )
        return _with_payload_hash(
            {
                **base,
                "matrix_construction": {"classification": "COMPLETE"},
                "matrix_evidence": evidence,
                "witness_evaluation": evaluation,
                "precision_containment": {
                    "all_higher_precision_components_contained": contained,
                    "first_failure": failed_entry,
                },
                "decision": decision,
            }
        )


def _negative_record(attempt: Mapping[str, Any]) -> dict[str, Any]:
    try:
        index = int(attempt["candidate_generation"]["negative_candidate_index"])
        record = attempt["candidate_generation"]["tested"][index]
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        raise WeilSearchCheckpointError(
            "negative discovery does not identify a tested witness"
        ) from exc
    evaluation = record.get("evaluation", {})
    if not (
        evaluation.get("classification") == "NEGATIVE"
        and evaluation.get("sign") == "NEGATIVE"
        and evaluation.get("upper_bound_is_negative") is True
        and evaluation.get("lower_bound_is_positive") is False
    ):
        raise WeilSearchCheckpointError(
            "declared negative candidate is not interval-negative"
        )
    return record


def _derived_search_decision(
    attempt: Mapping[str, Any], has_next_precision: bool
) -> str:
    construction = attempt.get("matrix_construction", {}).get("classification")
    terminal_inconclusive = (
        "ESCALATE_PRECISION"
        if has_next_precision
        else "INCONCLUSIVE_MAX_PRECISION"
    )
    if construction != "COMPLETE":
        return terminal_inconclusive

    ldlt = attempt.get("positive_definiteness", {})
    if ldlt.get("classification") == "POSITIVE_DEFINITE":
        pivots = ldlt.get("pivots")
        if not isinstance(pivots, list) or not pivots:
            raise WeilSearchCheckpointError("positive LDLT certificate has no pivots")
        for pivot in pivots:
            lower, _ = _record_bounds(pivot["value"])
            if not (lower > 0):
                raise WeilSearchCheckpointError(
                    "positive LDLT decision contains a nonpositive pivot"
                )
        generation = attempt.get("candidate_generation", {})
        if not (
            generation.get("classification") == "SKIPPED_LDLT_POSITIVE"
            and generation.get("tested") == []
            and generation.get("role") == "untrusted-search-hint-only"
        ):
            raise WeilSearchCheckpointError(
                "positive LDLT attempt has invalid candidate-generation state"
            )
        return "FINITE_POSITIVE_CERTIFIED"

    generation = attempt.get("candidate_generation")
    if not isinstance(generation, Mapping):
        raise WeilSearchCheckpointError(
            "nonpositive search attempt lacks candidate-generation evidence"
        )
    tested = generation.get("tested")
    if not isinstance(tested, list):
        raise WeilSearchCheckpointError("tested witnesses must be a list")
    generation_classification = generation.get("classification")
    if generation_classification not in {"COMPLETED", "INCONCLUSIVE"}:
        raise WeilSearchCheckpointError(
            "candidate-generation classification is invalid"
        )
    if generation_classification == "INCONCLUSIVE" and tested:
        raise WeilSearchCheckpointError(
            "inconclusive candidate generation cannot contain tested witnesses"
        )
    for record in tested:
        candidate = record.get("candidate", {})
        evaluation = record.get("evaluation", {})
        if candidate.get("witness") != evaluation.get("witness"):
            raise WeilSearchCheckpointError(
                "candidate and evaluated integer witnesses differ"
            )
        classification = evaluation.get("classification")
        sign = evaluation.get("sign")
        negative = evaluation.get("upper_bound_is_negative") is True
        positive = evaluation.get("lower_bound_is_positive") is True
        if classification == "NEGATIVE":
            valid = sign == "NEGATIVE" and negative and not positive
        elif classification == "POSITIVE_FOR_WITNESS":
            valid = sign == "POSITIVE" and positive and not negative
        elif classification == "INCONCLUSIVE":
            valid = sign == "INCONCLUSIVE" and not positive and not negative
        else:
            valid = False
        if not valid:
            raise WeilSearchCheckpointError(
                "integer witness sign fields are inconsistent"
            )
    negative_indices = [
        index
        for index, record in enumerate(tested)
        if record.get("evaluation", {}).get("classification") == "NEGATIVE"
        and record.get("evaluation", {}).get("upper_bound_is_negative") is True
    ]
    declared = generation.get("negative_candidate_index")
    if negative_indices:
        if (
            declared != str(negative_indices[0])
            or negative_indices[0] != len(tested) - 1
        ):
            raise WeilSearchCheckpointError(
                "candidate search did not stop at its first negative witness"
            )
        _negative_record(attempt)
        return "NEGATIVE_WITNESS_DISCOVERED"
    if declared is not None:
        raise WeilSearchCheckpointError(
            "negative candidate index exists without a negative witness"
        )
    return terminal_inconclusive


def _derived_confirmation_decision(
    discovery: Mapping[str, Any], confirmation: Mapping[str, Any]
) -> str:
    discovered = _negative_record(discovery)
    discovered_witness = discovered["evaluation"]["witness"]
    construction = confirmation.get("matrix_construction", {}).get(
        "classification"
    )
    if construction != "COMPLETE":
        if confirmation.get("witness") != discovered_witness:
            raise WeilSearchCheckpointError(
                "failed confirmation changed the discovered witness"
            )
        return "NEGATIVE_CONFIRMATION_FAILED"
    evaluation = confirmation.get("witness_evaluation")
    if not isinstance(evaluation, Mapping):
        raise WeilSearchCheckpointError("confirmation witness evaluation is missing")
    if evaluation.get("witness") != discovered_witness:
        raise WeilSearchCheckpointError("confirmation changed the discovered witness")
    containment = confirmation.get("precision_containment", {}).get(
        "all_higher_precision_components_contained"
    )
    negative = (
        evaluation.get("classification") == "NEGATIVE"
        and evaluation.get("sign") == "NEGATIVE"
        and evaluation.get("upper_bound_is_negative") is True
        and evaluation.get("lower_bound_is_positive") is False
    )
    return (
        "NEGATIVE_WITNESS_CONFIRMED"
        if negative and containment is True
        else "NEGATIVE_CONFIRMATION_FAILED"
    )


def _validate_attempt_hashes(attempts: Sequence[Mapping[str, Any]]) -> None:
    for expected_sequence, attempt in enumerate(attempts):
        if not isinstance(attempt, Mapping):
            raise WeilSearchCheckpointError("attempt must be an object")
        _require_payload_hash(attempt, WeilSearchCheckpointError)
        if attempt.get("schema") != ATTEMPT_SCHEMA:
            raise WeilSearchCheckpointError("unexpected attempt schema")
        if attempt.get("classification") != "EXPLORATORY":
            raise WeilSearchCheckpointError("attempt classification changed")
        if attempt.get("sequence") != str(expected_sequence):
            raise WeilSearchCheckpointError("attempt sequence is not contiguous")


AttemptSink = Callable[[dict[str, Any]], None]


def _append_sequence(attempt: Mapping[str, Any], sequence: int) -> dict[str, Any]:
    body = {key: value for key, value in attempt.items() if key != "payload_sha256"}
    body["sequence"] = str(sequence)
    return _with_payload_hash(body)


def _cell_payload(
    cell: WeilSearchCell,
    policy: WeilSearchPolicy,
    attempts: Sequence[Mapping[str, Any]],
    status: str,
    plan_sha256: str | None,
    candidate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": CELL_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "cell": cell.to_record(),
        "policy": policy.to_record(),
        "status": status,
        "attempts": list(attempts),
        "quarantined_candidate": candidate,
        "limitation": EXPLORATORY_LIMITATION,
    }
    if plan_sha256 is not None:
        payload["plan_sha256"] = plan_sha256
    return _with_payload_hash(payload)


def _search_weil_cell(
    cell: WeilSearchCell,
    policy: WeilSearchPolicy,
    *,
    existing_attempts: Sequence[Mapping[str, Any]] = (),
    attempt_sink: AttemptSink | None = None,
    plan_sha256: str | None = None,
) -> dict[str, Any]:
    attempts = [dict(attempt) for attempt in existing_attempts]
    _validate_attempt_hashes(attempts)
    cursor = 0

    def obtain(factory: Callable[[], dict[str, Any]]) -> dict[str, Any]:
        nonlocal cursor
        if cursor < len(attempts):
            stored = attempts[cursor]
            cursor += 1
            return stored
        new_attempt = factory()
        stored = _append_sequence(new_attempt, len(attempts))
        attempts.append(stored)
        cursor += 1
        if attempt_sink is not None:
            attempt_sink(stored)
        return stored

    for precision_index, precision_bits in enumerate(policy.attempt_bits):
        has_next = precision_index + 1 < len(policy.attempt_bits)
        attempt = obtain(
            lambda: _search_attempt(cell, policy, precision_bits, has_next)
        )
        if attempt.get("kind") != "SEARCH":
            raise WeilSearchCheckpointError("expected a search attempt")
        if attempt.get("cell") != cell.to_record():
            raise WeilSearchCheckpointError("attempt cell does not match plan")
        if attempt.get("precision_bits") != str(precision_bits):
            raise WeilSearchCheckpointError("attempt precision schedule changed")
        decision = attempt.get("decision")
        derived_decision = _derived_search_decision(attempt, has_next)
        if decision != derived_decision:
            raise WeilSearchCheckpointError(
                "search decision does not follow from its certificate evidence"
            )
        if decision == "FINITE_POSITIVE_CERTIFIED":
            if cursor != len(attempts):
                raise WeilSearchCheckpointError("attempts continue past terminal cell")
            return _cell_payload(
                cell,
                policy,
                attempts,
                "FINITE_POSITIVE_CERTIFIED",
                plan_sha256,
            )
        if decision == "NEGATIVE_WITNESS_DISCOVERED":
            negative = _negative_record(attempt)
            witness = tuple(int(value) for value in negative["evaluation"]["witness"])
            confirmation_bits = (
                policy.attempt_bits[precision_index + 1]
                if has_next
                else policy.confirmation_bits
            )
            confirmation = obtain(
                lambda: _confirmation_attempt(
                    cell, confirmation_bits, witness, attempt
                )
            )
            if confirmation.get("kind") != "CONFIRMATION":
                raise WeilSearchCheckpointError("expected a confirmation attempt")
            if confirmation.get("precision_bits") != str(confirmation_bits):
                raise WeilSearchCheckpointError(
                    "confirmation precision schedule changed"
                )
            if confirmation.get("cell") != cell.to_record():
                raise WeilSearchCheckpointError(
                    "confirmation cell does not match plan"
                )
            derived_confirmation = _derived_confirmation_decision(
                attempt, confirmation
            )
            if confirmation.get("decision") != derived_confirmation:
                raise WeilSearchCheckpointError(
                    "confirmation decision does not follow from its evidence"
                )
            if cursor != len(attempts):
                raise WeilSearchCheckpointError("attempts continue past confirmation")
            if confirmation.get("decision") == "NEGATIVE_WITNESS_CONFIRMED":
                candidate = {
                    "status": "NEGATIVE_CANDIDATE_QUARANTINED",
                    "promotion_status": (
                        "pending-independent-sign-normalization-source-equation-"
                        "and-backend-audit"
                    ),
                    "witness": negative["evaluation"]["witness"],
                    "discovery_precision_bits": attempt["precision_bits"],
                    "confirmation_precision_bits": confirmation["precision_bits"],
                    "discovery_evaluation": negative["evaluation"],
                    "confirmation_evaluation": confirmation["witness_evaluation"],
                    "same_backend_replay_only": True,
                }
                return _cell_payload(
                    cell,
                    policy,
                    attempts,
                    "NEGATIVE_CANDIDATE_QUARANTINED",
                    plan_sha256,
                    candidate,
                )
            return _cell_payload(
                cell,
                policy,
                attempts,
                "INCONCLUSIVE_MAX_PRECISION",
                plan_sha256,
                {
                    "status": "CONFLICTING_NEGATIVE_REPLAY_QUARANTINED",
                    "witness": negative["evaluation"]["witness"],
                    "discovery_precision_bits": attempt["precision_bits"],
                    "confirmation_precision_bits": confirmation["precision_bits"],
                },
            )
        if decision == "ESCALATE_PRECISION":
            if not has_next:
                raise WeilSearchCheckpointError(
                    "final precision cannot request escalation"
                )
            continue
        if decision == "INCONCLUSIVE_MAX_PRECISION":
            if has_next:
                raise WeilSearchCheckpointError(
                    "nonfinal precision stopped as max precision"
                )
            if cursor != len(attempts):
                raise WeilSearchCheckpointError("attempts continue past terminal cell")
            return _cell_payload(
                cell,
                policy,
                attempts,
                "INCONCLUSIVE_MAX_PRECISION",
                plan_sha256,
            )
        raise WeilSearchCheckpointError("unknown search attempt decision")
    raise AssertionError("precision policy exhausted without a terminal result")


def search_weil_cell(
    cell: WeilSearchCell,
    policy: WeilSearchPolicy = WeilSearchPolicy(),
) -> dict[str, Any]:
    """Search one finite cell without writing checkpoints."""

    return _search_weil_cell(cell, policy)


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _load_checkpoint(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise WeilSearchCheckpointError(f"cannot read checkpoint: {path}") from exc
    if not isinstance(value, dict):
        raise WeilSearchCheckpointError("checkpoint must contain an object")
    _require_payload_hash(value, WeilSearchCheckpointError)
    return value


def _run_manifest(plan: Mapping[str, Any]) -> dict[str, Any]:
    return _with_payload_hash(
        {
            "schema": RUN_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "plan": dict(plan),
            "plan_sha256": plan["plan_sha256"],
            "backend": _backend_record(),
            "engine_id": ENGINE_ID,
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def _validate_runtime_contract(plan: Mapping[str, Any]) -> None:
    if plan["backend_contract"] != _backend_record():
        raise WeilSearchPlanError(
            "installed python-flint/FLINT versions do not match the plan"
        )


def _load_attempts(attempt_directory: Path) -> list[dict[str, Any]]:
    if not attempt_directory.exists():
        return []
    paths = sorted(attempt_directory.glob("*.json"))
    attempts = [_load_checkpoint(path) for path in paths]
    _validate_attempt_hashes(attempts)
    return attempts


def _index_payload(
    plan: Mapping[str, Any],
    cells: Sequence[WeilSearchCell],
    completed: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    references: list[dict[str, Any]] = []
    for cell in cells:
        artifact = completed.get(cell.cell_id)
        if artifact is None:
            continue
        references.append(
            {
                "cell_id": cell.cell_id,
                "path": f"cells/{cell.cell_id}.json",
                "cell_payload_sha256": artifact["payload_sha256"],
                "status": artifact["status"],
            }
        )
    missing = [cell.cell_id for cell in cells if cell.cell_id not in completed]
    candidate_present = any(
        artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED"
        for artifact in completed.values()
    )
    if candidate_present:
        conclusion = "NEGATIVE_CANDIDATE_QUARANTINED"
    elif missing:
        conclusion = "INCOMPLETE"
    else:
        conclusion = "NO_NEGATIVE_WITNESS_IN_FINITE_GRID"
    statuses = {
        status: str(
            sum(artifact["status"] == status for artifact in completed.values())
        )
        for status in (
            "FINITE_POSITIVE_CERTIFIED",
            "NEGATIVE_CANDIDATE_QUARANTINED",
            "INCONCLUSIVE_MAX_PRECISION",
        )
    }
    return _with_payload_hash(
        {
            "schema": INDEX_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "conclusion": conclusion,
            "plan": dict(plan),
            "plan_sha256": plan["plan_sha256"],
            "backend": _backend_record(),
            "cells": references,
            "missing_cells": missing,
            "progress": {
                "planned_cells": str(len(cells)),
                "completed_cells": str(len(completed)),
                "status_counts": statuses,
            },
            "limitation": EXPLORATORY_LIMITATION,
        }
    )


def run_weil_search(
    plan: Mapping[str, Any],
    checkpoint_dir: Path,
    *,
    resume: bool = False,
    max_cells: int | None = None,
) -> dict[str, Any]:
    """Run or resume a deterministic serial grid search.

    ``max_cells`` limits cells newly completed by this invocation and exists to
    support bounded batches and interruption tests.  It does not change the
    plan hash or the final artifact produced after a later resume.
    """

    canonical_plan = canonicalize_plan(plan)
    _validate_runtime_contract(canonical_plan)
    checkpoint_dir = checkpoint_dir.resolve()
    if max_cells is not None:
        if isinstance(max_cells, bool) or not isinstance(max_cells, int):
            raise TypeError("max_cells must be an integer or None")
        if max_cells < 0:
            raise ValueError("max_cells must be nonnegative")

    run_path = checkpoint_dir / "run.json"
    expected_run = _run_manifest(canonical_plan)
    if run_path.exists():
        if not resume:
            raise WeilSearchCheckpointError(
                "checkpoint already exists; pass resume=True"
            )
        stored_run = _load_checkpoint(run_path)
        if stored_run != expected_run:
            raise WeilSearchCheckpointError(
                "checkpoint plan, backend, or engine does not match"
            )
    else:
        if resume:
            raise WeilSearchCheckpointError("resume requested without a checkpoint")
        if checkpoint_dir.exists() and any(checkpoint_dir.iterdir()):
            raise WeilSearchCheckpointError(
                "fresh checkpoint directory must be empty"
            )
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        _atomic_write_json(run_path, expected_run)

    cells = _cells_from_plan(canonical_plan)
    policy = _policy_from_canonical_plan(canonical_plan)
    completed: dict[str, dict[str, Any]] = {}
    for cell in cells:
        cell_path = checkpoint_dir / "cells" / f"{cell.cell_id}.json"
        if not cell_path.exists():
            continue
        artifact = _load_checkpoint(cell_path)
        if artifact.get("schema") != CELL_SCHEMA:
            raise WeilSearchCheckpointError("unexpected cell checkpoint schema")
        if artifact.get("cell") != cell.to_record():
            raise WeilSearchCheckpointError("cell checkpoint input changed")
        if artifact.get("plan_sha256") != canonical_plan["plan_sha256"]:
            raise WeilSearchCheckpointError("cell checkpoint plan changed")
        attempts = artifact.get("attempts")
        if not isinstance(attempts, list) or not attempts:
            raise WeilSearchCheckpointError(
                "completed cell checkpoint has no attempts"
            )
        if attempts[-1].get("decision") not in {
            "FINITE_POSITIVE_CERTIFIED",
            "INCONCLUSIVE_MAX_PRECISION",
            "NEGATIVE_WITNESS_CONFIRMED",
            "NEGATIVE_CONFIRMATION_FAILED",
        }:
            raise WeilSearchCheckpointError(
                "completed cell checkpoint is not terminal"
            )
        derived_artifact = _search_weil_cell(
            cell,
            policy,
            existing_attempts=attempts,
            plan_sha256=canonical_plan["plan_sha256"],
        )
        if artifact != derived_artifact:
            raise WeilSearchCheckpointError(
                "completed cell checkpoint does not follow from its attempts"
            )
        completed[cell.cell_id] = artifact

    newly_completed = 0
    candidate_already_present = any(
        artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED"
        for artifact in completed.values()
    )
    if not candidate_already_present:
        for cell in cells:
            if cell.cell_id in completed:
                continue
            if max_cells is not None and newly_completed >= max_cells:
                break
            attempt_directory = checkpoint_dir / "attempts" / cell.cell_id
            existing_attempts = _load_attempts(attempt_directory)

            def save_attempt(attempt: dict[str, Any]) -> None:
                sequence = int(attempt["sequence"])
                _atomic_write_json(
                    attempt_directory / f"{sequence:04d}.json", attempt
                )

            artifact = _search_weil_cell(
                cell,
                policy,
                existing_attempts=existing_attempts,
                attempt_sink=save_attempt,
                plan_sha256=canonical_plan["plan_sha256"],
            )
            _atomic_write_json(
                checkpoint_dir / "cells" / f"{cell.cell_id}.json", artifact
            )
            completed[cell.cell_id] = artifact
            newly_completed += 1
            if artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED":
                break

    index = _index_payload(canonical_plan, cells, completed)
    _atomic_write_json(checkpoint_dir / "index.json", index)
    return index


def _safe_artifact_path(root: Path, relative: str) -> Path:
    if not isinstance(relative, str) or not relative:
        raise WeilSearchVerificationError("cell artifact path is invalid")
    path = (root / relative).resolve()
    if not path.is_relative_to(root):
        raise WeilSearchVerificationError("cell artifact path escapes checkpoint")
    return path


def _replay_attempt(
    cell: WeilSearchCell, attempt: Mapping[str, Any]
) -> list[dict[str, Any]]:
    precision_bits = _parse_integer(attempt["precision_bits"], "precision")
    if attempt.get("matrix_construction", {}).get("classification") != "COMPLETE":
        with _working_precision(precision_bits):
            try:
                _build_components(cell)
            except ArithmeticError:
                return []
        raise WeilSearchVerificationError(
            "stored matrix-construction failure now succeeds"
        )

    with _working_precision(precision_bits):
        try:
            components = _build_components(cell)
        except ArithmeticError as exc:
            raise WeilSearchVerificationError(
                "stored complete matrix no longer regenerates"
            ) from exc
        evidence = _matrix_evidence(components)
        if attempt.get("matrix_evidence") != evidence:
            raise WeilSearchVerificationError(
                "matrix evidence does not match canonical regeneration"
            )
        matrix = components["q"]
        verified: list[dict[str, Any]] = []
        if attempt.get("kind") == "SEARCH":
            ldlt = _stable_ldlt(
                certify_positive_ldlt(matrix), precision_bits
            )
            if attempt.get("positive_definiteness") != ldlt:
                raise WeilSearchVerificationError("LDLT replay changed")
            stored_rump = attempt.get("secondary_rump_spectrum")
            replay_rump = _stabilize_ball_records(
                _rump_diagnostic(matrix), precision_bits
            )
            if stored_rump != replay_rump:
                raise WeilSearchVerificationError("Rump diagnostic replay changed")
            for tested in attempt.get("candidate_generation", {}).get("tested", []):
                witness = [int(value) for value in tested["evaluation"]["witness"]]
                evaluation = _stable_evaluation(
                    evaluate_integer_rayleigh(matrix, witness),
                    precision_bits,
                )
                if tested["evaluation"] != evaluation:
                    raise WeilSearchVerificationError(
                        "integer witness evaluation changed"
                    )
                verified.append(evaluation)
        elif attempt.get("kind") == "CONFIRMATION":
            stored = attempt.get("witness_evaluation")
            if stored is not None:
                witness = [int(value) for value in stored["witness"]]
                evaluation = _stable_evaluation(
                    evaluate_integer_rayleigh(matrix, witness),
                    precision_bits,
                )
                if stored != evaluation:
                    raise WeilSearchVerificationError(
                        "confirmation witness evaluation changed"
                    )
                verified.append(evaluation)
        else:
            raise WeilSearchVerificationError("unknown attempt kind")
        return verified


def verify_weil_search(
    index: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Replay all certificate-bearing evidence in a search index.

    Approximate eigenvector generation is intentionally not replayed as a
    trusted step.  Every recorded integer vector is checked against a
    canonically regenerated full matrix instead.
    """

    root = checkpoint_dir.resolve()
    if isinstance(index, Path):
        try:
            loaded = json.loads(index.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise WeilSearchVerificationError("cannot read search index") from exc
        if not isinstance(loaded, dict):
            raise WeilSearchVerificationError("search index must be an object")
        index_record = loaded
    elif isinstance(index, Mapping):
        index_record = dict(index)
    else:
        raise TypeError("index must be a mapping or Path")

    _require_payload_hash(index_record, WeilSearchVerificationError)
    if index_record.get("schema") != INDEX_SCHEMA:
        raise WeilSearchVerificationError("unexpected search index schema")
    if index_record.get("classification") != "EXPLORATORY":
        raise WeilSearchVerificationError("search index was improperly promoted")
    if index_record.get("hypothesis_status") != "UNRESOLVED":
        raise WeilSearchVerificationError("search index changed global RH status")
    try:
        plan = canonicalize_plan(index_record["plan"])
    except (KeyError, WeilSearchPlanError) as exc:
        raise WeilSearchVerificationError("search index plan is invalid") from exc
    if plan != index_record["plan"]:
        raise WeilSearchVerificationError("search index plan is not canonical")
    if index_record.get("plan_sha256") != plan["plan_sha256"]:
        raise WeilSearchVerificationError("search index plan hash changed")
    try:
        _validate_runtime_contract(plan)
    except WeilSearchPlanError as exc:
        raise WeilSearchVerificationError(str(exc)) from exc

    planned_cells = _cells_from_plan(plan)
    planned_by_id = {cell.cell_id: cell for cell in planned_cells}
    references = index_record.get("cells")
    if not isinstance(references, list):
        raise WeilSearchVerificationError("search index cells must be a list")
    seen: set[str] = set()
    completed: dict[str, dict[str, Any]] = {}
    negative_evaluations = 0
    for reference in references:
        if not isinstance(reference, Mapping):
            raise WeilSearchVerificationError("cell reference must be an object")
        cell_id = reference.get("cell_id")
        if cell_id not in planned_by_id or cell_id in seen:
            raise WeilSearchVerificationError("unknown or duplicate cell reference")
        seen.add(cell_id)
        artifact_path = _safe_artifact_path(root, reference.get("path"))
        try:
            artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise WeilSearchVerificationError(
                f"cannot read cell artifact: {cell_id}"
            ) from exc
        if not isinstance(artifact, dict):
            raise WeilSearchVerificationError("cell artifact must be an object")
        _require_payload_hash(artifact, WeilSearchVerificationError)
        if artifact.get("payload_sha256") != reference.get("cell_payload_sha256"):
            raise WeilSearchVerificationError("cell reference hash mismatch")
        if artifact.get("schema") != CELL_SCHEMA:
            raise WeilSearchVerificationError("unexpected cell artifact schema")
        if artifact.get("classification") != "EXPLORATORY":
            raise WeilSearchVerificationError("cell artifact was improperly promoted")
        cell = planned_by_id[cell_id]
        if artifact.get("cell") != cell.to_record():
            raise WeilSearchVerificationError("cell artifact input changed")
        if artifact.get("plan_sha256") != plan["plan_sha256"]:
            raise WeilSearchVerificationError("cell artifact plan changed")
        attempts = artifact.get("attempts")
        if not isinstance(attempts, list) or not attempts:
            raise WeilSearchVerificationError("cell attempts are missing")
        try:
            _validate_attempt_hashes(attempts)
        except WeilSearchCheckpointError as exc:
            raise WeilSearchVerificationError(str(exc)) from exc
        for attempt in attempts:
            try:
                evaluations = _replay_attempt(cell, attempt)
            except WeilSearchVerificationError:
                raise
            except (AttributeError, IndexError, KeyError, TypeError, ValueError) as exc:
                raise WeilSearchVerificationError(
                    "attempt evidence is malformed"
                ) from exc
            negative_evaluations += sum(
                value["classification"] == "NEGATIVE" for value in evaluations
            )

        policy = _policy_from_canonical_plan(plan)
        try:
            derived_artifact = _search_weil_cell(
                cell,
                policy,
                existing_attempts=attempts,
                plan_sha256=plan["plan_sha256"],
            )
        except WeilSearchCheckpointError as exc:
            raise WeilSearchVerificationError(
                "cell state machine is inconsistent with replayed evidence"
            ) from exc
        except (AttributeError, IndexError, KeyError, TypeError, ValueError) as exc:
            raise WeilSearchVerificationError(
                "cell state machine evidence is malformed"
            ) from exc
        if artifact != derived_artifact:
            raise WeilSearchVerificationError(
                "cell status, policy, or quarantined candidate is not derived "
                "from its attempts"
            )

        status = artifact.get("status")
        if status == "FINITE_POSITIVE_CERTIFIED":
            terminal = attempts[-1]
            if terminal.get("decision") != "FINITE_POSITIVE_CERTIFIED":
                raise WeilSearchVerificationError("positive cell lacks LDLT terminal")
        elif status == "NEGATIVE_CANDIDATE_QUARANTINED":
            candidate = artifact.get("quarantined_candidate")
            if not isinstance(candidate, Mapping):
                raise WeilSearchVerificationError("negative candidate record is missing")
            witness = candidate.get("witness")
            if not isinstance(witness, list):
                raise WeilSearchVerificationError("candidate witness is missing")
            discovery = attempts[-2]
            confirmation = attempts[-1]
            if discovery.get("decision") != "NEGATIVE_WITNESS_DISCOVERED":
                raise WeilSearchVerificationError("candidate discovery is missing")
            if confirmation.get("decision") != "NEGATIVE_WITNESS_CONFIRMED":
                raise WeilSearchVerificationError("candidate confirmation is missing")
            if discovery["matrix_construction"]["classification"] != "COMPLETE":
                raise WeilSearchVerificationError("candidate matrix is incomplete")
            if confirmation["matrix_construction"]["classification"] != "COMPLETE":
                raise WeilSearchVerificationError("confirmation matrix is incomplete")
            with _working_precision(int(confirmation["precision_bits"])):
                contained, _ = _precision_containment(
                    discovery["matrix_evidence"], confirmation["matrix_evidence"]
                )
            if not contained:
                raise WeilSearchVerificationError(
                    "candidate precision containment failed"
                )
            if candidate.get("discovery_evaluation", {}).get(
                "classification"
            ) != "NEGATIVE" or candidate.get("confirmation_evaluation", {}).get(
                "classification"
            ) != "NEGATIVE":
                raise WeilSearchVerificationError(
                    "candidate does not have two negative evaluations"
                )
        elif status != "INCONCLUSIVE_MAX_PRECISION":
            raise WeilSearchVerificationError("unknown cell status")
        if reference.get("status") != status:
            raise WeilSearchVerificationError("cell reference status mismatch")
        completed[cell_id] = artifact

    expected_index = _index_payload(plan, planned_cells, completed)
    if index_record != expected_index:
        raise WeilSearchVerificationError(
            "search index does not match its referenced cell artifacts"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_SEARCH",
        "conclusion": index_record["conclusion"],
        "verified_cells": str(len(completed)),
        "verified_negative_witness_evaluations": str(negative_evaluations),
        "same_backend_replay_only": True,
        "hypothesis_status": "UNRESOLVED",
    }


__all__ = [
    "ATTEMPT_SCHEMA",
    "CELL_SCHEMA",
    "ENGINE_ID",
    "INDEX_SCHEMA",
    "PLAN_SCHEMA",
    "WeilSearchCell",
    "WeilSearchCheckpointError",
    "WeilSearchError",
    "WeilSearchPlanError",
    "WeilSearchPolicy",
    "WeilSearchVerificationError",
    "canonicalize_plan",
    "evaluate_integer_rayleigh",
    "extract_integer_witnesses",
    "load_search_plan",
    "run_weil_search",
    "search_weil_cell",
    "verify_weil_search",
]
