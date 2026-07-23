"""Exploratory structural audits for finite Weil matrices.

The routines in this module are deliberately separate from the v1 finite
certificate and search paths.  They certify internal consistency properties
of *finite* matrices only: reversal parity, principal-submatrix nesting, and
exact-integer witness transport.  Every top-level record therefore remains
``EXPLORATORY`` and keeps the global Riemann Hypothesis status
``UNRESOLVED``.

Approximate eigensolvers and Rump eigenvalue enclosures are diagnostic aids.
They never decide the sign of a witness.  The authoritative sign calculation
is always the direct interval evaluation of ``x^T (P-R-S) x`` in the original
Fourier basis.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable, Mapping, Sequence, TypeVar

from flint import arb, arb_mat, ctx

from .artifacts import content_sha256
from .balls import arb_to_dyadic
from .weil import (
    _symmetric_components,
    certify_rump_spectrum,
    evaluate_integer_rayleigh,
)
from .weil_search import (
    WeilSearchCell,
    WeilSearchPolicy,
    _backend_record,
    _fraction_record,
    _midpoint_fraction,
    _primitive_witness,
    _record_bounds,
    _round_ties_even,
    _stable_arb_record,
    _stable_evaluation,
    _stabilize_ball_records,
    _with_payload_hash,
    _working_precision,
)


PARITY_AUDIT_SCHEMA = "rh-lab/weil-parity-audit/v1"
NESTING_AUDIT_SCHEMA = "rh-lab/weil-degree-nesting-audit/v1"
PARITY_ALGORITHM_ID = "finite-weil-reversal-parity-audit-v1"
NESTING_ALGORITHM_ID = "finite-weil-degree-nesting-audit-v1"

EXPLORATORY_LIMITATION = (
    "These are finite-dimensional structural consistency audits.  Passing "
    "parity or degree-nesting gates, or observing positive diagnostic "
    "spectra, does not prove RH.  A directly interval-negative integer "
    "witness remains quarantined for independent mathematical and backend "
    "reproduction; all other outcomes leave RH unresolved."
)

_COMPONENTS = (
    ("P", "pole"),
    ("R", "archimedean"),
    ("S", "prime_power"),
    ("A", "q"),
)


class WeilAuditVerificationError(ValueError):
    """Raised when a parity or nesting audit cannot be replayed safely."""


_OperationResult = TypeVar("_OperationResult")


def _cleanup_flint_backend() -> None:
    """Clear process-global FLINT caches between canonical audit replays."""

    ctx.cleanup()


def _isolated_flint_operation(
    operation: Callable[[], _OperationResult],
) -> _OperationResult:
    """Run one canonical operation from a clean FLINT cache and leave one behind.

    FLINT's cache state can change the final radius of an otherwise equivalent
    Arb enclosure by a few ulps after a higher-degree Rump eigensolve.  Audit
    artifacts intentionally compare canonical serialized enclosures exactly,
    so every public verifier must start from the same clean backend state.
    """

    previous_precision = ctx.prec
    _cleanup_flint_backend()
    try:
        return operation()
    finally:
        ctx.prec = previous_precision
        _cleanup_flint_backend()


def _validate_precision(value: int, name: str = "precision_bits") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 96:
        raise ValueError(f"{name} must be an integer >= 96")
    return value


def _parse_canonical_integer(
    value: Any, name: str, minimum: int | None = None
) -> int:
    if not isinstance(value, str):
        raise WeilAuditVerificationError(f"{name} must be a canonical integer string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise WeilAuditVerificationError(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise WeilAuditVerificationError(f"{name} is not canonically encoded")
    if minimum is not None and parsed < minimum:
        raise WeilAuditVerificationError(f"{name} must be at least {minimum}")
    return parsed


def _parse_cell_record(value: Any, name: str) -> WeilSearchCell:
    if not isinstance(value, Mapping):
        raise WeilAuditVerificationError(f"{name} must be an object")
    try:
        cutoff = value["cutoff_c"]
        if not isinstance(cutoff, Mapping):
            raise WeilAuditVerificationError(f"{name} cutoff must be an object")
        cell = WeilSearchCell(
            _parse_canonical_integer(
                cutoff["numerator"], f"{name} cutoff numerator", 1
            ),
            _parse_canonical_integer(
                cutoff["denominator"], f"{name} cutoff denominator", 1
            ),
            _parse_canonical_integer(value["degree"], f"{name} degree", 0),
        )
    except KeyError as exc:
        raise WeilAuditVerificationError(f"{name} is incomplete") from exc
    if dict(value) != cell.to_record():
        raise WeilAuditVerificationError(f"{name} is not canonical")
    return cell


def _parse_policy_record(value: Any) -> WeilSearchPolicy:
    if not isinstance(value, Mapping):
        raise WeilAuditVerificationError("parity audit policy must be an object")
    try:
        witness = value["witness_policy"]
        if not isinstance(witness, Mapping):
            raise WeilAuditVerificationError("witness policy must be an object")
        policy = WeilSearchPolicy(
            attempt_bits=tuple(
                _parse_canonical_integer(item, "attempt precision", 96)
                for item in value["attempt_bits"]
            ),
            confirmation_bits=_parse_canonical_integer(
                value["confirmation_bits"], "confirmation precision", 96
            ),
            eigenpair_count=_parse_canonical_integer(
                witness["eigenpair_count"], "eigenpair count", 1
            ),
            scale_bits=tuple(
                _parse_canonical_integer(item, "scale bits", 1)
                for item in witness["scale_bits"]
            ),
        )
    except (KeyError, TypeError) as exc:
        raise WeilAuditVerificationError("parity audit policy is malformed") from exc
    except ValueError as exc:
        raise WeilAuditVerificationError("parity audit policy is invalid") from exc
    if dict(value) != policy.to_record():
        raise WeilAuditVerificationError("parity audit policy is not canonical")
    return policy


def _parse_witness_record(
    value: Any, expected_size: int, name: str
) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise WeilAuditVerificationError(f"{name} must be a list")
    witness = tuple(
        _parse_canonical_integer(item, f"{name} entry") for item in value
    )
    try:
        return _validate_integer_witness(witness, expected_size)
    except (TypeError, ValueError) as exc:
        raise WeilAuditVerificationError(f"{name} is invalid") from exc


def _require_audit_payload_hash(artifact: Mapping[str, Any]) -> None:
    supplied = artifact.get("payload_sha256")
    body = {
        key: value for key, value in artifact.items() if key != "payload_sha256"
    }
    if not isinstance(supplied, str) or supplied != content_sha256(body):
        raise WeilAuditVerificationError("audit payload hash mismatch")


def _validate_integer_witness(
    witness: Sequence[int], expected_size: int
) -> tuple[int, ...]:
    values = tuple(witness)
    if len(values) != expected_size:
        raise ValueError("witness length must equal matrix dimension")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in values):
        raise TypeError("witness entries must be integers")
    if not any(values):
        raise ValueError("witness must be nonzero")
    return values


def _matrix_degree(matrix: Sequence[Sequence[Any]]) -> int:
    size = len(matrix)
    if size == 0 or size % 2 == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty, square, and odd-dimensional")
    for row in range(size):
        for column in range(row):
            if arb_to_dyadic(matrix[row][column]) != arb_to_dyadic(
                matrix[column][row]
            ):
                raise ValueError(
                    "matrix must use identical symmetric interval enclosures"
                )
    return (size - 1) // 2


def _quadratic_form(matrix: Sequence[Sequence[Any]], vector: Sequence[Any]) -> Any:
    total = arb(0)
    for row, coefficient in enumerate(vector):
        total += coefficient * coefficient * matrix[row][row]
        for column in range(row + 1, len(vector)):
            total += 2 * coefficient * vector[column] * matrix[row][column]
    return total


def _bilinear_form(
    left: Sequence[Any], matrix: Sequence[Sequence[Any]], right: Sequence[Any]
) -> Any:
    total = arb(0)
    for row, left_coefficient in enumerate(left):
        for column, right_coefficient in enumerate(right):
            total += left_coefficient * matrix[row][column] * right_coefficient
    return total


def _bounds_overlap(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    left_lower, left_upper = _record_bounds(left)
    right_lower, right_upper = _record_bounds(right)
    return left_lower <= right_upper and right_lower <= left_upper


def _bounds_contains(outer: Mapping[str, Any], inner: Mapping[str, Any]) -> bool:
    outer_lower, outer_upper = _record_bounds(outer)
    inner_lower, inner_upper = _record_bounds(inner)
    return outer_lower <= inner_lower and inner_upper <= outer_upper


def _serialize_symmetric_block(
    matrix: Sequence[Sequence[Any]], labels: Sequence[str], precision_bits: int
) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for row, row_label in enumerate(labels):
        for column in range(row, len(labels)):
            entries.append(
                {
                    "row": row_label,
                    "column": labels[column],
                    "value": _stable_arb_record(
                        matrix[row][column], precision_bits
                    ),
                }
            )
    return entries


def _serialize_cross_block(
    matrix: Sequence[Sequence[Any]],
    row_labels: Sequence[str],
    column_labels: Sequence[str],
    precision_bits: int,
) -> list[dict[str, Any]]:
    return [
        {
            "even_coordinate": row_labels[row],
            "odd_coordinate": column_labels[column],
            "value": _stable_arb_record(value, precision_bits),
            "contains_zero": value.contains(0),
        }
        for row, values in enumerate(matrix)
        for column, value in enumerate(values)
    ]


def _basis_columns(degree: int) -> tuple[list[list[Any]], list[list[Any]]]:
    size = 2 * degree + 1
    inverse_root_two = 1 / arb(2).sqrt()
    even = [[arb(0) for _ in range(degree + 1)] for _ in range(size)]
    odd = [[arb(0) for _ in range(degree)] for _ in range(size)]
    even[degree][0] = arb(1)
    for mode in range(1, degree + 1):
        negative_row = degree - mode
        positive_row = degree + mode
        even[negative_row][mode] = inverse_root_two
        even[positive_row][mode] = inverse_root_two
        odd[negative_row][mode - 1] = inverse_root_two
        odd[positive_row][mode - 1] = -inverse_root_two
    return even, odd


def _matrix_times_basis(
    matrix: Sequence[Sequence[Any]], basis: Sequence[Sequence[Any]]
) -> list[list[Any]]:
    column_count = len(basis[0]) if basis else 0
    return [
        [
            sum(
                (matrix[row][column] * basis[column][basis_column]
                 for column in range(len(matrix))),
                arb(0),
            )
            for basis_column in range(column_count)
        ]
        for row in range(len(matrix))
    ]


def _basis_transpose_times(
    basis: Sequence[Sequence[Any]], product: Sequence[Sequence[Any]]
) -> list[list[Any]]:
    column_count = len(basis[0]) if basis else 0
    product_columns = len(product[0]) if product else 0
    return [
        [
            sum(
                (basis[row][left_column] * product[row][right_column]
                 for row in range(len(basis))),
                arb(0),
            )
            for right_column in range(product_columns)
        ]
        for left_column in range(column_count)
    ]


def _identical_symmetric_copy(matrix: Sequence[Sequence[Any]]) -> list[list[Any]]:
    size = len(matrix)
    result = [[arb(0) for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for column in range(row, size):
            value = matrix[row][column]
            result[row][column] = value
            result[column][row] = value
    return result


def build_parity_blocks(matrix: Sequence[Sequence[Any]]) -> dict[str, Any]:
    """Construct explicit orthonormal reversal-even and reversal-odd blocks.

    The full basis is ordered by Fourier mode ``-N,...,N``.  Even coordinates
    are ``e_0`` followed by ``(e_-k + e_k)/sqrt(2)`` for ``k=1,...,N``;
    odd coordinates are ``(e_-k - e_k)/sqrt(2)``.  The returned object contains
    Arb matrices and is intended for in-process auditing, not JSON storage.
    """

    degree = _matrix_degree(matrix)
    even_basis, odd_basis = _basis_columns(degree)
    matrix_even = _matrix_times_basis(matrix, even_basis)
    matrix_odd = _matrix_times_basis(matrix, odd_basis)
    even = _identical_symmetric_copy(
        _basis_transpose_times(even_basis, matrix_even)
    )
    odd = _identical_symmetric_copy(
        _basis_transpose_times(odd_basis, matrix_odd)
    )
    cross = _basis_transpose_times(even_basis, matrix_odd)
    return {
        "degree": degree,
        "even_basis": even_basis,
        "odd_basis": odd_basis,
        "even": even,
        "odd": odd,
        "cross": cross,
    }


def _orthonormality_gate(blocks: Mapping[str, Any]) -> dict[str, Any]:
    even_basis = blocks["even_basis"]
    odd_basis = blocks["odd_basis"]
    checks: list[tuple[str, int, int, Any, int]] = []
    for label, left, right in (
        ("even_gram", even_basis, even_basis),
        ("odd_gram", odd_basis, odd_basis),
        ("even_odd_gram", even_basis, odd_basis),
    ):
        left_columns = len(left[0]) if left else 0
        right_columns = len(right[0]) if right else 0
        for row in range(left_columns):
            for column in range(right_columns):
                value = sum(
                    (left[index][row] * right[index][column]
                     for index in range(len(left))),
                    arb(0),
                )
                expected = 1 if label != "even_odd_gram" and row == column else 0
                checks.append((label, row, column, value, expected))
    failures = [item for item in checks if not item[3].contains(item[4])]
    first_failure = None
    if failures:
        label, row, column, _, expected = failures[0]
        first_failure = {
            "matrix": label,
            "row": str(row),
            "column": str(column),
            "expected": str(expected),
        }
    return {
        "classification": "PASS" if not failures else "FAIL",
        "criterion": "E^T E=I, O^T O=I, and E^T O=0 by interval containment",
        "checked_entries": str(len(checks)),
        "failed_entries": str(len(failures)),
        "first_failure": first_failure,
    }


def _centrosymmetry_gate(
    components: Mapping[str, Any], precision_bits: int
) -> dict[str, Any]:
    indices = components["indices"]
    size = len(indices)
    component_records: list[dict[str, Any]] = []
    all_pass = True
    for symbol, key in _COMPONENTS:
        matrix = components[key]
        failures: list[tuple[int, int, Any, Any]] = []
        for row in range(size):
            for column in range(size):
                reflected = matrix[size - 1 - row][size - 1 - column]
                if not matrix[row][column].overlaps(reflected):
                    failures.append((row, column, matrix[row][column], reflected))
        first_failure = None
        if failures:
            row, column, value, reflected = failures[0]
            first_failure = {
                "m": str(indices[row]),
                "n": str(indices[column]),
                "reflected_m": str(-indices[row]),
                "reflected_n": str(-indices[column]),
                "value": _stable_arb_record(value, precision_bits),
                "reflected_value": _stable_arb_record(
                    reflected, precision_bits
                ),
            }
        passed = not failures
        all_pass = all_pass and passed
        component_records.append(
            {
                "component": symbol,
                "source_matrix": key,
                "classification": "PASS" if passed else "FAIL",
                "checked_entries": str(size * size),
                "failed_entries": str(len(failures)),
                "first_failure": first_failure,
            }
        )
    return {
        "classification": "PASS" if all_pass else "FAIL",
        "criterion": "M[m,n] overlaps M[-m,-n] for P, R, S, and A=P-R-S",
        "components": component_records,
    }


def _cross_parity_gate(
    cross: Sequence[Sequence[Any]], precision_bits: int
) -> dict[str, Any]:
    entries = [value for row in cross for value in row]
    failures = [index for index, value in enumerate(entries) if not value.contains(0)]
    return {
        "classification": "PASS" if not failures else "FAIL",
        "criterion": "every entry of E^T A O contains zero",
        "checked_entries": str(len(entries)),
        "failed_entries": str(len(failures)),
        "entries": _serialize_cross_block(
            cross,
            [str(value) for value in range(len(cross))],
            [str(value + 1) for value in range(len(cross[0]))]
            if cross
            else [],
            precision_bits,
        ),
    }


def _rump_diagnostic(
    matrix: Sequence[Sequence[Any]], precision_bits: int
) -> dict[str, Any]:
    if not matrix:
        return {
            "role": "DIAGNOSTIC_ONLY",
            "classification": "EMPTY_BLOCK",
            "spectrum": None,
        }
    try:
        spectrum = certify_rump_spectrum(matrix)
    except (ArithmeticError, RuntimeError, ValueError):
        return {
            "role": "DIAGNOSTIC_ONLY",
            "classification": "INCONCLUSIVE",
            "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
            "spectrum": None,
        }
    return {
        "role": "DIAGNOSTIC_ONLY",
        "classification": "AVAILABLE",
        "spectrum": _stabilize_ball_records(spectrum, precision_bits),
    }


def _lift_parity_direction(
    direction: Sequence[Fraction], parity: str, degree: int
) -> tuple[Fraction, ...]:
    root_two = _midpoint_fraction(arb(2).sqrt())
    if root_two == 0:
        raise ArithmeticError("sqrt(2) midpoint vanished")
    inverse_root_two = Fraction(1, 1) / root_two
    full = [Fraction(0, 1) for _ in range(2 * degree + 1)]
    if parity == "even":
        full[degree] = direction[0]
        for mode in range(1, degree + 1):
            value = direction[mode] * inverse_root_two
            full[degree - mode] = value
            full[degree + mode] = value
    elif parity == "odd":
        for mode in range(1, degree + 1):
            value = direction[mode - 1] * inverse_root_two
            full[degree - mode] = value
            full[degree + mode] = -value
    else:
        raise ValueError("parity must be 'even' or 'odd'")
    return tuple(full)


def _candidate_hints(
    blocks: Mapping[str, Any], policy: WeilSearchPolicy
) -> list[dict[str, Any]]:
    degree = blocks["degree"]
    hints: list[dict[str, Any]] = []
    seen: set[tuple[int, ...]] = set()
    for parity in ("even", "odd"):
        block = blocks[parity]
        if not block:
            continue
        eigenvalues, right_vectors = arb_mat(block).eig(
            right=True, algorithm="approx"
        )
        order = sorted(
            range(len(eigenvalues)),
            key=lambda index: (
                _midpoint_fraction(eigenvalues[index].real),
                index,
            ),
        )[: min(policy.eigenpair_count, len(eigenvalues))]
        for rank, eigen_index in enumerate(order):
            eigenvalue_midpoint = _midpoint_fraction(
                eigenvalues[eigen_index].real
            )
            for part_name, attribute_name in (
                ("real", "real"),
                ("imaginary", "imag"),
            ):
                direction = tuple(
                    _midpoint_fraction(
                        getattr(right_vectors[row, eigen_index], attribute_name)
                    )
                    for row in range(len(block))
                )
                if not any(direction):
                    continue
                lifted = _lift_parity_direction(direction, parity, degree)
                pivot = max(
                    range(len(lifted)),
                    key=lambda index: (abs(lifted[index]), -index),
                )
                pivot_value = lifted[pivot]
                normalized = tuple(value / pivot_value for value in lifted)
                for scale_bits in policy.scale_bits:
                    scale = 1 << scale_bits
                    rounded = tuple(
                        _round_ties_even(value * scale) for value in normalized
                    )
                    witness = _primitive_witness(rounded)
                    if witness in seen:
                        continue
                    seen.add(witness)
                    hints.append(
                        {
                            "parity": parity.upper(),
                            "eigen_rank": rank,
                            "eigen_original_index": eigen_index,
                            "eigenvalue_real_midpoint": _fraction_record(
                                eigenvalue_midpoint
                            ),
                            "source_part": part_name,
                            "scale_bits": scale_bits,
                            "pivot_index": pivot,
                            "witness": witness,
                        }
                    )
    return hints


def _audit_full_vs_parity_quadratic(
    matrix: Sequence[Sequence[Any]],
    blocks: Mapping[str, Any],
    witness: Sequence[int],
    precision_bits: int,
) -> dict[str, Any]:
    degree = blocks["degree"]
    values = _validate_integer_witness(witness, 2 * degree + 1)
    primitive = _primitive_witness(values)
    root_two = arb(2).sqrt()
    even_coordinates = [arb(values[degree])]
    odd_coordinates: list[Any] = []
    for mode in range(1, degree + 1):
        negative = values[degree - mode]
        positive = values[degree + mode]
        even_coordinates.append(arb(negative + positive) / root_two)
        odd_coordinates.append(arb(negative - positive) / root_two)

    full_quadratic = _quadratic_form(matrix, [arb(value) for value in values])
    even_quadratic = _quadratic_form(blocks["even"], even_coordinates)
    odd_quadratic = (
        _quadratic_form(blocks["odd"], odd_coordinates)
        if odd_coordinates
        else arb(0)
    )
    cross_quadratic = (
        2 * _bilinear_form(even_coordinates, blocks["cross"], odd_coordinates)
        if odd_coordinates
        else arb(0)
    )
    parity_quadratic = even_quadratic + odd_quadratic + cross_quadratic
    norm_squared = sum(value * value for value in values)
    transformed_norm = sum(value * value for value in even_coordinates) + sum(
        value * value for value in odd_coordinates
    )
    quadratic_overlap = full_quadratic.overlaps(parity_quadratic)
    norm_contains = transformed_norm.contains(norm_squared)
    direct = _stable_evaluation(
        evaluate_integer_rayleigh(matrix, values), precision_bits
    )
    return {
        "classification": "PASS"
        if quadratic_overlap and norm_contains
        else "FAIL",
        "witness": [str(value) for value in values],
        "primitive_witness": [str(value) for value in primitive],
        "input_is_primitive_and_canonically_signed": values == primitive,
        "sign_authority": "direct_full_matrix_interval_evaluation",
        "authoritative_full_matrix_evaluation": direct,
        "full_quadratic": _stable_arb_record(full_quadratic, precision_bits),
        "parity_quadratic": _stable_arb_record(
            parity_quadratic, precision_bits
        ),
        "even_contribution": _stable_arb_record(
            even_quadratic, precision_bits
        ),
        "odd_contribution": _stable_arb_record(
            odd_quadratic, precision_bits
        ),
        "twice_cross_contribution": _stable_arb_record(
            cross_quadratic, precision_bits
        ),
        "full_and_parity_quadratics_overlap": quadratic_overlap,
        "exact_full_norm_squared": str(norm_squared),
        "transformed_norm": _stable_arb_record(
            transformed_norm, precision_bits
        ),
        "transformed_norm_contains_exact_norm": norm_contains,
    }


def audit_full_vs_parity_quadratic(
    matrix: Sequence[Sequence[Any]],
    witness: Sequence[int],
    precision_bits: int,
) -> dict[str, Any]:
    """Audit ``x^T A x`` against its explicit even/odd decomposition."""

    precision_bits = _validate_precision(precision_bits)
    with _working_precision(precision_bits):
        blocks = build_parity_blocks(matrix)
        return _audit_full_vs_parity_quadratic(
            matrix, blocks, witness, precision_bits
        )


def _evaluated_parity_candidates(
    matrix: Sequence[Sequence[Any]],
    blocks: Mapping[str, Any],
    policy: WeilSearchPolicy,
    precision_bits: int,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for hint in _candidate_hints(blocks, policy):
        witness = hint["witness"]
        identity = _audit_full_vs_parity_quadratic(
            matrix, blocks, witness, precision_bits
        )
        source = {
            "parity": hint["parity"],
            "eigen_rank": str(hint["eigen_rank"]),
            "eigen_original_index": str(hint["eigen_original_index"]),
            "eigenvalue_real_midpoint": hint["eigenvalue_real_midpoint"],
            "source_part": hint["source_part"],
            "scale_bits": str(hint["scale_bits"]),
            "pivot_index": str(hint["pivot_index"]),
            "role": "APPROXIMATE_UNTRUSTED_HINT_ONLY",
        }
        candidate_identity = {
            "source": source,
            "witness": [str(value) for value in witness],
        }
        records.append(
            {
                "candidate_id": content_sha256(candidate_identity),
                **candidate_identity,
                "sign_authority": "direct_full_matrix_interval_evaluation",
                "authoritative_full_matrix_evaluation": identity[
                    "authoritative_full_matrix_evaluation"
                ],
                "full_vs_parity_quadratic_audit": identity,
            }
        )
    return records


def extract_parity_integer_witnesses(
    matrix: Sequence[Sequence[Any]],
    precision_bits: int,
    policy: WeilSearchPolicy = WeilSearchPolicy(),
) -> list[dict[str, Any]]:
    """Generate parity-block hints and evaluate exact integer lifts directly.

    The approximate block eigenvectors only choose integer vectors.  Each
    returned record includes the authoritative direct full-matrix interval
    evaluation and a full-versus-parity quadratic overlap audit.
    """

    precision_bits = _validate_precision(precision_bits)
    if not isinstance(policy, WeilSearchPolicy):
        raise TypeError("policy must be a WeilSearchPolicy")
    with _working_precision(precision_bits):
        blocks = build_parity_blocks(matrix)
        return _evaluated_parity_candidates(
            matrix, blocks, policy, precision_bits
        )


def _certify_parity_audit(
    cell: WeilSearchCell,
    precision_bits: int = 192,
    policy: WeilSearchPolicy = WeilSearchPolicy(),
    witnesses: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    """Build and audit the reversal-parity decomposition of one Weil cell."""

    if not isinstance(cell, WeilSearchCell):
        raise TypeError("cell must be a WeilSearchCell")
    if not isinstance(policy, WeilSearchPolicy):
        raise TypeError("policy must be a WeilSearchPolicy")
    precision_bits = _validate_precision(precision_bits)
    base = {
        "schema": PARITY_AUDIT_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "algorithm_id": PARITY_ALGORITHM_ID,
        "cell": cell.to_record(),
        "precision_bits": str(precision_bits),
        "backend": _backend_record(),
        "policy": policy.to_record(),
        "limitation": EXPLORATORY_LIMITATION,
    }
    with _working_precision(precision_bits):
        try:
            components = _symmetric_components(
                cell.cutoff_numerator,
                cell.cutoff_denominator,
                cell.degree,
            )
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": (
                            "SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE"
                        ),
                    },
                    "decision": "AUDIT_INCONCLUSIVE",
                }
            )

        centrosymmetry = _centrosymmetry_gate(components, precision_bits)
        blocks = build_parity_blocks(components["q"])
        orthonormality = _orthonormality_gate(blocks)
        cross_gate = _cross_parity_gate(blocks["cross"], precision_bits)
        even_labels = ["0", *[str(value) for value in range(1, cell.degree + 1)]]
        odd_labels = [str(value) for value in range(1, cell.degree + 1)]
        spectra = {
            "role": "DIAGNOSTIC_ONLY_NOT_A_SIGN_CERTIFICATE",
            "even": _rump_diagnostic(blocks["even"], precision_bits),
            "odd": _rump_diagnostic(blocks["odd"], precision_bits),
        }
        candidate_failure = None
        try:
            candidates = _evaluated_parity_candidates(
                components["q"], blocks, policy, precision_bits
            )
        except (ArithmeticError, RuntimeError, ValueError) as exc:
            candidates = []
            candidate_failure = {
                "classification": "INCONCLUSIVE",
                "failure_code": "APPROXIMATE_EIGENVECTOR_HINT_FAILED",
                "error_type": type(exc).__name__,
            }
        requested_witnesses: list[dict[str, Any]] = []
        for witness in witnesses:
            supplied = _validate_integer_witness(
                witness, len(components["q"])
            )
            primitive = _primitive_witness(supplied)
            requested = _audit_full_vs_parity_quadratic(
                components["q"], blocks, primitive, precision_bits
            )
            requested_witnesses.append(
                {
                    **requested,
                    "supplied_witness": [str(value) for value in supplied],
                    "normalization": (
                        "primitive-gcd-and-first-nonzero-positive"
                    ),
                    "supplied_witness_was_canonical": supplied == primitive,
                }
            )

        structural_pass = all(
            gate["classification"] == "PASS"
            for gate in (centrosymmetry, orthonormality, cross_gate)
        ) and all(
            record["classification"] == "PASS"
            for record in requested_witnesses
        ) and all(
            record["full_vs_parity_quadratic_audit"]["classification"]
            == "PASS"
            for record in candidates
        )
        negative_ids = [
            record["candidate_id"]
            for record in candidates
            if record["authoritative_full_matrix_evaluation"]["sign"]
            == "NEGATIVE"
        ]
        requested_negative = [
            str(index)
            for index, record in enumerate(requested_witnesses)
            if record["authoritative_full_matrix_evaluation"]["sign"]
            == "NEGATIVE"
        ]
        decision = (
            "AUDIT_FAILED"
            if not structural_pass
            else "NEGATIVE_CANDIDATE_QUARANTINED"
            if negative_ids or requested_negative
            else "AUDIT_PASSED_EXPLORATORY"
        )
        payload = {
            **base,
            "matrix_construction": {
                "classification": "COMPLETE",
                "normalization": "A=P-R-S",
            },
            "basis": {
                "full_order": [
                    str(value) for value in range(-cell.degree, cell.degree + 1)
                ],
                "even_order": even_labels,
                "odd_order": odd_labels,
                "even_vectors": (
                    "e_0; then (e_-k+e_k)/sqrt(2), k=1,...,N"
                ),
                "odd_vectors": "(e_-k-e_k)/sqrt(2), k=1,...,N",
                "even_dimension": str(cell.degree + 1),
                "odd_dimension": str(cell.degree),
            },
            "gates": {
                "centrosymmetry_overlap": centrosymmetry,
                "orthonormal_basis": orthonormality,
                "cross_parity_contains_zero": cross_gate,
            },
            "blocks": {
                "even": {
                    "dimension": str(cell.degree + 1),
                    "entries_upper_triangle": _serialize_symmetric_block(
                        blocks["even"], even_labels, precision_bits
                    ),
                },
                "odd": {
                    "dimension": str(cell.degree),
                    "entries_upper_triangle": _serialize_symmetric_block(
                        blocks["odd"], odd_labels, precision_bits
                    ),
                },
            },
            "rump_block_spectra": spectra,
            "candidate_generation": {
                "classification": "AVAILABLE"
                if candidate_failure is None
                else "INCONCLUSIVE",
                "role": "APPROXIMATE_HINTS_ONLY",
                "policy": policy.to_record(),
                "failure": candidate_failure,
                "candidate_count": str(len(candidates)),
                "negative_candidate_ids": negative_ids,
                "candidates": candidates,
            },
            "requested_witness_audits": requested_witnesses,
            "requested_negative_witness_indices": requested_negative,
            "decision": decision,
        }
        return _with_payload_hash(payload)


def certify_parity_audit(
    cell: WeilSearchCell,
    precision_bits: int = 192,
    policy: WeilSearchPolicy = WeilSearchPolicy(),
    witnesses: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    """Build one parity artifact from the canonical clean FLINT state."""

    return _isolated_flint_operation(
        lambda: _certify_parity_audit(
            cell, precision_bits, policy, witnesses
        )
    )


def _component_snapshot(
    components: Mapping[str, Any], precision_bits: int
) -> dict[str, dict[tuple[int, int], dict[str, Any]]]:
    indices = components["indices"]
    snapshot: dict[str, dict[tuple[int, int], dict[str, Any]]] = {}
    for symbol, key in _COMPONENTS:
        entries: dict[tuple[int, int], dict[str, Any]] = {}
        for row, m in enumerate(indices):
            for column in range(row, len(indices)):
                n = indices[column]
                entries[(m, n)] = _stable_arb_record(
                    components[key][row][column], precision_bits
                )
        snapshot[symbol] = entries
    return snapshot


def _raw_component_entry(
    components: Mapping[str, Any], key: str, m: int, n: int
) -> Any:
    positions = {value: index for index, value in enumerate(components["indices"])}
    return components[key][positions[m]][positions[n]]


def _nesting_component_records(
    lower: Mapping[str, Any],
    higher: Mapping[str, Any],
    lower_replay: Mapping[str, Any],
    higher_replay: Mapping[str, Any],
    lower_snapshot: Mapping[str, Mapping[tuple[int, int], Mapping[str, Any]]],
    higher_snapshot: Mapping[str, Mapping[tuple[int, int], Mapping[str, Any]]],
    lower_replay_snapshot: Mapping[
        str, Mapping[tuple[int, int], Mapping[str, Any]]
    ],
    higher_replay_snapshot: Mapping[
        str, Mapping[tuple[int, int], Mapping[str, Any]]
    ],
) -> tuple[list[dict[str, Any]], bool]:
    indices = lower["indices"]
    records: list[dict[str, Any]] = []
    all_pass = True
    for symbol, key in _COMPONENTS:
        entries: list[dict[str, Any]] = []
        component_pass = True
        for row, m in enumerate(indices):
            for n in indices[row:]:
                lower_ball = lower_snapshot[symbol][(m, n)]
                higher_ball = higher_snapshot[symbol][(m, n)]
                lower_replay_ball = lower_replay_snapshot[symbol][(m, n)]
                higher_replay_ball = higher_replay_snapshot[symbol][(m, n)]
                base_overlap = _raw_component_entry(
                    lower, key, m, n
                ).overlaps(_raw_component_entry(higher, key, m, n))
                replay_overlap = _raw_component_entry(
                    lower_replay, key, m, n
                ).overlaps(_raw_component_entry(higher_replay, key, m, n))
                containment = {
                    "lower_contains_lower_replay": _bounds_contains(
                        lower_ball, lower_replay_ball
                    ),
                    "higher_contains_higher_replay": _bounds_contains(
                        higher_ball, higher_replay_ball
                    ),
                    "lower_contains_higher_central_replay": _bounds_contains(
                        lower_ball, higher_replay_ball
                    ),
                    "higher_central_contains_lower_replay": _bounds_contains(
                        higher_ball, lower_replay_ball
                    ),
                }
                passed = base_overlap and replay_overlap and all(containment.values())
                component_pass = component_pass and passed
                entries.append(
                    {
                        "m": str(m),
                        "n": str(n),
                        "lower_degree": lower_ball,
                        "higher_degree_central": higher_ball,
                        "lower_degree_replay": lower_replay_ball,
                        "higher_degree_central_replay": higher_replay_ball,
                        "same_precision_overlap": base_overlap,
                        "replay_precision_overlap": replay_overlap,
                        "containment": containment,
                        "classification": "PASS" if passed else "FAIL",
                    }
                )
        all_pass = all_pass and component_pass
        records.append(
            {
                "component": symbol,
                "source_matrix": key,
                "classification": "PASS" if component_pass else "FAIL",
                "compared_upper_triangle_entries": str(len(entries)),
                "entries_sha256": content_sha256(entries),
                "entries": entries,
            }
        )
    return records, all_pass


def _global_minimum_from_block_diagnostics(
    diagnostics: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    candidates: list[tuple[str, Mapping[str, Any]]] = []
    for parity in ("even", "odd"):
        diagnostic = diagnostics[parity]
        if diagnostic.get("classification") == "EMPTY_BLOCK":
            continue
        spectrum = diagnostic.get("spectrum")
        if not isinstance(spectrum, Mapping):
            return {
                "classification": "INCONCLUSIVE",
                "failure_code": f"{parity.upper()}_BLOCK_SPECTRUM_UNAVAILABLE",
            }
        smallest = spectrum.get("smallest_eigenvalue")
        if not isinstance(smallest, Mapping) or not smallest.get(
            "separated_from_all_others"
        ):
            return {
                "classification": "INCONCLUSIVE",
                "failure_code": (
                    f"{parity.upper()}_BLOCK_MINIMUM_NOT_SEPARATED"
                ),
            }
        imaginary = smallest.get("imaginary")
        if not isinstance(imaginary, Mapping):
            return {
                "classification": "INCONCLUSIVE",
                "failure_code": f"{parity.upper()}_BLOCK_IMAGINARY_MISSING",
            }
        imaginary_lower, imaginary_upper = _record_bounds(imaginary)
        if not (imaginary_lower <= 0 <= imaginary_upper):
            return {
                "classification": "INCONCLUSIVE",
                "failure_code": f"{parity.upper()}_BLOCK_MINIMUM_NOT_REAL",
            }
        candidates.append((parity, smallest["real"]))

    if not candidates:
        return {
            "classification": "INCONCLUSIVE",
            "failure_code": "NO_BLOCK_MINIMUM_AVAILABLE",
        }
    if len(candidates) == 1:
        parity, value = candidates[0]
        return {
            "classification": "SEPARATED_GLOBAL_MINIMUM",
            "parity": parity.upper(),
            "value": value,
        }
    (left_parity, left), (right_parity, right) = candidates
    left_lower, left_upper = _record_bounds(left)
    right_lower, right_upper = _record_bounds(right)
    if left_upper < right_lower:
        return {
            "classification": "SEPARATED_GLOBAL_MINIMUM",
            "parity": left_parity.upper(),
            "value": left,
        }
    if right_upper < left_lower:
        return {
            "classification": "SEPARATED_GLOBAL_MINIMUM",
            "parity": right_parity.upper(),
            "value": right,
        }
    return {
        "classification": "INCONCLUSIVE",
        "failure_code": "EVEN_AND_ODD_BLOCK_MINIMA_NOT_SEPARATED",
        "even_minimum": left if left_parity == "even" else right,
        "odd_minimum": right if right_parity == "odd" else left,
    }


def _rayleigh_ritz_diagnostic(
    lower_diagnostics: Mapping[str, Mapping[str, Any]],
    higher_diagnostics: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    lower = _global_minimum_from_block_diagnostics(lower_diagnostics)
    higher = _global_minimum_from_block_diagnostics(higher_diagnostics)
    base = {
        "role": "DIAGNOSTIC_GUARD_ONLY",
        "theorem": "lambda_min(A_M) <= lambda_min(A_N) for M>N",
        "lower_degree_minimum": lower,
        "higher_degree_minimum": higher,
    }
    if (
        lower["classification"] != "SEPARATED_GLOBAL_MINIMUM"
        or higher["classification"] != "SEPARATED_GLOBAL_MINIMUM"
    ):
        return {
            **base,
            "classification": "INCONCLUSIVE",
            "audit_failure": False,
        }
    lower_lower, lower_upper = _record_bounds(lower["value"])
    higher_lower, higher_upper = _record_bounds(higher["value"])
    if lower_upper < higher_lower:
        return {
            **base,
            "classification": "FAIL_SEPARATED_IMPOSSIBLE_ORDERING",
            "audit_failure": True,
        }
    if higher_upper < lower_lower:
        classification = "PASS_SEPARATED_EXPECTED_ORDERING"
    else:
        classification = "CONSISTENT_NOT_SEPARATED"
    return {
        **base,
        "classification": classification,
        "audit_failure": False,
    }


def _zero_padding_audit(
    lower_matrix: Sequence[Sequence[Any]],
    higher_matrix: Sequence[Sequence[Any]],
    witness: Sequence[int],
    precision_bits: int,
) -> dict[str, Any]:
    lower_degree = _matrix_degree(lower_matrix)
    higher_degree = _matrix_degree(higher_matrix)
    if higher_degree <= lower_degree:
        raise ValueError("higher matrix degree must be strictly greater")
    supplied = _validate_integer_witness(witness, len(lower_matrix))
    values = _primitive_witness(supplied)
    degree_gap = higher_degree - lower_degree
    padded = (0,) * degree_gap + values + (0,) * degree_gap
    lower_evaluation = _stable_evaluation(
        evaluate_integer_rayleigh(lower_matrix, values), precision_bits
    )
    higher_evaluation = _stable_evaluation(
        evaluate_integer_rayleigh(higher_matrix, padded), precision_bits
    )
    lower_quadratic = _quadratic_form(
        lower_matrix, [arb(value) for value in values]
    )
    higher_quadratic = _quadratic_form(
        higher_matrix, [arb(value) for value in padded]
    )
    numerator_overlap = lower_quadratic.overlaps(higher_quadratic)
    norm_squared = sum(value * value for value in values)
    lower_quotient = lower_quadratic / norm_squared
    higher_quotient = higher_quadratic / norm_squared
    quotient_overlap = lower_quotient.overlaps(higher_quotient)
    lower_negative = lower_evaluation["sign"] == "NEGATIVE"
    higher_negative = higher_evaluation["sign"] == "NEGATIVE"
    evidence_overlap = numerator_overlap and quotient_overlap
    if not evidence_overlap:
        classification = "FAIL_EVIDENCE_NONOVERLAP"
        audit_failure = True
        audit_inconclusive = False
    elif lower_negative and higher_negative:
        classification = "PASS"
        audit_failure = False
        audit_inconclusive = False
    elif lower_negative and higher_evaluation["sign"] == "POSITIVE":
        classification = "FAIL_INCOMPATIBLE_CERTIFIED_SIGNS"
        audit_failure = True
        audit_inconclusive = False
    elif lower_negative:
        classification = "INCONCLUSIVE_HIGHER_DEGREE_SIGN"
        audit_failure = False
        audit_inconclusive = True
    elif lower_evaluation["sign"] == "INCONCLUSIVE":
        classification = "INCONCLUSIVE_LOWER_DEGREE_SIGN"
        audit_failure = False
        audit_inconclusive = True
    else:
        classification = "NOT_APPLICABLE_LOWER_WITNESS_POSITIVE"
        audit_failure = False
        audit_inconclusive = False
    return {
        "classification": classification,
        "audit_failure": audit_failure,
        "audit_inconclusive": audit_inconclusive,
        "precision_bits": str(precision_bits),
        "degree_gap": str(degree_gap),
        "supplied_witness": [str(value) for value in supplied],
        "witness": [str(value) for value in values],
        "normalization": "primitive-gcd-and-first-nonzero-positive",
        "supplied_witness_was_canonical": supplied == values,
        "zero_padded_witness": [str(value) for value in padded],
        "sign_authority": "direct_full_matrix_interval_evaluation",
        "lower_degree_evaluation": lower_evaluation,
        "higher_degree_evaluation": higher_evaluation,
        "lower_and_zero_padded_numerators_overlap": numerator_overlap,
        "lower_and_zero_padded_rayleigh_quotients_overlap": quotient_overlap,
        "lower_numerator": _stable_arb_record(
            lower_quadratic, precision_bits
        ),
        "zero_padded_higher_numerator": _stable_arb_record(
            higher_quadratic, precision_bits
        ),
        "lower_rayleigh_quotient": _stable_arb_record(
            lower_quotient, precision_bits
        ),
        "zero_padded_higher_rayleigh_quotient": _stable_arb_record(
            higher_quotient, precision_bits
        ),
    }


def audit_zero_padding_negative_witness(
    lower_matrix: Sequence[Sequence[Any]],
    higher_matrix: Sequence[Sequence[Any]],
    witness: Sequence[int],
    precision_bits: int,
) -> dict[str, Any]:
    """Check zero-padded transport of a negative principal-matrix witness."""

    precision_bits = _validate_precision(precision_bits)
    with _working_precision(precision_bits):
        return _zero_padding_audit(
            lower_matrix, higher_matrix, witness, precision_bits
        )


def _replay_parity_diagnostics(
    components: Mapping[str, Any], precision_bits: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    centrosymmetry = _centrosymmetry_gate(components, precision_bits)
    blocks = build_parity_blocks(components["q"])
    cross = _cross_parity_gate(blocks["cross"], precision_bits)
    diagnostics = {
        "role": "DIAGNOSTIC_ONLY_NOT_A_SIGN_CERTIFICATE",
        "even": _rump_diagnostic(blocks["even"], precision_bits),
        "odd": _rump_diagnostic(blocks["odd"], precision_bits),
    }
    gates = {
        "centrosymmetry_overlap": centrosymmetry,
        "cross_parity_contains_zero": cross,
    }
    return gates, diagnostics


def certify_degree_nesting_audit(
    lower_cell: WeilSearchCell,
    higher_cell: WeilSearchCell,
    precision_bits: int = 192,
    replay_precision_bits: int = 384,
    negative_witnesses: Sequence[Sequence[int]] = (),
) -> dict[str, Any]:
    """Audit ordered-degree principal nesting at one exact rational cutoff."""

    if not isinstance(lower_cell, WeilSearchCell) or not isinstance(
        higher_cell, WeilSearchCell
    ):
        raise TypeError("lower_cell and higher_cell must be WeilSearchCell values")
    if (
        lower_cell.cutoff_numerator != higher_cell.cutoff_numerator
        or lower_cell.cutoff_denominator != higher_cell.cutoff_denominator
    ):
        raise ValueError("degree-nesting cells must use the same reduced cutoff")
    if higher_cell.degree <= lower_cell.degree:
        raise ValueError("higher nesting degree must be strictly greater")
    precision_bits = _validate_precision(precision_bits)
    replay_precision_bits = _validate_precision(
        replay_precision_bits, "replay_precision_bits"
    )
    if replay_precision_bits <= precision_bits:
        raise ValueError("replay_precision_bits must exceed precision_bits")

    base = {
        "schema": NESTING_AUDIT_SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "algorithm_id": NESTING_ALGORITHM_ID,
        "lower_cell": lower_cell.to_record(),
        "higher_cell": higher_cell.to_record(),
        "precision_bits": str(precision_bits),
        "replay_precision_bits": str(replay_precision_bits),
        "backend": _backend_record(),
        "limitation": EXPLORATORY_LIMITATION,
    }

    with _working_precision(precision_bits):
        try:
            lower = _symmetric_components(
                lower_cell.cutoff_numerator,
                lower_cell.cutoff_denominator,
                lower_cell.degree,
            )
            higher = _symmetric_components(
                higher_cell.cutoff_numerator,
                higher_cell.cutoff_denominator,
                higher_cell.degree,
            )
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": (
                            "BASE_SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE"
                        ),
                    },
                    "decision": "AUDIT_INCONCLUSIVE",
                }
            )
        lower_snapshot = _component_snapshot(lower, precision_bits)
        higher_snapshot = _component_snapshot(higher, precision_bits)
        base_zero_padding = [
            _zero_padding_audit(
                lower["q"], higher["q"], witness, precision_bits
            )
            for witness in negative_witnesses
        ]

    with _working_precision(replay_precision_bits):
        try:
            lower_replay = _symmetric_components(
                lower_cell.cutoff_numerator,
                lower_cell.cutoff_denominator,
                lower_cell.degree,
            )
            higher_replay = _symmetric_components(
                higher_cell.cutoff_numerator,
                higher_cell.cutoff_denominator,
                higher_cell.degree,
            )
        except ArithmeticError:
            return _with_payload_hash(
                {
                    **base,
                    "matrix_construction": {
                        "classification": "INCONCLUSIVE",
                        "failure_code": (
                            "REPLAY_SOURCE_ORACLE_OR_QUADRATURE_INCONCLUSIVE"
                        ),
                    },
                    "decision": "AUDIT_INCONCLUSIVE",
                }
            )
        lower_replay_snapshot = _component_snapshot(
            lower_replay, replay_precision_bits
        )
        higher_replay_snapshot = _component_snapshot(
            higher_replay, replay_precision_bits
        )
        lower_parity_gates, lower_spectra = _replay_parity_diagnostics(
            lower_replay, replay_precision_bits
        )
        higher_parity_gates, higher_spectra = _replay_parity_diagnostics(
            higher_replay, replay_precision_bits
        )
        replay_zero_padding = [
            _zero_padding_audit(
                lower_replay["q"],
                higher_replay["q"],
                witness,
                replay_precision_bits,
            )
            for witness in negative_witnesses
        ]

    component_records, nesting_pass = _nesting_component_records(
        lower,
        higher,
        lower_replay,
        higher_replay,
        lower_snapshot,
        higher_snapshot,
        lower_replay_snapshot,
        higher_replay_snapshot,
    )
    parity_gate_pass = all(
        gate["classification"] == "PASS"
        for gate in (
            lower_parity_gates["centrosymmetry_overlap"],
            lower_parity_gates["cross_parity_contains_zero"],
            higher_parity_gates["centrosymmetry_overlap"],
            higher_parity_gates["cross_parity_contains_zero"],
        )
    )
    rayleigh_ritz = _rayleigh_ritz_diagnostic(
        lower_spectra, higher_spectra
    )

    zero_padding_records: list[dict[str, Any]] = []
    zero_padding_failure = False
    zero_padding_inconclusive = False
    negative_quarantined = False
    for base_record, replay_record in zip(
        base_zero_padding, replay_zero_padding
    ):
        negative_at_either_precision = any(
            record["lower_degree_evaluation"]["sign"] == "NEGATIVE"
            for record in (base_record, replay_record)
        )
        certified_signs = {
            record["lower_degree_evaluation"]["sign"]
            for record in (base_record, replay_record)
            if record["lower_degree_evaluation"]["sign"]
            in {"NEGATIVE", "POSITIVE"}
        }
        incompatible_replay_signs = certified_signs == {
            "NEGATIVE",
            "POSITIVE",
        }
        negative_quarantined = (
            negative_quarantined or negative_at_either_precision
        )
        hard_failure = incompatible_replay_signs or any(
            record["audit_failure"]
            for record in (base_record, replay_record)
        )
        if hard_failure:
            classification = "FAIL"
            zero_padding_failure = True
        elif negative_at_either_precision and all(
            record["classification"] == "PASS"
            for record in (base_record, replay_record)
        ):
            classification = "PASS"
        elif negative_at_either_precision or any(
            record["audit_inconclusive"]
            for record in (base_record, replay_record)
        ):
            classification = "INCONCLUSIVE"
            zero_padding_inconclusive = True
        else:
            classification = "NOT_APPLICABLE_LOWER_WITNESS_POSITIVE"
        zero_padding_records.append(
            {
                "classification": classification,
                "incompatible_certified_replay_signs": (
                    incompatible_replay_signs
                ),
                "candidate_state": "QUARANTINED"
                if negative_at_either_precision
                else "NOT_A_CERTIFIED_NEGATIVE_CANDIDATE",
                "base_precision": base_record,
                "replay_precision": replay_record,
            }
        )

    structural_pass = (
        nesting_pass
        and parity_gate_pass
        and not zero_padding_failure
        and not rayleigh_ritz["audit_failure"]
    )
    decision = (
        "AUDIT_FAILED"
        if not structural_pass
        else "AUDIT_INCONCLUSIVE"
        if zero_padding_inconclusive
        else "NEGATIVE_CANDIDATE_QUARANTINED"
        if negative_quarantined
        else "AUDIT_PASSED_EXPLORATORY"
    )
    payload = {
        **base,
        "matrix_construction": {
            "classification": "COMPLETE_AT_BOTH_PRECISIONS",
            "normalization": "A=P-R-S",
        },
        "principal_submatrix_relation": {
            "description": (
                "the matrix for degree N is the central principal "
                "submatrix on modes -N,...,N of every higher-degree "
                "matrix A_M with M>N"
            ),
            "classification": "PASS" if nesting_pass else "FAIL",
            "components": component_records,
        },
        "replay_parity_gates": {
            "lower_degree": lower_parity_gates,
            "higher_degree": higher_parity_gates,
        },
        "rump_block_spectra": {
            "role": "DIAGNOSTIC_ONLY_NOT_A_SIGN_CERTIFICATE",
            "precision_bits": str(replay_precision_bits),
            "lower_degree": lower_spectra,
            "higher_degree": higher_spectra,
        },
        "rayleigh_ritz_minimum_check": rayleigh_ritz,
        "zero_padding_negative_witness_checks": zero_padding_records,
        "negative_candidate_state": "QUARANTINED"
        if negative_quarantined
        else "NONE_FOUND",
        "decision": decision,
    }
    return _with_payload_hash(payload)


def _validate_candidate_source(
    source: Any,
    witness: Sequence[int],
    cell: WeilSearchCell,
    policy: WeilSearchPolicy,
) -> dict[str, Any]:
    if not isinstance(source, Mapping):
        raise WeilAuditVerificationError("candidate source must be an object")
    expected_keys = {
        "parity",
        "eigen_rank",
        "eigen_original_index",
        "eigenvalue_real_midpoint",
        "source_part",
        "scale_bits",
        "pivot_index",
        "role",
    }
    if set(source) != expected_keys:
        raise WeilAuditVerificationError("candidate source fields are malformed")
    parity = source.get("parity")
    if parity not in {"EVEN", "ODD"}:
        raise WeilAuditVerificationError("candidate parity is invalid")
    if source.get("source_part") not in {"real", "imaginary"}:
        raise WeilAuditVerificationError("candidate source part is invalid")
    if source.get("role") != "APPROXIMATE_UNTRUSTED_HINT_ONLY":
        raise WeilAuditVerificationError("candidate source role was promoted")
    rank = _parse_canonical_integer(source["eigen_rank"], "eigen rank", 0)
    eigen_index = _parse_canonical_integer(
        source["eigen_original_index"], "eigen index", 0
    )
    scale_bits = _parse_canonical_integer(
        source["scale_bits"], "candidate scale bits", 1
    )
    pivot = _parse_canonical_integer(
        source["pivot_index"], "candidate pivot index", 0
    )
    block_dimension = cell.degree + 1 if parity == "EVEN" else cell.degree
    if rank >= min(policy.eigenpair_count, block_dimension) or (
        eigen_index >= block_dimension
    ):
        raise WeilAuditVerificationError("candidate eigenmode metadata is out of range")
    if scale_bits not in policy.scale_bits or pivot >= 2 * cell.degree + 1:
        raise WeilAuditVerificationError("candidate rationalization metadata is invalid")

    midpoint = source.get("eigenvalue_real_midpoint")
    if not isinstance(midpoint, Mapping) or set(midpoint) != {
        "numerator",
        "denominator",
    }:
        raise WeilAuditVerificationError("candidate eigenvalue midpoint is malformed")
    numerator = _parse_canonical_integer(
        midpoint["numerator"], "candidate eigenvalue numerator"
    )
    denominator = _parse_canonical_integer(
        midpoint["denominator"], "candidate eigenvalue denominator", 1
    )
    if dict(midpoint) != _fraction_record(Fraction(numerator, denominator)):
        raise WeilAuditVerificationError(
            "candidate eigenvalue midpoint is not reduced"
        )

    reverse = tuple(reversed(witness))
    if parity == "EVEN" and tuple(witness) != reverse:
        raise WeilAuditVerificationError("even candidate is not reversal-even")
    if parity == "ODD" and tuple(witness) != tuple(-value for value in reverse):
        raise WeilAuditVerificationError("odd candidate is not reversal-odd")
    return dict(source)


def _verify_parity_candidate_records(
    generation: Any,
    matrix: Sequence[Sequence[Any]],
    blocks: Mapping[str, Any],
    cell: WeilSearchCell,
    policy: WeilSearchPolicy,
    precision_bits: int,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not isinstance(generation, Mapping):
        raise WeilAuditVerificationError("candidate generation must be an object")
    expected_keys = {
        "classification",
        "role",
        "policy",
        "failure",
        "candidate_count",
        "negative_candidate_ids",
        "candidates",
    }
    if set(generation) != expected_keys:
        raise WeilAuditVerificationError("candidate-generation fields are malformed")
    if generation.get("role") != "APPROXIMATE_HINTS_ONLY":
        raise WeilAuditVerificationError("candidate-generation role was promoted")
    if generation.get("policy") != policy.to_record():
        raise WeilAuditVerificationError("candidate-generation policy changed")
    candidates = generation.get("candidates")
    if not isinstance(candidates, list):
        raise WeilAuditVerificationError("candidate list is malformed")
    count = _parse_canonical_integer(
        generation.get("candidate_count"), "candidate count", 0
    )
    if count != len(candidates):
        raise WeilAuditVerificationError("candidate count does not match its list")
    classification = generation.get("classification")
    failure = generation.get("failure")
    if classification == "AVAILABLE":
        if failure is not None:
            raise WeilAuditVerificationError(
                "available candidate generation contains a failure"
            )
    elif classification == "INCONCLUSIVE":
        if candidates or not isinstance(failure, Mapping):
            raise WeilAuditVerificationError(
                "inconclusive candidate generation is malformed"
            )
        if set(failure) != {"classification", "failure_code", "error_type"} or (
            failure.get("classification") != "INCONCLUSIVE"
            or failure.get("failure_code")
            != "APPROXIMATE_EIGENVECTOR_HINT_FAILED"
            or not isinstance(failure.get("error_type"), str)
            or not failure.get("error_type")
        ):
            raise WeilAuditVerificationError(
                "candidate-generation failure record is invalid"
            )
    else:
        raise WeilAuditVerificationError(
            "candidate-generation classification is invalid"
        )

    verified: list[dict[str, Any]] = []
    negative_ids: list[str] = []
    seen_ids: set[str] = set()
    seen_witnesses: set[tuple[int, ...]] = set()
    for candidate in candidates:
        if not isinstance(candidate, Mapping) or set(candidate) != {
            "candidate_id",
            "source",
            "witness",
            "sign_authority",
            "authoritative_full_matrix_evaluation",
            "full_vs_parity_quadratic_audit",
        }:
            raise WeilAuditVerificationError("candidate record is malformed")
        witness = _parse_witness_record(
            candidate["witness"], 2 * cell.degree + 1, "candidate witness"
        )
        if witness != _primitive_witness(witness):
            raise WeilAuditVerificationError(
                "candidate witness is not canonically primitive"
            )
        if witness in seen_witnesses:
            raise WeilAuditVerificationError("candidate witness is duplicated")
        seen_witnesses.add(witness)
        source = _validate_candidate_source(
            candidate["source"], witness, cell, policy
        )
        identity_record = {
            "source": source,
            "witness": [str(value) for value in witness],
        }
        candidate_id = content_sha256(identity_record)
        if candidate.get("candidate_id") != candidate_id or candidate_id in seen_ids:
            raise WeilAuditVerificationError("candidate id is invalid or duplicated")
        seen_ids.add(candidate_id)
        identity = _audit_full_vs_parity_quadratic(
            matrix, blocks, witness, precision_bits
        )
        expected = {
            "candidate_id": candidate_id,
            **identity_record,
            "sign_authority": "direct_full_matrix_interval_evaluation",
            "authoritative_full_matrix_evaluation": identity[
                "authoritative_full_matrix_evaluation"
            ],
            "full_vs_parity_quadratic_audit": identity,
        }
        if dict(candidate) != expected:
            raise WeilAuditVerificationError(
                "candidate direct evaluation or parity identity changed"
            )
        if identity["authoritative_full_matrix_evaluation"]["sign"] == "NEGATIVE":
            negative_ids.append(candidate_id)
        verified.append(expected)

    recorded_negative_ids = generation.get("negative_candidate_ids")
    if not isinstance(recorded_negative_ids, list) or recorded_negative_ids != negative_ids:
        raise WeilAuditVerificationError("negative candidate index is not derived")
    return verified, negative_ids


def _verify_requested_witness_records(
    records: Any,
    matrix: Sequence[Sequence[Any]],
    blocks: Mapping[str, Any],
    precision_bits: int,
) -> tuple[list[dict[str, Any]], list[str]]:
    if not isinstance(records, list):
        raise WeilAuditVerificationError("requested witness audits must be a list")
    verified: list[dict[str, Any]] = []
    negative_indices: list[str] = []
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise WeilAuditVerificationError("requested witness record is malformed")
        supplied = _parse_witness_record(
            record.get("supplied_witness"),
            len(matrix),
            "supplied requested witness",
        )
        primitive = _primitive_witness(supplied)
        expected = {
            **_audit_full_vs_parity_quadratic(
                matrix, blocks, primitive, precision_bits
            ),
            "supplied_witness": [str(value) for value in supplied],
            "normalization": "primitive-gcd-and-first-nonzero-positive",
            "supplied_witness_was_canonical": supplied == primitive,
        }
        if dict(record) != expected:
            raise WeilAuditVerificationError(
                "requested witness evaluation or normalization changed"
            )
        if expected["authoritative_full_matrix_evaluation"]["sign"] == "NEGATIVE":
            negative_indices.append(str(index))
        verified.append(expected)
    return verified, negative_indices


def _verify_parity_audit(artifact: Mapping[str, Any]) -> dict[str, Any]:
    record = dict(artifact)
    _require_audit_payload_hash(record)
    if record.get("schema") != PARITY_AUDIT_SCHEMA:
        raise WeilAuditVerificationError("unexpected parity audit schema")
    if record.get("classification") != "EXPLORATORY":
        raise WeilAuditVerificationError("parity audit was improperly promoted")
    if record.get("hypothesis_status") != "UNRESOLVED":
        raise WeilAuditVerificationError("parity audit changed global RH status")
    if record.get("algorithm_id") != PARITY_ALGORITHM_ID:
        raise WeilAuditVerificationError("parity audit algorithm changed")
    if record.get("backend") != _backend_record():
        raise WeilAuditVerificationError("parity audit backend does not match")
    if record.get("limitation") != EXPLORATORY_LIMITATION:
        raise WeilAuditVerificationError("parity audit limitation changed")
    cell = _parse_cell_record(record.get("cell"), "parity audit cell")
    precision_bits = _parse_canonical_integer(
        record.get("precision_bits"), "parity audit precision", 96
    )
    policy = _parse_policy_record(record.get("policy"))
    construction = record.get("matrix_construction")
    if not isinstance(construction, Mapping):
        raise WeilAuditVerificationError("matrix construction record is malformed")
    if construction.get("classification") == "INCONCLUSIVE":
        expected = _certify_parity_audit(cell, precision_bits, policy)
        if record != expected:
            raise WeilAuditVerificationError(
                "inconclusive parity audit does not canonically regenerate"
            )
        return {
            "classification": "REPRODUCED_EXPLORATORY_AUDIT",
            "audit_kind": "PARITY",
            "decision": "AUDIT_INCONCLUSIVE",
            "verified_integer_witnesses": "0",
            "same_backend_replay_only": True,
            "hypothesis_status": "UNRESOLVED",
        }
    if construction != {
        "classification": "COMPLETE",
        "normalization": "A=P-R-S",
    }:
        raise WeilAuditVerificationError("matrix construction record is invalid")

    expected_top_keys = {
        "schema",
        "classification",
        "hypothesis_status",
        "algorithm_id",
        "cell",
        "precision_bits",
        "backend",
        "policy",
        "limitation",
        "matrix_construction",
        "basis",
        "gates",
        "blocks",
        "rump_block_spectra",
        "candidate_generation",
        "requested_witness_audits",
        "requested_negative_witness_indices",
        "decision",
        "payload_sha256",
    }
    if set(record) != expected_top_keys:
        raise WeilAuditVerificationError("parity audit fields are malformed")

    with _working_precision(precision_bits):
        try:
            components = _symmetric_components(
                cell.cutoff_numerator,
                cell.cutoff_denominator,
                cell.degree,
            )
        except ArithmeticError as exc:
            raise WeilAuditVerificationError(
                "stored complete parity matrix no longer regenerates"
            ) from exc
        blocks = build_parity_blocks(components["q"])
        even_labels = ["0", *[str(value) for value in range(1, cell.degree + 1)]]
        odd_labels = [str(value) for value in range(1, cell.degree + 1)]
        expected_basis = {
            "full_order": [
                str(value) for value in range(-cell.degree, cell.degree + 1)
            ],
            "even_order": even_labels,
            "odd_order": odd_labels,
            "even_vectors": "e_0; then (e_-k+e_k)/sqrt(2), k=1,...,N",
            "odd_vectors": "(e_-k-e_k)/sqrt(2), k=1,...,N",
            "even_dimension": str(cell.degree + 1),
            "odd_dimension": str(cell.degree),
        }
        expected_gates = {
            "centrosymmetry_overlap": _centrosymmetry_gate(
                components, precision_bits
            ),
            "orthonormal_basis": _orthonormality_gate(blocks),
            "cross_parity_contains_zero": _cross_parity_gate(
                blocks["cross"], precision_bits
            ),
        }
        expected_blocks = {
            "even": {
                "dimension": str(cell.degree + 1),
                "entries_upper_triangle": _serialize_symmetric_block(
                    blocks["even"], even_labels, precision_bits
                ),
            },
            "odd": {
                "dimension": str(cell.degree),
                "entries_upper_triangle": _serialize_symmetric_block(
                    blocks["odd"], odd_labels, precision_bits
                ),
            },
        }
        expected_spectra = {
            "role": "DIAGNOSTIC_ONLY_NOT_A_SIGN_CERTIFICATE",
            "even": _rump_diagnostic(blocks["even"], precision_bits),
            "odd": _rump_diagnostic(blocks["odd"], precision_bits),
        }
        for name, expected in (
            ("basis", expected_basis),
            ("gates", expected_gates),
            ("blocks", expected_blocks),
            ("rump_block_spectra", expected_spectra),
        ):
            if record.get(name) != expected:
                raise WeilAuditVerificationError(
                    f"parity audit {name} changed on canonical replay"
                )

        candidates, negative_ids = _verify_parity_candidate_records(
            record["candidate_generation"],
            components["q"],
            blocks,
            cell,
            policy,
            precision_bits,
        )
        requested, requested_negative = _verify_requested_witness_records(
            record["requested_witness_audits"],
            components["q"],
            blocks,
            precision_bits,
        )

    if record.get("requested_negative_witness_indices") != requested_negative:
        raise WeilAuditVerificationError(
            "requested negative witness index is not derived"
        )
    structural_pass = all(
        gate["classification"] == "PASS" for gate in expected_gates.values()
    ) and all(
        item["full_vs_parity_quadratic_audit"]["classification"] == "PASS"
        for item in candidates
    ) and all(item["classification"] == "PASS" for item in requested)
    expected_decision = (
        "AUDIT_FAILED"
        if not structural_pass
        else "NEGATIVE_CANDIDATE_QUARANTINED"
        if negative_ids or requested_negative
        else "AUDIT_PASSED_EXPLORATORY"
    )
    if record.get("decision") != expected_decision:
        raise WeilAuditVerificationError("parity audit decision is not derived")
    return {
        "classification": "REPRODUCED_EXPLORATORY_AUDIT",
        "audit_kind": "PARITY",
        "decision": expected_decision,
        "verified_integer_witnesses": str(len(candidates) + len(requested)),
        "verified_negative_witnesses": str(
            len(negative_ids) + len(requested_negative)
        ),
        "same_backend_replay_only": True,
        "hypothesis_status": "UNRESOLVED",
    }


def verify_parity_audit(artifact: Mapping[str, Any]) -> dict[str, Any]:
    """Replay every certificate-bearing parity claim and exact witness sign."""

    if not isinstance(artifact, Mapping):
        raise WeilAuditVerificationError("parity audit must be an object")
    try:
        return _isolated_flint_operation(lambda: _verify_parity_audit(artifact))
    except WeilAuditVerificationError:
        raise
    except (
        ArithmeticError,
        AttributeError,
        IndexError,
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise WeilAuditVerificationError(
            "parity audit evidence is malformed"
        ) from exc


def _verify_degree_nesting_audit(artifact: Mapping[str, Any]) -> dict[str, Any]:
    record = dict(artifact)
    _require_audit_payload_hash(record)
    if record.get("schema") != NESTING_AUDIT_SCHEMA:
        raise WeilAuditVerificationError("unexpected degree-nesting audit schema")
    if record.get("classification") != "EXPLORATORY":
        raise WeilAuditVerificationError(
            "degree-nesting audit was improperly promoted"
        )
    if record.get("hypothesis_status") != "UNRESOLVED":
        raise WeilAuditVerificationError(
            "degree-nesting audit changed global RH status"
        )
    if record.get("algorithm_id") != NESTING_ALGORITHM_ID:
        raise WeilAuditVerificationError("degree-nesting algorithm changed")
    if record.get("backend") != _backend_record():
        raise WeilAuditVerificationError("degree-nesting backend does not match")
    if record.get("limitation") != EXPLORATORY_LIMITATION:
        raise WeilAuditVerificationError("degree-nesting limitation changed")
    lower_cell = _parse_cell_record(record.get("lower_cell"), "lower cell")
    higher_cell = _parse_cell_record(record.get("higher_cell"), "higher cell")
    precision_bits = _parse_canonical_integer(
        record.get("precision_bits"), "nesting precision", 96
    )
    replay_precision_bits = _parse_canonical_integer(
        record.get("replay_precision_bits"), "nesting replay precision", 96
    )

    witnesses: list[tuple[int, ...]] = []
    stored_checks = record.get("zero_padding_negative_witness_checks")
    if stored_checks is not None:
        if not isinstance(stored_checks, list):
            raise WeilAuditVerificationError(
                "zero-padding witness checks must be a list"
            )
        for check in stored_checks:
            if not isinstance(check, Mapping):
                raise WeilAuditVerificationError(
                    "zero-padding witness check is malformed"
                )
            base = check.get("base_precision")
            replay = check.get("replay_precision")
            if not isinstance(base, Mapping) or not isinstance(replay, Mapping):
                raise WeilAuditVerificationError(
                    "zero-padding precision records are malformed"
                )
            base_witness = _parse_witness_record(
                base.get("supplied_witness"),
                2 * lower_cell.degree + 1,
                "base supplied witness",
            )
            replay_witness = _parse_witness_record(
                replay.get("supplied_witness"),
                2 * lower_cell.degree + 1,
                "replay supplied witness",
            )
            if base_witness != replay_witness:
                raise WeilAuditVerificationError(
                    "base and replay supplied witnesses differ"
                )
            witnesses.append(base_witness)

    expected = certify_degree_nesting_audit(
        lower_cell,
        higher_cell,
        precision_bits,
        replay_precision_bits,
        witnesses,
    )
    if record != expected:
        raise WeilAuditVerificationError(
            "degree-nesting audit does not canonically regenerate"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_AUDIT",
        "audit_kind": "DEGREE_NESTING",
        "decision": expected["decision"],
        "verified_zero_padding_witnesses": str(len(witnesses)),
        "same_backend_replay_only": True,
        "hypothesis_status": "UNRESOLVED",
    }


def verify_degree_nesting_audit(
    artifact: Mapping[str, Any]
) -> dict[str, Any]:
    """Canonically regenerate a degree-nesting audit and its witness checks."""

    if not isinstance(artifact, Mapping):
        raise WeilAuditVerificationError("degree-nesting audit must be an object")
    try:
        return _isolated_flint_operation(
            lambda: _verify_degree_nesting_audit(artifact)
        )
    except WeilAuditVerificationError:
        raise
    except (
        ArithmeticError,
        AttributeError,
        IndexError,
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise WeilAuditVerificationError(
            "degree-nesting audit evidence is malformed"
        ) from exc


# Short, discoverable aliases for callers that do not use certificate naming.
audit_parity = certify_parity_audit
audit_degree_nesting = certify_degree_nesting_audit


__all__ = [
    "EXPLORATORY_LIMITATION",
    "NESTING_ALGORITHM_ID",
    "NESTING_AUDIT_SCHEMA",
    "PARITY_ALGORITHM_ID",
    "PARITY_AUDIT_SCHEMA",
    "WeilAuditVerificationError",
    "audit_degree_nesting",
    "audit_full_vs_parity_quadratic",
    "audit_parity",
    "audit_zero_padding_negative_witness",
    "build_parity_blocks",
    "certify_degree_nesting_audit",
    "certify_parity_audit",
    "extract_parity_integer_witnesses",
    "verify_degree_nesting_audit",
    "verify_parity_audit",
]
