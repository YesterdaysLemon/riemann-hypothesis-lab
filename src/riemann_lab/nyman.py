"""Rigorous finite natural-dilate calculations for the Nyman criterion.

The functions in this module use the modified ``L^2(0, infinity)``
normalization

``rho_a(x) = {1 / (a x)}``,  ``a = 1, 2, ...``.

For this normalization, with ``A(lambda)`` denoting the autocorrelation
kernel, the finite Gram system is

``G[a,b] = A(a / b) / a`` and
``target[a] = (log(a) + 1 - EulerGamma) / a``.

Every coefficient consumed by a certificate path is an exact dyadic
rational.  A numerical linear solve may be used to *propose* such a vector,
but the upper-bound decision is made only by direct interval evaluation of
``1 - 2 b^T c + c^T G c``.  A lower bound is decided only by fixed-order
interval ``LDL^T`` positivity of the corresponding augmented matrix.

All certificates here are finite-dimensional exploratory evidence.  They do
not prove or disprove the Riemann Hypothesis.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd
from typing import Any, Iterator, Sequence

from flint import arb, arb_mat, ctx

from .artifacts import content_sha256
from .balls import arb_record, arb_to_dyadic


NYMAN_UPPER_SCHEMA = "rh-lab/nyman-natural-distance-upper/v1"
NYMAN_LOWER_SCHEMA = "rh-lab/nyman-natural-distance-lower/v1"
NYMAN_BRACKET_SCHEMA = "rh-lab/nyman-natural-distance-bracket/v1"

NYMAN_NORMALIZATION = (
    "Modified L2(0,infinity) natural dilates rho_a(x)={1/(a*x)}; "
    "G_ab=A(a/b)/a and b_a=(log(a)+1-EulerGamma)/a."
)
NYMAN_LIMITATION = (
    "This is a finite-dimensional natural-dilate calculation.  A certified "
    "upper or lower bound at finitely many dilates does not establish the "
    "limiting Nyman criterion and therefore does not prove or disprove the "
    "Riemann Hypothesis."
)

# The nearest 2^-48 dyadics to the 256-bit interval solution of G c = b for
# the natural dilates 1,...,8.  These integers are frozen data, not a claim
# that the numerical solve is certificate evidence.
N8_FROZEN_DILATES = tuple(range(1, 9))
N8_FROZEN_DENOMINATOR_EXPONENT = 48
N8_FROZEN_NUMERATORS = (
    -264081585807853,
    251586556107910,
    252628013889294,
    48135412765604,
    188906090905926,
    -115437350602017,
    161913508068497,
    22696425045944,
)
N8_FROZEN_LOWER_BOUND = Fraction(3, 128)
N8_FROZEN_UPPER_BOUND = Fraction(25, 1024)


_CANONICAL_SYSTEM_PROVENANCE = object()


@dataclass(frozen=True, init=False)
class NaturalSystem:
    """A finite natural-dilate Gram system at the active Arb precision.

    Direct construction intentionally creates an untrusted value suitable for
    low-level algebra and tests.  Only :func:`build_natural_system` and
    :func:`prefix_natural_system` can create the commitment-bearing canonical
    systems accepted by public certificate functions.
    """

    dilates: tuple[int, ...]
    gram: tuple[tuple[Any, ...], ...]
    target: tuple[Any, ...]
    precision_bits: int | None = None
    _source_dilates: tuple[int, ...] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _canonical_commitment: str | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _canonical_provenance: object | None = field(
        default=None, init=False, repr=False, compare=False
    )

    def __init__(
        self,
        dilates: tuple[int, ...],
        gram: tuple[tuple[Any, ...], ...],
        target: tuple[Any, ...],
        precision_bits: int | None = None,
    ) -> None:
        object.__setattr__(self, "dilates", dilates)
        object.__setattr__(self, "gram", gram)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "precision_bits", precision_bits)
        object.__setattr__(self, "_source_dilates", None)
        object.__setattr__(self, "_canonical_commitment", None)
        object.__setattr__(self, "_canonical_provenance", None)


def _validate_integer(value: int, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _validate_precision(precision_bits: int) -> int:
    _validate_integer(precision_bits, "precision_bits", 96)
    return precision_bits


def _validate_fraction(
    value: Fraction,
    name: str,
    *,
    minimum: Fraction | None = None,
    maximum: Fraction | None = None,
) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError(f"{name} must be an exact Fraction")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    if maximum is not None and value >= maximum:
        raise ValueError(f"{name} must be below {maximum}")
    return value


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _arb_from_fraction(value: Fraction) -> Any:
    return arb(value.numerator) / value.denominator


def _is_power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


@contextmanager
def _working_precision(precision_bits: int) -> Iterator[None]:
    previous_precision = ctx.prec
    ctx.cleanup()
    ctx.prec = _validate_precision(precision_bits)
    try:
        yield
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()


def _with_payload_hash(payload: dict[str, Any]) -> dict[str, Any]:
    return {**payload, "payload_sha256": content_sha256(payload)}


def vasyunin_sum(p: int, q: int) -> Any:
    """Return the rigorous Vasyunin sum ``V(p,q)`` as an Arb ball.

    ``p`` and ``q`` must be positive and coprime.  The fractional parts are
    formed by exact integer remainders; binary floating point is never used.
    By convention ``V(p,1) = 0`` exactly.
    """

    _validate_integer(p, "p", 1)
    _validate_integer(q, "q", 1)
    if gcd(p, q) != 1:
        raise ValueError("p and q must be coprime")
    if q == 1:
        return arb(0)

    pi = arb.pi()
    total = arb(0)
    for k in range(1, q):
        fractional_part = arb((k * p) % q) / q
        total += fractional_part * (pi * k / q).cot()
    return total


def _autocorrelation_formula_with_sum(ratio: Fraction, vasyunin_pair_sum: Any) -> Any:
    p = ratio.numerator
    q = ratio.denominator
    lam = _arb_from_fraction(ratio)
    half = arb(1) / 2
    return (
        (1 - lam) * half * lam.log()
        + (1 + lam)
        * half
        * ((2 * arb.pi()).log() - arb.const_euler())
        - arb.pi()
        * vasyunin_pair_sum
        / (2 * q)
    )


def _autocorrelation_formula(ratio: Fraction) -> Any:
    """Evaluate the verified rational formula in its supplied orientation."""

    p = ratio.numerator
    q = ratio.denominator
    pair_sum = vasyunin_sum(p, q) + vasyunin_sum(q, p)
    return _autocorrelation_formula_with_sum(ratio, pair_sum)


def _reciprocal_autocorrelations(ratio: Fraction) -> tuple[Any, Any]:
    """Evaluate both rational orientations from one exact Vasyunin pair."""

    p = ratio.numerator
    q = ratio.denominator
    pair_sum = vasyunin_sum(p, q) + vasyunin_sum(q, p)
    return (
        _autocorrelation_formula_with_sum(ratio, pair_sum),
        _autocorrelation_formula_with_sum(1 / ratio, pair_sum),
    )


def autocorrelation_a(ratio: Fraction) -> Any:
    """Evaluate the natural-dilate autocorrelation ``A(ratio)`` rigorously.

    ``Fraction`` supplies an exact reduced ``p/q``.  The formula is evaluated
    in the supplied orientation; :func:`natural_gram_entry` independently
    checks the reciprocal orientation before storing a shared enclosure.
    """

    ratio = _validate_fraction(ratio, "ratio")
    if ratio <= 0:
        raise ValueError("ratio must be positive")
    return _autocorrelation_formula(ratio)


def autocorrelation_reciprocity_residual(ratio: Fraction) -> Any:
    """Return an interval enclosing ``A(lambda)-lambda*A(1/lambda)``."""

    ratio = _validate_fraction(ratio, "ratio")
    if ratio <= 0:
        raise ValueError("ratio must be positive")
    return autocorrelation_a(ratio) - _arb_from_fraction(ratio) * autocorrelation_a(
        1 / ratio
    )


def natural_gram_entry(a: int, b: int) -> Any:
    """Return ``<rho_a,rho_b> = A(a/b)/a`` in canonical symmetric form."""

    _validate_integer(a, "a", 1)
    _validate_integer(b, "b", 1)
    lower = min(a, b)
    upper = max(a, b)
    forward_a, reciprocal_a = _reciprocal_autocorrelations(
        Fraction(lower, upper)
    )
    forward = forward_a / lower
    reciprocal = reciprocal_a / upper
    if not forward.overlaps(reciprocal):
        raise ArithmeticError("autocorrelation reciprocity check failed")
    return forward.intersection(reciprocal)


def natural_target_entry(a: int) -> Any:
    """Return ``<chi_(0,1),rho_a>`` in the modified normalization."""

    _validate_integer(a, "a", 1)
    return (arb(a).log() + 1 - arb.const_euler()) / a


def _validate_dilates(dilates: Sequence[int]) -> tuple[int, ...]:
    if isinstance(dilates, (str, bytes)):
        raise TypeError("dilates must be a sequence of integers")
    result = tuple(dilates)
    if not result:
        raise ValueError("dilates must be nonempty")
    for index, value in enumerate(result):
        _validate_integer(value, f"dilates[{index}]", 1)
    if len(set(result)) != len(result):
        raise ValueError("dilates must be distinct")
    return result


def _system_commitment_payload(system: NaturalSystem) -> dict[str, Any]:
    source_dilates = system._source_dilates
    return {
        "schema": "rh-lab/nyman-natural-system/v1",
        "normalization": NYMAN_NORMALIZATION,
        "dilates": [str(value) for value in system.dilates],
        "precision_bits": str(system.precision_bits),
        "source_dilates": (
            [str(value) for value in source_dilates]
            if source_dilates is not None
            else None
        ),
        "gram": [
            [arb_record(value) for value in row] for row in system.gram
        ],
        "target": [arb_record(value) for value in system.target],
    }


def _make_canonical_system(
    dilates: tuple[int, ...],
    gram: tuple[tuple[Any, ...], ...],
    target: tuple[Any, ...],
    precision_bits: int,
    source_dilates: tuple[int, ...],
) -> NaturalSystem:
    system = object.__new__(NaturalSystem)
    object.__setattr__(system, "dilates", dilates)
    object.__setattr__(system, "gram", gram)
    object.__setattr__(system, "target", target)
    object.__setattr__(system, "precision_bits", precision_bits)
    object.__setattr__(system, "_source_dilates", source_dilates)
    object.__setattr__(system, "_canonical_commitment", None)
    object.__setattr__(
        system, "_canonical_provenance", _CANONICAL_SYSTEM_PROVENANCE
    )
    object.__setattr__(
        system,
        "_canonical_commitment",
        content_sha256(_system_commitment_payload(system)),
    )
    return system


def _build_natural_system(dilates: Sequence[int]) -> NaturalSystem:
    exact_dilates = _validate_dilates(dilates)
    size = len(exact_dilates)
    matrix = [[arb(0) for _ in range(size)] for _ in range(size)]
    autocorrelation_cache: dict[Fraction, tuple[Any, Any]] = {}

    for row, a in enumerate(exact_dilates):
        for column in range(row, size):
            b = exact_dilates[column]
            lower = min(a, b)
            upper = max(a, b)
            ratio = Fraction(lower, upper)
            if ratio not in autocorrelation_cache:
                autocorrelation_cache[ratio] = _reciprocal_autocorrelations(ratio)
            forward_a, reciprocal_a = autocorrelation_cache[ratio]
            forward = forward_a / lower
            reciprocal = reciprocal_a / upper
            if not forward.overlaps(reciprocal):
                raise ArithmeticError("autocorrelation reciprocity check failed")
            entry = forward.intersection(reciprocal)
            matrix[row][column] = entry
            matrix[column][row] = entry

    gram = tuple(tuple(row) for row in matrix)
    # Fail closed if a future refactor stops mirroring identical enclosures.
    for row in range(size):
        for column in range(row):
            if arb_to_dyadic(gram[row][column]) != arb_to_dyadic(
                gram[column][row]
            ):
                raise ArithmeticError("Gram construction lost exact symmetry")

    return _make_canonical_system(
        exact_dilates,
        gram,
        tuple(natural_target_entry(a) for a in exact_dilates),
        int(ctx.prec),
        exact_dilates,
    )


def build_natural_system(dilates: Sequence[int]) -> NaturalSystem:
    """Build a canonical finite Gram system from clean FLINT cache state."""

    previous_precision = ctx.prec
    ctx.cleanup()
    ctx.prec = previous_precision
    try:
        return _build_natural_system(dilates)
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()


def _validate_system(
    system: NaturalSystem, *, require_canonical: bool = False
) -> int:
    if not isinstance(system, NaturalSystem):
        raise TypeError("system must be a NaturalSystem")
    size = len(system.dilates)
    _validate_dilates(system.dilates)
    if (
        len(system.gram) != size
        or any(len(row) != size for row in system.gram)
        or len(system.target) != size
    ):
        raise ValueError("system dimensions are inconsistent")
    if system.precision_bits is not None:
        _validate_integer(system.precision_bits, "system precision_bits", 1)
    for row in range(size):
        target = system.target[row]
        if not hasattr(target, "is_finite") or not target.is_finite():
            raise TypeError("system target entries must be finite Arb enclosures")
        for column in range(size):
            entry = system.gram[row][column]
            if not hasattr(entry, "is_finite") or not entry.is_finite():
                raise TypeError("system Gram entries must be finite Arb enclosures")
        for column in range(row):
            if arb_to_dyadic(system.gram[row][column]) != arb_to_dyadic(
                system.gram[column][row]
            ):
                raise ValueError("system Gram matrix must be exactly symmetric")
    if require_canonical:
        if system._canonical_provenance is not _CANONICAL_SYSTEM_PROVENANCE:
            raise ValueError("supplied system lacks canonical builder provenance")
        source = system._source_dilates
        if source is None or source[:size] != system.dilates:
            raise ValueError("supplied system has invalid prefix provenance")
        expected = content_sha256(_system_commitment_payload(system))
        if system._canonical_commitment != expected:
            raise ValueError("supplied system content commitment changed")
    return size


def natural_system_content_sha256(system: NaturalSystem) -> str:
    """Return the canonical commitment binding all finite kernel evidence."""

    _validate_system(system, require_canonical=True)
    assert system._canonical_commitment is not None
    return system._canonical_commitment


def prefix_natural_system(system: NaturalSystem, size: int) -> NaturalSystem:
    """Return a provenance-preserving leading principal natural system."""

    full_size = _validate_system(system, require_canonical=True)
    _validate_integer(size, "size", 1)
    if size > full_size:
        raise ValueError("prefix size exceeds the supplied system dimension")
    assert system.precision_bits is not None
    assert system._source_dilates is not None
    return _make_canonical_system(
        tuple(system.dilates[:size]),
        tuple(tuple(row[:size]) for row in system.gram[:size]),
        tuple(system.target[:size]),
        system.precision_bits,
        system._source_dilates,
    )


def dyadic_coefficients(
    numerators: Sequence[int], denominator_exponent: int
) -> tuple[Fraction, ...]:
    """Decode signed integers over the common denominator ``2**exponent``."""

    _validate_integer(denominator_exponent, "denominator_exponent", 0)
    if isinstance(numerators, (str, bytes)):
        raise TypeError("numerators must be a sequence of integers")
    exact_numerators = tuple(numerators)
    if not exact_numerators:
        raise ValueError("numerators must be nonempty")
    for index, value in enumerate(exact_numerators):
        _validate_integer(value, f"numerators[{index}]")
    denominator = 1 << denominator_exponent
    return tuple(Fraction(value, denominator) for value in exact_numerators)


def _validate_dyadic_coefficients(
    coefficients: Sequence[Fraction], expected_size: int
) -> tuple[Fraction, ...]:
    if isinstance(coefficients, (str, bytes)):
        raise TypeError("coefficients must be a sequence of exact Fractions")
    result = tuple(coefficients)
    if len(result) != expected_size:
        raise ValueError("coefficient length must equal the system dimension")
    for index, value in enumerate(result):
        _validate_fraction(value, f"coefficients[{index}]")
        if not _is_power_of_two(value.denominator):
            raise ValueError("certificate coefficients must be dyadic Fractions")
    return result


def evaluate_natural_distance(
    system: NaturalSystem, coefficients: Sequence[Fraction]
) -> Any:
    """Directly enclose ``E(c)=1-2*b^T*c+c^T*G*c`` for exact dyadics."""

    size = _validate_system(system)
    exact_coefficients = _validate_dyadic_coefficients(coefficients, size)
    coefficient_balls = tuple(
        _arb_from_fraction(value) for value in exact_coefficients
    )

    energy = arb(1)
    for index in range(size):
        energy -= 2 * system.target[index] * coefficient_balls[index]
        energy += (
            system.gram[index][index]
            * coefficient_balls[index]
            * coefficient_balls[index]
        )
        for column in range(index + 1, size):
            energy += (
                2
                * system.gram[index][column]
                * coefficient_balls[index]
                * coefficient_balls[column]
            )
    return energy


def _fraction_from_exact_dyadic(mantissa: int, exponent: int) -> Fraction:
    if exponent >= 0:
        return Fraction(mantissa << exponent, 1)
    return Fraction(mantissa, 1 << -exponent)


def _round_fraction_ties_even(value: Fraction) -> int:
    floor_value = value.numerator // value.denominator
    remainder = value.numerator - floor_value * value.denominator
    doubled = 2 * remainder
    if doubled < value.denominator:
        return floor_value
    if doubled > value.denominator:
        return floor_value + 1
    return floor_value if floor_value % 2 == 0 else floor_value + 1


def round_arb_to_dyadic_numerator(value: Any, denominator_exponent: int) -> int:
    """Round an Arb enclosure to one provably determined ties-to-even integer.

    The returned integer represents the dyadic value ``integer / 2**exponent``.
    A nonzero-width enclosure that reaches a rounding boundary is rejected;
    an approximate solve is therefore never silently converted ambiguously.
    """

    _validate_integer(denominator_exponent, "denominator_exponent", 0)
    if not hasattr(value, "is_finite") or not value.is_finite():
        raise TypeError("value must be a finite Arb-compatible enclosure")
    scaled = value * (1 << denominator_exponent)
    midpoint_mantissa, midpoint_exponent = scaled.mid().man_exp()
    midpoint = _fraction_from_exact_dyadic(
        int(midpoint_mantissa), int(midpoint_exponent)
    )
    candidate = _round_fraction_ties_even(midpoint)
    if scaled.is_exact():
        return candidate

    lower_boundary = _arb_from_fraction(Fraction(2 * candidate - 1, 2))
    upper_boundary = _arb_from_fraction(Fraction(2 * candidate + 1, 2))
    if not (scaled > lower_boundary and scaled < upper_boundary):
        raise ArithmeticError("enclosure does not determine a unique dyadic rounding")
    return candidate


def round_arb_vector_to_dyadic(
    values: Sequence[Any], denominator_exponent: int
) -> tuple[int, ...]:
    """Round a proposed Arb vector to exact common-denominator dyadics."""

    if isinstance(values, (str, bytes)):
        raise TypeError("values must be a sequence of Arb enclosures")
    return tuple(
        round_arb_to_dyadic_numerator(value, denominator_exponent)
        for value in values
    )


def propose_dyadic_coefficients(
    system: NaturalSystem, denominator_exponent: int
) -> tuple[int, ...]:
    """Propose dyadic numerators from ``G c=b``; this is not a certificate.

    Only a subsequent direct call to :func:`certify_dyadic_upper_bound` can
    turn the frozen exact vector into finite upper-bound evidence.
    """

    _validate_integer(denominator_exponent, "denominator_exponent", 0)
    size = _validate_system(system)
    matrix = arb_mat(size, size)
    right_hand_side = arb_mat(size, 1)
    for row in range(size):
        right_hand_side[row, 0] = system.target[row]
        for column in range(size):
            matrix[row, column] = system.gram[row][column]
    solution = matrix.solve(right_hand_side)
    return round_arb_vector_to_dyadic(
        tuple(solution[index, 0] for index in range(size)),
        denominator_exponent,
    )


def build_augmented_lower_matrix(
    system: NaturalSystem, lower_bound: Fraction
) -> tuple[tuple[Any, ...], ...]:
    """Build ``[[G,-b],[-b^T,1-L]]`` for the lower-bound Schur argument."""

    lower_bound = _validate_fraction(
        lower_bound,
        "lower_bound",
        minimum=Fraction(0),
        maximum=Fraction(1),
    )
    size = _validate_system(system)
    matrix = [[arb(0) for _ in range(size + 1)] for _ in range(size + 1)]
    for row in range(size):
        for column in range(size):
            matrix[row][column] = system.gram[row][column]
        off_diagonal = -system.target[row]
        matrix[row][size] = off_diagonal
        matrix[size][row] = off_diagonal
    matrix[size][size] = 1 - _arb_from_fraction(lower_bound)
    return tuple(tuple(row) for row in matrix)


def fixed_order_interval_ldlt(matrix: Sequence[Sequence[Any]]) -> dict[str, Any]:
    """Prove positive definiteness using unpivoted interval ``LDL^T``.

    The order is exactly the supplied order.  No approximate eigenvalue or
    pivoting heuristic participates in the decision.
    """

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    for row in range(size):
        for column in range(size):
            value = matrix[row][column]
            if not hasattr(value, "is_finite") or not value.is_finite():
                raise TypeError("matrix entries must be finite Arb enclosures")
        for column in range(row):
            if arb_to_dyadic(matrix[row][column]) != arb_to_dyadic(
                matrix[column][row]
            ):
                raise ValueError("matrix must use identical symmetric enclosures")

    lower = [[arb(0) for _ in range(size)] for _ in range(size)]
    pivots: list[Any] = []
    for row in range(size):
        lower[row][row] = arb(1)
        for column in range(row):
            value = matrix[row][column]
            for prior in range(column):
                value -= lower[row][prior] * lower[column][prior] * pivots[prior]
            value /= pivots[column]
            lower[row][column] = value

        pivot = matrix[row][row]
        for prior in range(row):
            pivot -= lower[row][prior] ** 2 * pivots[prior]
        pivots.append(pivot)
        if not (pivot > 0):
            return {
                "classification": "INCONCLUSIVE",
                "fixed_order": [str(index) for index in range(size)],
                "failed_pivot_index": str(row),
                "pivots": [
                    {"index": str(index), "value": arb_record(value)}
                    for index, value in enumerate(pivots)
                ],
            }

    return {
        "classification": "POSITIVE_DEFINITE",
        "fixed_order": [str(index) for index in range(size)],
        "failed_pivot_index": None,
        "pivots": [
            {"index": str(index), "value": arb_record(value)}
            for index, value in enumerate(pivots)
        ],
    }


def _certificate_preamble(schema: str, precision_bits: int) -> dict[str, Any]:
    return {
        "schema": schema,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(precision_bits),
        "normalization": NYMAN_NORMALIZATION,
        "limitation": NYMAN_LIMITATION,
    }


def _certificate_system(
    dilates: tuple[int, ...],
    precision_bits: int,
    supplied: NaturalSystem | None,
) -> NaturalSystem:
    if supplied is None:
        generated = build_natural_system(dilates)
        _validate_system(generated, require_canonical=True)
        return generated
    _validate_system(supplied, require_canonical=True)
    if supplied.dilates != dilates:
        raise ValueError("supplied system dilates do not match the certificate scope")
    if supplied.precision_bits != precision_bits:
        raise ValueError("supplied system precision does not match precision_bits")
    return supplied


def _kernel_provenance_record(
    system: NaturalSystem, *, supplied_system_reused: bool
) -> dict[str, Any]:
    _validate_system(system, require_canonical=True)
    assert system._source_dilates is not None
    return {
        "construction": "canonical-rational-vasyunin-with-reciprocity-intersection",
        "normalization": NYMAN_NORMALIZATION,
        "system_content_sha256": natural_system_content_sha256(system),
        "precision_bits": str(system.precision_bits),
        "source_dilates": [str(value) for value in system._source_dilates],
        "source_dimension": str(len(system._source_dilates)),
        "scope_is_leading_principal_prefix": True,
        "supplied_system_reused": supplied_system_reused,
    }


def certify_dyadic_upper_bound(
    dilates: Sequence[int],
    numerators: Sequence[int],
    denominator_exponent: int,
    claimed_upper_bound: Fraction,
    *,
    precision_bits: int = 192,
    system: NaturalSystem | None = None,
) -> dict[str, Any]:
    """Certify ``E(c) <= U`` for one frozen exact common-dyadic vector."""

    exact_dilates = _validate_dilates(dilates)
    if isinstance(numerators, (str, bytes)):
        raise TypeError("numerators must be a sequence of integers")
    exact_numerators = tuple(numerators)
    exact_coefficients = dyadic_coefficients(
        exact_numerators, denominator_exponent
    )
    if len(exact_coefficients) != len(exact_dilates):
        raise ValueError("numerator length must equal the number of dilates")
    claimed_upper_bound = _validate_fraction(
        claimed_upper_bound,
        "claimed_upper_bound",
        minimum=Fraction(0),
        maximum=Fraction(1),
    )
    if claimed_upper_bound == 0:
        raise ValueError("claimed_upper_bound must be positive")
    precision_bits = _validate_precision(precision_bits)

    with _working_precision(precision_bits):
        certificate_system = _certificate_system(
            exact_dilates, precision_bits, system
        )
        energy = evaluate_natural_distance(
            certificate_system, exact_coefficients
        )
        margin = _arb_from_fraction(claimed_upper_bound) - energy
        certified = bool(margin >= 0)
        strictly_below = bool(margin > 0)
        payload = {
            **_certificate_preamble(NYMAN_UPPER_SCHEMA, precision_bits),
            "decision": (
                "UPPER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
            ),
            "scope": {
                "dilates": [str(value) for value in exact_dilates],
                "dimension": str(len(exact_dilates)),
                "finite_only": True,
            },
            "kernel_provenance": _kernel_provenance_record(
                certificate_system,
                supplied_system_reused=system is not None,
            ),
            "coefficients": {
                "encoding": "signed-integers-over-common-power-of-two",
                "numerators": [str(value) for value in exact_numerators],
                "denominator_exponent": str(denominator_exponent),
            },
            "claimed_upper_bound": _fraction_record(claimed_upper_bound),
            "energy": arb_record(energy),
            "margin": arb_record(margin),
            "checks": {
                "coefficients_are_exact_dyadics": True,
                "energy_at_most_claimed_upper_bound": certified,
                "energy_strictly_below_claimed_upper_bound": strictly_below,
                "direct_energy_evaluation_only": True,
                "approximate_solve_used_as_evidence": False,
            },
        }
        return _with_payload_hash(payload)


def certify_augmented_lower_bound(
    dilates: Sequence[int],
    lower_bound: Fraction,
    *,
    precision_bits: int = 192,
    system: NaturalSystem | None = None,
) -> dict[str, Any]:
    """Certify the finite minimum is strictly above ``L`` by interval LDL."""

    exact_dilates = _validate_dilates(dilates)
    lower_bound = _validate_fraction(
        lower_bound,
        "lower_bound",
        minimum=Fraction(0),
        maximum=Fraction(1),
    )
    precision_bits = _validate_precision(precision_bits)

    with _working_precision(precision_bits):
        certificate_system = _certificate_system(
            exact_dilates, precision_bits, system
        )
        augmented = build_augmented_lower_matrix(
            certificate_system, lower_bound
        )
        ldlt = fixed_order_interval_ldlt(augmented)
        certified = ldlt["classification"] == "POSITIVE_DEFINITE"
        payload = {
            **_certificate_preamble(NYMAN_LOWER_SCHEMA, precision_bits),
            "decision": (
                "LOWER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
            ),
            "scope": {
                "dilates": [str(value) for value in exact_dilates],
                "dimension": str(len(exact_dilates)),
                "finite_only": True,
            },
            "kernel_provenance": _kernel_provenance_record(
                certificate_system,
                supplied_system_reused=system is not None,
            ),
            "claimed_lower_bound": _fraction_record(lower_bound),
            "augmented_matrix": {
                "block_form": "[[G,-b],[-b^T,1-L]]",
                "quadratic_identity": "[c,1]^T M [c,1]=E(c)-L",
                "dimension": str(len(exact_dilates) + 1),
            },
            "augmented_ldlt": ldlt,
            "checks": {
                "fixed_order_interval_ldlt_positive": certified,
                "approximate_solve_used_as_evidence": False,
            },
        }
        return _with_payload_hash(payload)


def certify_frozen_n8_bracket(
    *, precision_bits: int = 192
) -> dict[str, Any]:
    """Reproduce the frozen finite N=8 bracket calibration."""

    precision_bits = _validate_precision(precision_bits)
    upper = certify_dyadic_upper_bound(
        N8_FROZEN_DILATES,
        N8_FROZEN_NUMERATORS,
        N8_FROZEN_DENOMINATOR_EXPONENT,
        N8_FROZEN_UPPER_BOUND,
        precision_bits=precision_bits,
    )
    lower = certify_augmented_lower_bound(
        N8_FROZEN_DILATES,
        N8_FROZEN_LOWER_BOUND,
        precision_bits=precision_bits,
    )
    certified = (
        upper["decision"] == "UPPER_BOUND_CERTIFIED"
        and lower["decision"] == "LOWER_BOUND_CERTIFIED"
    )
    payload = {
        **_certificate_preamble(NYMAN_BRACKET_SCHEMA, precision_bits),
        "decision": "FINITE_BRACKET_CERTIFIED" if certified else "INCONCLUSIVE",
        "scope": {
            "dilates": [str(value) for value in N8_FROZEN_DILATES],
            "dimension": "8",
            "finite_only": True,
        },
        "lower_bound": _fraction_record(N8_FROZEN_LOWER_BOUND),
        "upper_bound": _fraction_record(N8_FROZEN_UPPER_BOUND),
        "upper_certificate": upper,
        "lower_certificate": lower,
    }
    return _with_payload_hash(payload)
