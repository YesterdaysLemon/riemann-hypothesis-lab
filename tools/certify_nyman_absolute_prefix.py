"""Large-cutoff arbitrary-slope Nyman prefix certificates.

This module is a versioned successor to the legacy absolute-prefix helper
embedded in ``generate_nyman_rebased_schur_certificate.py``.  The legacy
implementation and its frozen p1--p4 artifacts intentionally remain
untouched.

For exact dyadic coefficients ``p_n = A_n / D``, put

    P = sum_n p_n / n,
    Q_M = D + sum_n A_n floor(M / n).

The exact energy through the integer interval ``M = T`` is

    (T + 1) P^2
    + sum_(M=1)^T Q_M^2 / (D^2 M (M + 1))
    - 2 P L_T,

where ``L_T`` is evaluated by exact Abel compression and Arb.  The ``Q_M``
recurrence is streamed in exact signed int64 chunks.  Each nonnegative
rational term is evaluated by scaling by the dyadic denominator and then
dividing sequentially by ``M`` and ``M + 1``.  This avoids forming the
binary64 product ``M(M+1)`` and supports cutoffs through ``2^53 - 1``.

The binary64 terms are reduced by an explicit, deterministic pairwise tree.
An exact rational forward-error radius covers integer conversion, squaring,
both divisions, every leaf reduction, and the final reduction of leaf roots.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
import operator
from typing import Any, Mapping

from flint import arb, ctx
import numpy as np

if __package__:
    from tools.certify_nyman_balanced_prefix import (
        BINARY64_UNIT_ROUNDOFF,
        INT64_MAX,
        common_dyadic_numerators,
        require_binary64_tcb,
    )
else:
    from certify_nyman_balanced_prefix import (
        BINARY64_UNIT_ROUNDOFF,
        INT64_MAX,
        common_dyadic_numerators,
        require_binary64_tcb,
    )


PREFIX_ALGORITHM_VERSION = "nyman-absolute-prefix-sequential-pairwise-v1"
TERM_EVALUATION_ALGORITHM = (
    "int64-to-binary64,exact-dyadic-scale,square,divide-M,divide-M-plus-1"
)
REDUCTION_ALGORITHM = "explicit-binary64-balanced-pairwise-leaves-v1"
BINARY64_TERM_ROUNDING_FACTOR_COUNT = 5
NORMAL_RANGE_LATTICE_EXPONENT = 900
MAX_EXACT_BINARY64_INTEGER = 1 << 53


@dataclass(frozen=True)
class PairwiseReduction:
    """The result and exact tree metadata for one binary64 reduction."""

    value: float
    input_count: int
    depth: int


@dataclass(frozen=True)
class AbsoluteEnergyPrefix:
    """A rigorous arbitrary-slope prefix enclosure and its proof metadata."""

    algorithm_version: str
    term_evaluation_algorithm: str
    reduction_algorithm: str
    cutoff: int
    chunk_size: int
    chunk_count: int
    reduction_leaf_size: int
    reduction_leaf_count: int
    maximum_leaf_length: int
    maximum_leaf_reduction_depth: int
    root_reduction_depth: int
    reduction_depth_bound: int
    coefficient_denominator: int
    coefficient_denominator_exponent: int
    recurrence_absolute_bound: int
    int64_headroom: int
    term_rounding_factor_count: int
    term_relative_error: Fraction
    reduction_gamma: Fraction
    rational_midpoint: Fraction
    rational_error: Fraction
    rational_lower: Fraction
    rational_upper: Fraction
    computed_absolute_term_sum_upper: Fraction
    term_evaluation_error: Fraction
    reduction_error: Fraction
    logarithmic_sum: Any
    logarithmic_cross_term: Any
    slope_square_term: Fraction
    prefix_energy: Any
    last_q_numerator: int
    maximum_q_numerator: int
    maximum_square_numerator: int


def _clean_coefficients(
    coefficients: Mapping[int, Fraction],
) -> dict[int, Fraction]:
    cleaned: dict[int, Fraction] = {}
    for raw_index, raw_value in coefficients.items():
        if isinstance(raw_index, (bool, np.bool_)):
            raise ValueError("coefficient indices must be exact integers")
        try:
            index = operator.index(raw_index)
        except TypeError as exc:
            raise ValueError(
                "coefficient indices must be exact integers"
            ) from exc
        if index < 1:
            raise ValueError("coefficient indices must be positive")
        if index in cleaned:
            raise ValueError("coefficient indices are not unique after conversion")
        value = Fraction(raw_value)
        if value:
            cleaned[index] = value
    return cleaned


def _power_of_two_exponent(value: int) -> int:
    if value < 1 or value & (value - 1):
        raise ValueError("coefficient denominator must be a power of two")
    return value.bit_length() - 1


def _exact_nonboolean_integer(value: object, label: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be an exact integer")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f"{label} must be an exact integer") from exc


def _arb_from_fraction(value: Fraction) -> Any:
    exact = Fraction(value)
    return arb(exact.numerator) / exact.denominator


def _radius_ball(value: Fraction, bits: int) -> Any:
    """Return an Arb zero ball with dyadic radius at least ``value``."""

    exact = Fraction(value)
    if exact < 0:
        raise ValueError("radius must be nonnegative")
    if not exact:
        return arb(0)
    scale = 1 << bits
    numerator = (
        exact.numerator * scale + exact.denominator - 1
    ) // exact.denominator
    return arb(0, (numerator, -bits))


def _gamma(operation_count: int) -> Fraction:
    if operation_count < 0:
        raise ValueError("operation count must be nonnegative")
    product = operation_count * BINARY64_UNIT_ROUNDOFF
    if product >= 1:
        raise ArithmeticError("binary64 pairwise-reduction bound exhausted")
    return product / (1 - product)


def _ceil_log2(value: int) -> int:
    if value < 1:
        raise ValueError("logarithm argument must be positive")
    return (value - 1).bit_length()


def explicit_pairwise_sum(values: np.ndarray) -> PairwiseReduction:
    """Reduce a nonempty binary64 vector in one declared balanced tree.

    Adjacent pairs are added by one elementwise ``np.add`` operation.  At an
    odd level the last value is copied forward without arithmetic.  Thus every
    input passes through at most ``ceil(log2(n))`` rounded additions.
    """

    work = np.asarray(values, dtype=np.float64)
    if work.ndim != 1:
        raise ValueError("pairwise reduction input must be one-dimensional")
    count = int(work.size)
    if count < 1:
        raise ValueError("pairwise reduction input must be nonempty")
    if not np.all(np.isfinite(work)):
        raise ArithmeticError("nonfinite pairwise-reduction input")

    depth = 0
    while work.size > 1:
        pair_count = int(work.size // 2)
        next_size = pair_count + int(work.size & 1)
        reduced = np.empty(next_size, dtype=np.float64)
        np.add(
            work[: 2 * pair_count : 2],
            work[1 : 2 * pair_count : 2],
            out=reduced[:pair_count],
        )
        if work.size & 1:
            reduced[-1] = work[-1]
        work = reduced
        depth += 1

    value = float(work[0])
    if not math.isfinite(value):
        raise ArithmeticError("nonfinite pairwise-reduction result")
    expected_depth = _ceil_log2(count)
    if depth != expected_depth:
        raise ArithmeticError("pairwise-reduction depth mismatch")
    return PairwiseReduction(
        value=value,
        input_count=count,
        depth=depth,
    )


def _harmonic_sum(
    coefficients: Mapping[int, Fraction],
) -> Fraction:
    return sum(
        (
            Fraction(coefficient) / int(index)
            for index, coefficient in coefficients.items()
        ),
        start=Fraction(),
    )


def _compressed_absolute_log_sum(
    coefficients: Mapping[int, Fraction],
    q_t: Fraction,
    cutoff: int,
) -> Any:
    result = _arb_from_fraction(q_t) * arb(cutoff + 1).log()
    for index, coefficient in coefficients.items():
        quotient = cutoff // int(index)
        if not quotient:
            continue
        compressed = (
            quotient * arb(index).log() + arb(quotient + 1).lgamma()
        )
        result -= _arb_from_fraction(Fraction(coefficient)) * compressed
    return result


def _recurrence_absolute_bound(
    denominator: int,
    numerators: Mapping[int, int],
    cutoff: int,
) -> int:
    return denominator + sum(
        abs(numerator) * (cutoff // index)
        for index, numerator in numerators.items()
    )


def _binary64_absolute_terms(
    q_values: np.ndarray,
    first: int,
    denominator_exponent: int,
) -> np.ndarray:
    """Evaluate one contiguous block of rational terms in binary64.

    The exact inputs are signed int64 numerators for intervals beginning at
    ``first``.  The dyadic scale is exact in the guarded normal range.  The
    square and sequential divisions by ``M`` and ``M + 1`` are the three
    arithmetic roundings beyond a possibly inexact integer conversion.
    """

    exact_q = np.asarray(q_values)
    if exact_q.ndim != 1 or exact_q.size < 1:
        raise ValueError("q values must be a nonempty one-dimensional array")
    if exact_q.dtype != np.dtype(np.int64):
        raise ValueError("q values must have signed int64 dtype")
    exact_first = _exact_nonboolean_integer(first, "first interval")
    if exact_first < 1:
        raise ValueError("first interval must be positive")
    exact_exponent = _exact_nonboolean_integer(
        denominator_exponent,
        "denominator exponent",
    )
    if exact_exponent < 0:
        raise ValueError("denominator exponent must be nonnegative")

    last = exact_first + int(exact_q.size) - 1
    if last + 1 > MAX_EXACT_BINARY64_INTEGER:
        raise ValueError(
            "interval block exceeds exact binary64 integer range"
        )

    terms = exact_q.astype(np.float64)
    np.ldexp(terms, -exact_exponent, out=terms)
    np.multiply(terms, terms, out=terms)
    intervals = np.arange(exact_q.size, dtype=np.float64)
    np.add(intervals, np.float64(exact_first), out=intervals)
    if (
        float(intervals[0]) != float(exact_first)
        or float(intervals[-1]) != float(last)
    ):
        raise ArithmeticError("binary64 interval construction mismatch")
    np.divide(terms, intervals, out=terms)
    np.add(intervals, np.float64(1), out=intervals)
    if (
        float(intervals[0]) != float(exact_first + 1)
        or float(intervals[-1]) != float(last + 1)
    ):
        raise ArithmeticError("binary64 successor interval mismatch")
    np.divide(terms, intervals, out=terms)
    if not np.all(np.isfinite(terms)) or np.any(terms < 0):
        raise ArithmeticError("invalid binary64 absolute-prefix term")
    if np.any((exact_q != 0) & (terms == 0)):
        raise ArithmeticError("binary64 absolute-prefix term underflow")
    return terms


def certified_absolute_energy_prefix(
    coefficients: Mapping[int, Fraction],
    cutoff: int,
    *,
    chunk_size: int = 1 << 24,
    reduction_leaf_size: int = 1 << 20,
    precision_bits: int = 768,
) -> AbsoluteEnergyPrefix:
    """Enclose the arbitrary-slope absolute energy through interval ``T``.

    The returned object contains both the Arb energy enclosure and every
    integer or rational scalar needed to audit the int64 and binary64 proof.
    Chunk boundaries affect only the exact recurrence schedule.  Reduction
    leaves and their final root reduction define the complete floating-point
    addition tree.
    """

    cutoff = _exact_nonboolean_integer(cutoff, "prefix cutoff")
    chunk_size = _exact_nonboolean_integer(chunk_size, "chunk size")
    reduction_leaf_size = _exact_nonboolean_integer(
        reduction_leaf_size,
        "reduction leaf size",
    )
    precision_bits = _exact_nonboolean_integer(
        precision_bits,
        "precision",
    )
    if cutoff < 1:
        raise ValueError("prefix cutoff must be positive")
    if chunk_size < 1:
        raise ValueError("chunk size must be positive")
    if reduction_leaf_size < 1:
        raise ValueError("reduction leaf size must be positive")
    if reduction_leaf_size > chunk_size:
        raise ValueError("reduction leaf size cannot exceed chunk size")
    if precision_bits < 64:
        raise ValueError("precision must be at least 64 bits")
    if cutoff + 1 > MAX_EXACT_BINARY64_INTEGER:
        raise ValueError(
            "cutoff is too large for exact binary64 interval integers"
        )
    require_binary64_tcb()

    cleaned = _clean_coefficients(coefficients)
    denominator, numerators = common_dyadic_numerators(cleaned)
    denominator_exponent = _power_of_two_exponent(denominator)
    if denominator > INT64_MAX:
        raise OverflowError("coefficient denominator exceeds signed int64")
    for numerator in numerators.values():
        if abs(numerator) > INT64_MAX:
            raise OverflowError("coefficient numerator exceeds signed int64")

    recurrence_bound = _recurrence_absolute_bound(
        denominator,
        numerators,
        cutoff,
    )
    if recurrence_bound > INT64_MAX:
        raise OverflowError("q divisor recurrence may overflow signed int64")
    if (
        denominator
        * denominator
        * cutoff
        * (cutoff + 1)
        > 1 << NORMAL_RANGE_LATTICE_EXPONENT
    ):
        raise ValueError("dyadic grid is too fine for the binary64 proof")

    active_indices = np.asarray(
        sorted(index for index in numerators if index <= cutoff),
        dtype=np.int64,
    )
    active_values = np.asarray(
        [numerators[int(index)] for index in active_indices],
        dtype=np.int64,
    )

    previous_q = denominator
    maximum_q = denominator
    leaf_roots: list[float] = []
    chunk_count = 0
    maximum_leaf_length = 0
    maximum_leaf_depth = 0

    for first in range(1, cutoff + 1, chunk_size):
        last = min(cutoff, first + chunk_size - 1)
        length = last - first + 1
        jumps = np.zeros(length, dtype=np.int64)
        for index, numerator in zip(
            active_indices,
            active_values,
            strict=True,
        ):
            divisor = int(index)
            initial = ((first + divisor - 1) // divisor) * divisor
            if initial <= last:
                jumps[initial - first :: divisor] += numerator

        np.cumsum(jumps, dtype=np.int64, out=jumps)
        np.add(jumps, np.int64(previous_q), out=jumps)
        chunk_minimum = int(np.min(jumps))
        chunk_maximum = int(np.max(jumps))
        maximum_q = max(
            maximum_q,
            abs(chunk_minimum),
            abs(chunk_maximum),
        )
        previous_q = int(jumps[-1])

        terms = _binary64_absolute_terms(
            jumps,
            first,
            denominator_exponent,
        )

        for offset in range(0, length, reduction_leaf_size):
            leaf = terms[offset : offset + reduction_leaf_size]
            reduction = explicit_pairwise_sum(leaf)
            leaf_roots.append(reduction.value)
            maximum_leaf_length = max(
                maximum_leaf_length,
                reduction.input_count,
            )
            maximum_leaf_depth = max(
                maximum_leaf_depth,
                reduction.depth,
            )
        chunk_count += 1

    root_reduction = explicit_pairwise_sum(
        np.asarray(leaf_roots, dtype=np.float64)
    )
    midpoint_float = root_reduction.value
    if midpoint_float < 0:
        raise ArithmeticError("negative absolute-prefix midpoint")
    midpoint = Fraction.from_float(midpoint_float)
    reduction_depth = maximum_leaf_depth + root_reduction.depth
    reduction_gamma = _gamma(reduction_depth)
    term_relative_error = (
        (1 + BINARY64_UNIT_ROUNDOFF)
        ** BINARY64_TERM_ROUNDING_FACTOR_COUNT
        - 1
    )
    computed_absolute_upper = midpoint / (1 - reduction_gamma)
    reduction_error = reduction_gamma * computed_absolute_upper
    term_evaluation_error = (
        term_relative_error
        / (1 - term_relative_error)
        * computed_absolute_upper
    )
    rational_error = reduction_error + term_evaluation_error

    direct_last = denominator + sum(
        numerator * (cutoff // index)
        for index, numerator in numerators.items()
    )
    if previous_q != direct_last:
        raise ArithmeticError("q divisor recurrence endpoint mismatch")
    if maximum_q > recurrence_bound:
        raise ArithmeticError("observed q exceeds its exact int64 bound")

    p_value = _harmonic_sum(cleaned)
    slope_square = (cutoff + 1) * p_value * p_value
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        logarithmic_sum = _compressed_absolute_log_sum(
            cleaned,
            Fraction(previous_q, denominator),
            cutoff,
        )
        logarithmic_cross = (
            -2 * _arb_from_fraction(p_value) * logarithmic_sum
        )
        rational_ball = _arb_from_fraction(midpoint) + _radius_ball(
            rational_error,
            precision_bits,
        )
        prefix_energy = (
            rational_ball
            + logarithmic_cross
            + _arb_from_fraction(slope_square)
        )
    finally:
        ctx.prec = previous_precision

    return AbsoluteEnergyPrefix(
        algorithm_version=PREFIX_ALGORITHM_VERSION,
        term_evaluation_algorithm=TERM_EVALUATION_ALGORITHM,
        reduction_algorithm=REDUCTION_ALGORITHM,
        cutoff=cutoff,
        chunk_size=chunk_size,
        chunk_count=chunk_count,
        reduction_leaf_size=reduction_leaf_size,
        reduction_leaf_count=len(leaf_roots),
        maximum_leaf_length=maximum_leaf_length,
        maximum_leaf_reduction_depth=maximum_leaf_depth,
        root_reduction_depth=root_reduction.depth,
        reduction_depth_bound=reduction_depth,
        coefficient_denominator=denominator,
        coefficient_denominator_exponent=denominator_exponent,
        recurrence_absolute_bound=recurrence_bound,
        int64_headroom=INT64_MAX - recurrence_bound,
        term_rounding_factor_count=BINARY64_TERM_ROUNDING_FACTOR_COUNT,
        term_relative_error=term_relative_error,
        reduction_gamma=reduction_gamma,
        rational_midpoint=midpoint,
        rational_error=rational_error,
        rational_lower=midpoint - rational_error,
        rational_upper=midpoint + rational_error,
        computed_absolute_term_sum_upper=computed_absolute_upper,
        term_evaluation_error=term_evaluation_error,
        reduction_error=reduction_error,
        logarithmic_sum=logarithmic_sum,
        logarithmic_cross_term=logarithmic_cross,
        slope_square_term=slope_square,
        prefix_energy=prefix_energy,
        last_q_numerator=previous_q,
        maximum_q_numerator=maximum_q,
        maximum_square_numerator=maximum_q * maximum_q,
    )
