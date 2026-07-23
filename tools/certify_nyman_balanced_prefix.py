"""Rigorous finite-prefix evaluator for one exact balanced Nyman direction.

The divisor recurrences and interval numerators are evaluated exactly in
signed 64-bit integers.  NumPy is used only to divide and reduce the resulting
exact integer terms in binary64 blocks.  A conservative, exact rational
forward-error bound encloses conversion, division, every within-block
reduction, and the final block reduction.

The logarithmic part is not summed over the cutoff.  Abel summation compresses
it to one log and one log-gamma evaluation per nonzero added coefficient, all
evaluated with Arb.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Any, Mapping

from flint import arb, ctx
import numpy as np

if __package__:
    from tools.certify_nyman_balanced_tail import harmonic_sum, require_balanced
else:
    from certify_nyman_balanced_tail import harmonic_sum, require_balanced


ZERO = Fraction(0)
INT64_MAX = (1 << 63) - 1
BINARY64_UNIT_ROUNDOFF = Fraction(1, 1 << 53)
CERTIFICATE_NUMPY_VERSION = "2.3.5"


def require_binary64_tcb() -> None:
    """Require the pinned NumPy/IEEE environment used by the error proof."""

    information = np.finfo(np.float64)
    if np.__version__ != CERTIFICATE_NUMPY_VERSION:
        raise RuntimeError(
            "the certificate requires NumPy " + CERTIFICATE_NUMPY_VERSION
        )
    if not (
        information.nmant == 52
        and information.maxexp == 1024
        and information.minexp == -1022
        and np.dtype(np.float64).itemsize == 8
    ):
        raise RuntimeError("the certificate requires IEEE-754 binary64")
    probes = np.asarray([1.0, 3.0, -54.0, -55.0, -53.0], dtype=np.float64)
    one = probes[0]
    three = probes[1]
    upward_offset = np.ldexp(three, int(probes[2]))
    downward_offset = np.ldexp(three, int(probes[3]))
    tie_offset = np.ldexp(one, int(probes[4]))
    next_up = np.nextafter(one, np.float64(np.inf))
    next_down = np.nextafter(one, np.float64(-np.inf))
    if np.add(one, upward_offset) != next_up:
        raise RuntimeError("binary64 round-to-nearest is required")
    if np.subtract(one, downward_offset) != next_down:
        raise RuntimeError("binary64 round-to-nearest is required")
    if np.add(one, tie_offset) != one:
        raise RuntimeError("binary64 round-to-nearest-ties-to-even is required")


def _power_of_two_exponent(value: int) -> int:
    if value < 1 or value & (value - 1):
        raise ValueError("coefficient denominators must be powers of two")
    return value.bit_length() - 1


def common_dyadic_numerators(
    values: Mapping[int, Fraction],
) -> tuple[int, dict[int, int]]:
    """Put exact dyadic coefficients over their least common denominator."""

    if not values:
        return 1, {}
    exponent = max(
        _power_of_two_exponent(Fraction(value).denominator)
        for value in values.values()
    )
    denominator = 1 << exponent
    numerators: dict[int, int] = {}
    for raw_index, raw_value in values.items():
        index = int(raw_index)
        value = Fraction(raw_value)
        if index < 1:
            raise ValueError("coefficient indices must be positive")
        numerator = value * denominator
        if numerator.denominator != 1:
            raise ArithmeticError("failed to lift dyadic coefficient")
        if numerator:
            numerators[index] = numerator.numerator
    return denominator, numerators


def _gamma(operation_count: int) -> Fraction:
    if operation_count < 0:
        raise ValueError("operation count must be nonnegative")
    product = operation_count * BINARY64_UNIT_ROUNDOFF
    if product >= 1:
        raise ArithmeticError("binary64 error bound exhausted")
    return product / (1 - product)


def _arb_from_fraction(value: Fraction) -> Any:
    exact = Fraction(value)
    return arb(exact.numerator) / exact.denominator


def _radius_ball(value: Fraction, bits: int = 160) -> Any:
    """Return an Arb zero ball whose dyadic radius is at least ``value``."""

    exact = Fraction(value)
    if exact < 0:
        raise ValueError("radius must be nonnegative")
    if not exact:
        return arb(0)
    scale = 1 << bits
    numerator = (exact.numerator * scale + exact.denominator - 1) // exact.denominator
    return arb(0, (numerator, -bits))


@dataclass(frozen=True)
class PrefixBound:
    cutoff: int
    block_size: int
    block_count: int
    added_denominator: int
    old_denominator: int
    rational_midpoint: Fraction
    rational_error: Fraction
    rational_lower: Fraction
    rational_upper: Fraction
    computed_absolute_term_sum_upper: Fraction
    logarithmic_gain: Any
    prefix_gain: Any
    last_added_step_numerator: int
    last_old_intercept_numerator: int
    maximum_added_step_numerator: int
    maximum_old_intercept_numerator: int
    maximum_interval_numerator: int


def compressed_logarithmic_sum(
    added_numerators: Mapping[int, int],
    added_denominator: int,
    cutoff: int,
) -> Any:
    """Evaluate ``sum_(M<=T) g_M log((M+1)/M)`` with Arb."""

    if cutoff < 1:
        raise ValueError("prefix cutoff must be positive")
    last_step = -sum(
        numerator * (cutoff // index)
        for index, numerator in added_numerators.items()
    )
    total = _arb_from_fraction(Fraction(last_step, added_denominator)) * arb(
        cutoff + 1
    ).log()
    for index, numerator in added_numerators.items():
        quotient = cutoff // index
        if not quotient:
            continue
        logarithm = quotient * arb(index).log() + arb(quotient + 1).lgamma()
        total += _arb_from_fraction(
            Fraction(numerator, added_denominator)
        ) * logarithm
    return total


def _recurrence_operation_bound(
    numerators: Mapping[int, int],
    cutoff: int,
) -> int:
    return sum(
        abs(numerator) * (cutoff // index)
        for index, numerator in numerators.items()
    )


def certified_prefix_gain(
    old: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    cutoff: int,
    *,
    block_size: int = 1 << 20,
    precision_bits: int = 192,
) -> PrefixBound:
    """Enclose the direct gain over the integer intervals ``1,...,T``.

    For ``g_M=-sum_n a_n floor(M/n)`` and
    ``q_M=1+sum_m p_m floor(M/m)``, the rational interval contribution is

        sum_M (2*g_M*q_M-g_M^2)/(M*(M+1)).

    The remaining contribution is

        -2*P1*sum_M g_M*log((M+1)/M),

    where ``P1=sum_m p_m/m``.  The first sum is enclosed by an explicit
    binary64 forward-error calculation; the second is evaluated by Arb after
    exact Abel compression.
    """

    if cutoff < 1:
        raise ValueError("prefix cutoff must be positive")
    if block_size < 1:
        raise ValueError("block size must be positive")
    if precision_bits < 64:
        raise ValueError("precision must be at least 64 bits")
    require_binary64_tcb()
    if cutoff * (cutoff + 1) >= 1 << 53:
        raise ValueError("cutoff is too large for exact binary64 denominators")
    require_balanced(added, "added vector")

    added_denominator, added_numerators = common_dyadic_numerators(added)
    old_denominator, old_numerators = common_dyadic_numerators(old)
    if added_denominator % old_denominator:
        raise ValueError(
            "added denominator must be a multiple of the old denominator"
        )
    # This conservative lattice guard rules out underflow not only in term
    # evaluation but also after arbitrary cancellation inside NumPy's
    # reduction tree.  Every nonzero computed term then has exponent at least
    # -901, so every exact partial sum is zero or a multiple of 2^-953 and is
    # therefore normal.  The frozen certificate is far inside this boundary.
    if (
        added_denominator
        * added_denominator
        * cutoff
        * (cutoff + 1)
        > 1 << 900
    ):
        raise ValueError("dyadic grid is too fine for the binary64 proof")
    cross_scale = 2 * (added_denominator // old_denominator)
    if cross_scale > INT64_MAX:
        raise OverflowError("cross-term scale exceeds signed int64")

    # These absolute recurrence bounds prove that every int64 cumsum below is
    # safe independently of cancellation.
    if _recurrence_operation_bound(added_numerators, cutoff) > INT64_MAX:
        raise OverflowError("added divisor recurrence may overflow int64")
    if (
        old_denominator
        + _recurrence_operation_bound(old_numerators, cutoff)
        > INT64_MAX
    ):
        raise OverflowError("old divisor recurrence may overflow int64")

    added_indices = np.asarray(sorted(added_numerators), dtype=np.int64)
    added_values = np.asarray(
        [added_numerators[int(index)] for index in added_indices],
        dtype=np.int64,
    )
    old_indices = np.asarray(sorted(old_numerators), dtype=np.int64)
    old_values = np.asarray(
        [old_numerators[int(index)] for index in old_indices],
        dtype=np.int64,
    )

    previous_added = np.int64(0)
    previous_old = np.int64(old_denominator)
    block_sums: list[float] = []
    block_absolute_sums: list[Fraction] = []
    block_lengths: list[int] = []
    maximum_added = 0
    maximum_old = old_denominator
    maximum_interval = 0

    for first in range(1, cutoff + 1, block_size):
        last = min(cutoff, first + block_size - 1)
        length = last - first + 1
        added_jumps = np.zeros(length, dtype=np.int64)
        old_jumps = np.zeros(length, dtype=np.int64)
        for index, numerator in zip(
            added_indices,
            added_values,
            strict=True,
        ):
            divisor = int(index)
            initial = ((first + divisor - 1) // divisor) * divisor
            if initial <= last:
                added_jumps[initial - first :: divisor] += numerator
        for index, numerator in zip(old_indices, old_values, strict=True):
            divisor = int(index)
            initial = ((first + divisor - 1) // divisor) * divisor
            if initial <= last:
                old_jumps[initial - first :: divisor] += numerator

        added_steps = previous_added - np.cumsum(
            added_jumps,
            dtype=np.int64,
        )
        old_intercepts = previous_old + np.cumsum(
            old_jumps,
            dtype=np.int64,
        )
        maximum_added_block = int(np.max(np.abs(added_steps)))
        maximum_old_block = int(np.max(np.abs(old_intercepts)))
        maximum_added = max(maximum_added, maximum_added_block)
        maximum_old = max(maximum_old, maximum_old_block)

        if maximum_added_block * maximum_old_block > INT64_MAX:
            raise OverflowError("interval cross product may overflow int64")
        products = added_steps * old_intercepts
        maximum_product = int(np.max(np.abs(products)))
        if (
            cross_scale * maximum_product
            + maximum_added_block * maximum_added_block
            > INT64_MAX
        ):
            raise OverflowError("interval numerator may overflow int64")
        interval_numerators = (
            cross_scale * products - added_steps * added_steps
        )
        maximum_interval = max(
            maximum_interval,
            int(np.max(np.abs(interval_numerators))),
        )

        intervals = np.arange(first, last + 1, dtype=np.float64)
        denominators = intervals * (intervals + 1.0)
        terms = (
            interval_numerators.astype(np.float64)
            / (added_denominator * added_denominator)
            / denominators
        )
        if not np.all(np.isfinite(terms)):
            raise ArithmeticError("nonfinite binary64 prefix term")
        block_sum = float(np.sum(terms, dtype=np.float64))
        block_absolute = float(
            np.sum(np.abs(terms), dtype=np.float64)
        )
        if not math.isfinite(block_sum) or not math.isfinite(block_absolute):
            raise ArithmeticError("nonfinite binary64 block reduction")
        block_sums.append(block_sum)
        block_absolute_sums.append(Fraction.from_float(block_absolute))
        block_lengths.append(length)
        previous_added = added_steps[-1]
        previous_old = old_intercepts[-1]

    rational_midpoint_float = 0.0
    for block_sum in block_sums:
        rational_midpoint_float += block_sum
    if not math.isfinite(rational_midpoint_float):
        raise ArithmeticError("nonfinite binary64 prefix midpoint")
    rational_midpoint = Fraction.from_float(rational_midpoint_float)

    direct_last_added = -sum(
        numerator * (cutoff // index)
        for index, numerator in added_numerators.items()
    )
    direct_last_old = old_denominator + sum(
        numerator * (cutoff // index)
        for index, numerator in old_numerators.items()
    )
    if int(previous_added) != direct_last_added:
        raise ArithmeticError("added divisor recurrence endpoint mismatch")
    if int(previous_old) != direct_last_old:
        raise ArithmeticError("old divisor recurrence endpoint mismatch")

    # One correctly rounded integer-to-binary64 conversion followed by one
    # correctly rounded division has relative error at most (1+u)^2-1.
    term_relative_error = (
        2 * BINARY64_UNIT_ROUNDOFF
        + BINARY64_UNIT_ROUNDOFF * BINARY64_UNIT_ROUNDOFF
    )
    computed_absolute_upper = ZERO
    block_reduction_error = ZERO
    for length, block_absolute in zip(
        block_lengths,
        block_absolute_sums,
        strict=True,
    ):
        reduction_gamma = _gamma(max(0, length - 1))
        exact_computed_absolute = block_absolute / (1 - reduction_gamma)
        computed_absolute_upper += exact_computed_absolute
        block_reduction_error += reduction_gamma * exact_computed_absolute
    term_evaluation_error = (
        term_relative_error
        / (1 - term_relative_error)
        * computed_absolute_upper
    )
    final_reduction_error = _gamma(len(block_sums)) * sum(
        (Fraction.from_float(abs(value)) for value in block_sums),
        start=ZERO,
    )
    rational_error = (
        term_evaluation_error
        + block_reduction_error
        + final_reduction_error
    )

    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        logarithmic_sum = compressed_logarithmic_sum(
            added_numerators,
            added_denominator,
            cutoff,
        )
        logarithmic_gain = (
            -2 * _arb_from_fraction(harmonic_sum(old)) * logarithmic_sum
        )
        rational_ball = _arb_from_fraction(rational_midpoint) + _radius_ball(
            rational_error
        )
        prefix_gain = rational_ball + logarithmic_gain
    finally:
        ctx.prec = previous_precision

    return PrefixBound(
        cutoff=cutoff,
        block_size=block_size,
        block_count=len(block_sums),
        added_denominator=added_denominator,
        old_denominator=old_denominator,
        rational_midpoint=rational_midpoint,
        rational_error=rational_error,
        rational_lower=rational_midpoint - rational_error,
        rational_upper=rational_midpoint + rational_error,
        computed_absolute_term_sum_upper=computed_absolute_upper,
        logarithmic_gain=logarithmic_gain,
        prefix_gain=prefix_gain,
        last_added_step_numerator=int(previous_added),
        last_old_intercept_numerator=int(previous_old),
        maximum_added_step_numerator=maximum_added,
        maximum_old_intercept_numerator=maximum_old,
        maximum_interval_numerator=maximum_interval,
    )
