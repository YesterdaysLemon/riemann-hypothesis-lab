"""Exact fixed-modulus spacing bounds for finite Nyman residual tails.

For a reduced active Fourier frequency ``a/d`` and another active reduced
frequency ``b/e``, circular Farey separation gives

    ||a/d-b/e|| >= gcd(d,e)/(d*e).

Fix a positive integer ``L`` and put ``m=gcd(d,L)``.  Since ``m`` divides
``d``, ``gcd(d,e) >= gcd(m,e)``.  Thus the reciprocal local gap at
denominator ``d`` is at most

    d * R_m,  R_m=max_(active e) e/gcd(m,e).

Montgomery--Vaughan's weighted periodic cosecant inequality then gives a
consecutive-interval discrepancy constant computable with only divisor sums,
Jordan ``J_2``, and one scan of the active denominators per divisor class of
``L``.  ``L=1`` recovers the earlier active-support ``(3/2) Q sigma`` bound.

The ``target`` parameter supports both the ordinary residual

    target*1_[1,infinity) - sum_n p_n {t/n}

and the target-zero increment needed for a rigorous Minkowski split.
Everything returned here is exact rational arithmetic.  This module supplies
finite-support tail bounds only; it makes no all-scale or RH claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys
from typing import Mapping


if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    from tools import generate_nyman_rebased_schur_certificate as parent
except ImportError:  # pragma: no cover - direct script import fallback
    import generate_nyman_rebased_schur_certificate as parent


@dataclass(frozen=True)
class FixedModulusTail:
    cutoff: int
    modulus: int
    support_limit: int
    active_denominator_count: int
    active_denominator_maximum: int
    group_reciprocal_bounds: tuple[tuple[int, int], ...]
    target: Fraction
    coefficient_sum: Fraction
    harmonic_sum: Fraction
    c0: Fraction
    rho: Fraction
    sigma: Fraction
    fixed_modulus_spacing_constant: Fraction
    active_q_spacing_constant: Fraction
    support_q_spacing_constant: Fraction
    active_q_to_fixed_modulus_improvement: Fraction
    absolute_residual_bound: Fraction
    mean_term: Fraction
    discrepancy_term: Fraction
    cross_term: Fraction
    slope_square_term: Fraction
    upper_bound: Fraction


def _validate_positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value


def _validate_target(value: Fraction) -> Fraction:
    if not isinstance(value, Fraction):
        raise TypeError("target must be an exact Fraction")
    return value


def fixed_modulus_tail_upper(
    coefficients: Mapping[int, Fraction],
    cutoff: int,
    modulus: int,
    *,
    target: Fraction = Fraction(1),
) -> FixedModulusTail:
    """Return an exact complete-tail upper bound with grouped local gaps."""

    _validate_positive_integer(cutoff, "cutoff")
    _validate_positive_integer(modulus, "modulus")
    target = _validate_target(target)
    cleaned = parent._clean(coefficients)
    support_limit = max(cleaned, default=0)
    coefficient_sum = parent._coefficient_sum(cleaned)
    harmonic_sum = parent._harmonic_sum(cleaned)
    c0 = target - coefficient_sum / 2

    if support_limit:
        divisor_sums = parent._divisor_harmonic_sums(
            cleaned,
            support_limit,
        )
        jordan = parent._jordan_j2_sieve(support_limit)
        active = tuple(
            divisor
            for divisor in range(2, support_limit + 1)
            if divisor_sums[divisor]
        )
        rho_sum = sum(
            (
                jordan[divisor]
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in active
            ),
            start=Fraction(),
        )
        sigma_sum = sum(
            (
                divisor
                * jordan[divisor]
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in active
            ),
            start=Fraction(),
        )
    else:
        divisor_sums = [Fraction()]
        jordan = [0]
        active = ()
        rho_sum = Fraction()
        sigma_sum = Fraction()

    rho = c0 * c0 + rho_sum / 12
    sigma = c0 * c0 + sigma_sum / 12
    if rho < 0 or sigma < 0:
        raise ArithmeticError("periodic residual moment became negative")

    if active:
        active_maximum = max(active)
        groups = tuple(sorted({gcd(divisor, modulus) for divisor in active}))
        reciprocal_by_group = {
            group: max(
                denominator // gcd(group, denominator)
                for denominator in active
            )
            for group in groups
        }
        weighted_mass = sum(
            (
                Fraction(
                    divisor
                    * reciprocal_by_group[gcd(divisor, modulus)]
                    * jordan[divisor],
                    12,
                )
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in active
            ),
            start=Fraction(),
        )
        fixed_constant = Fraction(3, 2) * (
            active_maximum * c0 * c0 + weighted_mass
        )
        active_q_constant = Fraction(3, 2) * active_maximum * sigma
        group_record = tuple(
            (group, reciprocal_by_group[group]) for group in groups
        )
    else:
        active_maximum = 0
        fixed_constant = Fraction()
        active_q_constant = Fraction()
        group_record = ()

    support_q_constant = (
        Fraction(3, 2) * support_limit * sigma
        if support_limit
        else Fraction()
    )
    if fixed_constant > active_q_constant:
        raise ArithmeticError("fixed-modulus spacing bound became weaker")
    improvement = (
        active_q_constant / fixed_constant
        if fixed_constant
        else Fraction()
    )

    absolute_residual_bound = abs(c0) + sum(
        (abs(value) for value in cleaned.values()),
        start=Fraction(),
    ) / 2
    mean_term = rho / (cutoff + 1)
    discrepancy_term = fixed_constant / (
        (cutoff + 1) * (cutoff + 2)
    )
    cross_term = (
        abs(harmonic_sum) * absolute_residual_bound / (cutoff + 1)
    )
    slope_square_term = (
        harmonic_sum * harmonic_sum / (4 * (cutoff + 1))
    )
    upper_bound = (
        mean_term
        + discrepancy_term
        + cross_term
        + slope_square_term
    )
    return FixedModulusTail(
        cutoff=cutoff,
        modulus=modulus,
        support_limit=support_limit,
        active_denominator_count=len(active),
        active_denominator_maximum=active_maximum,
        group_reciprocal_bounds=group_record,
        target=target,
        coefficient_sum=coefficient_sum,
        harmonic_sum=harmonic_sum,
        c0=c0,
        rho=rho,
        sigma=sigma,
        fixed_modulus_spacing_constant=fixed_constant,
        active_q_spacing_constant=active_q_constant,
        support_q_spacing_constant=support_q_constant,
        active_q_to_fixed_modulus_improvement=improvement,
        absolute_residual_bound=absolute_residual_bound,
        mean_term=mean_term,
        discrepancy_term=discrepancy_term,
        cross_term=cross_term,
        slope_square_term=slope_square_term,
        upper_bound=upper_bound,
    )


def minkowski_tail_upper(
    first_upper: Fraction,
    second_upper: Fraction,
    tradeoff: Fraction,
) -> Fraction:
    """Rationally bound ``(sqrt(first)+sqrt(second))**2`` from above.

    Young's inequality ``2*sqrt(a*b) <= t*a+b/t`` gives the returned
    expression.  A caller can choose an exact rational ``t`` near
    ``sqrt(second/first)`` without putting an irrational in an artifact.
    """

    for value, name in (
        (first_upper, "first_upper"),
        (second_upper, "second_upper"),
        (tradeoff, "tradeoff"),
    ):
        if not isinstance(value, Fraction):
            raise TypeError(f"{name} must be an exact Fraction")
    if first_upper < 0 or second_upper < 0:
        raise ValueError("tail upper bounds must be nonnegative")
    if tradeoff <= 0:
        raise ValueError("tradeoff must be positive")
    if not first_upper:
        return second_upper
    if not second_upper:
        return first_upper
    return (
        (1 + tradeoff) * first_upper
        + (1 + 1 / tradeoff) * second_upper
    )


__all__ = [
    "FixedModulusTail",
    "fixed_modulus_tail_upper",
    "minkowski_tail_upper",
]
