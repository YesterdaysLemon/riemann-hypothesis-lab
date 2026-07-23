from __future__ import annotations

from fractions import Fraction
from math import lcm
from pathlib import Path

import pytest

from tools import generate_nyman_rebased_third_step_certificate as third
from tools.certify_nyman_fixed_modulus_tail import (
    fixed_modulus_tail_upper,
    minkowski_tail_upper,
)


ROOT = Path(__file__).resolve().parents[1]


def _small_vector() -> dict[int, Fraction]:
    return {
        2: Fraction(1, 2),
        3: Fraction(-1, 4),
        6: Fraction(1, 8),
    }


def test_modulus_one_recovers_active_q_spacing_constant() -> None:
    tail = fixed_modulus_tail_upper(_small_vector(), 100, 1)
    assert tail.fixed_modulus_spacing_constant == (
        tail.active_q_spacing_constant
    )
    assert tail.active_q_to_fixed_modulus_improvement == 1
    assert tail.active_denominator_maximum <= tail.support_limit


def test_active_empty_n1_vector_has_zero_discrepancy_constant() -> None:
    tail = fixed_modulus_tail_upper({1: Fraction(1)}, 100, 840)
    assert tail.active_denominator_count == 0
    assert tail.active_denominator_maximum == 0
    assert tail.group_reciprocal_bounds == ()
    assert tail.fixed_modulus_spacing_constant == 0
    assert tail.discrepancy_term == 0


def test_fixed_modulus_constant_bounds_every_small_period_interval() -> None:
    coefficients = _small_vector()
    tail = fixed_modulus_tail_upper(coefficients, 100, 6)
    period = lcm(*coefficients)
    values: list[Fraction] = []
    for interval in range(period):
        centered = tail.c0 - sum(
            value
            * (
                Fraction(interval % index, index)
                - Fraction(index - 1, 2 * index)
            )
            for index, value in coefficients.items()
        )
        values.append(centered * centered - tail.rho)

    assert sum(values, start=Fraction()) == 0
    for first in range(period):
        running = Fraction()
        for length in range(1, period + 1):
            running += values[(first + length - 1) % period]
            assert abs(running) <= tail.fixed_modulus_spacing_constant


def test_target_zero_tail_is_exactly_quadratically_homogeneous() -> None:
    coefficients = _small_vector()
    scale = Fraction(-3, 7)
    scaled = {
        index: scale * value for index, value in coefficients.items()
    }
    first = fixed_modulus_tail_upper(
        coefficients,
        1 << 12,
        30,
        target=Fraction(),
    )
    second = fixed_modulus_tail_upper(
        scaled,
        1 << 12,
        30,
        target=Fraction(),
    )
    factor = scale * scale
    for name in (
        "rho",
        "sigma",
        "fixed_modulus_spacing_constant",
        "mean_term",
        "discrepancy_term",
        "cross_term",
        "slope_square_term",
        "upper_bound",
    ):
        assert getattr(second, name) == factor * getattr(first, name)


def test_minkowski_young_bound_is_exact_and_above_the_square_root_form() -> None:
    first = Fraction(5, 17)
    second = Fraction(7, 29)
    tradeoff = Fraction(9, 11)
    upper = minkowski_tail_upper(first, second, tradeoff)
    assert upper == (
        (1 + tradeoff) * first
        + (1 + 1 / tradeoff) * second
    )
    excess = upper - first - second
    assert excess >= 0
    assert excess * excess >= 4 * first * second
    assert minkowski_tail_upper(first, Fraction(), tradeoff) == first
    assert minkowski_tail_upper(Fraction(), second, tradeoff) == second


@pytest.fixture(scope="module")
def p3_tail():
    witness = third._construct_witness(ROOT)
    return fixed_modulus_tail_upper(witness.p3, 1 << 23, 840)


def test_p3_modulus_840_gives_the_frozen_exact_improvement(p3_tail) -> None:
    assert p3_tail.active_denominator_count == 17_157
    assert p3_tail.active_denominator_maximum == 32_768
    assert len(p3_tail.group_reciprocal_bounds) == 32
    assert p3_tail.group_reciprocal_bounds == (
        (1, 32768),
        (2, 28665),
        (3, 32768),
        (4, 28665),
        (5, 32768),
        (6, 28651),
        (7, 32768),
        (8, 28665),
        (10, 28651),
        (12, 28651),
        (14, 20465),
        (15, 32768),
        (20, 28651),
        (21, 32768),
        (24, 28651),
        (28, 20465),
        (30, 28651),
        (35, 32768),
        (40, 28651),
        (42, 20465),
        (56, 20465),
        (60, 28651),
        (70, 16384),
        (84, 20465),
        (105, 32768),
        (120, 28651),
        (140, 12279),
        (168, 20465),
        (210, 16384),
        (280, 12279),
        (420, 8192),
        (840, 4096),
    )
    assert Fraction(1816, 1) * 1_000_000 < (
        p3_tail.fixed_modulus_spacing_constant
    ) < Fraction(1817, 1) * 1_000_000
    assert p3_tail.active_q_to_fixed_modulus_improvement > Fraction(118, 100)
    assert p3_tail.upper_bound < Fraction(27, 1_000_000)


@pytest.mark.parametrize(
    ("cutoff", "modulus", "error"),
    [
        (0, 1, ValueError),
        (1, 0, ValueError),
        (True, 1, TypeError),
        (1, False, TypeError),
    ],
)
def test_rejects_invalid_integer_contracts(cutoff, modulus, error) -> None:
    with pytest.raises(error):
        fixed_modulus_tail_upper({}, cutoff, modulus)


def test_rejects_nonexact_target_and_invalid_minkowski_inputs() -> None:
    with pytest.raises(TypeError, match="target"):
        fixed_modulus_tail_upper({}, 1, 1, target=1)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="nonnegative"):
        minkowski_tail_upper(Fraction(-1), Fraction(), Fraction(1))
    with pytest.raises(ValueError, match="positive"):
        minkowski_tail_upper(Fraction(), Fraction(), Fraction())
