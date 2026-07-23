from __future__ import annotations

from fractions import Fraction
import inspect
from typing import Any

from flint import arb, ctx
import pytest

import riemann_lab.nyman_oracle as oracle


def _arb_fraction(value: Fraction) -> Any:
    return arb(value.numerator) / value.denominator


def _vasyunin_sum_for_test(numerator: int, denominator: int) -> Any:
    """Test-only closed-form side of the independent comparison."""

    total = arb(0)
    for index in range(1, denominator):
        fractional_numerator = (index * numerator) % denominator
        fractional_part = arb(fractional_numerator) / denominator
        total += fractional_part * (arb.pi() * index / denominator).cot()
    return total


def _closed_autocorrelation_for_test(
    ratio: Fraction, *, unreduced_pair: tuple[int, int] | None = None
) -> Any:
    if unreduced_pair is None:
        numerator, denominator = ratio.numerator, ratio.denominator
    else:
        numerator, denominator = unreduced_pair
    lambda_ball = arb(numerator) / denominator
    first = (arb(1) - lambda_ball) * lambda_ball.log() / 2
    second = (
        (arb(1) + lambda_ball)
        * ((arb(2) * arb.pi()).log() - arb.const_euler())
        / 2
    )
    correction = (
        arb.pi()
        * (
            _vasyunin_sum_for_test(numerator, denominator)
            + _vasyunin_sum_for_test(denominator, numerator)
        )
        / (2 * denominator)
    )
    return first + second - correction


def _mutated_closed_autocorrelation_for_test(
    ratio: Fraction, mutation: str
) -> Any:
    numerator = ratio.numerator
    denominator = ratio.denominator
    lambda_ball = arb(numerator) / denominator
    first = (arb(1) - lambda_ball) * lambda_ball.log() / 2
    if mutation == "reversed_log_term":
        first = -first
    second = (
        (arb(1) + lambda_ball)
        * ((arb(2) * arb.pi()).log() - arb.const_euler())
        / 2
    )
    forward = _vasyunin_sum_for_test(numerator, denominator)
    reciprocal = _vasyunin_sum_for_test(denominator, numerator)
    if mutation == "omit_forward_v":
        forward = arb(0)
    if mutation == "omit_reciprocal_v":
        reciprocal = arb(0)
    outer_denominator = numerator if mutation == "outer_p" else denominator
    correction = arb.pi() * (forward + reciprocal) / (2 * outer_denominator)
    return (
        first + second + correction
        if mutation == "reversed_v_sign"
        else first + second - correction
    )


def test_oracle_source_is_independent_of_formula_implementation() -> None:
    source = inspect.getsource(oracle)
    assert "from .nyman import" not in source
    assert "import riemann_lab.nyman" not in source
    assert "const_euler" not in source


@pytest.mark.parametrize(
    ("ratio", "truncation", "expected"),
    [
        (
            Fraction(2, 3),
            2,
            (Fraction(0), Fraction(1), Fraction(3, 2), Fraction(2)),
        ),
        (
            Fraction(3, 2),
            2,
            (
                Fraction(0),
                Fraction(2, 3),
                Fraction(1),
                Fraction(4, 3),
                Fraction(2),
            ),
        ),
    ],
)
def test_breakpoints_are_exact_floor_jumps(
    ratio: Fraction, truncation: int, expected: tuple[Fraction, ...]
) -> None:
    assert oracle.rational_breakpoints(ratio, truncation) == expected


def test_piece_formula_rejects_reversed_log_reciprocal_sign_and_lambda_factor_mutants() -> None:
    ratio = Fraction(2, 3)
    lower = Fraction(3, 2)
    upper = Fraction(2)
    floor_t = 1
    floor_ratio_t = 1

    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        trusted = oracle._integrate_floor_constant_piece(
            ratio, lower, upper, floor_t, floor_ratio_t
        )
        linear = _arb_fraction(ratio * (upper - lower))
        coefficient = _arb_fraction(Fraction(floor_ratio_t) + ratio * floor_t)
        log_forward = _arb_fraction(upper / lower).log()
        reciprocal = _arb_fraction(
            Fraction(floor_t * floor_ratio_t)
            * (Fraction(1) / lower - Fraction(1) / upper)
        )

        reversed_log = (
            linear - coefficient * _arb_fraction(lower / upper).log() + reciprocal
        )
        reversed_reciprocal_sign = linear - coefficient * log_forward - reciprocal
        missing_lambda_factor = (
            linear
            - arb(floor_ratio_t + floor_t) * log_forward
            + reciprocal
        )

        assert trusted > 0
        assert not trusted.overlaps(reversed_log)
        assert not trusted.overlaps(reversed_reciprocal_sign)
        assert not trusted.overlaps(missing_lambda_factor)
    finally:
        ctx.prec = previous_precision


def test_truncated_oracle_contains_known_a_of_one_and_restores_precision() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 192
        expected = (arb(2) * arb.pi()).log() - arb.const_euler()
        result = oracle.truncated_inner_product_oracle(
            Fraction(1), truncation=4096, precision_bits=128
        )
        assert ctx.prec == 192
        assert result.total.contains(expected)
        assert result.tail.contains(arb(0))
        assert result.tail.contains(arb(1) / 4096)
        assert result.breakpoint_count == 4097
    finally:
        ctx.prec = previous_precision


def test_supplied_ratio_reduction_is_explicit_and_numerically_benign() -> None:
    result = oracle.truncated_inner_product_oracle(
        (4, 6), truncation=64, precision_bits=96
    )
    record = result.to_record()
    assert result.ratio == Fraction(2, 3)
    assert record["ratio"] == {
        "supplied": {"numerator": "4", "denominator": "6"},
        "reduced": {"numerator": "2", "denominator": "3"},
        "supplied_was_reduced": False,
        "canonicalization": (
            "fractions.Fraction gcd reduction with positive denominator"
        ),
    }

    previous_precision = ctx.prec
    try:
        ctx.prec = 160
        reduced = _closed_autocorrelation_for_test(Fraction(2, 3))
        unreduced = _closed_autocorrelation_for_test(
            Fraction(2, 3), unreduced_pair=(4, 6)
        )
        assert reduced.overlaps(unreduced)
    finally:
        ctx.prec = previous_precision


def test_harmonic_sum_encloses_one_minus_euler_gamma_without_builtin_constant() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 192
        expected = arb(1) - arb.const_euler()
        result = oracle.one_minus_euler_gamma_oracle(4096, 128)
        assert ctx.prec == 192
        assert result.enclosure.contains(expected)
        record = result.to_record()
        assert record["classification"] == "EXPLORATORY"
        assert record["hypothesis_status"] == "UNRESOLVED"
        assert record["bound"]["correction_lower"] == {
            "numerator": "1",
            "denominator": "8194",
        }
        assert record["bound"]["correction_upper"] == {
            "numerator": "1",
            "denominator": "8192",
        }
        assert record == result.to_record()
    finally:
        ctx.prec = previous_precision


def test_frozen_ratio_table_is_exactly_ordered_reproducible_and_twelve_bit_scoped() -> None:
    assert oracle.FROZEN_RATIO_PAIRS == (
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 3),
        (3, 2),
        (4, 6),
        (5, 7),
        (7, 5),
        (31, 32),
        (32, 31),
        (127, 128),
        (128, 127),
        (251, 256),
        (256, 251),
    )
    first = oracle.generate_frozen_ratio_oracle_table()
    second = oracle.generate_frozen_ratio_oracle_table()
    assert first == second
    assert first["classification"] == "EXPLORATORY"
    assert first["hypothesis_status"] == "UNRESOLVED"
    assert first["generic_tail_upper_bound"] == {
        "numerator": "1",
        "denominator": "4096",
    }
    assert "about 12 bits" in first["resolution_scope"]
    assert len(first["entries"]) == 14
    assert first["entries"][5]["ratio"]["supplied"] == {
        "numerator": "4",
        "denominator": "6",
    }
    assert first["entries"][5]["ratio"]["reduced"] == {
        "numerator": "2",
        "denominator": "3",
    }


def test_frozen_audit_accepts_independent_closed_formula_at_tail_resolution() -> None:
    artifact = oracle.audit_frozen_ratio_normalization(
        _closed_autocorrelation_for_test,
        evaluator_label="test-only reduced Vasyunin formula",
    )
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["audit_outcome"] == (
        "ALL_14_CONSISTENT_WITH_TRUNCATION_ORACLE"
    )
    assert len(artifact["comparisons"]) == 14
    assert {item["outcome"] for item in artifact["comparisons"]} == {
        "CONSISTENT_WITH_TRUNCATION_ORACLE"
    }
    assert artifact["comparisons"][5]["ratio"]["supplied_was_reduced"] is False


@pytest.mark.parametrize(
    "mutation",
    (
        "reversed_v_sign",
        "outer_p",
        "omit_forward_v",
        "omit_reciprocal_v",
        "reversed_log_term",
    ),
)
def test_frozen_audit_rejects_closed_formula_normalization_mutants(
    mutation: str,
) -> None:
    artifact = oracle.audit_frozen_ratio_normalization(
        lambda ratio: _mutated_closed_autocorrelation_for_test(ratio, mutation),
        evaluator_label=f"test-only {mutation} mutant",
    )
    assert artifact["audit_outcome"] == (
        "AT_LEAST_ONE_RATIO_DISJOINT_FROM_TRUNCATION_ORACLE"
    )
    assert any(
        comparison["outcome"] == "DISJOINT_FROM_TRUNCATION_ORACLE"
        for comparison in artifact["comparisons"]
    )


@pytest.mark.parametrize(
    ("ratio", "truncation", "precision", "exception"),
    [
        ((0, 1), 8, 96, ValueError),
        ((1, 0), 8, 96, ValueError),
        ((1, 1), 0, 96, ValueError),
        ((1, 1), 8, 63, ValueError),
        ((True, 1), 8, 96, TypeError),
    ],
)
def test_invalid_oracle_inputs_are_rejected(
    ratio: oracle.RatioInput,
    truncation: int,
    precision: int,
    exception: type[Exception],
) -> None:
    with pytest.raises(exception):
        oracle.truncated_inner_product_oracle(ratio, truncation, precision)
