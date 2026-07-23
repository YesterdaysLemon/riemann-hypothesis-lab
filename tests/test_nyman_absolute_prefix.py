from __future__ import annotations

from fractions import Fraction
import math

from flint import arb, ctx
import numpy as np
import pytest

from tools import certify_nyman_absolute_prefix as prefix_v2
from tools import generate_nyman_rebased_schur_certificate as legacy


ZERO = Fraction()


def _exact_rational_sum(
    coefficients: dict[int, Fraction],
    cutoff: int,
) -> Fraction:
    return sum(
        (
            (
                1
                + sum(
                    coefficient * (interval // index)
                    for index, coefficient in coefficients.items()
                )
            )
            ** 2
            / (interval * (interval + 1))
            for interval in range(1, cutoff + 1)
        ),
        start=ZERO,
    )


def _direct_arb_prefix(
    coefficients: dict[int, Fraction],
    cutoff: int,
) -> object:
    slope = sum(
        (
            coefficient / index
            for index, coefficient in coefficients.items()
        ),
        start=ZERO,
    )
    total = prefix_v2._arb_from_fraction((cutoff + 1) * slope * slope)
    for interval in range(1, cutoff + 1):
        q_value = 1 + sum(
            coefficient * (interval // index)
            for index, coefficient in coefficients.items()
        )
        total += (
            prefix_v2._arb_from_fraction(q_value * q_value)
            / interval
            / (interval + 1)
            - 2
            * prefix_v2._arb_from_fraction(slope * q_value)
            * (arb(interval + 1) / interval).log()
        )
    return total


@pytest.mark.parametrize(
    ("coefficients", "cutoff", "chunk_size", "leaf_size"),
    [
        ({}, 19, 7, 3),
        ({1: Fraction(1)}, 16, 7, 3),
        (
            {
                1: Fraction(3, 2),
                2: Fraction(-5, 4),
                5: Fraction(7, 8),
                29: Fraction(-9, 16),
            },
            23,
            7,
            3,
        ),
        (
            {
                2: Fraction(-7, 8),
                3: Fraction(5, 16),
                11: Fraction(13, 32),
            },
            31,
            11,
            5,
        ),
    ],
)
def test_exact_small_prefix_is_contained_for_odd_chunk_and_leaf_layouts(
    coefficients: dict[int, Fraction],
    cutoff: int,
    chunk_size: int,
    leaf_size: int,
) -> None:
    bound = prefix_v2.certified_absolute_energy_prefix(
        coefficients,
        cutoff,
        chunk_size=chunk_size,
        reduction_leaf_size=leaf_size,
        precision_bits=256,
    )
    exact_rational = _exact_rational_sum(coefficients, cutoff)

    assert bound.rational_lower <= exact_rational <= bound.rational_upper
    assert (
        bound.rational_error
        == bound.term_evaluation_error + bound.reduction_error
    )
    assert bound.rational_lower == (
        bound.rational_midpoint - bound.rational_error
    )
    assert bound.rational_upper == (
        bound.rational_midpoint + bound.rational_error
    )

    previous_precision = ctx.prec
    ctx.prec = 512
    try:
        direct = _direct_arb_prefix(coefficients, cutoff)
        assert bound.prefix_energy.contains(direct)
    finally:
        ctx.prec = previous_precision


def test_odd_layout_metadata_describes_the_complete_declared_tree() -> None:
    bound = prefix_v2.certified_absolute_energy_prefix(
        {1: Fraction(3, 2), 2: Fraction(-5, 4)},
        20,
        chunk_size=7,
        reduction_leaf_size=3,
        precision_bits=192,
    )

    # Chunk lengths are 7, 7, 6.  The first two yield three leaves and the
    # last yields two; the final leaves of the first two chunks are short.
    assert bound.chunk_count == 3
    assert bound.reduction_leaf_count == 8
    assert bound.maximum_leaf_length == 3
    assert bound.maximum_leaf_reduction_depth == 2
    assert bound.root_reduction_depth == 3
    assert bound.reduction_depth_bound == 5
    assert bound.reduction_gamma == prefix_v2._gamma(5)


@pytest.mark.parametrize("length", [1, 2, 3, 5, 7, 9, 17])
def test_pairwise_reduction_depth_and_odd_carry_rule(length: int) -> None:
    values = np.arange(1, length + 1, dtype=np.float64)
    reduced = prefix_v2.explicit_pairwise_sum(values)

    assert reduced.value == float(length * (length + 1) // 2)
    assert reduced.input_count == length
    assert reduced.depth == math.ceil(math.log2(length))


def test_pairwise_reduction_handles_zeros_and_disparate_magnitudes() -> None:
    zeros = prefix_v2.explicit_pairwise_sum(
        np.asarray([0.0, -0.0, 0.0, 0.0, -0.0])
    )
    assert zeros.value == 0.0
    assert zeros.depth == 3

    # This also pins the declared adjacent-pair tree.  Its rounded result is
    # 0.0 even though the exact sum of the four binary64 inputs is 2.
    disparate = np.asarray([1.0e16, 1.0, -1.0e16, 1.0])
    assert sum(
        (Fraction.from_float(float(value)) for value in disparate),
        start=ZERO,
    ) == 2
    reduced = prefix_v2.explicit_pairwise_sum(disparate)
    assert reduced.value == 0.0
    assert reduced.depth == 2


@pytest.mark.parametrize(
    "values",
    [
        np.asarray([], dtype=np.float64),
        np.zeros((2, 2), dtype=np.float64),
    ],
)
def test_pairwise_reduction_rejects_invalid_shapes(values: np.ndarray) -> None:
    with pytest.raises(ValueError):
        prefix_v2.explicit_pairwise_sum(values)


@pytest.mark.parametrize(
    "values",
    [
        np.asarray([1.0, np.inf]),
        np.asarray([np.nan]),
        np.asarray([np.finfo(np.float64).max] * 2),
    ],
)
def test_pairwise_reduction_rejects_nonfinite_input_or_result(
    values: np.ndarray,
) -> None:
    with np.errstate(over="ignore"), pytest.raises(
        ArithmeticError,
        match="nonfinite",
    ):
        prefix_v2.explicit_pairwise_sum(values)


def test_sequential_divisions_cross_the_former_product_exactness_boundary() -> None:
    old_last_accepted = (
        math.isqrt(1 + 4 * (1 << 53)) - 1
    ) // 2
    first = old_last_accepted
    q_values = np.asarray(
        [1, (1 << 53) + 1, -prefix_v2.INT64_MAX],
        dtype=np.int64,
    )
    assert first * (first + 1) < 1 << 53
    assert (first + 1) * (first + 2) >= 1 << 53
    with pytest.raises(
        ValueError,
        match="too large for exact binary64 denominators",
    ):
        legacy.absolute_energy_prefix(
            {},
            first + 1,
            block_size=1,
            precision_bits=64,
        )

    terms = prefix_v2._binary64_absolute_terms(
        q_values,
        first,
        denominator_exponent=17,
    )
    relative_bound = (
        (1 + prefix_v2.BINARY64_UNIT_ROUNDOFF) ** 5 - 1
    )
    for offset, (q_value, term) in enumerate(
        zip(q_values, terms, strict=True)
    ):
        interval = first + offset
        exact = Fraction(
            int(q_value) * int(q_value),
            (1 << 34) * interval * (interval + 1),
        )
        computed = Fraction.from_float(float(term))
        assert abs(computed - exact) <= relative_bound * exact


def test_term_kernel_accepts_the_largest_exact_successor_interval() -> None:
    first = (1 << 53) - 3
    q_values = np.asarray([0, 1, -7], dtype=np.int64)
    terms = prefix_v2._binary64_absolute_terms(q_values, first, 0)

    assert terms[0] == 0
    for offset in (1, 2):
        interval = first + offset
        exact = Fraction(
            int(q_values[offset]) ** 2,
            interval * (interval + 1),
        )
        computed = Fraction.from_float(float(terms[offset]))
        relative_bound = (
            (1 + prefix_v2.BINARY64_UNIT_ROUNDOFF) ** 5 - 1
        )
        assert abs(computed - exact) <= relative_bound * exact


@pytest.mark.parametrize(
    ("q_values", "first", "exponent", "message"),
    [
        (np.asarray([], dtype=np.int64), 1, 0, "nonempty"),
        (np.zeros((1, 1), dtype=np.int64), 1, 0, "one-dimensional"),
        (np.asarray([1], dtype=np.int32), 1, 0, "signed int64 dtype"),
        (np.asarray([1], dtype=np.uint64), 1, 0, "signed int64 dtype"),
        (np.asarray([1], dtype=np.int64), 0, 0, "must be positive"),
        (np.asarray([1], dtype=np.int64), 1.0, 0, "exact integer"),
        (np.asarray([1], dtype=np.int64), 1, -1, "nonnegative"),
        (np.asarray([1], dtype=np.int64), 1, 0.0, "exact integer"),
        (
            np.asarray([1, 2], dtype=np.int64),
            (1 << 53) - 1,
            0,
            "exceeds exact binary64 integer range",
        ),
    ],
)
def test_term_kernel_guards_exact_integer_domain(
    q_values: np.ndarray,
    first: object,
    exponent: object,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        prefix_v2._binary64_absolute_terms(q_values, first, exponent)


def test_interval_integer_guard_rejects_only_after_the_new_limit() -> None:
    with pytest.raises(
        ValueError,
        match="too large for exact binary64 interval integers",
    ):
        prefix_v2.certified_absolute_energy_prefix(
            {},
            1 << 53,
            chunk_size=1,
            reduction_leaf_size=1,
        )


@pytest.mark.parametrize(
    "invalid_index",
    [True, np.bool_(False), 1.0, "1", Fraction(3, 2)],
)
def test_coefficient_indices_must_be_exact_nonboolean_integers(
    invalid_index: object,
) -> None:
    with pytest.raises(ValueError, match="must be exact integers"):
        prefix_v2.certified_absolute_energy_prefix(
            {invalid_index: Fraction(1)},
            1,
            chunk_size=1,
            reduction_leaf_size=1,
        )


def test_numpy_integer_coefficient_index_is_accepted_exactly() -> None:
    bound = prefix_v2.certified_absolute_energy_prefix(
        {np.int64(1): Fraction(1)},
        3,
        chunk_size=2,
        reduction_leaf_size=1,
        precision_bits=128,
    )
    assert bound.last_q_numerator == 4


@pytest.mark.parametrize(
    ("argument", "value", "message"),
    [
        ("cutoff", True, "prefix cutoff must be an exact integer"),
        ("cutoff", 1.0, "prefix cutoff must be an exact integer"),
        ("chunk_size", np.bool_(True), "chunk size must be an exact integer"),
        ("chunk_size", "1", "chunk size must be an exact integer"),
        (
            "reduction_leaf_size",
            False,
            "reduction leaf size must be an exact integer",
        ),
        (
            "reduction_leaf_size",
            Fraction(1),
            "reduction leaf size must be an exact integer",
        ),
        ("precision_bits", True, "precision must be an exact integer"),
        ("precision_bits", 128.0, "precision must be an exact integer"),
    ],
)
def test_public_schedule_inputs_must_be_exact_nonboolean_integers(
    argument: str,
    value: object,
    message: str,
) -> None:
    arguments: dict[str, object] = {
        "cutoff": 1,
        "chunk_size": 1,
        "reduction_leaf_size": 1,
        "precision_bits": 128,
    }
    arguments[argument] = value
    cutoff = arguments.pop("cutoff")
    with pytest.raises(ValueError, match=message):
        prefix_v2.certified_absolute_energy_prefix(
            {},
            cutoff,
            **arguments,
        )


def test_public_schedule_accepts_numpy_integers_exactly() -> None:
    bound = prefix_v2.certified_absolute_energy_prefix(
        {},
        np.int64(3),
        chunk_size=np.int64(2),
        reduction_leaf_size=np.int64(1),
        precision_bits=np.int64(128),
    )
    assert bound.cutoff == 3
    assert bound.chunk_size == 2
    assert bound.reduction_leaf_size == 1


@pytest.mark.parametrize(
    ("coefficients", "message"),
    [
        ({1: Fraction(1, 3)}, "powers of two"),
        ({1: Fraction(1, 1 << 63)}, "denominator exceeds signed int64"),
        ({1: Fraction(1 << 63)}, "numerator exceeds signed int64"),
        ({1: Fraction(prefix_v2.INT64_MAX)}, "recurrence may overflow"),
    ],
)
def test_dyadic_and_int64_guards_fail_closed(
    coefficients: dict[int, Fraction],
    message: str,
) -> None:
    with pytest.raises((ValueError, OverflowError), match=message):
        prefix_v2.certified_absolute_energy_prefix(
            coefficients,
            1,
            chunk_size=1,
            reduction_leaf_size=1,
        )


@pytest.mark.parametrize(
    "keyword_arguments",
    [
        {"chunk_size": 0},
        {"reduction_leaf_size": 0},
        {"chunk_size": 2, "reduction_leaf_size": 3},
        {"precision_bits": 63},
    ],
)
def test_schedule_and_precision_guards_fail_closed(
    keyword_arguments: dict[str, int],
) -> None:
    with pytest.raises(ValueError):
        prefix_v2.certified_absolute_energy_prefix(
            {},
            1,
            **keyword_arguments,
        )


def test_large_q_square_is_python_exact_and_not_an_int64_product() -> None:
    coefficient = Fraction(4_000_000_000)
    coefficients = {1: coefficient}
    bound = prefix_v2.certified_absolute_energy_prefix(
        coefficients,
        1,
        chunk_size=1,
        reduction_leaf_size=1,
        precision_bits=256,
    )
    q_numerator = 4_000_000_001
    exact_rational = Fraction(q_numerator * q_numerator, 2)

    assert q_numerator < prefix_v2.INT64_MAX
    assert q_numerator * q_numerator > prefix_v2.INT64_MAX
    assert bound.rational_lower <= exact_rational <= bound.rational_upper
    assert bound.maximum_q_numerator == q_numerator
    assert bound.maximum_square_numerator == q_numerator * q_numerator


def test_endpoint_and_headroom_records_match_direct_exact_recurrence() -> None:
    coefficients = {
        1: Fraction(3, 2),
        2: Fraction(-5, 4),
        5: Fraction(7, 8),
        19: Fraction(-11, 16),
    }
    cutoff = 13
    bound = prefix_v2.certified_absolute_energy_prefix(
        coefficients,
        cutoff,
        chunk_size=5,
        reduction_leaf_size=2,
        precision_bits=192,
    )
    denominator = 16
    numerators = {1: 24, 2: -20, 5: 14, 19: -11}
    q_values = [
        denominator
        + sum(
            numerator * (interval // index)
            for index, numerator in numerators.items()
        )
        for interval in range(1, cutoff + 1)
    ]
    recurrence_bound = denominator + sum(
        abs(numerator) * (cutoff // index)
        for index, numerator in numerators.items()
    )

    assert bound.coefficient_denominator == denominator
    assert bound.coefficient_denominator_exponent == 4
    assert bound.last_q_numerator == q_values[-1]
    assert bound.maximum_q_numerator == max(
        denominator,
        *(abs(value) for value in q_values),
    )
    assert bound.maximum_square_numerator == (
        bound.maximum_q_numerator * bound.maximum_q_numerator
    )
    assert bound.recurrence_absolute_bound == recurrence_bound
    assert bound.int64_headroom == prefix_v2.INT64_MAX - recurrence_bound


def test_independent_generation_and_replay_layouts_overlap() -> None:
    coefficients = {
        1: Fraction(3, 2),
        2: Fraction(-5, 4),
        3: Fraction(11, 16),
        7: Fraction(-13, 32),
        17: Fraction(9, 64),
        263: Fraction(-5, 128),
    }
    cutoff = 257
    generation = prefix_v2.certified_absolute_energy_prefix(
        coefficients,
        cutoff,
        chunk_size=37,
        reduction_leaf_size=11,
        precision_bits=192,
    )
    replay = prefix_v2.certified_absolute_energy_prefix(
        coefficients,
        cutoff,
        chunk_size=29,
        reduction_leaf_size=7,
        precision_bits=384,
    )
    exact_rational = _exact_rational_sum(coefficients, cutoff)

    assert generation.algorithm_version == prefix_v2.PREFIX_ALGORITHM_VERSION
    assert (
        generation.term_evaluation_algorithm
        == prefix_v2.TERM_EVALUATION_ALGORITHM
    )
    assert generation.reduction_algorithm == prefix_v2.REDUCTION_ALGORITHM
    assert generation.rational_lower <= exact_rational
    assert exact_rational <= generation.rational_upper
    assert replay.rational_lower <= exact_rational <= replay.rational_upper
    assert max(generation.rational_lower, replay.rational_lower) <= min(
        generation.rational_upper,
        replay.rational_upper,
    )
    assert generation.prefix_energy.overlaps(replay.prefix_energy)
    assert generation.last_q_numerator == replay.last_q_numerator
    assert generation.maximum_q_numerator == replay.maximum_q_numerator
