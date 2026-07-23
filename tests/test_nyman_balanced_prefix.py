from __future__ import annotations

from fractions import Fraction
from pathlib import Path

from flint import arb, ctx
import pytest

from tools.certify_nyman_balanced_prefix import (
    _arb_from_fraction,
    _radius_ball,
    certified_prefix_gain,
    common_dyadic_numerators,
    compressed_logarithmic_sum,
)
from tools.certify_nyman_balanced_tail import (
    build_exact_vectors,
    convolve_coefficients,
    tail_center_radius,
)


ROOT = Path(__file__).resolve().parents[1]


def _small_vectors() -> tuple[dict[int, Fraction], dict[int, Fraction]]:
    old = {
        1: Fraction(1, 2),
        2: Fraction(-1, 4),
        5: Fraction(3, 4),
    }
    shell = {2: Fraction(3, 4), 3: Fraction(-3, 4)}
    multiplier = {1: Fraction(1), 2: Fraction(-2)}
    return old, convolve_coefficients(shell, multiplier)


def _exact_rational_prefix(
    old: dict[int, Fraction],
    added: dict[int, Fraction],
    cutoff: int,
) -> Fraction:
    total = Fraction(0)
    for interval in range(1, cutoff + 1):
        step = -sum(
            value * (interval // index)
            for index, value in added.items()
        )
        intercept = 1 + sum(
            value * (interval // index)
            for index, value in old.items()
        )
        total += (2 * step * intercept - step * step) / (
            interval * (interval + 1)
        )
    return total


def test_binary64_rational_radius_contains_exact_small_prefix() -> None:
    old, added = _small_vectors()
    cutoff = 200
    bound = certified_prefix_gain(
        old,
        added,
        cutoff,
        block_size=37,
        precision_bits=192,
    )
    exact = _exact_rational_prefix(old, added, cutoff)
    assert bound.rational_lower <= exact <= bound.rational_upper
    assert bound.last_added_step_numerator == -sum(
        numerator * (cutoff // index)
        for index, numerator in common_dyadic_numerators(added)[1].items()
    )
    old_denominator, old_numerators = common_dyadic_numerators(old)
    assert bound.last_old_intercept_numerator == old_denominator + sum(
        numerator * (cutoff // index)
        for index, numerator in old_numerators.items()
    )


def test_compressed_logarithmic_identity_overlaps_direct_arb_sum() -> None:
    _, added = _small_vectors()
    denominator, numerators = common_dyadic_numerators(added)
    cutoff = 137
    previous_precision = ctx.prec
    ctx.prec = 256
    try:
        compressed = compressed_logarithmic_sum(
            numerators,
            denominator,
            cutoff,
        )
        direct = arb(0)
        for interval in range(1, cutoff + 1):
            step = -sum(
                value * (interval // index)
                for index, value in added.items()
            )
            direct += _arb_from_fraction(step) * (
                arb(interval + 1) / interval
            ).log()
        assert compressed.overlaps(direct)
    finally:
        ctx.prec = previous_precision


def test_binary64_path_rejects_a_grid_too_fine_for_underflow_proof() -> None:
    old, added = _small_vectors()
    tiny_added = {
        index: value / (1 << 500)
        for index, value in added.items()
    }
    with pytest.raises(ValueError, match="grid is too fine"):
        certified_prefix_gain(old, tiny_added, 1)


@pytest.fixture(scope="module")
def frozen_vectors():
    return build_exact_vectors(
        ROOT,
        n=256,
        multiplier_limit=16,
        shell_bits=9,
        z_bits=9,
        old_bits=16,
    )


def test_frozen_full_tail_gain_has_a_strictly_positive_lower_endpoint(
    frozen_vectors,
) -> None:
    cutoff = 1 << 26
    first = certified_prefix_gain(
        frozen_vectors.old_coefficients,
        frozen_vectors.added_coefficients,
        cutoff,
        block_size=1 << 20,
        precision_bits=192,
    )
    replay = certified_prefix_gain(
        frozen_vectors.old_coefficients,
        frozen_vectors.added_coefficients,
        cutoff,
        block_size=1 << 19,
        precision_bits=384,
    )
    tail = tail_center_radius(
        frozen_vectors.old_coefficients,
        frozen_vectors.added_coefficients,
        cutoff,
    )

    assert first.added_denominator == 1 << 18
    assert first.old_denominator == 1 << 16
    assert first.maximum_interval_numerator == 421_449_858_668_025_872
    assert first.last_added_step_numerator == -1_086_805
    assert first.last_old_intercept_numerator == 12_321_937_299
    assert first.prefix_gain.overlaps(replay.prefix_gain)
    assert max(first.rational_lower, replay.rational_lower) <= min(
        first.rational_upper,
        replay.rational_upper,
    )

    full_gain = (
        replay.prefix_gain
        + _arb_from_fraction(tail.center)
        + _radius_ball(tail.radius, 384)
    )
    assert full_gain.lower() > arb("0.000142")
    assert full_gain.upper() < arb("0.000368")
