from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil

import pytest
from flint import arb, ctx

import tools.certify_nyman_balanced_tail as tail_core
from tools.certify_nyman_balanced_tail import (
    _canonical_sha256,
    _ideal_shell_arb,
    absolute_moments,
    alias_step,
    build_exact_vectors,
    centered_fractional_part,
    coefficient_sum,
    convolve_coefficients,
    harmonic_sum,
    jordan_gcd_means,
    lcm_absolute_cross,
    lcm_absolute_self,
    old_periodic_center,
    pairwise_lcm_absolute,
    phi_covariance,
    require_balanced,
    rounded_balanced_multiplier,
    tail_center_radius,
)


ROOT = Path(__file__).resolve().parents[1]


def _balanced_example() -> dict[int, Fraction]:
    shell = {2: Fraction(3, 4), 3: Fraction(-3, 4)}
    multiplier = {1: Fraction(1), 2: Fraction(-2)}
    result = convolve_coefficients(shell, multiplier)
    require_balanced(result)
    return result


@pytest.mark.parametrize(
    ("left", "right"),
    [(1, 1), (2, 3), (4, 6), (5, 15), (8, 12)],
)
def test_centered_fractional_part_covariance_is_exact(
    left: int,
    right: int,
) -> None:
    period = math.lcm(left, right)
    brute = sum(
        (
            centered_fractional_part(interval, left)
            * centered_fractional_part(interval, right)
            for interval in range(1, period + 1)
        ),
        start=Fraction(0),
    ) / period
    assert brute == phi_covariance(left, right)


def test_jordan_means_match_one_complete_period() -> None:
    added = _balanced_example()
    old = {
        1: Fraction(2, 5),
        2: Fraction(-1, 3),
        5: Fraction(7, 11),
    }
    period = math.lcm(*old, *added)
    mu, nu = jordan_gcd_means(old, added)
    brute_mu = sum(
        (alias_step(added, interval) ** 2 for interval in range(1, period + 1)),
        start=Fraction(0),
    ) / period
    brute_nu = sum(
        (
            alias_step(added, interval)
            * old_periodic_center(old, interval)
            for interval in range(1, period + 1)
        ),
        start=Fraction(0),
    ) / period
    assert mu == brute_mu
    assert nu == brute_nu


def test_lcm_divisor_aggregation_matches_pairwise_baseline() -> None:
    left = {
        2: Fraction(-3, 8),
        4: Fraction(5, 16),
        9: Fraction(-7, 32),
    }
    right = {
        3: Fraction(11, 16),
        6: Fraction(-1, 4),
        10: Fraction(9, 32),
    }
    assert lcm_absolute_self(left) == pairwise_lcm_absolute(left, left)
    assert lcm_absolute_cross(left, right) == pairwise_lcm_absolute(
        left,
        right,
    )


def test_exact_builder_rounds_then_enforces_both_balances() -> None:
    vectors = build_exact_vectors(
        ROOT,
        n=256,
        multiplier_limit=16,
        shell_bits=10,
        z_bits=10,
        old_bits=16,
    )
    assert vectors.n == 256
    assert vectors.multiplier_limit == 16
    assert vectors.old_bits == 16
    assert coefficient_sum(vectors.shell) == 0
    assert vectors.multiplier[1] == 1
    assert harmonic_sum(vectors.multiplier) == 0
    assert coefficient_sum(vectors.added_coefficients) == 0
    assert harmonic_sum(vectors.added_coefficients) == 0
    assert max(vectors.added_coefficients) <= 2 * 256 * 16
    assert all(
        value.denominator <= 1 << 16
        for value in vectors.old_coefficients.values()
    )


def test_every_frozen_shell_rounding_bin_is_proved_by_arb() -> None:
    bits = 9
    vectors = build_exact_vectors(
        ROOT,
        n=256,
        multiplier_limit=16,
        shell_bits=bits,
        z_bits=bits,
        old_bits=16,
    )
    previous_precision = ctx.prec
    ctx.prec = 384
    try:
        ideal = _ideal_shell_arb(vectors.source_candidate)
        half = arb(1) / 2
        for offset, value in enumerate(ideal[:-1], start=257):
            rounded = vectors.shell.get(offset, Fraction(0))
            rounded_integer = rounded * (1 << bits)
            assert rounded_integer.denominator == 1
            difference = value * (1 << bits) - rounded_integer.numerator
            assert difference.lower() > -half
            assert difference.upper() < half
    finally:
        ctx.prec = previous_precision


def _copy_builder_inputs(tmp_path: Path) -> tuple[Path, Path]:
    scout_relative = Path("results/nyman-balanced-multiplier-scout-v1.json")
    candidate_relative = Path("results/nyman-natural-v1/cells/n-0256.json")
    for relative in (scout_relative, candidate_relative):
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, destination)
    return tmp_path / scout_relative, tmp_path / candidate_relative


def _selected_cell(scout: dict[str, object]) -> dict[str, object]:
    cells = scout["n256_width_sweep"]["cells"]
    return next(
        cell
        for cell in cells
        if int(cell["n"]) == 256
        and int(cell["multiplier_limit"]) == 16
        and cell["objective"] == "direct-gain"
    )


def _rehash(record: dict[str, object]) -> None:
    body = {key: value for key, value in record.items() if key != "payload_sha256"}
    record["payload_sha256"] = _canonical_sha256(body)


def test_builder_rejects_a_mutated_scout_cell(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scout_path, _ = _copy_builder_inputs(tmp_path)
    scout = json.loads(scout_path.read_text(encoding="utf-8"))
    cell = _selected_cell(scout)
    cell["coefficients"][3] += 1e-9
    _rehash(scout)
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_PAYLOAD_SHA256",
        scout["payload_sha256"],
    )
    scout_path.write_text(json.dumps(scout), encoding="utf-8")
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_RAW_SHA256",
        hashlib.sha256(scout_path.read_bytes()).hexdigest(),
    )
    with pytest.raises(ValueError, match="scout cell canonical payload"):
        build_exact_vectors(tmp_path)


def test_builder_rejects_decimal_token_mutation_hidden_by_binary64(
    tmp_path: Path,
) -> None:
    scout_path, _ = _copy_builder_inputs(tmp_path)
    raw = scout_path.read_text(encoding="utf-8")
    changed = raw.replace(
        "-1.2530062330140541",
        "-1.25300623301405409",
        1,
    )
    assert changed != raw
    scout_path.write_text(changed, encoding="utf-8")
    with pytest.raises(ValueError, match="scout raw SHA-256 mismatch"):
        build_exact_vectors(tmp_path)


def test_builder_rejects_a_mutated_candidate_payload(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scout_path, candidate_path = _copy_builder_inputs(tmp_path)
    candidate = json.loads(candidate_path.read_text(encoding="utf-8"))
    candidate["candidate"]["coefficients"]["numerators"][0] = "0"
    _rehash(candidate)
    candidate_path.write_text(json.dumps(candidate), encoding="utf-8")

    scout = json.loads(scout_path.read_text(encoding="utf-8"))
    cell = _selected_cell(scout)
    cell["source_binding"]["raw_sha256"] = hashlib.sha256(
        candidate_path.read_bytes()
    ).hexdigest()
    _rehash(cell)
    _rehash(scout)
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_PAYLOAD_SHA256",
        scout["payload_sha256"],
    )
    scout_path.write_text(json.dumps(scout), encoding="utf-8")
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_RAW_SHA256",
        hashlib.sha256(scout_path.read_bytes()).hexdigest(),
    )
    with pytest.raises(ValueError, match="candidate payload SHA-256 mismatch"):
        build_exact_vectors(tmp_path)


def test_builder_rejects_wrong_multiplier_coefficient_count(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scout_path, _ = _copy_builder_inputs(tmp_path)
    scout = json.loads(scout_path.read_text(encoding="utf-8"))
    cell = _selected_cell(scout)
    cell["coefficients"].pop()
    _rehash(cell)
    _rehash(scout)
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_PAYLOAD_SHA256",
        scout["payload_sha256"],
    )
    scout_path.write_text(json.dumps(scout), encoding="utf-8")
    monkeypatch.setattr(
        tail_core,
        "FROZEN_SCOUT_RAW_SHA256",
        hashlib.sha256(scout_path.read_bytes()).hexdigest(),
    )
    with pytest.raises(ValueError, match="coefficient count mismatch"):
        build_exact_vectors(tmp_path)


def test_balanced_multiplier_rounds_z_not_c_and_replaces_c2() -> None:
    supplied = [
        Decimal("1"),
        Decimal("-1.234567"),
        Decimal("0.333333"),
        Decimal("-0.271828"),
    ]
    multiplier = rounded_balanced_multiplier(supplied, bits=7)
    assert multiplier[3] / 3 == Fraction(
        round(Fraction(Decimal("0.333333")) / 3 * (1 << 7)),
        1 << 7,
    )
    assert multiplier[4] / 4 == Fraction(
        round(Fraction(Decimal("-0.271828")) / 4 * (1 << 7)),
        1 << 7,
    )
    assert harmonic_sum(multiplier) == 0
    assert multiplier[2] != Fraction(Decimal("-1.234567"))


def test_tail_center_radius_contains_brute_force_tail() -> None:
    added = _balanced_example()
    old = {
        1: Fraction(2, 5),
        2: Fraction(-1, 3),
        5: Fraction(7, 11),
    }
    cutoff = 19
    last = 200_000
    bound = tail_center_radius(old, added, cutoff)
    assert bound.radius == bound.periodic_radius + bound.delta_radius
    assert bound.lower == bound.center - bound.radius
    assert bound.upper == bound.center + bound.radius

    p1 = harmonic_sum(old)
    with localcontext() as context:
        context.prec = 50
        total = Decimal(0)
        maximum_periodic_term = Fraction(0)
        for interval in range(cutoff + 1, last + 1):
            d_value = alias_step(added, interval)
            r_bar = old_periodic_center(old, interval)
            periodic_term = 2 * d_value * r_bar - d_value * d_value
            maximum_periodic_term = max(
                maximum_periodic_term,
                abs(periodic_term),
            )
            weight = Fraction(1, interval * (interval + 1))
            main = periodic_term * weight
            m = Decimal(interval)
            delta = (
                (Decimal(1) + Decimal(1) / m).ln()
                - Decimal(1) / (m + 1)
                - Decimal(1) / (2 * m * (m + 1))
            )
            total += (
                Decimal(main.numerator) / Decimal(main.denominator)
                - Decimal(2)
                * Decimal(p1.numerator)
                / Decimal(p1.denominator)
                * Decimal(d_value.numerator)
                / Decimal(d_value.denominator)
                * delta
            )

        # The uncomputed periodic remainder is at most sup|h|/(last+1).
        # The same delta estimate used by the certificate, restarted at
        # ``last``, bounds the remaining logarithmic correction.
        a0, _ = absolute_moments(added)
        remainder = (
            Decimal(maximum_periodic_term.numerator)
            / Decimal(maximum_periodic_term.denominator)
            / Decimal(last + 1)
            + Decimal(abs(p1).numerator)
            / Decimal(abs(p1).denominator)
            * Decimal(a0.numerator)
            / Decimal(a0.denominator)
            / Decimal(12 * last * last)
        )
        lower = Decimal(bound.lower.numerator) / Decimal(
            bound.lower.denominator
        )
        upper = Decimal(bound.upper.numerator) / Decimal(
            bound.upper.denominator
        )
        assert total - remainder >= lower
        assert total + remainder <= upper


@pytest.mark.parametrize("interval", [1, 2, 3, 7, 31, 100])
def test_logarithmic_delta_bound(interval: int) -> None:
    with localcontext() as context:
        context.prec = 70
        m = Decimal(interval)
        delta = (
            (Decimal(1) + Decimal(1) / m).ln()
            - Decimal(1) / (m + 1)
            - Decimal(1) / (2 * m * (m + 1))
        )
        upper = Decimal(1) / (6 * m * m * m)
        assert delta <= 0
        assert abs(delta) <= upper
