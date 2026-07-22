from __future__ import annotations

import copy

import pytest
from flint import acb, arb, ctx

from riemann_lab.artifacts import content_sha256
from riemann_lab.weil import (
    CONTROL_WITNESS,
    WeilCertificateError,
    archimedean_entry,
    archimedean_entry_integral,
    certify_positive_ldlt,
    certify_weil_matrix,
    pole_entry,
    pole_entry_integral,
    prime_power_entry,
    prime_powers_leq,
    verify_weil_certificate,
    weil_kernel,
)


def test_prime_power_cutoff_uses_exact_rational_boundaries() -> None:
    assert prime_powers_leq(399, 100) == [(2, 1, 2), (3, 1, 3)]
    assert prime_powers_leq(4) == [(2, 1, 2), (3, 1, 3), (2, 2, 4)]
    assert prime_powers_leq(8) == [
        (2, 1, 2),
        (3, 1, 3),
        (2, 2, 4),
        (5, 1, 5),
        (7, 1, 7),
        (2, 3, 8),
    ]
    assert prime_powers_leq(9)[-1] == (3, 2, 9)


def test_endpoint_prime_power_is_retained_but_contributes_exact_zero() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        length = arb(4).log()
        powers = prime_powers_leq(4)
        value = prime_power_entry(2, 2, length, powers, 4, 1)
        endpoint_only = prime_power_entry(2, 2, length, [(2, 2, 4)], 4, 1)
        without_endpoint = prime_power_entry(
            2,
            2,
            length,
            [record for record in powers if record[2] != 4],
            4,
            1,
        )
        assert endpoint_only.is_zero()
        assert value.overlaps(without_endpoint)
    finally:
        ctx.prec = previous_precision


def test_kernel_symmetries_and_endpoints() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 192
        length = (arb(5) / 2).log()
        zero = arb(0)
        endpoint = length
        for m, n in [(-4, -4), (-4, 3), (-1, 2), (0, 0), (3, 4)]:
            assert weil_kernel(m, n, zero, length).contains(2 if m == n else 0)
            assert weil_kernel(m, n, endpoint, length).contains(0)
            assert weil_kernel(m, n, length / 3, length).overlaps(
                weil_kernel(n, m, length / 3, length)
            )
            assert weil_kernel(m, n, length / 3, length).overlaps(
                weil_kernel(-m, -n, length / 3, length)
            )
    finally:
        ctx.prec = previous_precision


def test_kernel_matches_direct_basis_correlation() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        imaginary_unit = acb(0, 1)
        pi = acb.pi()
        t = arb(1) / 3
        length = (arb(5) / 2).log()
        for m, n in [(-4, 3), (-1, 2), (0, 0), (3, 4)]:
            correlation = 2 * acb.integral(
                lambda s, analytic: (
                    (-2 * pi * imaginary_unit * m * s).exp()
                    * (2 * pi * imaginary_unit * n * (s + t)).exp()
                ),
                0,
                1 - t,
            ).real
            assert weil_kernel(m, n, length * t, length).overlaps(correlation)
    finally:
        ctx.prec = previous_precision


def test_closed_forms_contain_independently_frozen_reference_values() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 192
        cutoff = arb(5) / 2
        length = cutoff.log()
        assert pole_entry(0, 0, length).overlaps(
            arb("1.86486105064915965618926101350174929409 +/- 1e-39")
        )
        assert pole_entry(-4, 3, length).overlaps(
            arb("0.000825854793887869019291101884468264249 +/- 1e-39")
        )
        assert archimedean_entry(0, 0, cutoff, length).overlaps(
            arb("1.54708302433451456561402240769105220519 +/- 1e-35")
        )
        assert archimedean_entry(-4, 3, cutoff, length).overlaps(
            arb("-0.0694927163399683319143157058240112224872 +/- 1e-40")
        )
        for m, n in [(-4, 3), (0, 0), (2, 4)]:
            assert pole_entry(m, n, length).overlaps(
                pole_entry_integral(m, n, length)
            )
            assert archimedean_entry(m, n, cutoff, length).overlaps(
                archimedean_entry_integral(m, n, cutoff, length)
            )
    finally:
        ctx.prec = previous_precision


def test_prime_term_has_exact_weight_and_no_extra_factor_two() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 192
        length = (arb(5) / 2).log()
        powers = prime_powers_leq(5, 2)
        log_two = arb(2).log()
        for m, n in [(-4, 3), (0, 0), (2, 4)]:
            expected = (
                log_two
                / arb(2).sqrt()
                * weil_kernel(m, n, log_two, length)
            )
            actual = prime_power_entry(m, n, length, powers, 5, 2)
            assert actual.overlaps(expected)
    finally:
        ctx.prec = previous_precision


def test_frozen_matrix_and_mutation_control_are_certified() -> None:
    artifact = certify_weil_matrix(precision_bits=128)
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["scope"]["dimension"] == "9"
    assert [
        (entry["prime"], entry["exponent"], entry["power"])
        for entry in artifact["prime_power_transcript"]
    ] == [("2", "1", "2")]
    assert artifact["positive_definiteness"]["classification"] == "POSITIVE_DEFINITE"
    assert len(artifact["positive_definiteness"]["pivots"]) == 9
    assert artifact["secondary_rump_spectrum"]["classification"] == "POSITIVE_SPECTRUM"
    assert artifact["secondary_rump_spectrum"]["smallest_eigenvalue"][
        "separated_from_all_others"
    ]
    assert artifact["checks"]["archimedean_closed_form_matches_integral_oracle"]
    assert artifact["checks"]["pole_closed_form_matches_integral_oracle"]
    assert artifact["negative_control"]["classification"] == "CONTROL_NEGATIVE"
    assert artifact["negative_control"]["witness"] == [
        str(value) for value in CONTROL_WITNESS
    ]
    assert artifact["negative_control"]["norm_squared"] == "2114603971100"
    assert artifact["negative_control"]["upper_bound_is_negative"] is True


def test_weil_certificate_replays_at_higher_precision() -> None:
    artifact = certify_weil_matrix(precision_bits=128)
    replay = verify_weil_certificate(artifact, replay_precision_bits=256)
    assert replay["classification"] == "REPRODUCED"
    assert replay["dimension"] == "9"
    assert replay["prime_power_count"] == "1"
    with pytest.raises(WeilCertificateError, match="replay_precision_bits"):
        verify_weil_certificate(artifact, replay_precision_bits=64)


def test_weil_replay_rejects_semantic_forgery_with_valid_hash() -> None:
    artifact = certify_weil_matrix(precision_bits=128)
    forged = copy.deepcopy(artifact)
    forged["normalization"]["matrix"] = "A=P-R"
    payload = {key: value for key, value in forged.items() if key != "payload_sha256"}
    forged["payload_sha256"] = content_sha256(payload)
    with pytest.raises(WeilCertificateError, match="canonical regeneration"):
        verify_weil_certificate(forged, replay_precision_bits=256)


def test_ldlt_rejects_merely_overlapping_asymmetric_entries() -> None:
    matrix = [
        [arb(3), arb("1 +/- 0.1")],
        [arb("1.05 +/- 0.1"), arb(3)],
    ]
    with pytest.raises(ValueError, match="identical symmetric"):
        certify_positive_ldlt(matrix)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cutoff_numerator": 1, "cutoff_denominator": 1},
        {"cutoff_numerator": 5, "cutoff_denominator": 0},
        {"degree": -1},
        {"precision_bits": 64},
    ],
)
def test_weil_certificate_rejects_bad_parameters(kwargs: dict[str, int]) -> None:
    with pytest.raises((TypeError, ValueError)):
        certify_weil_matrix(**kwargs)
