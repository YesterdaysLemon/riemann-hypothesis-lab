from __future__ import annotations

import copy

import pytest
from flint import arb

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_from_dyadic, arb_record
from riemann_lab.lagarias import (
    LagariasCertificateError,
    certify_lagarias_range,
    verify_lagarias_certificate,
)
from riemann_lab.zeros import (
    ZeroCertificateError,
    certify_critical_line_zeros,
    verify_zero_certificate,
)


def test_lagarias_small_range_is_rigorously_positive() -> None:
    result = certify_lagarias_range(250, precision_bits=128)
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["checks"] == {
        "equality_at_one": True,
        "inconclusive_n": [],
        "violating_n": [],
    }
    assert "does not prove RH" in result["limitation"]


def test_first_ten_zeros_and_count_are_certified() -> None:
    result = certify_critical_line_zeros(
        first_index=1,
        count=10,
        precision_bits=128,
        count_block_size=5,
    )
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["checks"]["hardy_z_roots_embedded_on_critical_line"] is True
    assert result["checks"]["complete_prefix_by_turing_counts"] is True
    assert [block["total_zero_count"] for block in result["checks"]["blocks"]] == [
        "5",
        "10",
    ]
    assert "cannot prove" in result["limitation"]
    replay = verify_zero_certificate(result, replay_precision_bits=256)
    assert replay["classification"] == "REPRODUCED"
    assert replay["reproduced_block_counts"] == ["5", "10"]
    assert replay["reproduced_zero_count"] == "10"


def test_lagarias_certificate_replays_at_higher_precision() -> None:
    result = certify_lagarias_range(250, precision_bits=96)
    replay = verify_lagarias_certificate(result)
    assert replay["classification"] == "REPRODUCED"
    assert replay["replay_precision_bits"] == "192"


def test_singleton_lagarias_certificate_replays() -> None:
    result = certify_lagarias_range(1, precision_bits=96)
    assert result["smallest_recorded_margin_lower_bound"] is None
    assert verify_lagarias_certificate(result)["classification"] == "REPRODUCED"


def test_lagarias_replay_rejects_semantic_forgery_with_valid_hash() -> None:
    result = certify_lagarias_range(25, precision_bits=96)
    forged = copy.deepcopy(result)
    forged["claim"] = "RH proved"
    unhashed = {key: value for key, value in forged.items() if key != "payload_sha256"}
    forged["payload_sha256"] = content_sha256(unhashed)
    with pytest.raises(LagariasCertificateError, match="canonical regeneration"):
        verify_lagarias_certificate(forged)


def test_zero_replay_rejects_reordered_transcript_even_with_rehashed_payload() -> None:
    result = certify_critical_line_zeros(count=10, precision_bits=128)
    forged = copy.deepcopy(result)
    forged["zeros"][0], forged["zeros"][1] = forged["zeros"][1], forged["zeros"][0]
    forged["zeros_sha256"] = content_sha256(forged["zeros"])
    unhashed = {key: value for key, value in forged.items() if key != "payload_sha256"}
    forged["payload_sha256"] = content_sha256(unhashed)
    with pytest.raises(ZeroCertificateError, match="indices are not consecutive"):
        verify_zero_certificate(forged, replay_precision_bits=256)


def test_zero_replay_rejects_forged_interior_ball_with_valid_hashes() -> None:
    result = certify_critical_line_zeros(count=10, precision_bits=128)
    forged = copy.deepcopy(result)
    forged["zeros"][3]["imaginary"] = arb_record(arb(29))
    forged["zeros_sha256"] = content_sha256(forged["zeros"])
    unhashed = {key: value for key, value in forged.items() if key != "payload_sha256"}
    forged["payload_sha256"] = content_sha256(unhashed)
    with pytest.raises(ZeroCertificateError, match="does not contain replayed root"):
        verify_zero_certificate(forged, replay_precision_bits=256)


def test_zero_replay_rejects_deleted_count_schedule() -> None:
    result = certify_critical_line_zeros(count=10, precision_bits=128)
    forged = copy.deepcopy(result)
    forged["checks"]["blocks"] = []
    unhashed = {key: value for key, value in forged.items() if key != "payload_sha256"}
    forged["payload_sha256"] = content_sha256(unhashed)
    with pytest.raises(ZeroCertificateError, match="no count blocks"):
        verify_zero_certificate(forged, replay_precision_bits=256)


def test_dyadic_reconstruction_preserves_enclosure() -> None:
    original = arb("1.23456789 +/- 1e-20")
    reconstructed = arb_from_dyadic(arb_record(original)["dyadic"])
    assert reconstructed.contains(original)


def test_zero_certificate_rejects_bad_parameters() -> None:
    for kwargs in ({"count": 0}, {"first_index": 0}, {"precision_bits": 32}):
        try:
            certify_critical_line_zeros(**kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid parameters accepted: {kwargs}")
