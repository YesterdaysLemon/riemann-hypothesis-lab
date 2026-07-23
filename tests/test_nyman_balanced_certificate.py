from __future__ import annotations

import json
from pathlib import Path

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
import tools.generate_nyman_balanced_tail_certificate as certificate_tool
from tools.generate_nyman_balanced_tail_certificate import (
    CLAIMED_GAIN_LOWER_BOUND,
    SCHEMA,
    _strict_lower_proved,
    _wide_ball_record,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-balanced-full-tail-v1.json"


def test_checked_in_balanced_full_tail_certificate_replays() -> None:
    result = verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    assert result["claimed_gain_lower_bound"] == "7/50000"


def test_certificate_metadata_mutation_is_rejected(
    tmp_path: Path,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert artifact["schema"] == SCHEMA
    assert CLAIMED_GAIN_LOWER_BOUND.numerator == 7
    artifact["classification"] = "PROVED"
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)
    mutated = tmp_path / "mutated.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="not the frozen v1 payload"):
        verify(ROOT, mutated)


def test_semantic_gate_still_rejects_a_rehashed_false_classification(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["classification"] = "PROVED"
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "semantic-mutation.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="classification changed"):
        verify(ROOT, mutated)


def test_certificate_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    duplicated = raw.replace("{", '{"schema":"false",', 1)
    mutated = tmp_path / "duplicate.json"
    mutated.write_text(duplicated, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        verify(ROOT, mutated)


def test_arb_overlap_cannot_pass_a_strict_lower_gate() -> None:
    candidate = arb(1, (1, -10))
    endpoint = candidate.lower()
    threshold = arb(endpoint.mid(), (1, -80))
    assert not (endpoint <= threshold)
    assert not (endpoint > threshold)
    assert not _strict_lower_proved(candidate, threshold)


def test_fabricated_narrow_full_gain_cannot_replace_the_enclosure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["full_gain"] = _wide_ball_record(
        arb("0.0002 +/- 1e-30")
    )
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "narrow-full-gain.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="does not contain replay"):
        verify(ROOT, mutated)
