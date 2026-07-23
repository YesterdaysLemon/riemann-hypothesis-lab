from __future__ import annotations

import json
from pathlib import Path

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
import tools.generate_nyman_large_sieve_tail_certificate as certificate_tool
from tools.generate_nyman_large_sieve_tail_certificate import (
    CLAIMED_GAIN_LOWER_BOUND,
    FROZEN_ARTIFACT_PAYLOAD_SHA256,
    SCHEMA,
    _wide_ball_record,
    generate,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-large-sieve-tail-v1.json"


def _rehash(artifact: dict[str, object]) -> None:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)


def test_checked_in_large_sieve_certificate_replays() -> None:
    result = verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    assert result["claimed_gain_lower_bound"] == "1/5000"
    assert int(result["periodic_bound_reduction_factor_lower"]) >= 632


def test_generation_reproduces_frozen_payload() -> None:
    artifact = generate(ROOT)
    assert artifact["schema"] == SCHEMA
    assert artifact["payload_sha256"] == FROZEN_ARTIFACT_PAYLOAD_SHA256
    assert artifact == json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert CLAIMED_GAIN_LOWER_BOUND.numerator == 1
    assert CLAIMED_GAIN_LOWER_BOUND.denominator == 5_000


def test_large_sieve_certificate_rejects_duplicate_json_keys(
    tmp_path: Path,
) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    duplicated = raw.replace("{", '{"schema":"false",', 1)
    mutated = tmp_path / "duplicate.json"
    mutated.write_text(duplicated, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        verify(ROOT, mutated)


def test_rehashed_spacing_mutation_fails_semantic_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["tail"]["farey_spacing_reciprocal"] = "1"
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "spacing.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="tail farey_spacing_reciprocal changed"):
        verify(ROOT, mutated)


def test_fabricated_large_sieve_full_gain_cannot_replace_enclosure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["full_gain"] = _wide_ball_record(arb("0.00021 +/- 1e-30"))
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "narrow-full-gain.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="not canonical"):
        verify(ROOT, mutated)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("display", "Arb display changed"),
        ("is_exact", "Arb exactness flag changed"),
        ("integer_spelling", "mid_exponent is not canonical"),
    ],
)
def test_rehashed_arb_metadata_and_encoding_mutations_are_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    message: str,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    enclosure = artifact["full_gain"]["enclosure"]
    if mutation == "display":
        enclosure["display"] = "PROVED RH"
    elif mutation == "is_exact":
        enclosure["is_exact"] = not enclosure["is_exact"]
    else:
        exponent = enclosure["dyadic"]["mid_exponent"]
        enclosure["dyadic"]["mid_exponent"] = exponent.replace("-", "-0", 1)
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / f"arb-{mutation}.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        verify(ROOT, mutated)
