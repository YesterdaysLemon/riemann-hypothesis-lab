from __future__ import annotations

import json
from pathlib import Path

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
import tools.generate_nyman_large_sieve_scaling_certificate as certificate_tool
from tools.generate_nyman_large_sieve_scaling_certificate import (
    CLAIMED_GAIN_LOWER_BOUNDS,
    FROZEN_ARTIFACT_PAYLOAD_SHA256,
    N_VALUES,
    SCHEMA,
    generate,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-large-sieve-scaling-v1.json"


def _rehash(artifact: dict[str, object]) -> None:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)


def test_checked_in_large_sieve_scaling_certificate_replays() -> None:
    result = verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    assert [int(row["n"]) for row in result["replay_cells"]] == list(N_VALUES)


def test_scaling_generation_reproduces_frozen_payload() -> None:
    artifact = generate(ROOT)
    assert artifact["schema"] == SCHEMA
    assert artifact["payload_sha256"] == FROZEN_ARTIFACT_PAYLOAD_SHA256
    assert artifact == json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert set(CLAIMED_GAIN_LOWER_BOUNDS) == set(N_VALUES)


def test_scaling_certificate_rejects_duplicate_json_keys(
    tmp_path: Path,
) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    duplicated = raw.replace("{", '{"schema":"false",', 1)
    mutated = tmp_path / "duplicate.json"
    mutated.write_text(duplicated, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        verify(ROOT, mutated)


def test_rehashed_scaling_classification_mutation_fails_semantics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["classification"] = "PROVED"
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "classification.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="classification changed"):
        verify(ROOT, mutated)


def test_rehashed_scaling_tail_mutation_fails_exact_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["cells"][0]["tail"]["farey_spacing_reciprocal"] = "1"
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "tail.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="tail farey_spacing_reciprocal changed"):
        verify(ROOT, mutated)


def test_rehashed_scaling_extra_cell_claim_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["cells"][0]["unverified_claim"] = "RH_PROVED"
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "extra-cell-field.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="cell fields changed"):
        verify(ROOT, mutated)


def test_rehashed_scaling_fabricated_enclosure_cannot_replace_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["cells"][0]["full_gain"] = {
        "enclosure": arb_record(arb("0.007 +/- 1e-30"))
    }
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    mutated = tmp_path / "fabricated-enclosure.json"
    mutated.write_text(json.dumps(artifact), encoding="utf-8")
    with pytest.raises(ValueError, match="does not contain replay"):
        verify(ROOT, mutated)
