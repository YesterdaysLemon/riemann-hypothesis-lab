from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from tools import generate_nyman_fixed_modulus_refinement_certificate as cert


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-fixed-modulus-refinement-v1.json"
RAW_LF_SHA256 = (
    "45571b59b62c613fe5a107db76a1c80a732d223254065d0d23ab08230894db2b"
)
CANONICAL_SHA256 = (
    "5666bfbd314da21f0d286924f2a4b93c3ff97fc025d269dc769f396e627efe54"
)


def _artifact() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def _fraction(record: object) -> Fraction:
    assert isinstance(record, dict)
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def _rehash(artifact: dict[str, object]) -> str:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)
    return str(artifact["payload_sha256"])


def _write(tmp_path: Path, artifact: dict[str, object]) -> Path:
    path = tmp_path / "certificate.json"
    path.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return path


def test_frozen_artifact_has_exact_scope_hashes_and_boolean_ledger() -> None:
    raw = ARTIFACT.read_bytes()
    artifact = _artifact()
    assert b"\r" not in raw and raw.endswith(b"\n")
    assert hashlib.sha256(raw).hexdigest() == RAW_LF_SHA256
    assert content_sha256(artifact) == CANONICAL_SHA256
    assert artifact["payload_sha256"] == cert.FROZEN_ARTIFACT_PAYLOAD_SHA256
    assert artifact["schema"] == cert.SCHEMA
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["claimed_gain_lower_bound"] == {
        "numerator": "1",
        "denominator": "17500",
    }
    assert artifact["checks"] == cert.EXPECTED_CHECKS
    assert all(type(value) is bool for value in artifact["checks"].values())
    assert artifact["checks"]["optimizer_used_as_evidence"] is False
    assert artifact["checks"]["resolves_rh"] is False
    assert "no fourth step" in artifact["limitation"]


def test_parent_lineage_and_fixed_modulus_tail_are_exactly_bound() -> None:
    artifact = _artifact()
    assert artifact["source"]["parent_artifact"] == {
        "path": "results/nyman-rebased-third-step-v1.json",
        "raw_lf_sha256": cert.PARENT_ARTIFACT_RAW_LF_SHA256,
        "canonical_sha256": cert.PARENT_ARTIFACT_CANONICAL_SHA256,
        "payload_sha256": cert.PARENT_ARTIFACT_PAYLOAD_SHA256,
        "schema": cert.third.SCHEMA,
        "p3_payload_sha256": cert.PARENT_P3_PAYLOAD_SHA256,
    }
    tail = artifact["proof"]["new_fixed_modulus_tail_upper"]
    assert tail["modulus"] == "840"
    assert tail["cutoff"] == str(1 << 23)
    assert tail["active_denominator_count"] == "17157"
    assert tail["active_denominator_maximum"] == "32768"
    assert len(tail["group_reciprocal_bounds"]) == 32
    fixed = _fraction(tail["fixed_modulus_spacing_constant"])
    active_q = _fraction(tail["active_q_spacing_constant"])
    assert fixed < active_q
    assert active_q / fixed > Fraction(118, 100)
    assert _fraction(tail["upper_bound"]) < Fraction(27, 1_000_000)


def test_independent_replay_proves_the_stronger_finite_claim() -> None:
    result = cert.verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    lower = arb(result["replay_gain_lower_bound_computation"]["lower"])
    assert lower.lower() > arb(1) / 17_500
    assert result["replay_environment"]["block_size"] == "250003"
    assert result["replay_environment"]["precision_bits"] == "512"


def test_generation_and_replay_enclosures_must_overlap() -> None:
    generation = (
        object(),
        SimpleNamespace(prefix_energy=arb(1)),
        SimpleNamespace(prefix_energy=arb(2)),
        object(),
        arb(3),
        arb(4),
    )
    cert._require_overlap(generation, generation)
    replay = (*generation[:-1], arb(5))
    with pytest.raises(ArithmeticError, match="independent replay disagree"):
        cert._require_overlap(generation, replay)


def test_rejects_unhashed_semantic_mutation(tmp_path: Path) -> None:
    artifact = _artifact()
    artifact["classification"] = "PROVED_RH"
    with pytest.raises(ValueError, match="payload SHA-256 mismatch"):
        cert.verify(ROOT, _write(tmp_path, artifact))


def test_rejects_rehashed_unknown_field(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["unexpected"] = "fail closed"
    frozen = _rehash(artifact)
    monkeypatch.setattr(cert, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    with pytest.raises(ValueError, match="top-level fields changed"):
        cert.verify(ROOT, _write(tmp_path, artifact))


def test_rejects_rehashed_parent_binding_before_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["source"]["parent_artifact"]["p3_payload_sha256"] = "0" * 64
    frozen = _rehash(artifact)
    monkeypatch.setattr(cert, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    with pytest.raises(ValueError, match="parent or source binding changed"):
        cert.verify(ROOT, _write(tmp_path, artifact))


def test_rejects_rehashed_numeric_boolean_ledger(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["checks"] = {
        name: int(value) for name, value in artifact["checks"].items()
    }
    frozen = _rehash(artifact)
    monkeypatch.setattr(cert, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    with pytest.raises(ValueError, match="ledger values are not booleans"):
        cert.verify(ROOT, _write(tmp_path, artifact))


def test_rejects_rehashed_numeric_alias_for_proof_boolean(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    enclosure = artifact["proof"]["gain_lower_bound_computation"]["enclosure"]
    enclosure["is_exact"] = int(enclosure["is_exact"])
    frozen = _rehash(artifact)
    monkeypatch.setattr(cert, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    with pytest.raises(ValueError, match="invalid Arb exactness flag"):
        cert.verify(ROOT, _write(tmp_path, artifact))


def test_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    mutated = raw.replace(
        '  "schema":',
        '  "schema": "duplicate must fail",\n  "schema":',
        1,
    )
    path = tmp_path / "duplicate.json"
    path.write_text(mutated, encoding="utf-8", newline="\n")
    with pytest.raises(ValueError, match="duplicate JSON key 'schema'"):
        cert.verify(ROOT, path)
