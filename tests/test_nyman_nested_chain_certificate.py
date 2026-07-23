from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
from tools.certify_nyman_balanced_tail import harmonic_sum
import tools.generate_nyman_nested_chain_certificate as certificate_tool
from tools.generate_nyman_nested_chain_certificate import (
    FROZEN_ARTIFACT_PAYLOAD_SHA256,
    SCHEMA,
    STEP2_CLAIM,
    _construct_chain,
    generate,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-nested-chain-v1.json"
SCALING_ARTIFACT = (
    ROOT / "results" / "nyman-large-sieve-scaling-v1.json"
)


def _rehash(artifact: dict[str, object]) -> None:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)


def _write_rehashed_mutation(
    artifact: dict[str, object],
    path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _rehash(artifact)
    monkeypatch.setattr(
        certificate_tool,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        artifact["payload_sha256"],
    )
    path.write_text(json.dumps(artifact), encoding="utf-8")


def test_checked_in_nested_chain_certificate_replays() -> None:
    result = verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    assert [row["name"] for row in result["replay_steps"]] == [
        "step_1",
        "step_2",
    ]


def test_nested_chain_generation_reproduces_frozen_payload() -> None:
    artifact = generate(ROOT)
    assert artifact["schema"] == SCHEMA
    assert artifact["payload_sha256"] == FROZEN_ARTIFACT_PAYLOAD_SHA256
    assert artifact == json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert STEP2_CLAIM == Fraction(1, 9_000)


def test_nested_chain_commits_to_actual_step_1_updated_vector() -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    step1 = artifact["chain"][0]["vector_commitments"]
    step2 = artifact["chain"][1]["vector_commitments"]
    assert step1["updated"] == step2["old"]
    assert step2["old"]["payload_sha256"] == (
        "ca182ac99ba08892bf227d806293028bab2825db3327f7cb1a7a215e42f6f84b"
    )
    assert step2["shell"]["payload_sha256"] == (
        "4a5c707de249d9ae497c0f6cb4638de35717f193fa1da4f29df387c7d296509a"
    )
    assert step2["added"]["payload_sha256"] == (
        "5fd1364cfea00316991f6d71e96df8bc67bb6036498fc3a6e049a2e8d9a0282e"
    )
    assert step2["updated"]["payload_sha256"] == (
        "b3ed188142359e2ba26a45b7c2a0b8f187c6236b52b3d095539598a1ae1e3afb"
    )


def test_nested_chain_preserves_the_exact_frozen_prefix_moment() -> None:
    chain = _construct_chain(ROOT)
    moment = Fraction(241_057, 27_525_120)
    assert min(chain.a1) > max(chain.p0)
    assert min(chain.a2) > max(chain.p1)
    assert harmonic_sum(chain.a1) == 0
    assert harmonic_sum(chain.a2) == 0
    assert harmonic_sum(chain.p0) == moment
    assert harmonic_sum(chain.p1) == moment
    assert harmonic_sum(chain.p2) == moment
    assert moment * moment == Fraction(
        58_108_477_249,
        757_632_231_014_400,
    )


def test_nested_chain_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    duplicated = raw.replace("{", '{"schema":"false",', 1)
    mutated = tmp_path / "duplicate.json"
    mutated.write_text(duplicated, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        verify(ROOT, mutated)


def test_rehashed_nested_chain_rh_status_mutation_fails_semantics(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["hypothesis_status"] = "PROVED"
    mutated = tmp_path / "proved.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="hypothesis_status changed"):
        verify(ROOT, mutated)


def test_rehashed_independent_candidate_substitution_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    scaling = json.loads(SCALING_ARTIFACT.read_text(encoding="utf-8"))
    independent_n16 = next(
        cell for cell in scaling["cells"] if cell["n"] == "16"
    )
    artifact["chain"][0]["vector_commitments"]["old"] = independent_n16[
        "vector_commitments"
    ]["old"]
    mutated = tmp_path / "independent-candidate.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="step_1 vector commitments changed"):
        verify(ROOT, mutated)


def test_rehashed_step_2_old_substitution_breaks_nesting(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["chain"][1]["vector_commitments"]["old"] = artifact["chain"][0][
        "vector_commitments"
    ]["old"]
    mutated = tmp_path / "non-nested-old.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="step_2 vector commitments changed"):
        verify(ROOT, mutated)


def test_rehashed_step_2_multiplier_substitution_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["configuration"]["step2_multiplier"][1] = {
        "numerator": "-1",
        "denominator": "1",
    }
    mutated = tmp_path / "multiplier.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="configuration changed"):
        verify(ROOT, mutated)


def test_rehashed_fabricated_step_2_gain_cannot_replace_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["chain"][1]["full_gain"] = {
        "enclosure": arb_record(arb("0.001 +/- 1e-30"))
    }
    mutated = tmp_path / "fabricated-gain.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="does not contain replay"):
        verify(ROOT, mutated)


def test_rehashed_extra_step_claim_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    artifact["chain"][0]["unverified_claim"] = "RH_PROVED"
    mutated = tmp_path / "extra-step-field.json"
    _write_rehashed_mutation(artifact, mutated, monkeypatch)
    with pytest.raises(ValueError, match="step fields changed"):
        verify(ROOT, mutated)
