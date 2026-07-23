from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
from types import SimpleNamespace

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from tools import generate_nyman_rebased_fourth_step_certificate as certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-rebased-fourth-step-v1.json"


def _artifact() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def _rehash(artifact: dict[str, object]) -> str:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    artifact["payload_sha256"] = content_sha256(body)
    return str(artifact["payload_sha256"])


def _write_artifact(
    tmp_path: Path,
    artifact: dict[str, object],
) -> Path:
    path = tmp_path / "certificate.json"
    path.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return path


@pytest.fixture(scope="module")
def witness() -> certificate.FourthStepWitness:
    certificate._require_environment()
    return certificate._construct_witness(ROOT)


@pytest.fixture(scope="module")
def tail(
    witness: certificate.FourthStepWitness,
) -> certificate.third.LocalSpacingTail:
    return certificate.third.local_spacing_tail_upper(
        witness.p4,
        certificate.CUTOFF,
    )


def test_frozen_artifact_has_honest_scope_and_exact_commitments() -> None:
    artifact = _artifact()
    assert artifact["schema"] == certificate.SCHEMA
    assert artifact["payload_sha256"] == (
        certificate.FROZEN_ARTIFACT_PAYLOAD_SHA256
    )
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["claimed_gain_lower_bound"] == {
        "numerator": "1",
        "denominator": "100000",
    }
    assert artifact["configuration"]["cutoff"] == str(1 << 26)
    assert artifact["configuration"]["shell_4_multiplier_limit"] == "1"
    assert artifact["configuration"]["generation_block_size"] == str(1 << 20)
    assert artifact["configuration"]["generation_precision_bits"] == "512"
    assert artifact["configuration"]["replay_block_size"] == "250003"
    assert artifact["configuration"]["replay_precision_bits"] == "640"
    assert artifact["checks"] == certificate.EXPECTED_CHECKS
    assert all(type(value) is bool for value in artifact["checks"].values())
    assert artifact["checks"]["optimizer_used_as_evidence"] is False
    assert artifact["checks"]["resolves_rh"] is False
    assert artifact["reference"]["doi"] == "10.1112/jlms/s2-8.1.73"
    assert "remains unresolved" in artifact["limitation"]
    assert "No fifth step" in artifact["limitation"]

    source = artifact["source"]["parent_artifact"]
    assert source == {
        "path": "results/nyman-rebased-third-step-v1.json",
        "raw_lf_sha256": certificate.PARENT_ARTIFACT_RAW_LF_SHA256,
        "canonical_sha256": certificate.PARENT_ARTIFACT_CANONICAL_SHA256,
        "payload_sha256": certificate.PARENT_ARTIFACT_PAYLOAD_SHA256,
        "schema": certificate.third.SCHEMA,
        "p3_payload_sha256": certificate.PARENT_P3_PAYLOAD_SHA256,
    }
    vectors = artifact["vectors"]
    assert vectors["p3"]["payload_sha256"] == (
        certificate.PARENT_P3_PAYLOAD_SHA256
    )
    assert vectors["shell_4"] == {
        "denominator_exponent": "9",
        "support_count": "32685",
        "maximum_index": "65536",
        "payload_sha256": certificate.SHELL_4_PAYLOAD_SHA256,
    }
    assert vectors["weights_4"] == {
        "denominator_exponent": "16",
        "support_count": "33",
        "maximum_index": "33",
        "payload_sha256": certificate.WEIGHTS_4_PAYLOAD_SHA256,
    }
    assert vectors["p4"] == {
        "denominator_exponent": "25",
        "support_count": "47345",
        "maximum_index": "65536",
        "payload_sha256": certificate.P4_PAYLOAD_SHA256,
    }
    assert vectors["shell_4_parent_p3_payload_sha256"] == (
        certificate.PARENT_P3_PAYLOAD_SHA256
    )


def test_integer_serialization_limit_covers_exact_tail_record() -> None:
    original_limit = sys.get_int_max_str_digits()
    certificate._require_environment()
    assert sys.get_int_max_str_digits() == original_limit
    with certificate._integer_serialization_context():
        raised_limit = sys.get_int_max_str_digits()
        assert raised_limit == 0 or (
            raised_limit >= certificate.MAX_INTEGER_DECIMAL_DIGITS
        )
        artifact = _artifact()
        largest_declared_denominator = artifact["proof"]["new_tail_upper"][
            "upper_bound"
        ]["denominator"]
        assert len(largest_declared_denominator) > 4_300
        assert raised_limit == 0 or (
            len(largest_declared_denominator) < raised_limit
        )
    assert sys.get_int_max_str_digits() == original_limit


def test_reconstructed_fourth_shell_and_p4_are_exactly_bound(
    witness: certificate.FourthStepWitness,
) -> None:
    assert certificate.third.parent._vector_commitment(
        witness.parent_witness.p3
    )["payload_sha256"] == certificate.PARENT_P3_PAYLOAD_SHA256
    assert witness.shell4 == certificate.third.parent.rounded_ideal_shell(
        witness.parent_witness.p3,
        certificate.SHELL_BITS,
        support_cutoff=32_768,
    )
    assert len(witness.shell4) == 32_685
    assert min(witness.shell4) == 32_769
    assert max(witness.shell4) == 65_536
    assert certificate.third.parent._coefficient_sum(witness.shell4) == 0
    assert sum(
        (abs(value) for value in witness.shell4.values()),
        start=Fraction(),
    ) == Fraction(3_159_119, 256)
    assert all(
        (value * (1 << certificate.SHELL_BITS)).denominator == 1
        for value in witness.shell4.values()
    )
    assert tuple(
        witness.weights4[index] for index in range(1, 34)
    ) == certificate.WEIGHTS_4
    assert certificate.third.parent._vector_commitment(
        witness.weights4
    )["payload_sha256"] == certificate.WEIGHTS_4_PAYLOAD_SHA256
    assert len(witness.p4) == 47_345
    assert min(witness.p4) == 1
    assert max(witness.p4) == 65_536
    assert certificate.third.parent._vector_commitment(
        witness.p4
    )["payload_sha256"] == certificate.P4_PAYLOAD_SHA256
    assert certificate.third.parent._coefficient_sum(
        witness.p4
    ) == Fraction(63_751, 32_768)
    assert sum(
        (abs(value) for value in witness.p4.values()),
        start=Fraction(),
    ) == Fraction(19_391_840_491, 8_388_608)


def test_changed_and_unchanged_old_coordinates_are_pinned(
    witness: certificate.FourthStepWitness,
) -> None:
    p3 = witness.parent_witness.p3
    changed = sum(
        witness.p4.get(index, Fraction()) != value
        for index, value in p3.items()
    )
    unchanged = sum(
        witness.p4.get(index, Fraction()) == value
        for index, value in p3.items()
    )
    assert changed == 13_537
    assert unchanged == 1_123
    assert changed + unchanged == len(p3) == 14_660
    assert {
        index for index in witness.p4 if index <= 32_768
    } == set(p3)

    delta = certificate.third.parent._clean(
        {
            index: witness.p4.get(index, Fraction())
            - p3.get(index, Fraction())
            for index in set(witness.p4) | set(p3)
        }
    )
    assert certificate.third.parent._coefficient_sum(delta) == Fraction(
        -1,
        65_536,
    )
    assert certificate.third.parent._harmonic_sum(delta) != 0


def test_candidate_specific_tail_constants_and_stored_record(
    witness: certificate.FourthStepWitness,
    tail: certificate.third.LocalSpacingTail,
) -> None:
    assert tail.cutoff == 1 << 26
    assert tail.support_limit == 65_536
    assert tail.coefficient_sum == Fraction(63_751, 32_768)
    assert tail.c0 == Fraction(1_785, 65_536)
    assert tail.harmonic_sum == certificate.third.parent._harmonic_sum(
        witness.p4
    )
    assert tail.absolute_residual_bound == Fraction(
        19_392_297_451,
        16_777_216,
    )
    assert tail.upper_bound == (
        tail.mean_term
        + tail.discrepancy_term
        + tail.cross_term
        + tail.slope_square_term
    )
    assert tail.local_spacing_constant == (
        Fraction(3, 2) * tail.support_limit * tail.sigma
    )
    assert tail.global_spacing_constant == (
        tail.support_limit * (tail.support_limit - 1) * tail.rho
    )
    assert 3.4776 < float(tail.global_to_local_improvement) < 3.4777
    assert 0.0000052313 < float(tail.upper_bound) < 0.0000052314
    with certificate._integer_serialization_context():
        record = certificate.third._tail_record(tail)
    assert content_sha256(record) == certificate.TAIL_RECORD_SHA256
    assert _artifact()["proof"]["new_tail_upper"] == record


def test_parent_artifact_matches_fresh_exact_reconstruction(
    witness: certificate.FourthStepWitness,
) -> None:
    parent_artifact = certificate._load_parent_artifact(ROOT)
    assert parent_artifact["payload_sha256"] == (
        certificate.PARENT_ARTIFACT_PAYLOAD_SHA256
    )
    assert content_sha256(parent_artifact) == (
        certificate.PARENT_ARTIFACT_CANONICAL_SHA256
    )
    assert parent_artifact["source"] == witness.parent_witness.source
    assert parent_artifact["vectors"] == certificate.third._vectors_record(
        witness.parent_witness
    )


@pytest.mark.slow
def test_independent_replay_closes_above_conservative_claim() -> None:
    result = certificate.verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    replay_lower = arb(
        result["replay_gain_lower_bound_computation"]["lower"]
    )
    assert replay_lower.lower() > arb(1) / 100_000
    assert result["replay_environment"]["block_size"] == "250003"
    assert result["replay_environment"]["precision_bits"] == "640"


def test_generation_and_replay_enclosures_must_overlap() -> None:
    generation = (
        object(),
        SimpleNamespace(prefix_energy=arb(1)),
        SimpleNamespace(prefix_energy=arb(2)),
        object(),
        arb(3),
        arb(4),
    )
    replay = generation
    certificate._require_overlap(generation, replay)

    disagreeing_replay = (*replay[:-1], arb(5))
    with pytest.raises(
        ArithmeticError,
        match="generation and independent replay disagree",
    ):
        certificate._require_overlap(generation, disagreeing_replay)


@pytest.mark.parametrize("record_name", ["source", "vectors"])
def test_construct_witness_rejects_parent_reconstruction_mismatch(
    monkeypatch: pytest.MonkeyPatch,
    record_name: str,
) -> None:
    parent_artifact = deepcopy(certificate._load_parent_artifact(ROOT))
    if record_name == "source":
        parent_artifact["source"]["parent_artifact"]["p2_payload_sha256"] = (
            "0" * 64
        )
        match = "parent source lineage differs"
    else:
        parent_artifact["vectors"]["p3"]["payload_sha256"] = "0" * 64
        match = "parent vectors differ"
    monkeypatch.setattr(
        certificate,
        "_load_parent_artifact",
        lambda root: parent_artifact,
    )
    with pytest.raises(ArithmeticError, match=match):
        certificate._construct_witness(ROOT)


def test_rejects_unhashed_artifact_mutation(tmp_path: Path) -> None:
    artifact = _artifact()
    artifact["classification"] = "PROVED_RH"
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="payload SHA-256 mismatch"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_hypothesis_status_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["hypothesis_status"] = "RESOLVED"
    frozen = _rehash(artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="hypothesis_status changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_unknown_field(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["unexpected"] = "must fail closed"
    frozen = _rehash(artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="top-level fields changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_weight_mutation_before_energy_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["configuration"]["weights_4"][32] = {
        "numerator": "2391",
        "denominator": "65536",
    }
    frozen = _rehash(artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="configuration changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_parent_binding_mutation_before_energy_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["source"]["parent_artifact"]["p3_payload_sha256"] = "0" * 64
    frozen = _rehash(artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="parent or source binding changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_numeric_boolean_ledger_before_energy_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["checks"] = {
        name: int(value)
        for name, value in artifact["checks"].items()
    }
    frozen = _rehash(artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="ledger values are not booleans"):
        certificate.verify(ROOT, path)


def test_generated_proof_rejects_numeric_alias_for_arb_boolean() -> None:
    generated = deepcopy(_artifact()["proof"])
    stored = deepcopy(generated)
    exactness = stored["old_complete_energy_lower_from_prefix"][
        "prefix_lower_endpoint"
    ]["is_exact"]
    assert type(exactness) is bool
    stored["old_complete_energy_lower_from_prefix"][
        "prefix_lower_endpoint"
    ]["is_exact"] = int(exactness)
    assert stored == generated
    with pytest.raises(ValueError, match="invalid Arb exactness flag"):
        certificate.third._verify_generated_proof(stored, generated)


def test_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    raw = ARTIFACT.read_text(encoding="utf-8")
    mutated = raw.replace(
        '  "schema":',
        '  "schema": "duplicate must fail",\n  "schema":',
        1,
    )
    assert mutated != raw
    path = tmp_path / "duplicate.json"
    path.write_text(mutated, encoding="utf-8", newline="\n")
    with pytest.raises(ValueError, match="duplicate JSON key 'schema'"):
        certificate.verify(ROOT, path)
