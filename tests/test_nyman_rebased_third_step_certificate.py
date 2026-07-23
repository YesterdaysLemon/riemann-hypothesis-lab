from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
from types import SimpleNamespace

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from tools import generate_nyman_rebased_third_step_certificate as certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-rebased-third-step-v1.json"


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
        "denominator": "20000",
    }
    assert artifact["configuration"]["cutoff"] == str(1 << 23)
    assert artifact["checks"] == certificate.EXPECTED_CHECKS
    assert all(type(value) is bool for value in artifact["checks"].values())
    assert artifact["checks"]["optimizer_used_as_evidence"] is False
    assert artifact["checks"]["resolves_rh"] is False
    assert artifact["reference"]["doi"] == "10.1112/jlms/s2-8.1.73"
    assert "periodic weighted cosecant" in (
        artifact["reference"]["used_result"]
    )
    assert "remains unresolved" in artifact["limitation"]

    source = artifact["source"]["parent_artifact"]
    assert source == {
        "path": "results/nyman-rebased-schur-v1.json",
        "raw_lf_sha256": certificate.PARENT_ARTIFACT_RAW_LF_SHA256,
        "canonical_sha256": certificate.PARENT_ARTIFACT_CANONICAL_SHA256,
        "payload_sha256": certificate.PARENT_ARTIFACT_PAYLOAD_SHA256,
        "schema": certificate.parent.SCHEMA,
        "p2_payload_sha256": certificate.PARENT_P2_PAYLOAD_SHA256,
    }
    vectors = artifact["vectors"]
    assert vectors["p2"]["payload_sha256"] == (
        certificate.PARENT_P2_PAYLOAD_SHA256
    )
    assert vectors["shell_3"] == {
        "denominator_exponent": "9",
        "support_count": "2048",
        "maximum_index": "4096",
        "payload_sha256": certificate.SHELL_3_PAYLOAD_SHA256,
    }
    assert vectors["weights_3"] == {
        "denominator_exponent": "16",
        "support_count": "32",
        "maximum_index": "32",
        "payload_sha256": certificate.WEIGHTS_3_PAYLOAD_SHA256,
    }
    assert vectors["p3"] == {
        "denominator_exponent": "25",
        "support_count": "14660",
        "maximum_index": "32768",
        "payload_sha256": certificate.P3_PAYLOAD_SHA256,
    }
    assert vectors["shell_3_parent_p2_payload_sha256"] == (
        certificate.PARENT_P2_PAYLOAD_SHA256
    )


def test_independent_replay_closes_above_conservative_claim() -> None:
    result = certificate.verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["classification"] == "CERTIFIED_FINITE"
    assert result["hypothesis_status"] == "UNRESOLVED"
    replay_lower = arb(
        result["replay_gain_lower_bound_computation"]["lower"]
    )
    assert replay_lower.lower() > arb(1) / 20_000
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
    replay = generation
    certificate._require_replay_overlap(generation, replay)

    disagreeing_replay = (*replay[:-1], arb(5))
    with pytest.raises(
        ArithmeticError,
        match="generation and independent replay disagree",
    ):
        certificate._require_replay_overlap(generation, disagreeing_replay)


def test_reconstructed_third_shell_and_p3_are_exactly_bound() -> None:
    witness = certificate._construct_witness(ROOT)
    assert certificate.parent._vector_commitment(
        witness.parent_witness.p2
    )["payload_sha256"] == certificate.PARENT_P2_PAYLOAD_SHA256
    assert witness.shell3 == certificate.parent.rounded_ideal_shell(
        witness.parent_witness.p2,
        certificate.SHELL_BITS,
        support_cutoff=2048,
    )
    assert len(witness.shell3) == 2048
    assert min(witness.shell3) == 2049
    assert max(witness.shell3) == 4096
    assert certificate.parent._coefficient_sum(witness.shell3) == 0
    assert all(
        (value * (1 << certificate.SHELL_BITS)).denominator == 1
        for value in witness.shell3.values()
    )
    assert tuple(
        witness.weights3[index] for index in range(1, 33)
    ) == certificate.WEIGHTS_3
    assert len(witness.p3) == 14_660
    assert max(witness.p3) == 32_768


def test_candidate_specific_absolute_tail_constants_and_cutoff() -> None:
    witness = certificate._construct_witness(ROOT)
    tail = certificate.local_spacing_tail_upper(
        witness.p3,
        certificate.CUTOFF,
    )
    assert tail.cutoff == 1 << 23
    assert tail.support_limit == 32_768
    assert tail.coefficient_sum == Fraction(127_503, 65_536)
    assert tail.c0 == Fraction(3_569, 131_072)
    assert tail.absolute_residual_bound == Fraction(
        8_875_081_453,
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
    divisor_sums = certificate.parent._divisor_harmonic_sums(
        witness.p3,
        tail.support_limit,
    )
    jordan = certificate.parent._jordan_j2_sieve(tail.support_limit)
    assert tail.sigma == tail.c0 * tail.c0 + sum(
        (
            divisor
            * jordan[divisor]
            * divisor_sums[divisor]
            * divisor_sums[divisor]
            for divisor in range(2, tail.support_limit + 1)
        ),
        start=Fraction(),
    ) / 12
    assert 3.1395 < float(tail.global_to_local_improvement) < 3.1396
    assert 0.0000316177 < float(tail.upper_bound) < 0.0000316178


def test_local_spacing_constant_bounds_small_exact_periods() -> None:
    coefficients = {
        1: Fraction(1, 2),
        2: Fraction(-1, 4),
        3: Fraction(3, 8),
        5: Fraction(-1, 8),
    }
    tail = certificate.local_spacing_tail_upper(coefficients, 64)
    period = math.lcm(*coefficients)
    values: list[Fraction] = []
    for interval in range(period):
        centered = tail.c0 - sum(
            (
                coefficient
                * (
                    Fraction(interval % index, index)
                    - Fraction(index - 1, 2 * index)
                )
                for index, coefficient in coefficients.items()
            ),
            start=Fraction(),
        )
        values.append(centered * centered - tail.rho)

    assert sum(values, start=Fraction()) == 0
    for first in range(period):
        running = Fraction()
        for length in range(1, period + 1):
            running += values[(first + length - 1) % period]
            assert abs(running) <= tail.local_spacing_constant


def test_parent_artifact_is_strictly_pinned() -> None:
    parent_artifact = certificate._load_parent_artifact(ROOT)
    assert parent_artifact["payload_sha256"] == (
        certificate.PARENT_ARTIFACT_PAYLOAD_SHA256
    )
    assert content_sha256(parent_artifact) == (
        certificate.PARENT_ARTIFACT_CANONICAL_SHA256
    )
    assert parent_artifact["vectors"]["p2"]["payload_sha256"] == (
        certificate.PARENT_P2_PAYLOAD_SHA256
    )


def test_rejects_unhashed_artifact_mutation(tmp_path: Path) -> None:
    artifact = _artifact()
    artifact["classification"] = "PROVED_RH"
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="payload SHA-256 mismatch"):
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
    artifact["configuration"]["weights_3"][31] = {
        "numerator": "20",
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
    artifact["source"]["parent_artifact"]["p2_payload_sha256"] = "0" * 64
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
