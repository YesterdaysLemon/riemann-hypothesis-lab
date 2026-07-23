from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from flint import arb, ctx
import pytest

from riemann_lab.artifacts import content_sha256
from tools import generate_nyman_rebased_schur_certificate as certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "nyman-rebased-schur-v1.json"


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


def _write_artifact(tmp_path: Path, artifact: dict[str, object]) -> Path:
    path = tmp_path / "certificate.json"
    path.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
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
        "denominator": "5000",
    }
    assert artifact["checks"]["optimizer_used_as_evidence"] is False
    assert artifact["checks"]["resolves_rh"] is False
    assert "remains unresolved" in artifact["limitation"]

    assert artifact["vectors"] == {
        "p0_source": {
            "denominator_exponent": "256",
            "support_count": "8",
            "maximum_index": "8",
            "payload_sha256": (
                "4200517a85a84817b64be599aa5680abb"
                "0a7e45ad35b6ada53692410b56ee2d6"
            ),
        },
        "shell_1": {
            "denominator_exponent": "9",
            "support_count": "7",
            "maximum_index": "16",
            "payload_sha256": (
                "376513fa3ce5dfe30d6b1b738ea0ddf11"
                "63e643da2600426bc5c7ca68683eeb8"
            ),
        },
        "weights_1": {
            "denominator_exponent": "16",
            "support_count": "16",
            "maximum_index": "16",
            "payload_sha256": (
                "c93e6ad2f5d90ee3e19dc86249aa2195"
                "d335b763c70c02e5bbfa0c97a74c21f1"
            ),
        },
        "p1": {
            "denominator_exponent": "25",
            "support_count": "55",
            "maximum_index": "128",
            "payload_sha256": (
                "4d1f3bb4359542d5007f560b4d57b404"
                "26d858b5831377b0fd4f6d8757ae377b"
            ),
        },
        "shell_2": {
            "denominator_exponent": "9",
            "support_count": "128",
            "maximum_index": "256",
            "payload_sha256": (
                "100d0408424f33cf659c64cdd60895a28"
                "dc9c82aedefc40e703d2358a0cadbb6"
            ),
        },
        "weights_2": {
            "denominator_exponent": "16",
            "support_count": "24",
            "maximum_index": "24",
            "payload_sha256": (
                "65aa8df808e491f02db2e9c20573c430"
                "ff3b32200089f849eab6c7a5ebd6110d"
            ),
        },
        "p2": {
            "denominator_exponent": "25",
            "support_count": "915",
            "maximum_index": "2048",
            "payload_sha256": (
                "b392a0c8d9ea2bf09b5fdf98575906bb"
                "8c17b75044c4a01586bdc12fa299520a"
            ),
        },
        "shell_2_parent_p1_payload_sha256": (
            "4d1f3bb4359542d5007f560b4d57b404"
            "26d858b5831377b0fd4f6d8757ae377b"
        ),
    }


def test_independent_replay_closes_above_conservative_claim() -> None:
    result = certificate.verify(ROOT, ARTIFACT)
    assert result["verified"] is True
    assert result["hypothesis_status"] == "UNRESOLVED"
    replay_lower = result["replay_gain_lower_bound_computation"]["lower"]
    assert replay_lower.startswith("[0.000317646")
    assert result["replay_environment"]["block_size"] == "32749"
    assert result["replay_environment"]["precision_bits"] == "448"


def test_reconstructed_second_shell_is_bound_to_exact_p1() -> None:
    witness = certificate._construct_witness(ROOT)
    assert certificate._vector_commitment(witness.p1)["payload_sha256"] == (
        "4d1f3bb4359542d5007f560b4d57b404"
        "26d858b5831377b0fd4f6d8757ae377b"
    )
    assert witness.shell2 == certificate.rounded_ideal_shell(
        witness.p1,
        certificate.SHELL_BITS,
        support_cutoff=128,
    )
    assert len(witness.shell2) == 128
    assert max(witness.shell2) == 256
    assert certificate._coefficient_sum(witness.shell2) == 0


def test_candidate_specific_absolute_tail_constants() -> None:
    witness = certificate._construct_witness(ROOT)
    tail = certificate.absolute_energy_tail_upper(
        witness.p2,
        certificate.CUTOFF,
    )
    assert tail.support_limit == 2048
    assert tail.farey_spacing_reciprocal == 2048 * 2047
    assert tail.c0 == Fraction(223, 8192)
    assert tail.absolute_residual_bound == Fraction(
        2_234_890_213,
        33_554_432,
    )
    assert tail.upper_bound == (
        tail.mean_term
        + tail.discrepancy_term
        + tail.cross_term
        + tail.slope_square_term
    )
    assert 0.00034690 < float(tail.upper_bound) < 0.00034691


def test_tail_formula_keeps_constant_term_for_empty_vector() -> None:
    tail = certificate.absolute_energy_tail_upper({}, 7)
    assert tail.support_limit == 0
    assert tail.c0 == 1
    assert tail.rho == 1
    assert tail.mean_term == Fraction(1, 8)
    assert tail.discrepancy_term == 0
    assert tail.cross_term == 0
    assert tail.slope_square_term == 0
    assert tail.upper_bound == Fraction(1, 8)


def test_absolute_prefix_contains_direct_small_arb_sum() -> None:
    coefficients = {1: Fraction(1)}
    cutoff = 16
    prefix = certificate.absolute_energy_prefix(
        coefficients,
        cutoff,
        block_size=7,
        precision_bits=256,
    )
    previous_precision = ctx.prec
    ctx.prec = 256
    try:
        p_value = Fraction(1)
        direct = arb(1)
        for interval in range(1, cutoff + 1):
            q_value = 1 + interval
            direct += (
                arb(q_value * q_value) / (interval * (interval + 1))
                - 2
                * certificate._arb_from_fraction(p_value)
                * q_value
                * (arb(interval + 1) / interval).log()
                + 1
            )
        assert prefix.prefix_energy.contains(direct)
    finally:
        ctx.prec = previous_precision


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
    monkeypatch.setattr(certificate, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="top-level fields changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_weight_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["configuration"]["weights_2"][23] = {
        "numerator": "40",
        "denominator": "2048",
    }
    frozen = _rehash(artifact)
    monkeypatch.setattr(certificate, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="configuration changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_parent_binding_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["vectors"]["shell_2_parent_p1_payload_sha256"] = "0" * 64
    frozen = _rehash(artifact)
    monkeypatch.setattr(certificate, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="vector commitments changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_fabricated_gain_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    record = artifact["proof"]["gain_lower_bound_computation"]["enclosure"]
    record["display"] = "PROVED RH"
    frozen = _rehash(artifact)
    monkeypatch.setattr(certificate, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
    path = _write_artifact(tmp_path, artifact)
    with pytest.raises(ValueError, match="generation proof record changed"):
        certificate.verify(ROOT, path)


def test_rejects_rehashed_numeric_boolean_ledger(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _artifact()
    artifact["checks"] = {
        name: int(value) for name, value in artifact["checks"].items()
    }
    frozen = _rehash(artifact)
    monkeypatch.setattr(certificate, "FROZEN_ARTIFACT_PAYLOAD_SHA256", frozen)
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
    path.write_text(mutated, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key 'schema'"):
        certificate.verify(ROOT, path)
