from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from tools import generate_nyman_rebased_fifth_step_certificate as certificate


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def witness() -> certificate.FifthStepWitness:
    """Perform the expensive exact p5 reconstruction only once."""

    certificate._require_environment()
    return certificate._construct_witness(ROOT)


def _fake_tail() -> certificate.fourth.third.LocalSpacingTail:
    support_limit = certificate.P5_MAXIMUM_INDEX
    sigma = Fraction(2)
    rho = Fraction(3)
    harmonic_sum = Fraction(1, 7)
    absolute_residual_bound = Fraction(4)
    local_spacing_constant = (
        Fraction(3, 2) * support_limit * sigma
    )
    global_spacing_constant = (
        support_limit * (support_limit - 1) * rho
    )
    mean_term = rho / (certificate.CUTOFF + 1)
    discrepancy_term = local_spacing_constant / (
        (certificate.CUTOFF + 1) * (certificate.CUTOFF + 2)
    )
    cross_term = (
        abs(harmonic_sum)
        * absolute_residual_bound
        / (certificate.CUTOFF + 1)
    )
    slope_square_term = (
        harmonic_sum
        * harmonic_sum
        / (4 * (certificate.CUTOFF + 1))
    )
    return certificate.fourth.third.LocalSpacingTail(
        cutoff=certificate.CUTOFF,
        support_limit=support_limit,
        coefficient_sum=certificate.P5_COEFFICIENT_SUM,
        harmonic_sum=harmonic_sum,
        c0=certificate.P5_C0,
        rho=rho,
        sigma=sigma,
        local_spacing_constant=local_spacing_constant,
        global_spacing_constant=global_spacing_constant,
        global_to_local_improvement=(
            global_spacing_constant / local_spacing_constant
        ),
        absolute_residual_bound=absolute_residual_bound,
        mean_term=mean_term,
        discrepancy_term=discrepancy_term,
        cross_term=cross_term,
        slope_square_term=slope_square_term,
        upper_bound=(
            mean_term
            + discrepancy_term
            + cross_term
            + slope_square_term
        ),
    )


def _fake_prefix(
    *,
    chunk_size: int,
    reduction_leaf_size: int,
    prefix_energy: Any,
) -> certificate.absolute_prefix.AbsoluteEnergyPrefix:
    maximum_q = certificate.P5_MAXIMUM_Q_NUMERATOR
    schedule = certificate._expected_schedule_signatures(
        chunk_size,
        reduction_leaf_size,
        (
            certificate.GENERATION_PRECISION_BITS
            if chunk_size == certificate.GENERATION_CHUNK_SIZE
            else certificate.REPLAY_PRECISION_BITS
        ),
    )
    if schedule is None:
        layout = (
            chunk_size,
            (certificate.CUTOFF + chunk_size - 1) // chunk_size,
            reduction_leaf_size,
            3,
            reduction_leaf_size,
        )
        depths = (2, 2, 4)
    else:
        layout, depths = schedule
    term_relative_error = (
        (1 + certificate.absolute_prefix.BINARY64_UNIT_ROUNDOFF)
        ** certificate.absolute_prefix.BINARY64_TERM_ROUNDING_FACTOR_COUNT
        - 1
    )
    reduction_gamma = certificate.absolute_prefix._gamma(depths[2])
    rational_midpoint = Fraction(1)
    computed_upper = rational_midpoint / (1 - reduction_gamma)
    reduction_error = reduction_gamma * computed_upper
    term_evaluation_error = (
        term_relative_error
        / (1 - term_relative_error)
        * computed_upper
    )
    rational_error = reduction_error + term_evaluation_error
    return certificate.absolute_prefix.AbsoluteEnergyPrefix(
        algorithm_version=(
            certificate.absolute_prefix.PREFIX_ALGORITHM_VERSION
        ),
        term_evaluation_algorithm=(
            certificate.absolute_prefix.TERM_EVALUATION_ALGORITHM
        ),
        reduction_algorithm=(
            certificate.absolute_prefix.REDUCTION_ALGORITHM
        ),
        cutoff=certificate.CUTOFF,
        chunk_size=chunk_size,
        chunk_count=layout[1],
        reduction_leaf_size=reduction_leaf_size,
        reduction_leaf_count=layout[3],
        maximum_leaf_length=layout[4],
        maximum_leaf_reduction_depth=depths[0],
        root_reduction_depth=depths[1],
        reduction_depth_bound=depths[2],
        coefficient_denominator=1
        << certificate.P5_DENOMINATOR_EXPONENT,
        coefficient_denominator_exponent=(
            certificate.P5_DENOMINATOR_EXPONENT
        ),
        recurrence_absolute_bound=(
            certificate.P5_RECURRENCE_ABSOLUTE_BOUND
        ),
        int64_headroom=certificate.P5_INT64_HEADROOM,
        term_rounding_factor_count=(
            certificate.absolute_prefix
            .BINARY64_TERM_ROUNDING_FACTOR_COUNT
        ),
        term_relative_error=term_relative_error,
        reduction_gamma=reduction_gamma,
        rational_midpoint=rational_midpoint,
        rational_error=rational_error,
        rational_lower=rational_midpoint - rational_error,
        rational_upper=rational_midpoint + rational_error,
        computed_absolute_term_sum_upper=computed_upper,
        term_evaluation_error=term_evaluation_error,
        reduction_error=reduction_error,
        logarithmic_sum=arb(0),
        logarithmic_cross_term=arb(0),
        slope_square_term=Fraction(0),
        prefix_energy=prefix_energy,
        last_q_numerator=certificate.P5_LAST_Q_NUMERATOR,
        maximum_q_numerator=maximum_q,
        maximum_square_numerator=maximum_q * maximum_q,
    )


def _fake_components(
    witness: certificate.FifthStepWitness,
    *,
    chunk_size: int = certificate.GENERATION_CHUNK_SIZE,
    reduction_leaf_size: int = (
        certificate.GENERATION_REDUCTION_LEAF_SIZE
    ),
    precision_bits: int = certificate.GENERATION_PRECISION_BITS,
) -> tuple[Any, ...]:
    old_prefix = _fake_prefix(
        chunk_size=chunk_size,
        reduction_leaf_size=reduction_leaf_size,
        prefix_energy=arb(2),
    )
    new_prefix = _fake_prefix(
        chunk_size=chunk_size,
        reduction_leaf_size=reduction_leaf_size,
        prefix_energy=arb(1),
    )
    tail = _fake_tail()
    new_upper, gain = certificate._combine_proof(
        old_prefix,
        new_prefix,
        tail,
        precision_bits,
    )
    return witness, old_prefix, new_prefix, tail, new_upper, gain


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


def _generated_fake_artifact(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> dict[str, object]:
    components = _fake_components(witness)
    monkeypatch.setattr(
        certificate,
        "_prepare_exact_data",
        lambda root: (witness, components[3]),
    )
    monkeypatch.setattr(
        certificate,
        "_generate_components",
        lambda *args, **kwargs: components,
    )
    monkeypatch.setattr(
        certificate,
        "_checks",
        lambda *args, **kwargs: dict(certificate.EXPECTED_CHECKS),
    )
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        None,
    )
    return certificate.generate(ROOT)


def test_configuration_pins_large_cutoff_and_distinct_replay_schedule() -> None:
    configuration = certificate._configuration_record()

    assert certificate.CUTOFF == 1 << 30
    assert certificate.CLAIMED_GAIN_LOWER_BOUND == Fraction(1, 300_000)
    assert configuration["cutoff"] == str(1 << 30)
    assert configuration["shell_5_multiplier_limit"] == "1"
    assert configuration["prefix_algorithm_version"] == (
        certificate.absolute_prefix.PREFIX_ALGORITHM_VERSION
    )
    assert configuration["generation_chunk_size"] == str(1 << 24)
    assert configuration["generation_reduction_leaf_size"] == str(1 << 20)
    assert configuration["generation_precision_bits"] == "768"
    assert configuration["replay_chunk_size"] == "16000003"
    assert configuration["replay_reduction_leaf_size"] == "1000003"
    assert configuration["replay_precision_bits"] == "896"
    assert certificate.GENERATION_CHUNK_SIZE != certificate.REPLAY_CHUNK_SIZE
    assert (
        certificate.GENERATION_REDUCTION_LEAF_SIZE
        != certificate.REPLAY_REDUCTION_LEAF_SIZE
    )
    assert (
        certificate.GENERATION_PRECISION_BITS
        != certificate.REPLAY_PRECISION_BITS
    )
    assert len(configuration["weights_5"]) == 34


def test_exact_witness_reconstruction_is_bound_to_all_frozen_commitments(
    witness: certificate.FifthStepWitness,
) -> None:
    parent_commitment = certificate.fourth.third.parent._vector_commitment(
        witness.parent_witness.p4
    )
    shell_commitment = certificate.fourth.third.parent._vector_commitment(
        witness.shell5
    )
    weight_commitment = certificate.fourth.third.parent._vector_commitment(
        witness.weights5
    )
    p5_commitment = certificate.fourth.third.parent._vector_commitment(
        witness.p5
    )

    assert parent_commitment["payload_sha256"] == (
        certificate.PARENT_P4_PAYLOAD_SHA256
    )
    assert shell_commitment == {
        "denominator_exponent": str(certificate.SHELL_BITS),
        "support_count": str(certificate.SHELL_5_SUPPORT_COUNT),
        "maximum_index": str(certificate.SHELL_5_MAXIMUM_INDEX),
        "payload_sha256": certificate.SHELL_5_PAYLOAD_SHA256,
    }
    assert weight_commitment == {
        "denominator_exponent": str(certificate.WEIGHT_BITS),
        "support_count": "34",
        "maximum_index": "34",
        "payload_sha256": certificate.WEIGHTS_5_PAYLOAD_SHA256,
    }
    assert p5_commitment == {
        "denominator_exponent": str(
            certificate.P5_DENOMINATOR_EXPONENT
        ),
        "support_count": str(certificate.P5_SUPPORT_COUNT),
        "maximum_index": str(certificate.P5_MAXIMUM_INDEX),
        "payload_sha256": certificate.P5_PAYLOAD_SHA256,
    }
    assert min(witness.shell5) == certificate.SHELL_5_MINIMUM_INDEX
    assert certificate._l1_norm(witness.shell5) == (
        certificate.SHELL_5_L1_NORM
    )
    assert certificate.fourth.third.parent._coefficient_sum(
        witness.shell5
    ) == 0
    assert tuple(
        witness.weights5[index] for index in range(1, 35)
    ) == certificate.WEIGHTS_5
    assert certificate.fourth.third.parent._coefficient_sum(
        witness.p5
    ) == certificate.P5_COEFFICIENT_SUM
    assert certificate._l1_norm(witness.p5) == certificate.P5_L1_NORM
    assert (
        1
        - certificate.P5_COEFFICIENT_SUM / 2
        == certificate.P5_C0
    )
    parent = witness.parent_witness.p4
    delta = certificate.fourth.third.parent._clean(
        {
            index: witness.p5.get(index, Fraction())
            - parent.get(index, Fraction())
            for index in set(witness.p5) | set(parent)
        }
    )
    assert {
        index for index in witness.p5 if index <= 65_536
    } == set(parent)
    changed = sum(witness.p5[index] != value for index, value in parent.items())
    assert changed == certificate.PARENT_CHANGED_COUNT
    assert len(parent) - changed == certificate.PARENT_UNCHANGED_COUNT
    assert len(delta) == certificate.DELTA_SUPPORT_COUNT
    assert certificate.fourth.third.parent._coefficient_sum(delta) == (
        certificate.DELTA_COEFFICIENT_SUM
    )


def test_integer_serialization_context_restores_limit_after_exception() -> None:
    original_limit = sys.get_int_max_str_digits()

    with pytest.raises(RuntimeError, match="sentinel"):
        with certificate._integer_serialization_context():
            raised_limit = sys.get_int_max_str_digits()
            assert raised_limit == 0 or (
                raised_limit >= certificate.MAX_INTEGER_DECIMAL_DIGITS
            )
            raise RuntimeError("sentinel")

    assert sys.get_int_max_str_digits() == original_limit


def test_self_check_reports_prefight_not_certificate_verification(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    tail = _fake_tail()
    tail_hash = content_sha256(certificate.fourth.third._tail_record(tail))
    monkeypatch.setattr(certificate, "_construct_witness", lambda root: witness)
    monkeypatch.setattr(
        certificate.fourth.third,
        "local_spacing_tail_upper",
        lambda coefficients, cutoff: tail,
    )
    monkeypatch.setattr(certificate, "TAIL_RECORD_SHA256", tail_hash)

    result = certificate.self_check(ROOT)

    assert result["mode"] == "EXACT_DATA_PREFLIGHT"
    assert result["self_check_passed"] is True
    assert result["certificate_verified"] is False
    assert result["prefix_evaluated"] is False
    assert result["hypothesis_status"] == "UNRESOLVED"
    assert "verified" not in result


def test_generate_components_wires_both_prefixes_to_one_declared_schedule(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    tail = _fake_tail()
    tail_hash = content_sha256(certificate.fourth.third._tail_record(tail))
    calls: list[tuple[object, int, int, int, int]] = []

    def fake_prefix(
        coefficients: object,
        cutoff: int,
        *,
        chunk_size: int,
        reduction_leaf_size: int,
        precision_bits: int,
    ) -> certificate.absolute_prefix.AbsoluteEnergyPrefix:
        calls.append(
            (
                coefficients,
                cutoff,
                chunk_size,
                reduction_leaf_size,
                precision_bits,
            )
        )
        return _fake_prefix(
            chunk_size=chunk_size,
            reduction_leaf_size=reduction_leaf_size,
            prefix_energy=arb(2 if len(calls) == 1 else 1),
        )

    monkeypatch.setattr(certificate, "_construct_witness", lambda root: witness)
    monkeypatch.setattr(
        certificate.absolute_prefix,
        "certified_absolute_energy_prefix",
        fake_prefix,
    )
    monkeypatch.setattr(
        certificate.fourth.third,
        "local_spacing_tail_upper",
        lambda coefficients, cutoff: tail,
    )
    monkeypatch.setattr(certificate, "TAIL_RECORD_SHA256", tail_hash)

    components = certificate._generate_components(
        witness,
        tail,
        chunk_size=17,
        reduction_leaf_size=7,
        precision_bits=192,
    )

    assert calls == [
        (witness.parent_witness.p4, certificate.CUTOFF, 17, 7, 192),
        (witness.p5, certificate.CUTOFF, 17, 7, 192),
    ]
    assert components[0] is witness
    assert components[3] is tail
    assert components[5].lower() > arb(1) / 300_000


def test_generate_builds_honest_schema_and_boolean_ledger_without_huge_prefix(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    artifact = _generated_fake_artifact(monkeypatch, witness)

    assert set(artifact) == certificate.TOP_LEVEL_FIELDS
    assert artifact["schema"] == certificate.SCHEMA
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["checks"] == certificate.EXPECTED_CHECKS
    assert all(type(value) is bool for value in artifact["checks"].values())
    assert artifact["checks"]["optimizer_used_as_evidence"] is False
    assert artifact["checks"]["resolves_rh"] is False
    assert artifact["payload_sha256"] == content_sha256(
        {
            key: value
            for key, value in artifact.items()
            if key != "payload_sha256"
        }
    )
    assert "one finite fifth" in artifact["limitation"].lower()
    assert "remains unresolved" in artifact["limitation"].lower()


def test_generate_fails_closed_when_any_check_does_not_close(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    components = _fake_components(witness)
    bad_checks = dict(certificate.EXPECTED_CHECKS)
    bad_checks["complete_gain_strictly_above_claim"] = False
    monkeypatch.setattr(
        certificate,
        "_prepare_exact_data",
        lambda root: (witness, components[3]),
    )
    monkeypatch.setattr(
        certificate,
        "_generate_components",
        lambda *args, **kwargs: components,
    )
    monkeypatch.setattr(
        certificate,
        "_checks",
        lambda *args, **kwargs: bad_checks,
    )

    with pytest.raises(
        ArithmeticError,
        match="certificate checks did not close",
    ):
        certificate.generate(ROOT)


def test_prefix_record_schema_and_arb_proof_records_fail_closed(
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    components = _fake_components(witness)
    proof = certificate._proof_record(*components[1:])
    monkeypatch.setattr(
        certificate,
        "TAIL_RECORD_SHA256",
        content_sha256(proof["new_tail_upper"]),
    )
    certificate._verify_stored_proof_shape(proof)

    stored = deepcopy(proof)
    exactness = stored["gain_lower_bound_computation"]["enclosure"][
        "is_exact"
    ]
    assert type(exactness) is bool
    stored["gain_lower_bound_computation"]["enclosure"]["is_exact"] = int(
        exactness
    )
    generated = deepcopy(stored)
    with pytest.raises(ValueError, match="invalid Arb exactness flag"):
        certificate._verify_stored_proof_shape(stored)


def test_generation_and_replay_enclosures_must_overlap(
    witness: certificate.FifthStepWitness,
) -> None:
    generation = (
        witness,
        _fake_prefix(
            chunk_size=certificate.GENERATION_CHUNK_SIZE,
            reduction_leaf_size=(
                certificate.GENERATION_REDUCTION_LEAF_SIZE
            ),
            prefix_energy=arb(2),
        ),
        _fake_prefix(
            chunk_size=certificate.GENERATION_CHUNK_SIZE,
            reduction_leaf_size=(
                certificate.GENERATION_REDUCTION_LEAF_SIZE
            ),
            prefix_energy=arb(1),
        ),
        _fake_tail(),
        arb(1),
        arb(1),
    )
    replay = (
        generation[0],
        _fake_prefix(
            chunk_size=certificate.REPLAY_CHUNK_SIZE,
            reduction_leaf_size=certificate.REPLAY_REDUCTION_LEAF_SIZE,
            prefix_energy=arb(2),
        ),
        _fake_prefix(
            chunk_size=certificate.REPLAY_CHUNK_SIZE,
            reduction_leaf_size=certificate.REPLAY_REDUCTION_LEAF_SIZE,
            prefix_energy=arb(1),
        ),
        generation[3],
        arb(1),
        arb(1),
    )
    certificate._require_distinct_layout_consistency(generation, replay)

    disagreeing_replay = (*replay[:-1], arb(5))
    with pytest.raises(
        ArithmeticError,
        match="generation and distinct-layout replay disagree",
    ):
        certificate._require_distinct_layout_consistency(
            generation,
            disagreeing_replay,
        )


def test_verify_refuses_before_reading_when_payload_is_not_frozen(
    tmp_path: Path,
) -> None:
    missing = tmp_path / "does-not-exist.json"
    assert certificate.FROZEN_ARTIFACT_PAYLOAD_SHA256 is None

    with pytest.raises(RuntimeError, match="payload is not frozen"):
        certificate.verify(ROOT, missing)


def test_verify_rejects_unhashed_mutation_before_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    artifact = _generated_fake_artifact(monkeypatch, witness)
    frozen = str(artifact["payload_sha256"])
    artifact["hypothesis_status"] = "RESOLVED"
    path = _write_artifact(tmp_path, artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    monkeypatch.setattr(
        certificate,
        "_generate_components",
        lambda *args, **kwargs: pytest.fail("huge prefix was reached"),
    )

    with pytest.raises(ValueError, match="payload SHA-256 mismatch"):
        certificate.verify(ROOT, path)


def test_verify_rejects_rehashed_semantic_mutation_before_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    artifact = _generated_fake_artifact(monkeypatch, witness)
    artifact["hypothesis_status"] = "RESOLVED"
    frozen = _rehash(artifact)
    path = _write_artifact(tmp_path, artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    monkeypatch.setattr(
        certificate,
        "_generate_components",
        lambda *args, **kwargs: pytest.fail("huge prefix was reached"),
    )

    with pytest.raises(ValueError, match="hypothesis_status changed"):
        certificate.verify(ROOT, path)


def test_verify_rejects_numeric_boolean_ledger_before_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    artifact = _generated_fake_artifact(monkeypatch, witness)
    artifact["checks"] = {
        name: int(value)
        for name, value in artifact["checks"].items()
    }
    frozen = _rehash(artifact)
    path = _write_artifact(tmp_path, artifact)
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )
    monkeypatch.setattr(
        certificate,
        "_generate_components",
        lambda *args, **kwargs: pytest.fail("huge prefix was reached"),
    )

    with pytest.raises(ValueError, match="ledger values are not booleans"):
        certificate.verify(ROOT, path)


def test_verify_rejects_duplicate_json_keys_before_replay(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    witness: certificate.FifthStepWitness,
) -> None:
    artifact = _generated_fake_artifact(monkeypatch, witness)
    frozen = str(artifact["payload_sha256"])
    raw = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    mutated = raw.replace(
        '  "schema":',
        '  "schema": "duplicate must fail",\n  "schema":',
        1,
    )
    path = tmp_path / "duplicate.json"
    path.write_text(mutated, encoding="utf-8", newline="\n")
    monkeypatch.setattr(
        certificate,
        "FROZEN_ARTIFACT_PAYLOAD_SHA256",
        frozen,
    )

    with pytest.raises(ValueError, match="duplicate JSON key 'schema'"):
        certificate.verify(ROOT, path)


def test_cli_self_check_prints_preflight_semantics_without_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    artifact = tmp_path / "must-not-be-written.json"
    result = {
        "mode": "EXACT_DATA_PREFLIGHT",
        "self_check_passed": True,
        "certificate_verified": False,
        "prefix_evaluated": False,
        "hypothesis_status": "UNRESOLVED",
    }
    monkeypatch.setattr(certificate, "self_check", lambda root: result)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_nyman_rebased_fifth_step_certificate.py",
            "--root",
            str(ROOT),
            "--artifact",
            str(artifact),
            "--self-check",
        ],
    )

    certificate.main()

    assert json.loads(capsys.readouterr().out) == result
    assert not artifact.exists()
