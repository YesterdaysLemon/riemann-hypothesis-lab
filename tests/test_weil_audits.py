from __future__ import annotations

import copy
from fractions import Fraction
import json
from typing import Any, Callable, Mapping

from flint import arb
import pytest

import riemann_lab.weil_audits as weil_audits
from riemann_lab.artifacts import content_sha256
from riemann_lab.weil_audits import (
    NESTING_AUDIT_SCHEMA,
    PARITY_AUDIT_SCHEMA,
    WeilAuditVerificationError,
    audit_zero_padding_negative_witness,
    certify_degree_nesting_audit,
    certify_parity_audit,
    verify_degree_nesting_audit,
    verify_parity_audit,
)
from riemann_lab.weil_search import WeilSearchCell, WeilSearchPolicy


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {key: value for key, value in record.items() if key != "payload_sha256"}
    record["payload_sha256"] = content_sha256(body)
    return record


def _diagonal(values: list[Any]) -> list[list[Any]]:
    return [
        [value if row == column else arb(0) for column in range(len(values))]
        for row, value in enumerate(values)
    ]


@pytest.fixture(scope="module")
def audit_policy() -> WeilSearchPolicy:
    return WeilSearchPolicy(
        attempt_bits=(96,),
        confirmation_bits=192,
        eigenpair_count=1,
        scale_bits=(4,),
    )


@pytest.fixture(scope="module")
def parity_artifact(audit_policy: WeilSearchPolicy) -> dict[str, Any]:
    return certify_parity_audit(
        WeilSearchCell(3, 2, 1),
        96,
        audit_policy,
        witnesses=((2, 4, 2), (1, 0, -1)),
    )


@pytest.fixture(scope="module")
def nesting_artifact() -> dict[str, Any]:
    return certify_degree_nesting_audit(
        WeilSearchCell(3, 2, 0),
        WeilSearchCell(3, 2, 4),
        96,
        160,
        negative_witnesses=((2,),),
    )


def test_parity_generation_has_explicit_orthonormal_blocks_and_zero_cross_gate(
    parity_artifact: dict[str, Any],
    audit_policy: WeilSearchPolicy,
) -> None:
    artifact = parity_artifact
    assert artifact["schema"] == PARITY_AUDIT_SCHEMA
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["policy"] == audit_policy.to_record()
    assert artifact["candidate_generation"]["policy"] == audit_policy.to_record()
    assert artifact["basis"]["even_dimension"] == "2"
    assert artifact["basis"]["odd_dimension"] == "1"
    assert artifact["blocks"]["even"]["dimension"] == "2"
    assert artifact["blocks"]["odd"]["dimension"] == "1"
    assert len(artifact["blocks"]["even"]["entries_upper_triangle"]) == 3
    assert len(artifact["blocks"]["odd"]["entries_upper_triangle"]) == 1

    gates = artifact["gates"]
    assert gates["centrosymmetry_overlap"]["classification"] == "PASS"
    assert gates["orthonormal_basis"]["classification"] == "PASS"
    cross = gates["cross_parity_contains_zero"]
    assert cross["classification"] == "PASS"
    assert cross["checked_entries"] == "2"
    assert all(entry["contains_zero"] for entry in cross["entries"])

    requested = artifact["requested_witness_audits"]
    assert [item["classification"] for item in requested] == ["PASS", "PASS"]
    assert requested[0]["supplied_witness_was_canonical"] is False
    assert requested[0]["witness"] == ["1", "2", "1"]
    assert requested[1]["witness"] == ["1", "0", "-1"]
    assert all(
        item["authoritative_full_matrix_evaluation"]["sign"]
        in {"POSITIVE", "INCONCLUSIVE"}
        for item in requested
    )
    assert artifact["payload_sha256"] == content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )


def test_parity_replay_regenerates_structural_evidence_and_exact_witnesses(
    parity_artifact: dict[str, Any],
) -> None:
    replay = verify_parity_audit(parity_artifact)
    candidate_count = int(parity_artifact["candidate_generation"]["candidate_count"])
    assert replay == {
        "classification": "REPRODUCED_EXPLORATORY_AUDIT",
        "audit_kind": "PARITY",
        "decision": "AUDIT_PASSED_EXPLORATORY",
        "verified_integer_witnesses": str(candidate_count + 2),
        "verified_negative_witnesses": "0",
        "same_backend_replay_only": True,
        "hypothesis_status": "UNRESOLVED",
    }


def test_degree_nesting_supports_arbitrary_gap_and_canonical_replay(
    nesting_artifact: dict[str, Any],
) -> None:
    artifact = nesting_artifact
    assert artifact["schema"] == NESTING_AUDIT_SCHEMA
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["lower_cell"]["degree"] == "0"
    assert artifact["higher_cell"]["degree"] == "4"
    relation = artifact["principal_submatrix_relation"]
    assert relation["classification"] == "PASS"
    assert "M>N" in relation["description"]
    assert all(item["classification"] == "PASS" for item in relation["components"])
    assert {item["component"] for item in relation["components"]} == {
        "P",
        "R",
        "S",
        "A",
    }

    witness = artifact["zero_padding_negative_witness_checks"][0]
    assert witness["classification"] == "NOT_APPLICABLE_LOWER_WITNESS_POSITIVE"
    assert witness["base_precision"]["degree_gap"] == "4"
    assert witness["base_precision"]["supplied_witness"] == ["2"]
    assert witness["base_precision"]["witness"] == ["1"]
    assert witness["base_precision"]["zero_padded_witness"] == [
        "0",
        "0",
        "0",
        "0",
        "1",
        "0",
        "0",
        "0",
        "0",
    ]

    replay = verify_degree_nesting_audit(artifact)
    assert replay["classification"] == "REPRODUCED_EXPLORATORY_AUDIT"
    assert replay["audit_kind"] == "DEGREE_NESTING"
    assert replay["decision"] == "AUDIT_PASSED_EXPLORATORY"
    assert replay["verified_zero_padding_witnesses"] == "1"
    assert replay["hypothesis_status"] == "UNRESOLVED"


def test_sequential_nesting_then_parity_replay_is_backend_isolated(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend_state = {"dirty": True}
    cleanup_calls: list[bool] = []

    def cleanup() -> None:
        cleanup_calls.append(backend_state["dirty"])
        backend_state["dirty"] = False

    def nesting_replay(artifact: Mapping[str, Any]) -> dict[str, Any]:
        assert artifact == {"kind": "nesting"}
        assert backend_state["dirty"] is False
        backend_state["dirty"] = True
        return {"audit_kind": "DEGREE_NESTING"}

    def parity_replay(artifact: Mapping[str, Any]) -> dict[str, Any]:
        assert artifact == {"kind": "parity"}
        assert backend_state["dirty"] is False
        backend_state["dirty"] = True
        return {"audit_kind": "PARITY"}

    monkeypatch.setattr(weil_audits, "_cleanup_flint_backend", cleanup)
    monkeypatch.setattr(
        weil_audits, "_verify_degree_nesting_audit", nesting_replay
    )
    monkeypatch.setattr(weil_audits, "_verify_parity_audit", parity_replay)

    assert verify_degree_nesting_audit({"kind": "nesting"}) == {
        "audit_kind": "DEGREE_NESTING"
    }
    assert verify_parity_audit({"kind": "parity"}) == {
        "audit_kind": "PARITY"
    }
    assert cleanup_calls == [True, True, False, True]
    assert backend_state["dirty"] is False


def test_parity_generation_uses_the_same_clean_backend_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    backend_state = {"dirty": True}
    cleanup_calls: list[bool] = []
    cell = WeilSearchCell(3, 2, 1)
    policy = WeilSearchPolicy(
        attempt_bits=(96,),
        confirmation_bits=192,
        eigenpair_count=1,
        scale_bits=(4,),
    )

    def cleanup() -> None:
        cleanup_calls.append(backend_state["dirty"])
        backend_state["dirty"] = False

    def generate(
        supplied_cell: WeilSearchCell,
        precision_bits: int,
        supplied_policy: WeilSearchPolicy,
        witnesses: tuple[tuple[int, ...], ...],
    ) -> dict[str, Any]:
        assert (supplied_cell, precision_bits, supplied_policy, witnesses) == (
            cell,
            96,
            policy,
            (),
        )
        assert backend_state["dirty"] is False
        backend_state["dirty"] = True
        return {"decision": "AUDIT_PASSED_EXPLORATORY"}

    monkeypatch.setattr(weil_audits, "_cleanup_flint_backend", cleanup)
    monkeypatch.setattr(weil_audits, "_certify_parity_audit", generate)

    assert certify_parity_audit(cell, 96, policy) == {
        "decision": "AUDIT_PASSED_EXPLORATORY"
    }
    assert cleanup_calls == [True, True]
    assert backend_state["dirty"] is False


def test_small_fixture_nesting_then_parity_replays_canonically(
    nesting_artifact: dict[str, Any], parity_artifact: dict[str, Any]
) -> None:
    assert verify_degree_nesting_audit(nesting_artifact)["decision"] == (
        "AUDIT_PASSED_EXPLORATORY"
    )
    assert verify_parity_audit(parity_artifact)["decision"] == (
        "AUDIT_PASSED_EXPLORATORY"
    )


def test_zero_padding_positive_and_inconclusive_states_are_not_false_failures() -> None:
    positive_lower = _diagonal([arb(2), arb(1), arb(2)])
    positive_higher = _diagonal(
        [arb(5), arb(5), arb(2), arb(1), arb(2), arb(5), arb(5)]
    )
    positive = audit_zero_padding_negative_witness(
        positive_lower, positive_higher, (0, 2, 0), 96
    )
    assert positive["classification"] == "NOT_APPLICABLE_LOWER_WITNESS_POSITIVE"
    assert positive["audit_failure"] is False
    assert positive["audit_inconclusive"] is False
    assert positive["witness"] == ["0", "1", "0"]

    negative_lower = _diagonal([arb(2), arb(-1), arb(2)])
    inconclusive_higher = _diagonal(
        [arb(5), arb(5), arb(2), arb("0 +/- 2"), arb(2), arb(5), arb(5)]
    )
    inconclusive = audit_zero_padding_negative_witness(
        negative_lower, inconclusive_higher, (0, 1, 0), 96
    )
    assert inconclusive["classification"] == "INCONCLUSIVE_HIGHER_DEGREE_SIGN"
    assert inconclusive["lower_degree_evaluation"]["sign"] == "NEGATIVE"
    assert inconclusive["higher_degree_evaluation"]["sign"] == "INCONCLUSIVE"
    assert inconclusive["audit_failure"] is False
    assert inconclusive["audit_inconclusive"] is True
    assert inconclusive["lower_and_zero_padded_numerators_overlap"] is True


def test_unmodified_hash_is_required_before_semantic_replay(
    parity_artifact: dict[str, Any],
) -> None:
    forged = copy.deepcopy(parity_artifact)
    forged["decision"] = "AUDIT_FAILED"
    with pytest.raises(WeilAuditVerificationError, match="payload hash mismatch"):
        verify_parity_audit(forged)


def _mutate_schema(record: dict[str, Any]) -> None:
    record["schema"] = "rh-lab/weil-parity-audit/v999"


def _mutate_status(record: dict[str, Any]) -> None:
    record["hypothesis_status"] = "PROVED"


def _mutate_classification(record: dict[str, Any]) -> None:
    record["classification"] = "CERTIFIED"


def _mutate_cell(record: dict[str, Any]) -> None:
    record["cell"]["dimension"] = "999"


def _mutate_policy(record: dict[str, Any]) -> None:
    record["policy"]["witness_policy"]["scale_bits"] = ["4", "4"]


def _mutate_precision(record: dict[str, Any]) -> None:
    record["precision_bits"] = "096"


@pytest.mark.parametrize(
    "mutation",
    [
        _mutate_schema,
        _mutate_status,
        _mutate_classification,
        _mutate_cell,
        _mutate_policy,
        _mutate_precision,
    ],
    ids=["schema", "status", "classification", "cell", "policy", "precision"],
)
def test_rehashed_contract_tampering_is_rejected(
    parity_artifact: dict[str, Any],
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    forged = copy.deepcopy(parity_artifact)
    mutation(forged)
    _rehash(forged)
    with pytest.raises(WeilAuditVerificationError):
        verify_parity_audit(forged)


def test_rehashed_matrix_and_direct_witness_evidence_tampering_is_rejected(
    parity_artifact: dict[str, Any],
) -> None:
    block_forgery = copy.deepcopy(parity_artifact)
    block_forgery["blocks"]["even"]["entries_upper_triangle"][0]["row"] = "9"
    _rehash(block_forgery)
    with pytest.raises(WeilAuditVerificationError, match="blocks changed"):
        verify_parity_audit(block_forgery)

    witness_forgery = copy.deepcopy(parity_artifact)
    candidate = witness_forgery["candidate_generation"]["candidates"][0]
    candidate["authoritative_full_matrix_evaluation"]["sign"] = "NEGATIVE"
    _rehash(witness_forgery)
    with pytest.raises(WeilAuditVerificationError, match="direct evaluation"):
        verify_parity_audit(witness_forgery)


def test_approximate_candidate_origin_is_not_a_replay_dependency(
    parity_artifact: dict[str, Any],
) -> None:
    altered = copy.deepcopy(parity_artifact)
    candidate = altered["candidate_generation"]["candidates"][0]
    midpoint = candidate["source"]["eigenvalue_real_midpoint"]
    changed = Fraction(int(midpoint["numerator"]), int(midpoint["denominator"])) + 1
    midpoint["numerator"] = str(changed.numerator)
    midpoint["denominator"] = str(changed.denominator)
    identity = {
        "source": candidate["source"],
        "witness": candidate["witness"],
    }
    candidate["candidate_id"] = content_sha256(identity)
    altered["candidate_generation"]["negative_candidate_ids"] = []
    _rehash(altered)

    replay = verify_parity_audit(altered)
    assert replay["classification"] == "REPRODUCED_EXPLORATORY_AUDIT"
    assert replay["verified_negative_witnesses"] == "0"


def test_rehashed_nesting_evidence_tampering_is_rejected(
    nesting_artifact: dict[str, Any],
) -> None:
    forged = copy.deepcopy(nesting_artifact)
    forged["principal_submatrix_relation"]["components"][0][
        "classification"
    ] = "FAIL"
    _rehash(forged)
    with pytest.raises(
        WeilAuditVerificationError, match="does not canonically regenerate"
    ):
        verify_degree_nesting_audit(forged)


@pytest.mark.parametrize(
    "verifier",
    [verify_parity_audit, verify_degree_nesting_audit],
)
def test_verifiers_wrap_all_malformed_inputs_in_dedicated_error(
    verifier: Callable[[Any], dict[str, Any]],
) -> None:
    with pytest.raises(WeilAuditVerificationError):
        verifier([])
    malformed = _rehash({"schema": None})
    with pytest.raises(WeilAuditVerificationError):
        verifier(malformed)


def test_audits_never_emit_top_level_resolution_claims(
    parity_artifact: dict[str, Any], nesting_artifact: dict[str, Any]
) -> None:
    for artifact in (parity_artifact, nesting_artifact):
        assert artifact["classification"] == "EXPLORATORY"
        assert artifact["hypothesis_status"] == "UNRESOLVED"
        assert artifact["decision"] in {
            "AUDIT_PASSED_EXPLORATORY",
            "AUDIT_INCONCLUSIVE",
            "AUDIT_FAILED",
            "NEGATIVE_CANDIDATE_QUARANTINED",
        }
        encoded = json.dumps(artifact, sort_keys=True)
        assert '"PROVED"' not in encoded
        assert '"DISPROVED"' not in encoded
