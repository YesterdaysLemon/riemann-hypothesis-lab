from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import pytest

from riemann_lab.artifacts import content_sha256, write_json
import riemann_lab.weil_search as search
import riemann_lab.weil_transition as transition
import riemann_lab.weil_transition_summary as summary


FROZEN_PLAN_PATH = Path("plans/weil-transition-q7-q9-v2.json")


def _with_hash(body: dict[str, Any]) -> dict[str, Any]:
    return {**body, "payload_sha256": content_sha256(body)}


def _ball(midpoint: int, *, exponent: int = -10, radius: int = 1) -> dict[str, Any]:
    return {
        "display": f"mock-dyadic-{midpoint}-times-2^{exponent}",
        "dyadic": {
            "mid_mantissa": str(midpoint),
            "mid_exponent": str(exponent),
            "radius_mantissa": str(radius),
            "radius_exponent": str(exponent - 8),
        },
        "is_exact": radius == 0,
    }


def _rump(midpoint: int, dimension: int) -> dict[str, Any]:
    eigenvalues = [
        {
            "index": str(index),
            "real": _ball(midpoint + 10_000 * index),
            "imaginary": _ball(0),
            "imaginary_contains_zero": True,
            "real_part_is_positive": midpoint + 10_000 * index > 0,
        }
        for index in range(dimension)
    ]
    positive = all(item["real_part_is_positive"] for item in eigenvalues)
    return {
        "algorithm": "flint-arb_mat-eig-rump",
        "classification": "POSITIVE_SPECTRUM" if positive else "INCONCLUSIVE",
        "ordering": "ascending by pairwise-separated real enclosures",
        "eigenvalues": eigenvalues,
        "smallest_eigenvalue": {
            "ordered_index": "0",
            "real": eigenvalues[0]["real"],
            "imaginary": eigenvalues[0]["imaginary"],
            "separated_from_all_others": True,
        },
    }


def _positive_artifact(
    contract: transition.WeilTransitionCellContract,
    rump_midpoint: int,
) -> dict[str, Any]:
    attempt = _with_hash(
        {
            "schema": search.ATTEMPT_SCHEMA,
            "classification": "EXPLORATORY",
            "kind": "SEARCH",
            "cell": contract.v1_cell.to_record(),
            "precision_bits": str(contract.policy.attempt_bits[0]),
            "limitation": search.EXPLORATORY_LIMITATION,
            "matrix_construction": {"classification": "COMPLETE"},
            "matrix_evidence": {"mocked": True},
            "positive_definiteness": {
                "classification": "POSITIVE_DEFINITE",
                "pivots": [{"value": _ball(1, exponent=0, radius=0)}],
            },
            "secondary_rump_spectrum": _rump(
                rump_midpoint, 2 * contract.degree + 1
            ),
            "candidate_generation": {
                "classification": "SKIPPED_LDLT_POSITIVE",
                "role": "untrusted-search-hint-only",
                "tested": [],
            },
            "decision": "FINITE_POSITIVE_CERTIFIED",
            "sequence": "0",
        }
    )
    return search._search_weil_cell(
        contract.v1_cell,
        contract.policy,
        existing_attempts=[attempt],
        plan_sha256=contract.to_record()["cell_contract_sha256"],
    )


def _negative_evaluation(witness: list[str], negative: bool) -> dict[str, Any]:
    if negative:
        return {
            "witness": witness,
            "classification": "NEGATIVE",
            "sign": "NEGATIVE",
            "upper_bound_is_negative": True,
            "lower_bound_is_positive": False,
        }
    return {
        "witness": witness,
        "classification": "INCONCLUSIVE",
        "sign": "INCONCLUSIVE",
        "upper_bound_is_negative": False,
        "lower_bound_is_positive": False,
    }


def _candidate_artifact(
    contract: transition.WeilTransitionCellContract,
    rump_midpoint: int,
) -> dict[str, Any]:
    witness = ["1"] + ["0"] * (2 * contract.degree)
    evaluation = _negative_evaluation(witness, True)
    discovery = _with_hash(
        {
            "schema": search.ATTEMPT_SCHEMA,
            "classification": "EXPLORATORY",
            "kind": "SEARCH",
            "cell": contract.v1_cell.to_record(),
            "precision_bits": str(contract.policy.attempt_bits[0]),
            "limitation": search.EXPLORATORY_LIMITATION,
            "matrix_construction": {"classification": "COMPLETE"},
            "matrix_evidence": {"mocked": True},
            "positive_definiteness": {"classification": "INCONCLUSIVE"},
            "secondary_rump_spectrum": _rump(
                rump_midpoint, 2 * contract.degree + 1
            ),
            "candidate_generation": {
                "classification": "COMPLETED",
                "role": "untrusted-search-hint-only",
                "tested": [
                    {
                        "candidate": {"witness": witness},
                        "evaluation": evaluation,
                    }
                ],
                "negative_candidate_index": "0",
            },
            "decision": "NEGATIVE_WITNESS_DISCOVERED",
            "sequence": "0",
        }
    )
    confirmation = _with_hash(
        {
            "schema": search.ATTEMPT_SCHEMA,
            "classification": "EXPLORATORY",
            "kind": "CONFIRMATION",
            "cell": contract.v1_cell.to_record(),
            "precision_bits": str(contract.policy.attempt_bits[1]),
            "limitation": search.EXPLORATORY_LIMITATION,
            "matrix_construction": {"classification": "COMPLETE"},
            "matrix_evidence": {"mocked": True},
            "witness_evaluation": evaluation,
            "precision_containment": {
                "all_higher_precision_components_contained": True,
                "first_failure": None,
            },
            "decision": "NEGATIVE_WITNESS_CONFIRMED",
            "sequence": "1",
        }
    )
    return search._search_weil_cell(
        contract.v1_cell,
        contract.policy,
        existing_attempts=[discovery, confirmation],
        plan_sha256=contract.to_record()["cell_contract_sha256"],
    )


def _dedicated_confirmation(
    contract: transition.WeilTransitionCellContract,
    artifact: dict[str, Any],
    *,
    confirmed: bool,
) -> dict[str, Any]:
    witness = artifact["quarantined_candidate"]["witness"]
    body = {
        "schema": "rh-lab/weil-transition-dedicated-confirmation/v2",
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "contract_id": contract.contract_id,
        "cell_contract_sha256": contract.to_record()["cell_contract_sha256"],
        "matrix_input_sha256": contract.matrix_input_sha256,
        "discovery_attempt_payload_sha256": artifact["attempts"][0][
            "payload_sha256"
        ],
        "witness": witness,
        "precision_bits": str(contract.policy.confirmation_bits),
        "limitation": transition.EXPLORATORY_LIMITATION,
        "mode": "DEDICATED_REGENERATION",
        "matrix_construction": {"classification": "COMPLETE"},
        "matrix_evidence": {"mocked": True},
        "witness_evaluation": _negative_evaluation(witness, confirmed),
        "precision_containment": {
            "all_higher_precision_components_contained": True,
            "first_failure": None,
        },
        "decision": (
            "V2_NEGATIVE_CANDIDATE_CONFIRMED"
            if confirmed
            else "V2_DEDICATED_CONFIRMATION_INCONCLUSIVE"
        ),
    }
    return _with_hash(body)


def _mocked_checkpoint(
    root: Path,
    *,
    candidate_indices: dict[int, bool] | None = None,
    rump_failure_indices: set[int] | None = None,
) -> tuple[dict[str, Any], list[transition.WeilTransitionCellContract]]:
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    contracts = transition._contracts_from_plan(plan)
    candidate_indices = candidate_indices or {}
    rump_failure_indices = rump_failure_indices or set()
    write_json(root / "plan.json", plan)
    write_json(root / "run.json", transition._run_manifest(plan))

    completed: dict[str, tuple[dict[str, Any], dict[str, Any] | None]] = {}
    for ordinal, contract in enumerate(contracts):
        rump_midpoint = -500 if ordinal == 17 else 1000 + ordinal
        if ordinal in candidate_indices:
            artifact = _candidate_artifact(contract, rump_midpoint)
            confirmation = _dedicated_confirmation(
                contract, artifact, confirmed=candidate_indices[ordinal]
            )
            write_json(
                root / "confirmations" / f"{contract.contract_id}.json",
                confirmation,
            )
        else:
            artifact = _positive_artifact(contract, rump_midpoint)
            confirmation = None
        if ordinal in rump_failure_indices:
            attempt = artifact["attempts"][0]
            attempt["secondary_rump_spectrum"] = {
                "classification": "INCONCLUSIVE",
                "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
            }
            attempt["payload_sha256"] = content_sha256(
                {
                    key: value
                    for key, value in attempt.items()
                    if key != "payload_sha256"
                }
            )
            artifact["payload_sha256"] = content_sha256(
                {
                    key: value
                    for key, value in artifact.items()
                    if key != "payload_sha256"
                }
            )
        write_json(root / "cells" / f"{contract.contract_id}.json", artifact)
        completed[contract.contract_id] = (artifact, confirmation)

    index = transition._index_payload(plan, contracts, completed)
    write_json(root / "index.json", index)
    return index, contracts


def test_summary_is_compact_deterministic_and_does_no_numerical_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    index, contracts = _mocked_checkpoint(
        tmp_path, candidate_indices={0: True, 1: False}
    )

    def forbidden(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("numerical replay was attempted")

    monkeypatch.setattr(search, "_search_attempt", forbidden)
    monkeypatch.setattr(search, "_confirmation_attempt", forbidden)
    monkeypatch.setattr(transition, "_dedicated_confirmation", forbidden)

    artifact = summary.generate_weil_transition_summary(tmp_path)
    repeated = summary.generate_weil_transition_summary(tmp_path)
    assert artifact == repeated
    assert artifact["schema"] == summary.SUMMARY_SCHEMA
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["source"]["index"]["payload_sha256"] == index[
        "payload_sha256"
    ]
    assert artifact["source"]["plan"]["plan_sha256"] == transition.load_transition_plan(
        FROZEN_PLAN_PATH
    )["plan_sha256"]
    assert artifact["payload_sha256"] == content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )

    counts = artifact["counts"]
    assert counts["cells"] == "81"
    assert counts["by_status"] == {
        "FINITE_POSITIVE_CERTIFIED": "79",
        "INCONCLUSIVE_MAX_PRECISION": "0",
        "NEGATIVE_CANDIDATE_QUARANTINED": "1",
        "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE": "1",
    }
    assert counts["by_q"] == {"7": "27", "8": "27", "9": "27"}
    assert counts["by_position"] == {"above": "36", "below": "36", "exact": "9"}
    assert counts["by_degree"] == {"12": "9", "16": "27", "20": "27", "24": "18"}
    assert counts["by_terminal_precision_bits"] == {
        "192": "34",
        "384": "45",
        "768": "1",
        "1536": "1",
    }
    assert counts["total_v1_attempts"] == "83"
    assert counts["stored_candidate_evaluation_records"] == {
        "definition": (
            "hash-bound stored records counted from candidate-generation tested "
            "arrays or confirmation fields; this summary neither replays their "
            "numerical values nor validates that recorded witness entries are "
            "canonical integers"
        ),
        "v1_search": "2",
        "v1_confirmation": "2",
        "v2_dedicated_confirmation": "2",
        "total": "6",
    }

    quarantines = artifact["quarantine_references"]
    assert [item["contract_id"] for item in quarantines["negative_candidates"]] == [
        contracts[0].contract_id
    ]
    assert [
        item["contract_id"]
        for item in quarantines["conflicting_negative_observations"]
    ] == [contracts[1].contract_id]
    rump = artifact["rump_diagnostic"]
    assert rump["role"] == "DIAGNOSTIC_ONLY"
    assert rump["coverage"] == {
        "search_attempts": "81",
        "usable_spectra": "81",
        "usable_eigenvalue_enclosures": "3105",
        "pairwise_separated_usable_spectra": "81",
        "unseparated_usable_spectra": "0",
        "failure_sentinels": "0",
        "matrix_construction_inconclusive_without_spectrum": "0",
        "failure_references": [],
        "unseparated_references": [],
        "matrix_construction_inconclusive_references": [],
    }
    least = rump["least_lower_endpoint_among_usable_stored_enclosures"]
    assert least["contract_id"] == contracts[17].contract_id
    assert least["stored_eigenvalue_index"] == "0"
    assert least["real_enclosure"] == summary._normalized_enclosure(_ball(-500))
    assert "display" not in json.dumps(least["real_enclosure"])
    assert artifact["summary_method"]["numerical_replay_performed"] is False

    verification = summary.verify_weil_transition_summary(artifact, tmp_path)
    assert verification["classification"] == "REPRODUCED_EXPLORATORY_SUMMARY"
    assert verification["verified_cells"] == "81"
    assert verification["numerical_replay_performed"] is False


def test_verifier_rejects_rehashed_summary_tampering(tmp_path: Path) -> None:
    _mocked_checkpoint(tmp_path)
    artifact = summary.generate_weil_transition_summary(tmp_path)
    tampered = copy.deepcopy(artifact)
    tampered["counts"]["cells"] = "80"
    tampered_body = {
        key: value for key, value in tampered.items() if key != "payload_sha256"
    }
    tampered["payload_sha256"] = content_sha256(tampered_body)

    with pytest.raises(
        summary.WeilTransitionSummaryVerificationError,
        match="does not canonically regenerate",
    ):
        summary.verify_weil_transition_summary(tampered, tmp_path)


def test_generator_rejects_an_incomplete_rehashed_index(tmp_path: Path) -> None:
    index, _ = _mocked_checkpoint(tmp_path)
    damaged = copy.deepcopy(index)
    damaged["cells"].pop()
    damaged["payload_sha256"] = content_sha256(
        {key: value for key, value in damaged.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", damaged)

    with pytest.raises(summary.WeilTransitionSummaryError, match="81 index references"):
        summary.generate_weil_transition_summary(tmp_path)


def test_generator_rejects_rehashed_noncanonical_cell_state(tmp_path: Path) -> None:
    index, _ = _mocked_checkpoint(tmp_path)
    reference = index["cells"][0]
    cell_path = tmp_path / reference["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))
    cell["status"] = "INCONCLUSIVE_MAX_PRECISION"
    cell["payload_sha256"] = content_sha256(
        {key: value for key, value in cell.items() if key != "payload_sha256"}
    )
    write_json(cell_path, cell)
    reference["v1_cell_payload_sha256"] = cell["payload_sha256"]
    reference["v1_status"] = cell["status"]
    reference["status"] = cell["status"]
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    with pytest.raises(summary.WeilTransitionSummaryError, match="not canonical"):
        summary.generate_weil_transition_summary(tmp_path)


def test_rump_failure_is_visible_in_coverage_but_malformed_failure_is_rejected(
    tmp_path: Path,
) -> None:
    frozen_contracts = transition._contracts_from_plan(
        transition.load_transition_plan(FROZEN_PLAN_PATH)
    )
    live_failure_index = next(
        index
        for index, contract in enumerate(frozen_contracts)
        if contract.contract_id == "q-9-j-4-above-n-16"
    )
    index, contracts = _mocked_checkpoint(
        tmp_path, rump_failure_indices={live_failure_index}
    )

    artifact = summary.generate_weil_transition_summary(tmp_path)
    diagnostic = artifact["rump_diagnostic"]
    coverage = diagnostic["coverage"]
    assert coverage["search_attempts"] == "81"
    assert coverage["usable_spectra"] == "80"
    assert coverage["usable_eigenvalue_enclosures"] == str(
        sum(
            2 * contract.degree + 1
            for ordinal, contract in enumerate(contracts)
            if ordinal != live_failure_index
        )
    )
    assert coverage["failure_sentinels"] == "1"
    assert coverage["unseparated_usable_spectra"] == "0"
    assert coverage["failure_references"] == [
        {
            "contract_id": contracts[live_failure_index].contract_id,
            "cell_payload_sha256": index["cells"][live_failure_index][
                "v1_cell_payload_sha256"
            ],
            "attempt_payload_sha256": json.loads(
                (tmp_path / index["cells"][live_failure_index]["path"]).read_text(
                    encoding="utf-8"
                )
            )["attempts"][0]["payload_sha256"],
            "attempt_sequence": "0",
            "precision_bits": str(
                contracts[live_failure_index].policy.attempt_bits[0]
            ),
            "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
        }
    ]
    assert diagnostic[
        "least_lower_endpoint_among_usable_stored_enclosures"
    ]["contract_id"] == contracts[17].contract_id

    reference = index["cells"][live_failure_index]
    cell_path = tmp_path / reference["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))
    attempt = cell["attempts"][0]
    attempt["secondary_rump_spectrum"]["failure_code"] = "UNRECOGNIZED_FAILURE"
    attempt["payload_sha256"] = content_sha256(
        {key: value for key, value in attempt.items() if key != "payload_sha256"}
    )
    cell["payload_sha256"] = content_sha256(
        {key: value for key, value in cell.items() if key != "payload_sha256"}
    )
    write_json(cell_path, cell)
    reference["v1_cell_payload_sha256"] = cell["payload_sha256"]
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    with pytest.raises(summary.WeilTransitionSummaryError, match="failure record"):
        summary.generate_weil_transition_summary(tmp_path)


def test_rump_selection_scans_all_enclosures_in_an_unseparated_spectrum(
    tmp_path: Path,
) -> None:
    index, contracts = _mocked_checkpoint(tmp_path)
    contract = contracts[0]
    reference = index["cells"][0]
    cell_path = tmp_path / reference["path"]
    artifact = json.loads(cell_path.read_text(encoding="utf-8"))
    attempt = artifact["attempts"][0]
    spectrum = attempt["secondary_rump_spectrum"]
    eigenvalues = spectrum["eigenvalues"]

    eigenvalues[0]["real"] = _ball(5000)
    eigenvalues[1]["real"] = _ball(-500)
    eigenvalues[1]["real_part_is_positive"] = False
    spectrum["classification"] = "INCONCLUSIVE"
    spectrum["smallest_eigenvalue"]["real"] = eigenvalues[0]["real"]
    spectrum["smallest_eigenvalue"]["separated_from_all_others"] = False
    attempt["payload_sha256"] = content_sha256(
        {key: value for key, value in attempt.items() if key != "payload_sha256"}
    )
    artifact["payload_sha256"] = content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )
    write_json(cell_path, artifact)
    reference["v1_cell_payload_sha256"] = artifact["payload_sha256"]
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    generated = summary.generate_weil_transition_summary(tmp_path)
    diagnostic = generated["rump_diagnostic"]
    selected = diagnostic[
        "least_lower_endpoint_among_usable_stored_enclosures"
    ]
    assert selected["contract_id"] == contract.contract_id
    assert selected["stored_eigenvalue_index"] == "1"
    assert selected["real_enclosure"] == summary._normalized_enclosure(
        eigenvalues[1]["real"]
    )
    coverage = diagnostic["coverage"]
    assert coverage["usable_spectra"] == "81"
    assert coverage["pairwise_separated_usable_spectra"] == "80"
    assert coverage["unseparated_usable_spectra"] == "1"
    assert coverage["unseparated_references"] == [
        {
            "contract_id": contract.contract_id,
            "cell_payload_sha256": artifact["payload_sha256"],
            "attempt_payload_sha256": attempt["payload_sha256"],
            "attempt_sequence": "0",
            "precision_bits": str(contract.policy.attempt_bits[0]),
            "spectrum_classification": "INCONCLUSIVE",
        }
    ]


def test_generator_rejects_rehashed_unknown_attempt_field(tmp_path: Path) -> None:
    index, _ = _mocked_checkpoint(tmp_path)
    reference = index["cells"][0]
    cell_path = tmp_path / reference["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))
    attempt = cell["attempts"][0]
    attempt["unknown_attempt_field"] = "must-not-be-accepted"
    attempt["payload_sha256"] = content_sha256(
        {key: value for key, value in attempt.items() if key != "payload_sha256"}
    )
    cell["payload_sha256"] = content_sha256(
        {key: value for key, value in cell.items() if key != "payload_sha256"}
    )
    write_json(cell_path, cell)
    reference["v1_cell_payload_sha256"] = cell["payload_sha256"]
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    with pytest.raises(
        summary.WeilTransitionSummaryError,
        match="attempt has noncanonical fields",
    ):
        summary.generate_weil_transition_summary(tmp_path)


@pytest.mark.parametrize("mutation", ("missing", "malformed"))
def test_generator_reports_dedicated_confirmation_reference_errors(
    tmp_path: Path, mutation: str
) -> None:
    index, _ = _mocked_checkpoint(tmp_path)
    reference = index["cells"][0]
    if mutation == "missing":
        del reference["dedicated_confirmation"]
    else:
        reference["dedicated_confirmation"] = "not-an-object"
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    with pytest.raises(
        summary.WeilTransitionSummaryError,
        match="dedicated_confirmation|confirmation reference",
    ):
        summary.generate_weil_transition_summary(tmp_path)


@pytest.mark.parametrize(
    "mutation",
    (
        "missing_smallest",
        "nonobject_spectrum",
        "false_exact_flag",
        "false_positive_spectrum",
    ),
)
def test_rump_adversarial_records_are_rejected(mutation: str) -> None:
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    contract = transition._contracts_from_plan(plan)[0]
    artifact = _positive_artifact(contract, 1000)
    attempt = copy.deepcopy(artifact["attempts"][0])
    spectrum = attempt["secondary_rump_spectrum"]

    if mutation == "missing_smallest":
        del spectrum["smallest_eigenvalue"]
    elif mutation == "nonobject_spectrum":
        attempt["secondary_rump_spectrum"] = None
    elif mutation == "false_exact_flag":
        spectrum["smallest_eigenvalue"]["real"]["is_exact"] = True
        spectrum["eigenvalues"][0]["real"]["is_exact"] = True
    else:
        nonpositive = _ball(0)
        spectrum["smallest_eigenvalue"]["real"] = nonpositive
        spectrum["eigenvalues"][0]["real"] = nonpositive
        spectrum["eigenvalues"][0]["real_part_is_positive"] = False
        assert spectrum["classification"] == "POSITIVE_SPECTRUM"

    with pytest.raises(summary.WeilTransitionSummaryError, match="Rump"):
        summary._consider_rump(
            None,
            contract=contract,
            artifact=artifact,
            attempt=attempt,
        )


def test_rump_display_text_is_not_propagated() -> None:
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    contract = transition._contracts_from_plan(plan)[0]
    artifact = _positive_artifact(contract, 1000)
    attempt = copy.deepcopy(artifact["attempts"][0])
    spectrum = attempt["secondary_rump_spectrum"]
    spectrum["smallest_eigenvalue"]["real"]["display"] = "unchecked text"
    spectrum["eigenvalues"][0]["real"]["display"] = "unchecked text"

    selected, observation = summary._consider_rump(
        None,
        contract=contract,
        artifact=artifact,
        attempt=attempt,
    )
    assert selected is not None
    assert observation["kind"] == "USABLE_SPECTRUM"
    assert "unchecked text" not in json.dumps(selected[1])
    assert "display" not in selected[1]["real_enclosure"]


def test_generator_rejects_rehashed_frozen_plan_with_bogus_backend(
    tmp_path: Path,
) -> None:
    index, _ = _mocked_checkpoint(tmp_path)
    raw_plan = json.loads((tmp_path / "plan.json").read_text(encoding="utf-8"))
    raw_plan.pop("plan_sha256")
    raw_plan["backend_contract"] = {
        "python_flint": "999.0.0",
        "flint": "999.0.0",
    }
    bogus_plan = transition.canonicalize_transition_plan(raw_plan)
    write_json(tmp_path / "plan.json", bogus_plan)
    write_json(tmp_path / "run.json", transition._run_manifest(bogus_plan))
    index["plan"]["plan_sha256"] = bogus_plan["plan_sha256"]
    index["payload_sha256"] = content_sha256(
        {key: value for key, value in index.items() if key != "payload_sha256"}
    )
    write_json(tmp_path / "index.json", index)

    with pytest.raises(summary.WeilTransitionSummaryError, match="plan hash changed"):
        summary.generate_weil_transition_summary(tmp_path)
