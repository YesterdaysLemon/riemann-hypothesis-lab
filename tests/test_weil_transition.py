from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx
import pytest

from riemann_lab.artifacts import canonical_json_file_sha256, content_sha256
from riemann_lab.weil import prime_power_entry, prime_powers_leq
import riemann_lab.weil_search as search
import riemann_lab.weil_transition as transition


FROZEN_PLAN_PATH = Path("plans/weil-transition-q7-q9-v2.json")
FROZEN_PLAN_SHA256 = (
    "33185881cc619fda99b7835f5e71422b0c5cf745cb5411162638a8e08dd05336"
)
FROZEN_PLAN_FILE_SHA256 = (
    "87b1d8b5ff10b77107b2bced652e6c9149f6e034b8b51acfc7a0d01ee63907b1"
)

EXPECTED_CUTOFFS_BY_Q = {
    7: (
        (112, 17),
        (448, 65),
        (1792, 257),
        (7168, 1025),
        (7, 1),
        (7175, 1024),
        (1799, 256),
        (455, 64),
        (119, 16),
    ),
    8: (
        (128, 17),
        (512, 65),
        (2048, 257),
        (8192, 1025),
        (8, 1),
        (1025, 128),
        (257, 32),
        (65, 8),
        (17, 2),
    ),
    9: (
        (144, 17),
        (576, 65),
        (2304, 257),
        (9216, 1025),
        (9, 1),
        (9225, 1024),
        (2313, 256),
        (585, 64),
        (153, 16),
    ),
}

EXPECTED_DEGREES_BY_Q = {
    7: (12, 16, 20),
    8: (16, 20, 24),
    9: (16, 20, 24),
}

EXPECTED_PRECISIONS = {
    (7, 12): ((192, 384), 768),
    (7, 16): ((192, 384, 768), 1536),
    (7, 20): ((384, 768), 1536),
    (8, 16): ((192, 384, 768), 1536),
    (8, 20): ((384, 768), 1536),
    (8, 24): ((384, 768), 1536),
    (9, 16): ((192, 384, 768), 1536),
    (9, 20): ((384, 768), 1536),
    (9, 24): ((384, 768), 1536),
}

EXPECTED_WITNESS_POLICY = {
    "eigenpair_count": "3",
    "eigensolver_role": "approximate-untrusted-hint-only",
    "normalization": "primitive-gcd-and-first-nonzero-positive",
    "parts": ["real", "imaginary"],
    "projections": ["even", "odd", "raw"],
    "rounding": "exact-rational-round-ties-to-even",
    "scale_bits": ["8", "12", "16", "24", "32", "48"],
}


def _raw_frozen_plan() -> dict[str, Any]:
    return json.loads(FROZEN_PLAN_PATH.read_text(encoding="utf-8"))


def _cutoff(cell: dict[str, Any]) -> Fraction:
    cutoff = cell["matrix_input"]["cutoff_c"]
    return Fraction(int(cutoff["numerator"]), int(cutoff["denominator"]))


def _cutoff_pair(cell: dict[str, Any]) -> tuple[int, int]:
    value = _cutoff(cell)
    return value.numerator, value.denominator


def _cells_for(
    plan: dict[str, Any],
    *,
    q: int,
    degree: int | None = None,
    position: str | None = None,
) -> list[dict[str, Any]]:
    cells = [cell for cell in plan["cells"] if int(cell["transition"]["q"]) == q]
    if degree is not None:
        cells = [cell for cell in cells if int(cell["matrix_input"]["degree"]) == degree]
    if position is not None:
        cells = [cell for cell in cells if cell["transition"]["position"] == position]
    return cells


def _tiny_plan() -> dict[str, Any]:
    policy = search.WeilSearchPolicy(
        attempt_bits=(96,),
        confirmation_bits=128,
        eigenpair_count=1,
        scale_bits=(4,),
    )
    cells = [
        transition.WeilTransitionCellContract(
            q=2,
            prime=2,
            exponent=1,
            position=position,
            stencil_exponent=None if position == "exact" else 4,
            cutoff_numerator=numerator,
            cutoff_denominator=denominator,
            degree=0,
            policy=policy,
        ).to_record()
        for position, numerator, denominator in (
            ("below", 32, 17),
            ("exact", 2, 1),
            ("above", 17, 8),
        )
    ]
    return transition.canonicalize_transition_plan(
        {
            "schema": transition.PLAN_SCHEMA,
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "cells": cells,
            "backend_contract": search._backend_record(),
        }
    )


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def test_frozen_transition_plan_has_the_exact_81_cell_contract() -> None:
    raw = _raw_frozen_plan()
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)

    assert plan == raw
    assert plan["schema"] == transition.PLAN_SCHEMA
    assert plan["classification"] == "EXPLORATORY"
    assert plan["hypothesis_status"] == "UNRESOLVED"
    assert plan["frozen_batch"] == transition.FROZEN_Q7_Q9_BATCH
    assert plan["plan_sha256"] == FROZEN_PLAN_SHA256
    assert canonical_json_file_sha256(FROZEN_PLAN_PATH) == FROZEN_PLAN_FILE_SHA256
    assert plan["plan_sha256"] == content_sha256(
        {key: value for key, value in plan.items() if key != "plan_sha256"}
    )
    assert len(plan["cells"]) == 81

    actual_order = [
        (
            _cutoff(cell),
            int(cell["matrix_input"]["degree"]),
            int(cell["matrix_input"]["cutoff_c"]["numerator"]),
            int(cell["matrix_input"]["cutoff_c"]["denominator"]),
        )
        for cell in plan["cells"]
    ]
    assert actual_order == sorted(actual_order)

    expected_global_cutoffs = sorted(
        Fraction(numerator, denominator)
        for cutoffs in EXPECTED_CUTOFFS_BY_Q.values()
        for numerator, denominator in cutoffs
    )
    actual_global_cutoffs = sorted({_cutoff(cell) for cell in plan["cells"]})
    assert actual_global_cutoffs == expected_global_cutoffs
    assert len(actual_global_cutoffs) == 27

    # The q=8 and q=9 neighborhoods overlap: global rational ordering, rather
    # than transition-group ordering, must put 144/17 before 17/2.
    assert actual_global_cutoffs.index(Fraction(144, 17)) < actual_global_cutoffs.index(
        Fraction(17, 2)
    )

    for q, expected_cutoffs in EXPECTED_CUTOFFS_BY_Q.items():
        q_cells = _cells_for(plan, q=q)
        assert len(q_cells) == 27
        assert sorted({_cutoff_pair(cell) for cell in q_cells}, key=lambda pair: Fraction(*pair)) == list(
            expected_cutoffs
        )
        assert sorted({int(cell["matrix_input"]["degree"]) for cell in q_cells}) == list(
            EXPECTED_DEGREES_BY_Q[q]
        )
        for degree in EXPECTED_DEGREES_BY_Q[q]:
            assert len(_cells_for(plan, q=q, degree=degree)) == 9


def test_frozen_contract_hashes_and_per_cell_policies_are_exact_and_unique() -> None:
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    contract_ids: set[str] = set()
    matrix_hashes: set[str] = set()
    contract_hashes: set[str] = set()

    for cell in plan["cells"]:
        q = int(cell["transition"]["q"])
        degree = int(cell["matrix_input"]["degree"])
        attempt_bits, confirmation_bits = EXPECTED_PRECISIONS[(q, degree)]
        assert cell["policy"] == {
            "attempt_bits": [str(value) for value in attempt_bits],
            "candidate_confirmation": (
                "same primitive integer witness is interval-negative at a "
                "strictly higher precision"
            ),
            "confirmation_bits": str(confirmation_bits),
            "positive_terminal_gate": "fixed-order-interval-ldlt",
            "witness_policy": EXPECTED_WITNESS_POLICY,
        }

        matrix_input = cell["matrix_input"]
        assert matrix_input["dimension"] == str(2 * degree + 1)
        expected_matrix_hash = content_sha256(
            {
                "cutoff_c": matrix_input["cutoff_c"],
                "degree": matrix_input["degree"],
                "dimension": matrix_input["dimension"],
            }
        )
        assert matrix_input["input_sha256"] == expected_matrix_hash
        assert cell["matrix_input_sha256"] == expected_matrix_hash
        assert cell["cell_contract_sha256"] == content_sha256(
            {
                key: value
                for key, value in cell.items()
                if key != "cell_contract_sha256"
            }
        )

        contract_ids.add(cell["contract_id"])
        matrix_hashes.add(cell["matrix_input_sha256"])
        contract_hashes.add(cell["cell_contract_sha256"])

    assert len(contract_ids) == 81
    assert len(matrix_hashes) == 81
    assert len(contract_hashes) == 81


def test_each_stencil_is_exact_reciprocal_and_changes_one_prime_power() -> None:
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    expected_transition = {7: (7, 1, 7), 8: (2, 3, 8), 9: (3, 2, 9)}

    for q, endpoint in expected_transition.items():
        exact_cutoffs = {_cutoff(cell) for cell in _cells_for(plan, q=q, position="exact")}
        assert exact_cutoffs == {Fraction(q, 1)}
        for j in (4, 6, 8, 10):
            below = {
                _cutoff(cell)
                for cell in _cells_for(plan, q=q, position="below")
                if int(cell["transition"]["stencil_exponent"]) == j
            }
            above = {
                _cutoff(cell)
                for cell in _cells_for(plan, q=q, position="above")
                if int(cell["transition"]["stencil_exponent"]) == j
            }
            assert len(below) == len(above) == 1
            c_minus = next(iter(below))
            c_plus = next(iter(above))
            assert c_minus * c_plus == q * q

            below_transcript = set(
                prime_powers_leq(c_minus.numerator, c_minus.denominator)
            )
            above_transcript = set(
                prime_powers_leq(c_plus.numerator, c_plus.denominator)
            )
            assert above_transcript - below_transcript == {endpoint}
            assert below_transcript - above_transcript == set()

    previous_precision = ctx.prec
    try:
        ctx.prec = 128
        for q, (prime, exponent, _) in expected_transition.items():
            length = arb(q).log()
            endpoint_only = prime_power_entry(
                2,
                2,
                length,
                [(prime, exponent, q)],
                q,
                1,
            )
            assert endpoint_only.is_zero()
    finally:
        ctx.prec = previous_precision


@pytest.mark.parametrize(
    ("label", "mutate"),
    [
        (
            "composite prime",
            lambda raw: raw["cells"][0]["transition"].__setitem__("prime", "4"),
        ),
        (
            "wrong exponent",
            lambda raw: raw["cells"][0]["transition"].__setitem__("exponent", "2"),
        ),
        (
            "wrong q",
            lambda raw: raw["cells"][0]["transition"].__setitem__("q", "8"),
        ),
        (
            "wrong side",
            lambda raw: raw["cells"][0]["transition"].__setitem__("position", "above"),
        ),
        (
            "wrong stencil exponent",
            lambda raw: raw["cells"][0]["transition"].__setitem__(
                "stencil_exponent", "6"
            ),
        ),
        (
            "wrong formula label",
            lambda raw: raw["cells"][0]["transition"].__setitem__(
                "formula", "not-the-stencil"
            ),
        ),
        (
            "wrong matrix hash",
            lambda raw: raw["cells"][0].__setitem__(
                "matrix_input_sha256", "0" * 64
            ),
        ),
        (
            "wrong contract hash",
            lambda raw: raw["cells"][0].__setitem__(
                "cell_contract_sha256", "0" * 64
            ),
        ),
        (
            "wrong cell schema",
            lambda raw: raw["cells"][0].__setitem__("schema", "v0"),
        ),
    ],
)
def test_frozen_plan_cell_mutations_fail_closed(label: str, mutate: Any) -> None:
    raw = _raw_frozen_plan()
    mutate(raw)
    with pytest.raises(transition.WeilTransitionPlanError):
        transition.canonicalize_transition_plan(raw)


def test_duplicate_inputs_wrong_plan_hash_and_intervening_transition_fail_closed() -> None:
    duplicate = _raw_frozen_plan()
    duplicate["cells"].append(copy.deepcopy(duplicate["cells"][0]))
    with pytest.raises(transition.WeilTransitionPlanError, match="duplicate"):
        transition.canonicalize_transition_plan(duplicate)

    wrong_hash = _raw_frozen_plan()
    wrong_hash["plan_sha256"] = "0" * 64
    with pytest.raises(transition.WeilTransitionPlanError, match="plan hash"):
        transition.canonicalize_transition_plan(wrong_hash)

    crossing = _raw_frozen_plan()
    cell = crossing["cells"][0]
    degree = int(cell["matrix_input"]["degree"])
    cell["transition"] = {
        "q": "2",
        "prime": "2",
        "exponent": "1",
        "position": "above",
        "stencil_exponent": "1",
        "formula": "q*(2^j+1)/2^j",
    }
    cell["matrix_input"] = search.WeilSearchCell(3, 1, degree).to_record()
    for key in (
        "contract_id",
        "matrix_input_sha256",
        "cell_contract_sha256",
        "numerical_engine_id",
        "v1_plan_binding",
    ):
        cell.pop(key, None)
    with pytest.raises(transition.WeilTransitionPlanError, match="reaches or crosses"):
        transition.canonicalize_transition_plan(crossing)


@pytest.mark.parametrize(
    "drifted_policy",
    [
        search.WeilSearchPolicy(
            attempt_bits=(192, 384),
            confirmation_bits=1024,
            eigenpair_count=3,
            scale_bits=(8, 12, 16, 24, 32, 48),
        ),
        search.WeilSearchPolicy(
            attempt_bits=(192, 384),
            confirmation_bits=768,
            eigenpair_count=4,
            scale_bits=(8, 12, 16, 24, 32, 48, 64),
        ),
    ],
    ids=("precision-schedule", "witness-policy"),
)
def test_frozen_batch_rejects_rehashed_internally_valid_policy_drift(
    drifted_policy: search.WeilSearchPolicy,
) -> None:
    raw = _raw_frozen_plan()
    original = raw["cells"][0]
    metadata = original["transition"]
    matrix_input = original["matrix_input"]
    cutoff = matrix_input["cutoff_c"]
    drifted_contract = transition.WeilTransitionCellContract(
        q=int(metadata["q"]),
        prime=int(metadata["prime"]),
        exponent=int(metadata["exponent"]),
        position=metadata["position"],
        stencil_exponent=int(metadata["stencil_exponent"]),
        cutoff_numerator=int(cutoff["numerator"]),
        cutoff_denominator=int(cutoff["denominator"]),
        degree=int(matrix_input["degree"]),
        policy=drifted_policy,
    ).to_record()
    raw["cells"][0] = drifted_contract

    assert drifted_contract["cell_contract_sha256"] == content_sha256(
        {
            key: value
            for key, value in drifted_contract.items()
            if key != "cell_contract_sha256"
        }
    )
    raw["plan_sha256"] = content_sha256(
        {key: value for key, value in raw.items() if key != "plan_sha256"}
    )
    assert raw["frozen_batch"] == transition.FROZEN_Q7_Q9_BATCH

    # The drift is a valid general v2 contract; only the frozen release policy
    # makes it invalid. This prevents a fully rehashed plan from weakening the
    # declared precision or candidate-generation schedule.
    unfrozen = copy.deepcopy(raw)
    unfrozen.pop("frozen_batch")
    unfrozen.pop("plan_sha256")
    canonical_unfrozen = transition.canonicalize_transition_plan(unfrozen)
    assert canonical_unfrozen["cells"][0]["policy"] == drifted_policy.to_record()

    with pytest.raises(
        transition.WeilTransitionPlanError,
        match="per-cell precision or witness policy changed",
    ):
        transition.canonicalize_transition_plan(raw)


def test_tiny_transition_run_checkpoints_resumes_replays_and_rejects_tampering(
    tmp_path: Path,
) -> None:
    plan = _tiny_plan()
    checkpoint = tmp_path / "tiny-transition"

    partial = transition.run_weil_transition_search(
        plan,
        checkpoint,
        max_cells=1,
    )
    assert partial["progress"]["planned_cells"] == "3"
    assert partial["progress"]["completed_cells"] == "1"
    assert len(partial["missing_contract_ids"]) == 2
    assert partial["conclusion"] == "INCOMPLETE_BOUNDED_TRANSITION_SEARCH"
    with pytest.raises(
        transition.WeilTransitionCheckpointError,
        match="pass resume=True",
    ):
        transition.run_weil_transition_search(plan, checkpoint)

    complete = transition.run_weil_transition_search(
        plan,
        checkpoint,
        resume=True,
    )
    assert complete["progress"] == {
        "planned_cells": "3",
        "completed_cells": "3",
        "status_counts": {
            "FINITE_POSITIVE_CERTIFIED": "3",
            "NEGATIVE_CANDIDATE_QUARANTINED": "0",
            "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE": "0",
            "INCONCLUSIVE_MAX_PRECISION": "0",
        },
    }
    assert complete["missing_contract_ids"] == []
    assert complete["classification"] == "EXPLORATORY"
    assert complete["hypothesis_status"] == "UNRESOLVED"

    reproduced = transition.verify_weil_transition_search(
        checkpoint / "index.json",
        checkpoint,
    )
    assert reproduced["classification"] == (
        "REPRODUCED_EXPLORATORY_TRANSITION_SEARCH"
    )
    assert reproduced["hypothesis_status"] == "UNRESOLVED"
    assert reproduced["verified_cells"] == "3"

    # A no-op resume must replay the existing checkpoint and reproduce the
    # same compact index without scheduling any additional cells.
    resumed_again = transition.run_weil_transition_search(
        plan,
        checkpoint,
        resume=True,
    )
    assert resumed_again == complete

    tampered_index = copy.deepcopy(complete)
    tampered_index["conclusion"] = "NO_NEGATIVE_WITNESS_IN_FINITE_TRANSITION_BATCH"
    tampered_index = _rehash(tampered_index)
    with pytest.raises(
        transition.WeilTransitionVerificationError,
        match="does not match referenced artifacts",
    ):
        transition.verify_weil_transition_search(tampered_index, checkpoint)

    stored_plan = json.loads((checkpoint / "plan.json").read_text(encoding="utf-8"))
    stored_plan["cells"][0]["transition"]["formula"] = "tampered"
    (checkpoint / "plan.json").write_text(
        json.dumps(stored_plan, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(
        transition.WeilTransitionCheckpointError,
        match="stored transition plan is invalid",
    ):
        transition.run_weil_transition_search(plan, checkpoint, resume=True)


def test_quarantined_candidate_does_not_stop_later_cells_and_requires_confirmation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = _tiny_plan()
    checkpoint = tmp_path / "candidate-transition"
    searched: list[str] = []
    confirmations: list[str] = []

    def fake_search_cell(
        cell: search.WeilSearchCell,
        policy: search.WeilSearchPolicy,
        **kwargs: Any,
    ) -> dict[str, Any]:
        del policy, kwargs
        searched.append(cell.cell_id)
        status = (
            "NEGATIVE_CANDIDATE_QUARANTINED"
            if len(searched) == 1
            else "FINITE_POSITIVE_CERTIFIED"
        )
        return _rehash(
            {
                "schema": search.CELL_SCHEMA,
                "classification": "EXPLORATORY",
                "hypothesis_status": "UNRESOLVED",
                "cell": cell.to_record(),
                "status": status,
                "attempts": [],
            }
        )

    def fake_dedicated_confirmation(
        contract: transition.WeilTransitionCellContract,
        artifact: dict[str, Any],
    ) -> dict[str, Any]:
        assert artifact["status"] == "NEGATIVE_CANDIDATE_QUARANTINED"
        confirmations.append(contract.contract_id)
        return _rehash(
            {
                "schema": "rh-lab/weil-transition-dedicated-confirmation/v2",
                "classification": "EXPLORATORY",
                "hypothesis_status": "UNRESOLVED",
                "contract_id": contract.contract_id,
                "decision": "V2_NEGATIVE_CANDIDATE_CONFIRMED",
            }
        )

    monkeypatch.setattr(search, "_search_weil_cell", fake_search_cell)
    monkeypatch.setattr(
        transition,
        "_dedicated_confirmation",
        fake_dedicated_confirmation,
    )
    index = transition.run_weil_transition_search(plan, checkpoint)

    assert len(searched) == 3
    assert len(confirmations) == 1
    assert confirmations == [index["cells"][0]["contract_id"]]
    assert index["progress"]["completed_cells"] == "3"
    assert index["progress"]["status_counts"] == {
        "FINITE_POSITIVE_CERTIFIED": "2",
        "NEGATIVE_CANDIDATE_QUARANTINED": "1",
        "NEGATIVE_OBSERVATION_QUARANTINED_INCONCLUSIVE": "0",
        "INCONCLUSIVE_MAX_PRECISION": "0",
    }
    assert index["conclusion"] == (
        "QUARANTINED_NEGATIVE_CANDIDATE_IN_FINITE_BATCH"
    )
    assert index["cells"][0]["dedicated_confirmation"] is not None

    with pytest.raises(
        transition.WeilTransitionCheckpointError,
        match="lacks dedicated v2 confirmation",
    ):
        transition._v2_status(
            {"status": "NEGATIVE_CANDIDATE_QUARANTINED"},
            None,
        )

    confirmation_path = (
        checkpoint
        / index["cells"][0]["dedicated_confirmation"]["path"]
    )
    confirmation_path.unlink()
    monkeypatch.setattr(
        transition,
        "_derive_v1_artifact",
        lambda contract, artifact, *, replay: artifact,
    )
    with pytest.raises(
        transition.WeilTransitionVerificationError,
        match="cannot read checkpoint",
    ):
        transition.verify_weil_transition_search(index, checkpoint)
