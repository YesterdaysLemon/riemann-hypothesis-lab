from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Callable

from flint import arb
import pytest

from riemann_lab.artifacts import canonical_json_file_sha256, content_sha256
import riemann_lab.nyman_search as search


FROZEN_PLAN_PATH = Path("plans/nyman-natural-v1.json")
FROZEN_PLAN_SHA256 = (
    "27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a"
)
FROZEN_PLAN_FILE_SHA256 = (
    "563b136f43bb6129ce182d0f29806dc5bbb70c30983cee8593806417aacfa9ba"
)


def _raw_plan() -> dict[str, Any]:
    return json.loads(FROZEN_PLAN_PATH.read_text(encoding="utf-8"))


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _install_fake_core(
    monkeypatch: pytest.MonkeyPatch,
    *,
    lower_inconclusive_n: set[int] | None = None,
    constant_energy: bool = False,
    replay_upper_failure: Callable[[], bool] | None = None,
) -> dict[str, list[Any]]:
    lower_inconclusive_n = lower_inconclusive_n or set()
    replay_upper_failure = replay_upper_failure or (lambda: False)
    calls: dict[str, list[Any]] = {
        "build": [],
        "upper": [],
        "lower": [],
        "round": [],
    }

    def fake_build(dilates: Any) -> Any:
        exact = tuple(dilates)
        calls["build"].append(exact)
        size = len(exact)
        gram = tuple(
            tuple(arb(1 if row == column else 0) for column in range(size))
            for row in range(size)
        )
        target = tuple(arb(0) for _ in range(size))
        return search.core.NaturalSystem(exact, gram, target, search.ctx.prec)

    def fake_prefix(system: Any, n: int) -> Any:
        return search.core.NaturalSystem(
            tuple(system.dilates[:n]),
            tuple(tuple(row[:n]) for row in system.gram[:n]),
            tuple(system.target[:n]),
            system.precision_bits,
        )

    def fake_system_hash(system: Any) -> str:
        return content_sha256(
            {
                "dilates": [str(value) for value in system.dilates],
                "precision_bits": str(system.precision_bits),
            }
        )

    def fake_round(values: Any, exponent: int) -> tuple[int, ...]:
        exact = tuple(values)
        calls["round"].append((len(exact), exponent))
        return tuple(0 for _ in exact)

    def fake_energy(system: Any, coefficients: Any) -> Any:
        assert len(tuple(coefficients)) == len(system.dilates)
        denominator = 64 if constant_energy else len(system.dilates) + 32
        return arb(1) / denominator

    def fake_upper(
        dilates: Any,
        numerators: Any,
        denominator_exponent: int,
        claimed_upper_bound: Fraction,
        *,
        precision_bits: int,
        system: Any,
    ) -> dict[str, Any]:
        exact_dilates = tuple(dilates)
        exact_numerators = tuple(numerators)
        assert tuple(system.dilates) == exact_dilates
        assert system.precision_bits == precision_bits
        calls["upper"].append((len(exact_dilates), precision_bits))
        certified = not (
            precision_bits == search.FROZEN_REPLAY_BITS
            and replay_upper_failure()
        )
        return search._with_payload_hash(
            {
                "schema": "fake-upper/v1",
                "classification": "EXPLORATORY",
                "hypothesis_status": "UNRESOLVED",
                "precision_bits": str(precision_bits),
                "decision": (
                    "UPPER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
                ),
                "scope": [str(value) for value in exact_dilates],
                "coefficients": {
                    "numerators": [str(value) for value in exact_numerators],
                    "denominator_exponent": str(denominator_exponent),
                },
                "claimed_upper_bound": search._fraction_record(
                    claimed_upper_bound
                ),
                "evidence": {"direct_energy_marker": "unchanged"},
                "checks": {
                    "energy_at_most_claimed_upper_bound": certified,
                    "energy_strictly_below_claimed_upper_bound": certified,
                },
            }
        )

    def fake_lower(
        dilates: Any,
        lower_bound: Fraction,
        *,
        precision_bits: int,
        system: Any,
    ) -> dict[str, Any]:
        exact_dilates = tuple(dilates)
        size = len(exact_dilates)
        assert tuple(system.dilates) == exact_dilates
        assert system.precision_bits == precision_bits
        calls["lower"].append((size, precision_bits))
        certified = size not in lower_inconclusive_n
        return search._with_payload_hash(
            {
                "schema": "fake-lower/v1",
                "classification": "EXPLORATORY",
                "hypothesis_status": "UNRESOLVED",
                "precision_bits": str(precision_bits),
                "decision": (
                    "LOWER_BOUND_CERTIFIED" if certified else "INCONCLUSIVE"
                ),
                "scope": [str(value) for value in exact_dilates],
                "claimed_lower_bound": search._fraction_record(lower_bound),
                "augmented_ldlt": {
                    "classification": (
                        "POSITIVE_DEFINITE" if certified else "INCONCLUSIVE"
                    ),
                    "pivots": [{"index": "0", "marker": "unchanged"}],
                },
                "checks": {
                    "fixed_order_interval_ldlt_positive": certified,
                },
            }
        )

    monkeypatch.setattr(search.core, "build_natural_system", fake_build)
    monkeypatch.setattr(search.core, "prefix_natural_system", fake_prefix)
    monkeypatch.setattr(
        search.core, "natural_system_content_sha256", fake_system_hash
    )
    monkeypatch.setattr(search.core, "round_arb_vector_to_dyadic", fake_round)
    monkeypatch.setattr(search.core, "evaluate_natural_distance", fake_energy)
    monkeypatch.setattr(search.core, "certify_dyadic_upper_bound", fake_upper)
    monkeypatch.setattr(search.core, "certify_augmented_lower_bound", fake_lower)
    return calls


def test_frozen_plan_is_exact_canonical_and_hash_bound() -> None:
    raw = _raw_plan()
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)

    assert plan == raw
    assert plan["plan_sha256"] == FROZEN_PLAN_SHA256
    assert canonical_json_file_sha256(FROZEN_PLAN_PATH) == FROZEN_PLAN_FILE_SHA256
    assert plan["plan_sha256"] == content_sha256(
        {key: value for key, value in plan.items() if key != "plan_sha256"}
    )
    assert plan["classification"] == "EXPLORATORY"
    assert plan["hypothesis_status"] == "UNRESOLVED"
    assert plan["frozen_experiment"] == search.FROZEN_EXPERIMENT
    assert [int(cell["n"]) for cell in plan["cells"]] == list(
        search.FROZEN_N_VALUES
    )
    assert plan["policy"] == search._policy_record(
        768, 1536, 256, 128, 120
    )


def test_all_cells_bind_exact_prefixes_of_one_shared_n256_kernel() -> None:
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    kernel = plan["engine"]["shared_kernel"]
    assert kernel == search._shared_kernel_contract(256)

    contract_hashes: set[str] = set()
    for cell, n in zip(plan["cells"], search.FROZEN_N_VALUES, strict=True):
        assert cell == search.NymanSearchCellContract(n, 256).to_record()
        assert cell["dimension"] == str(n)
        assert cell["dilates"] == [str(value) for value in range(1, n + 1)]
        assert cell["shared_kernel_binding"] == {
            "kernel_contract_sha256": kernel["kernel_contract_sha256"],
            "shared_kernel_n": "256",
            "prefix_length": str(n),
        }
        contract_hashes.add(cell["cell_contract_sha256"])
    assert len(contract_hashes) == 6


@pytest.mark.parametrize(
    ("label", "mutate", "message"),
    [
        (
            "classification",
            lambda raw: raw.__setitem__("classification", "PROVED"),
            "EXPLORATORY",
        ),
        (
            "hypothesis status",
            lambda raw: raw.__setitem__("hypothesis_status", "PROVED"),
            "resolve RH",
        ),
        (
            "N sequence",
            lambda raw: raw["cells"].pop(),
            "N sequence",
        ),
        (
            "generation precision",
            lambda raw: raw["policy"].__setitem__(
                "generation_precision_bits", "767"
            ),
            "arithmetic policy",
        ),
        (
            "lower slack",
            lambda raw: raw["policy"]["target_lower_slack"].__setitem__(
                "bits", "119"
            ),
            "policy is not canonical",
        ),
        (
            "prefix hash",
            lambda raw: raw["cells"][0]["shared_kernel_binding"].__setitem__(
                "kernel_contract_sha256", "0" * 64
            ),
            "cell contract is not canonical",
        ),
        (
            "contract hash",
            lambda raw: raw["cells"][0].__setitem__(
                "cell_contract_sha256", "0" * 64
            ),
            "cell contract is not canonical",
        ),
        (
            "engine",
            lambda raw: raw["engine"].__setitem__("objective", "mutated"),
            "engine contract changed",
        ),
        (
            "plan hash",
            lambda raw: raw.__setitem__("plan_sha256", "0" * 64),
            "plan hash mismatch",
        ),
    ],
)
def test_frozen_plan_mutations_fail_closed(
    label: str,
    mutate: Callable[[dict[str, Any]], None],
    message: str,
) -> None:
    raw = _raw_plan()
    mutate(raw)
    with pytest.raises(search.NymanSearchPlanError, match=message):
        search.canonicalize_nyman_plan(raw)


def test_search_runs_all_six_checkpoints_resumes_and_replays_stored_candidates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = _install_fake_core(monkeypatch)
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    checkpoint = tmp_path / "nyman"

    partial = search.run_nyman_search(plan, checkpoint, max_cells=2)
    assert partial["conclusion"] == "INCOMPLETE_BOUNDED_NYMAN_SEARCH"
    assert partial["progress"] == {
        "planned_cells": "6",
        "completed_cells": "2",
        "status_counts": {
            "FINITE_DISTANCE_BRACKET_CERTIFIED": "2",
            "INCONCLUSIVE": "0",
        },
    }
    assert partial["missing_cell_ids"] == [
        "n-0032",
        "n-0064",
        "n-0128",
        "n-0256",
    ]
    with pytest.raises(search.NymanSearchCheckpointError, match="pass resume=True"):
        search.run_nyman_search(plan, checkpoint)

    complete = search.run_nyman_search(plan, checkpoint, resume=True)
    assert complete["conclusion"] == "ALL_FINITE_DISTANCE_BRACKETS_CERTIFIED"
    assert complete["progress"]["completed_cells"] == "6"
    assert complete["progress"]["status_counts"] == {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": "6",
        "INCONCLUSIVE": "0",
    }
    assert complete["missing_cell_ids"] == []
    assert calls["round"] == [(8, 256), (16, 256), (32, 256), (64, 256), (128, 256), (256, 256)]
    assert calls["build"] == [tuple(range(1, 257)), tuple(range(1, 257))]

    for n in search.FROZEN_N_VALUES:
        artifact = json.loads(
            (checkpoint / "cells" / f"n-{n:04d}.json").read_text(
                encoding="utf-8"
            )
        )
        assert artifact["status"] == "FINITE_DISTANCE_BRACKET_CERTIFIED"
        assert artifact["terminal_statement"] == "L_N < d_N^2 <= U_N"
        candidate = artifact["candidate"]
        assert len(candidate["coefficients"]["numerators"]) == n
        assert candidate["coefficients"]["denominator_exponent"] == "256"
        _, upper = search._parse_dyadic_record(
            candidate["upper_bound"],
            "U",
            expected_exponent=128,
            error_type=search.NymanSearchVerificationError,
        )
        _, lower = search._parse_dyadic_record(
            candidate["lower_bound"],
            "L",
            expected_exponent=128,
            error_type=search.NymanSearchVerificationError,
        )
        assert upper - lower == Fraction(1, 1 << 120)

    no_op = search.run_nyman_search(plan, checkpoint, resume=True)
    assert no_op == complete
    assert len(calls["build"]) == 3

    # Candidate generation is not part of replay.  The verifier does rebuild
    # one canonical max-N kernel at each precision, but it never solves again
    # and never demands that a regenerated coefficient vector equal stored c.
    monkeypatch.setattr(
        search.core,
        "round_arb_vector_to_dyadic",
        lambda values, exponent: (_ for _ in ()).throw(
            AssertionError("regenerated")
        ),
    )
    reproduced = search.verify_nyman_search(
        checkpoint / "index.json", checkpoint
    )
    assert reproduced["classification"] == "REPRODUCED_EXPLORATORY_NYMAN_SEARCH"
    assert reproduced["hypothesis_status"] == "UNRESOLVED"
    assert reproduced["conclusion"] == "ALL_FINITE_DISTANCE_BRACKETS_CERTIFIED"
    assert reproduced["verified_cells"] == "6"
    assert reproduced["replay_precision_bits"] == "1536"
    assert reproduced["replay_status_counts"] == {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": "6",
        "INCONCLUSIVE": "0",
    }
    assert reproduced["generation_kernel_payload_sha256"] == complete[
        "generation_kernel"
    ]["payload_sha256"]
    assert isinstance(reproduced["replay_kernel_payload_sha256"], str)
    assert reproduced["stored_candidates_replayed_without_regeneration"] is True
    assert reproduced["same_backend_replay_only"] is True
    assert {bits for _, bits in calls["upper"]} == {768, 1536}
    assert {bits for _, bits in calls["lower"]} == {768, 1536}
    assert len(calls["build"]) == 5


def test_zero_cell_dry_run_is_checkpoint_only_and_kernel_lazy(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = _install_fake_core(monkeypatch)
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    checkpoint = tmp_path / "dry-run"

    dry = search.run_nyman_search(plan, checkpoint, max_cells=0)
    assert dry["progress"]["completed_cells"] == "0"
    assert dry["generation_kernel"] is None
    assert calls["build"] == []
    assert not (checkpoint / "kernels" / "generation.json").exists()

    reproduced = search.verify_nyman_search(
        checkpoint / "index.json", checkpoint
    )
    assert reproduced["verified_cells"] == "0"
    assert reproduced["generation_kernel_payload_sha256"] is None
    assert reproduced["replay_kernel_payload_sha256"] is None
    assert calls["build"] == []

    one = search.run_nyman_search(
        plan, checkpoint, resume=True, max_cells=1
    )
    assert one["progress"]["completed_cells"] == "1"
    assert one["generation_kernel"] is not None
    assert len(calls["build"]) == 1


def test_failed_frozen_lower_target_is_inconclusive_without_adaptation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_core(monkeypatch, lower_inconclusive_n={32})
    index = search.run_nyman_search(
        search.load_nyman_plan(FROZEN_PLAN_PATH), tmp_path / "inconclusive"
    )

    assert index["progress"]["status_counts"] == {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": "5",
        "INCONCLUSIVE": "1",
    }
    assert index["conclusion"] == (
        "FINITE_DISTANCE_BATCH_CONTAINS_INCONCLUSIVE_CELLS"
    )
    failed = next(cell for cell in index["cells"] if cell["cell_id"] == "n-0032")
    assert failed["status"] == "INCONCLUSIVE"
    artifact = json.loads(
        (tmp_path / "inconclusive" / failed["path"]).read_text(encoding="utf-8")
    )
    assert artifact["terminal_statement"] is None
    assert artifact["generation"]["lower_certificate"]["decision"] == "INCONCLUSIVE"
    assert artifact["candidate"]["upper_bound"]["denominator_exponent"] == "128"
    assert artifact["candidate"]["lower_bound"]["denominator_exponent"] == "128"
    by_pair = {
        (item["smaller_n"], item["larger_n"]): item
        for item in index["strict_improvement_diagnostics"]
    }
    assert by_pair[("16", "32")]["decision"] == (
        "UNAVAILABLE_UNCERTIFIED_BRACKET"
    )
    assert by_pair[("32", "64")]["decision"] == (
        "UNAVAILABLE_UNCERTIFIED_BRACKET"
    )


def test_strict_improvement_is_a_separate_non_gate_diagnostic(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_core(monkeypatch, constant_energy=True)
    index = search.run_nyman_search(
        search.load_nyman_plan(FROZEN_PLAN_PATH), tmp_path / "flat"
    )

    assert index["conclusion"] == "ALL_FINITE_DISTANCE_BRACKETS_CERTIFIED"
    assert len(index["strict_improvement_diagnostics"]) == 5
    assert all(
        item["decision"] == "STRICT_IMPROVEMENT_NOT_CERTIFIED"
        and item["cell_validity_gate"] is False
        for item in index["strict_improvement_diagnostics"]
    )


@pytest.mark.parametrize("mutation", ["coefficient", "pivot", "classification"])
def test_rehashed_cell_evidence_mutations_fail_reproduction(
    mutation: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_core(monkeypatch)
    checkpoint = tmp_path / mutation
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    search.run_nyman_search(plan, checkpoint, max_cells=1)
    path = checkpoint / "cells" / "n-0008.json"
    artifact = json.loads(path.read_text(encoding="utf-8"))

    if mutation == "coefficient":
        artifact["candidate"]["coefficients"]["numerators"][0] = "1"
        candidate_body = {
            key: value
            for key, value in artifact["candidate"].items()
            if key != "candidate_sha256"
        }
        artifact["candidate"]["candidate_sha256"] = content_sha256(candidate_body)
    elif mutation == "pivot":
        artifact["generation"]["lower_certificate"]["augmented_ldlt"]["pivots"][0][
            "marker"
        ] = "mutated"
        artifact["generation"]["lower_certificate"] = _rehash(
            artifact["generation"]["lower_certificate"]
        )
    else:
        artifact["classification"] = "PROVED"
    artifact = _rehash(artifact)
    _write_json(path, artifact)

    expected = {
        "coefficient": "upper certificate does not reproduce",
        "pivot": "lower certificate does not reproduce",
        "classification": "improperly promoted",
    }[mutation]
    with pytest.raises(search.NymanSearchVerificationError, match=expected):
        search.verify_nyman_search(checkpoint / "index.json", checkpoint)


def test_terminal_cell_must_remain_certified_at_frozen_replay_precision(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fail_replay = {"value": False}
    _install_fake_core(
        monkeypatch,
        replay_upper_failure=lambda: fail_replay["value"],
    )
    checkpoint = tmp_path / "replay-failure"
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    search.run_nyman_search(plan, checkpoint, max_cells=1)
    fail_replay["value"] = True

    with pytest.raises(
        search.NymanSearchVerificationError,
        match="failed requested-precision replay",
    ):
        search.verify_nyman_search(checkpoint / "index.json", checkpoint)

    with pytest.raises(search.NymanSearchVerificationError, match="frozen requested"):
        search.verify_nyman_search(
            checkpoint / "index.json",
            checkpoint,
            replay_precision_bits=1024,
        )


def test_index_and_checkpoint_mutations_fail_closed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _install_fake_core(monkeypatch)
    checkpoint = tmp_path / "mutations"
    plan = search.load_nyman_plan(FROZEN_PLAN_PATH)
    index = search.run_nyman_search(plan, checkpoint, max_cells=1)

    tampered = copy.deepcopy(index)
    tampered["cells"][0]["path"] = "../outside.json"
    tampered = _rehash(tampered)
    with pytest.raises(search.NymanSearchVerificationError, match="path changed"):
        search.verify_nyman_search(tampered, checkpoint)

    tampered = copy.deepcopy(index)
    tampered["conclusion"] = "PROVED_RH"
    tampered = _rehash(tampered)
    with pytest.raises(
        search.NymanSearchVerificationError,
        match="does not match referenced artifacts",
    ):
        search.verify_nyman_search(tampered, checkpoint)

    stored_plan = json.loads((checkpoint / "plan.json").read_text(encoding="utf-8"))
    stored_plan["classification"] = "PROVED"
    _write_json(checkpoint / "plan.json", stored_plan)
    with pytest.raises(
        search.NymanSearchCheckpointError,
        match="stored Nyman plan is invalid",
    ):
        search.run_nyman_search(plan, checkpoint, resume=True)
