from __future__ import annotations

import copy
from fractions import Fraction
from functools import reduce
import json
from math import gcd
from pathlib import Path
from typing import Any

import flint
from flint import arb, ctx
import pytest

from riemann_lab.artifacts import canonical_json_file_sha256, content_sha256
import riemann_lab.weil_search as search


def _raw_plan(
    *,
    cutoffs: tuple[tuple[int, int], ...] = ((3, 2),),
    degrees: tuple[int, ...] = (0,),
    attempt_bits: tuple[int, ...] = (96,),
    confirmation_bits: int = 128,
    eigenpair_count: int = 1,
    scale_bits: tuple[int, ...] = (4,),
) -> dict[str, Any]:
    return {
        "schema": search.PLAN_SCHEMA,
        "classification": "EXPLORATORY",
        "cutoffs": [
            {"numerator": str(numerator), "denominator": str(denominator)}
            for numerator, denominator in cutoffs
        ],
        "degrees": [str(degree) for degree in degrees],
        "precision_policy": {
            "attempt_bits": [str(value) for value in attempt_bits],
            "confirmation_bits": str(confirmation_bits),
        },
        "witness_policy": {
            "eigenpair_count": str(eigenpair_count),
            "scale_bits": [str(value) for value in scale_bits],
        },
        "backend_contract": {
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
    }


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _omitted_p2_components(cell: search.WeilSearchCell) -> dict[str, Any]:
    """Return the known corrupted P-R control while preserving ball rigor."""

    components = search._symmetric_components(
        cell.cutoff_numerator,
        cell.cutoff_denominator,
        cell.degree,
    )
    size = len(components["indices"])
    prime_power = [[arb(0) for _ in range(size)] for _ in range(size)]
    mutated = [[arb(0) for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for column in range(row, size):
            value = (
                components["pole"][row][column]
                - components["archimedean"][row][column]
            )
            mutated[row][column] = value
            mutated[column][row] = value
    components["prime_powers"] = []
    components["prime_power"] = prime_power
    components["q"] = mutated
    return components


def test_plan_canonicalization_reduces_sorts_hashes_and_loads(
    tmp_path: Path,
) -> None:
    raw = _raw_plan(
        cutoffs=((5, 2), (6, 4)),
        degrees=(2, 0),
        attempt_bits=(96, 192),
        confirmation_bits=256,
        scale_bits=(4, 8),
    )
    canonical = search.canonicalize_plan(raw)

    assert canonical["cutoffs"] == [
        {"numerator": "3", "denominator": "2"},
        {"numerator": "5", "denominator": "2"},
    ]
    assert canonical["degrees"] == ["0", "2"]
    assert canonical["classification"] == "EXPLORATORY"
    assert canonical["hypothesis_status"] == "UNRESOLVED"
    body = {
        key: value for key, value in canonical.items() if key != "plan_sha256"
    }
    assert canonical["plan_sha256"] == content_sha256(body)
    assert search.canonicalize_plan(canonical) == canonical

    equivalent = _raw_plan(
        cutoffs=((3, 2), (5, 2)),
        degrees=(0, 2),
        attempt_bits=(96, 192),
        confirmation_bits=256,
        scale_bits=(4, 8),
    )
    assert search.canonicalize_plan(equivalent) == canonical

    plan_path = tmp_path / "plan.json"
    _write_json(plan_path, canonical)
    assert search.load_search_plan(plan_path) == canonical


def test_frozen_grid_release_is_hash_bound_and_exploratory() -> None:
    plan_path = Path("plans/weil-grid-v1.json")
    index_path = Path("results/weil-search-grid-v1/index.json")
    plan = search.load_search_plan(plan_path)
    index = json.loads(index_path.read_text(encoding="utf-8"))

    assert plan["plan_sha256"] == (
        "68afe80851df02f1db4eb0b8119818bf45e9afef2d1ade7de4f73200f32827fc"
    )
    assert canonical_json_file_sha256(plan_path) == (
        "83d39d8dadae8fc2375589155863bca61b38619736c3c20b27a28db26b0eab42"
    )
    assert canonical_json_file_sha256(index_path) == (
        "e7e9e2cf1ae05b228d8e91380eeaf5e22ce525f9e849398b9b95a32ec2e3f762"
    )
    assert index["payload_sha256"] == (
        "f0595f9630b8738d518dfabfeb2d4d9921f4e0fcb0d4d09a9c201fe6f47e7d87"
    )
    assert index["classification"] == "EXPLORATORY"
    assert index["hypothesis_status"] == "UNRESOLVED"
    assert index["conclusion"] == "NO_NEGATIVE_WITNESS_IN_FINITE_GRID"
    assert index["missing_cells"] == []
    assert index["progress"] == {
        "planned_cells": "32",
        "completed_cells": "32",
        "status_counts": {
            "FINITE_POSITIVE_CERTIFIED": "32",
            "NEGATIVE_CANDIDATE_QUARANTINED": "0",
            "INCONCLUSIVE_MAX_PRECISION": "0",
        },
    }
    assert len(index["cells"]) == 32
    for reference in index["cells"]:
        cell_path = index_path.parent / reference["path"]
        cell = json.loads(cell_path.read_text(encoding="utf-8"))
        assert cell["payload_sha256"] == reference["cell_payload_sha256"]
        assert cell["payload_sha256"] == content_sha256(
            {
                key: value
                for key, value in cell.items()
                if key != "payload_sha256"
            }
        )
        assert cell["classification"] == "EXPLORATORY"
        assert cell["hypothesis_status"] == "UNRESOLVED"
        assert cell["status"] == "FINITE_POSITIVE_CERTIFIED"


def test_malformed_or_hash_inconsistent_plans_are_rejected(tmp_path: Path) -> None:
    cases: list[dict[str, Any]] = []

    bad_schema = _raw_plan()
    bad_schema["schema"] = "rh-lab/weil-search-plan/v0"
    cases.append(bad_schema)

    promoted = _raw_plan()
    promoted["classification"] = "PROVED"
    cases.append(promoted)

    empty_grid = _raw_plan()
    empty_grid["cutoffs"] = []
    cases.append(empty_grid)

    duplicate_cutoff = _raw_plan()
    duplicate_cutoff["cutoffs"] = [
        {"numerator": "3", "denominator": "2"},
        {"numerator": "6", "denominator": "4"},
    ]
    cases.append(duplicate_cutoff)

    noncanonical_integer = _raw_plan()
    noncanonical_integer["degrees"] = ["00"]
    cases.append(noncanonical_integer)

    boolean_integer = _raw_plan()
    boolean_integer["cutoffs"][0]["numerator"] = True
    cases.append(boolean_integer)

    repeated_precision = _raw_plan(attempt_bits=(96, 192))
    repeated_precision["precision_policy"]["attempt_bits"] = ["96", "96"]
    cases.append(repeated_precision)

    weak_confirmation = _raw_plan(attempt_bits=(96, 192))
    weak_confirmation["precision_policy"]["confirmation_bits"] = "192"
    cases.append(weak_confirmation)

    repeated_scale = _raw_plan(scale_bits=(4, 8))
    repeated_scale["witness_policy"]["scale_bits"] = ["4", "4"]
    cases.append(repeated_scale)

    missing_backend = _raw_plan()
    missing_backend["backend_contract"]["flint"] = ""
    cases.append(missing_backend)

    wrong_hash = _raw_plan()
    wrong_hash["plan_sha256"] = "0" * 64
    cases.append(wrong_hash)

    for malformed in cases:
        with pytest.raises(search.WeilSearchPlanError):
            search.canonicalize_plan(malformed)

    malformed_path = tmp_path / "malformed.json"
    malformed_path.write_text("{not-json", encoding="utf-8")
    with pytest.raises(search.WeilSearchPlanError, match="cannot read search plan"):
        search.load_search_plan(malformed_path)

    wrong_backend = _raw_plan()
    wrong_backend["backend_contract"]["flint"] = "0.0.0"
    with pytest.raises(search.WeilSearchPlanError, match="do not match"):
        search.run_weil_search(wrong_backend, tmp_path / "wrong-backend")


def test_exact_ties_to_even_and_primitive_witness_normalization() -> None:
    expected = {
        Fraction(1, 2): 0,
        Fraction(3, 2): 2,
        Fraction(5, 2): 2,
        Fraction(7, 2): 4,
        Fraction(-1, 2): 0,
        Fraction(-3, 2): -2,
        Fraction(-5, 2): -2,
        Fraction(-7, 2): -4,
        Fraction(8, 3): 3,
    }
    assert {
        value: search._round_ties_even(value) for value in expected
    } == expected

    assert search._primitive_witness((0, -6, 9)) == (0, 2, -3)
    assert search._primitive_witness((4, 0, -8)) == (1, 0, -2)
    assert search._primitive_witness((-10, 20)) == (1, -2)
    with pytest.raises(ValueError, match="nonzero"):
        search._primitive_witness((0, 0, 0))


def test_candidate_extraction_is_deterministic_and_primitive() -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 96
        zero = arb(0)
        matrix = [
            [arb(-2), zero, zero],
            [zero, arb(3), zero],
            [zero, zero, arb(5)],
        ]
        policy = search.WeilSearchPolicy(
            attempt_bits=(96,),
            confirmation_bits=128,
            eigenpair_count=1,
            scale_bits=(2, 4),
        )
        first = search.extract_integer_witnesses(matrix, policy)
        second = search.extract_integer_witnesses(matrix, policy)
    finally:
        ctx.prec = previous_precision

    assert first == second
    assert first
    witnesses = [tuple(candidate["witness"]) for candidate in first]
    assert len(witnesses) == len(set(witnesses))
    assert (1, 0, 0) in witnesses
    for witness in witnesses:
        common = reduce(gcd, (abs(value) for value in witness if value), 0)
        assert common == 1
        assert next(value for value in witness if value) > 0


def test_integer_rayleigh_signs_and_exact_symmetry_guard() -> None:
    zero = arb(0)
    positive = [[arb(2), zero], [zero, arb(3)]]
    positive_result = search.evaluate_integer_rayleigh(positive, [-2, 0])
    assert positive_result["classification"] == "POSITIVE_FOR_WITNESS"
    assert positive_result["sign"] == "POSITIVE"
    assert positive_result["witness"] == ["1", "0"]
    assert positive_result["norm_squared"] == "1"
    assert positive_result["lower_bound_is_positive"] is True

    negative = [[arb(-2), zero], [zero, arb(3)]]
    negative_result = search.evaluate_integer_rayleigh(negative, [6, 0])
    assert negative_result["classification"] == "NEGATIVE"
    assert negative_result["sign"] == "NEGATIVE"
    assert negative_result["witness"] == ["1", "0"]
    assert negative_result["upper_bound_is_negative"] is True

    inconclusive = search.evaluate_integer_rayleigh([[arb("0 +/- 1")]], [1])
    assert inconclusive["classification"] == "INCONCLUSIVE"
    assert inconclusive["sign"] == "INCONCLUSIVE"
    assert inconclusive["upper_bound_is_negative"] is False
    assert inconclusive["lower_bound_is_positive"] is False

    asymmetric = [
        [arb(1), arb("0 +/- 0.1")],
        [arb("0 +/- 0.2"), arb(1)],
    ]
    with pytest.raises(ValueError, match="identical symmetric enclosures"):
        search.evaluate_integer_rayleigh(asymmetric, [1, 1])
    with pytest.raises(TypeError, match="integers"):
        search.evaluate_integer_rayleigh([[arb(1)]], [True])
    with pytest.raises(ValueError, match="nonzero"):
        search.evaluate_integer_rayleigh([[arb(1)]], [0])


def test_genuine_one_cell_search_is_finite_positive_and_exploratory() -> None:
    policy = search.WeilSearchPolicy(
        attempt_bits=(96,),
        confirmation_bits=128,
        eigenpair_count=1,
        scale_bits=(4,),
    )
    artifact = search.search_weil_cell(search.WeilSearchCell(3, 2, 0), policy)

    assert artifact["schema"] == search.CELL_SCHEMA
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["status"] == "FINITE_POSITIVE_CERTIFIED"
    assert artifact["quarantined_candidate"] is None
    assert len(artifact["attempts"]) == 1
    attempt = artifact["attempts"][0]
    assert attempt["decision"] == "FINITE_POSITIVE_CERTIFIED"
    assert attempt["positive_definiteness"]["classification"] == (
        "POSITIVE_DEFINITE"
    )
    assert attempt["candidate_generation"] == {
        "classification": "SKIPPED_LDLT_POSITIVE",
        "role": "untrusted-search-hint-only",
        "tested": [],
    }
    checks = attempt["matrix_evidence"]["checks"]
    assert checks["all_pole_oracle_pairs_checked"] is True
    assert checks["all_archimedean_oracle_pairs_checked"] is True
    assert artifact["payload_sha256"] == content_sha256(
        {
            key: value
            for key, value in artifact.items()
            if key != "payload_sha256"
        }
    )


def test_corrupted_omit_p2_path_extracts_and_quarantines_candidate(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = _raw_plan(
        cutoffs=((5, 2),),
        degrees=(4,),
        attempt_bits=(96,),
        confirmation_bits=192,
        eigenpair_count=3,
        scale_bits=(8, 12, 16, 24),
    )
    checkpoint = tmp_path / "omitted-p2"

    with monkeypatch.context() as context:
        context.setattr(search, "_build_components", _omitted_p2_components)
        index = search.run_weil_search(plan, checkpoint)
        assert index["classification"] == "EXPLORATORY"
        assert index["hypothesis_status"] == "UNRESOLVED"
        assert index["conclusion"] == "NEGATIVE_CANDIDATE_QUARANTINED"

        cell_path = checkpoint / index["cells"][0]["path"]
        cell = json.loads(cell_path.read_text(encoding="utf-8"))
        assert cell["status"] == "NEGATIVE_CANDIDATE_QUARANTINED"
        assert [attempt["decision"] for attempt in cell["attempts"]] == [
            "NEGATIVE_WITNESS_DISCOVERED",
            "NEGATIVE_WITNESS_CONFIRMED",
        ]
        candidate = cell["quarantined_candidate"]
        assert candidate["promotion_status"].startswith("pending-independent")
        assert candidate["same_backend_replay_only"] is True
        assert candidate["discovery_evaluation"]["classification"] == "NEGATIVE"
        assert candidate["confirmation_evaluation"]["classification"] == "NEGATIVE"
        assert cell["attempts"][-1]["precision_containment"][
            "all_higher_precision_components_contained"
        ] is True

        witness = tuple(int(value) for value in candidate["witness"])
        assert reduce(gcd, (abs(value) for value in witness if value), 0) == 1
        assert next(value for value in witness if value) > 0

        replay = search.verify_weil_search(index, checkpoint)
        assert replay["classification"] == "REPRODUCED_EXPLORATORY_SEARCH"
        assert replay["verified_negative_witness_evaluations"] == "2"
        assert replay["hypothesis_status"] == "UNRESOLVED"

    with pytest.raises(
        search.WeilSearchVerificationError,
        match="matrix evidence does not match canonical regeneration",
    ):
        search.verify_weil_search(index, checkpoint)


def test_mid_cell_checkpoint_resume_is_byte_identical_to_full_run(
    tmp_path: Path,
) -> None:
    plan = _raw_plan(
        cutoffs=((5, 1),),
        degrees=(8,),
        attempt_bits=(96, 192),
        confirmation_bits=384,
        eigenpair_count=1,
        scale_bits=(4,),
    )
    full_dir = tmp_path / "full"
    resumed_dir = tmp_path / "resumed"
    full_index = search.run_weil_search(plan, full_dir)

    partial = search.run_weil_search(plan, resumed_dir, max_cells=0)
    assert partial["conclusion"] == "INCOMPLETE"
    canonical = search.canonicalize_plan(plan)
    cell = search._cells_from_plan(canonical)[0]
    policy = search._policy_from_canonical_plan(canonical)
    first_attempt = search._search_attempt(
        cell,
        policy,
        precision_bits=96,
        has_next_precision=True,
    )
    assert first_attempt["decision"] == "ESCALATE_PRECISION"
    first_attempt = search._append_sequence(first_attempt, 0)
    search._atomic_write_json(
        resumed_dir / "attempts" / cell.cell_id / "0000.json",
        first_attempt,
    )

    resumed_index = search.run_weil_search(plan, resumed_dir, resume=True)
    assert resumed_index == full_index
    assert resumed_index["conclusion"] == "NO_NEGATIVE_WITNESS_IN_FINITE_GRID"
    for relative in (
        "run.json",
        "index.json",
        f"attempts/{cell.cell_id}/0000.json",
        f"attempts/{cell.cell_id}/0001.json",
        f"cells/{cell.cell_id}.json",
    ):
        assert (resumed_dir / relative).read_bytes() == (
            full_dir / relative
        ).read_bytes()


def test_genuine_index_verifier_and_reference_path_safety(tmp_path: Path) -> None:
    checkpoint = tmp_path / "verified"
    index = search.run_weil_search(_raw_plan(), checkpoint)
    replay = search.verify_weil_search(checkpoint / "index.json", checkpoint)
    assert replay == {
        "classification": "REPRODUCED_EXPLORATORY_SEARCH",
        "conclusion": "NO_NEGATIVE_WITNESS_IN_FINITE_GRID",
        "verified_cells": "1",
        "verified_negative_witness_evaluations": "0",
        "same_backend_replay_only": True,
        "hypothesis_status": "UNRESOLVED",
    }

    escaped = copy.deepcopy(index)
    escaped["cells"][0]["path"] = "../outside.json"
    escaped = _rehash(escaped)
    with pytest.raises(search.WeilSearchVerificationError, match="escapes checkpoint"):
        search.verify_weil_search(escaped, checkpoint)

    promoted = copy.deepcopy(index)
    promoted["classification"] = "PROVED"
    promoted = _rehash(promoted)
    with pytest.raises(search.WeilSearchVerificationError, match="improperly promoted"):
        search.verify_weil_search(promoted, checkpoint)


def test_rehashed_semantic_cell_tamper_is_rejected(tmp_path: Path) -> None:
    checkpoint = tmp_path / "tampered"
    index = search.run_weil_search(_raw_plan(), checkpoint)
    cell_path = checkpoint / index["cells"][0]["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))

    attempt = copy.deepcopy(cell["attempts"][0])
    attempt["matrix_evidence"]["normalization"] = "A=P-R"
    cell["attempts"][0] = _rehash(attempt)
    cell = _rehash(cell)
    _write_json(cell_path, cell)

    tampered_index = copy.deepcopy(index)
    tampered_index["cells"][0]["cell_payload_sha256"] = cell["payload_sha256"]
    tampered_index = _rehash(tampered_index)
    with pytest.raises(
        search.WeilSearchVerificationError,
        match="matrix evidence does not match canonical regeneration",
    ):
        search.verify_weil_search(tampered_index, checkpoint)


def test_rehashed_malformed_attempt_is_cleanly_rejected(tmp_path: Path) -> None:
    checkpoint = tmp_path / "malformed-attempt"
    index = search.run_weil_search(_raw_plan(), checkpoint)
    cell_path = checkpoint / index["cells"][0]["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))

    attempt = copy.deepcopy(cell["attempts"][0])
    del attempt["precision_bits"]
    cell["attempts"][0] = _rehash(attempt)
    cell = _rehash(cell)
    _write_json(cell_path, cell)

    malformed_index = copy.deepcopy(index)
    malformed_index["cells"][0]["cell_payload_sha256"] = cell["payload_sha256"]
    malformed_index = _rehash(malformed_index)
    with pytest.raises(
        search.WeilSearchVerificationError,
        match="attempt evidence is malformed",
    ):
        search.verify_weil_search(malformed_index, checkpoint)


def test_inconclusive_ldlt_cannot_be_rehashed_into_a_positive_terminal(
    tmp_path: Path,
) -> None:
    plan = _raw_plan(
        cutoffs=((8, 1),),
        degrees=(16,),
        attempt_bits=(96,),
        confirmation_bits=192,
        eigenpair_count=1,
        scale_bits=(4,),
    )
    checkpoint = tmp_path / "forged-positive"
    index = search.run_weil_search(plan, checkpoint)
    cell_path = checkpoint / index["cells"][0]["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))
    assert cell["status"] == "INCONCLUSIVE_MAX_PRECISION"
    assert cell["attempts"][0]["positive_definiteness"]["classification"] == (
        "INCONCLUSIVE"
    )

    attempt = copy.deepcopy(cell["attempts"][0])
    attempt["decision"] = "FINITE_POSITIVE_CERTIFIED"
    cell["attempts"][0] = _rehash(attempt)
    cell["status"] = "FINITE_POSITIVE_CERTIFIED"
    cell["quarantined_candidate"] = None
    cell = _rehash(cell)
    _write_json(cell_path, cell)

    forged_index = copy.deepcopy(index)
    forged_index["cells"][0]["cell_payload_sha256"] = cell["payload_sha256"]
    forged_index["cells"][0]["status"] = "FINITE_POSITIVE_CERTIFIED"
    counts = forged_index["progress"]["status_counts"]
    counts["FINITE_POSITIVE_CERTIFIED"] = "1"
    counts["INCONCLUSIVE_MAX_PRECISION"] = "0"
    forged_index = _rehash(forged_index)

    with pytest.raises(
        search.WeilSearchVerificationError,
        match="cell state machine is inconsistent with replayed evidence",
    ):
        search.verify_weil_search(forged_index, checkpoint)


def test_checkpoint_rejects_resume_tampering_and_invalid_modes(tmp_path: Path) -> None:
    plan = _raw_plan()
    checkpoint = tmp_path / "checkpoint"
    search.run_weil_search(plan, checkpoint, max_cells=0)

    run_path = checkpoint / "run.json"
    run_record = json.loads(run_path.read_text(encoding="utf-8"))
    run_record["engine_id"] = "tampered-engine"
    _write_json(run_path, run_record)
    with pytest.raises(search.WeilSearchCheckpointError, match="payload hash mismatch"):
        search.run_weil_search(plan, checkpoint, resume=True)

    with pytest.raises(search.WeilSearchCheckpointError, match="without a checkpoint"):
        search.run_weil_search(plan, tmp_path / "missing", resume=True)
    with pytest.raises(TypeError, match="max_cells"):
        search.run_weil_search(plan, tmp_path / "bool-limit", max_cells=True)
    with pytest.raises(ValueError, match="nonnegative"):
        search.run_weil_search(plan, tmp_path / "negative-limit", max_cells=-1)

    semantic_dir = tmp_path / "semantic-cell"
    semantic_index = search.run_weil_search(plan, semantic_dir)
    cell_path = semantic_dir / semantic_index["cells"][0]["path"]
    cell = json.loads(cell_path.read_text(encoding="utf-8"))
    attempt = copy.deepcopy(cell["attempts"][0])
    attempt["decision"] = "INCONCLUSIVE_MAX_PRECISION"
    cell["attempts"][0] = _rehash(attempt)
    cell["status"] = "INCONCLUSIVE_MAX_PRECISION"
    cell = _rehash(cell)
    _write_json(cell_path, cell)
    with pytest.raises(
        search.WeilSearchCheckpointError,
        match="search decision does not follow",
    ):
        search.run_weil_search(plan, semantic_dir, resume=True)
