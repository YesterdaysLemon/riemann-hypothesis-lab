from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
import tarfile
from typing import Any

import pytest

import riemann_lab.weil_search as search
import riemann_lab.weil_transition as transition
from tools import package_weil_transition_v2 as packager


REPOSITORY = "example/research"
RELEASE_TAG = "weil-transition-q7-q9-v2"
PRODUCER_COMMIT = "1" * 40
FROZEN_PLAN_PATH = Path("plans/weil-transition-q7-q9-v2.json")


def _hashed(body: dict[str, Any]) -> dict[str, Any]:
    return {**body, "payload_sha256": packager._content_sha256(body)}


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _ball(midpoint: int) -> dict[str, Any]:
    return {
        "display": f"mock-dyadic-{midpoint}",
        "dyadic": {
            "mid_mantissa": str(midpoint),
            "mid_exponent": "-10",
            "radius_mantissa": "1",
            "radius_exponent": "-18",
        },
        "is_exact": False,
    }


def _rump(degree: int) -> dict[str, Any]:
    eigenvalues = [
        {
            "index": str(index),
            "real": _ball(1000 * (index + 1)),
            "imaginary": _ball(0),
            "imaginary_contains_zero": True,
            "real_part_is_positive": True,
        }
        for index in range(2 * degree + 1)
    ]
    return {
        "algorithm": "flint-arb_mat-eig-rump",
        "classification": "POSITIVE_SPECTRUM",
        "eigenvalues": eigenvalues,
        "ordering": "ascending by pairwise-separated real enclosures",
        "smallest_eigenvalue": {
            "ordered_index": "0",
            "real": eigenvalues[0]["real"],
            "imaginary": eigenvalues[0]["imaginary"],
            "separated_from_all_others": True,
        },
    }


def _positive_artifact(
    contract: transition.WeilTransitionCellContract,
    ordinal: int,
) -> dict[str, Any]:
    attempt = _hashed(
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
                "pivots": [{"value": _ball(1)}],
            },
            "secondary_rump_spectrum": (
                _rump(contract.degree)
                if ordinal == 0
                else {
                    "classification": "INCONCLUSIVE",
                    "failure_code": "RUMP_EIGENVALUE_ISOLATION_FAILED",
                }
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


def _parity_artifact(contract: transition.WeilTransitionCellContract) -> dict[str, Any]:
    return _hashed(
        {
            "schema": "rh-lab/weil-parity-audit/v1",
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "cell": contract.v1_cell.to_record(),
            "decision": "AUDIT_PASSED_EXPLORATORY",
        }
    )


def _nesting_artifact(
    lower: transition.WeilTransitionCellContract,
    higher: transition.WeilTransitionCellContract,
) -> dict[str, Any]:
    return _hashed(
        {
            "schema": "rh-lab/weil-degree-nesting-audit/v1",
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "lower_cell": lower.v1_cell.to_record(),
            "higher_cell": higher.v1_cell.to_record(),
            "decision": "AUDIT_PASSED_EXPLORATORY",
        }
    )


def _frozen_complete_tree(root: Path) -> Path:
    """Write a small-byte, semantically canonical version of all 81 cells."""

    evidence = root / "weil-transition-q7-q9-v2"
    plan = transition.load_transition_plan(FROZEN_PLAN_PATH)
    contracts = transition._contracts_from_plan(plan)
    _write_json(evidence / "plan.json", plan)
    _write_json(evidence / "run.json", transition._run_manifest(plan))

    completed: dict[str, tuple[dict[str, Any], None]] = {}
    for ordinal, contract in enumerate(contracts):
        artifact = _positive_artifact(contract, ordinal)
        _write_json(evidence / "cells" / f"{contract.contract_id}.json", artifact)
        completed[contract.contract_id] = (artifact, None)
    _write_json(
        evidence / "index.json",
        transition._index_payload(plan, contracts, completed),
    )

    exact = [contract for contract in contracts if contract.position == "exact"]
    assert len(exact) == 9
    for contract in exact:
        _write_json(
            evidence / "audits" / "parity" / f"{contract.contract_id}.json",
            _parity_artifact(contract),
        )
    for q in (7, 8, 9):
        ladder = sorted(
            (contract for contract in exact if contract.q == q),
            key=lambda contract: contract.degree,
        )
        for lower, higher in zip(ladder, ladder[1:]):
            _write_json(
                evidence
                / "audits"
                / "nesting"
                / f"q-{q}-n-{lower.degree}-to-{higher.degree}.json",
                _nesting_artifact(lower, higher),
            )

    # Attempts are deliberately malformed and ignored. Terminal cell files
    # already embed every ordered attempt used by the summary generator.
    attempt = evidence / "attempts" / contracts[0].contract_id / "0000.json"
    attempt.parent.mkdir(parents=True)
    attempt.write_text("not primary evidence\n", encoding="utf-8")
    return evidence


@pytest.fixture(autouse=True)
def cheap_audit_replay(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Keep package tests fast while proving both verifier paths are called."""

    calls: list[str] = []

    def parity(artifact: dict[str, Any]) -> dict[str, Any]:
        calls.append("PARITY")
        return {
            "classification": "REPRODUCED_EXPLORATORY_AUDIT",
            "audit_kind": "PARITY",
            "decision": artifact["decision"],
            "hypothesis_status": "UNRESOLVED",
        }

    def nesting(artifact: dict[str, Any]) -> dict[str, Any]:
        calls.append("DEGREE_NESTING")
        return {
            "classification": "REPRODUCED_EXPLORATORY_AUDIT",
            "audit_kind": "DEGREE_NESTING",
            "decision": artifact["decision"],
            "hypothesis_status": "UNRESOLVED",
        }

    monkeypatch.setattr(packager.weil_audits, "verify_parity_audit", parity)
    monkeypatch.setattr(
        packager.weil_audits,
        "verify_degree_nesting_audit",
        nesting,
    )
    return calls


def _package(evidence: Path, archive: Path) -> dict[str, Any]:
    return packager.package_release(
        evidence,
        archive,
        repository=REPOSITORY,
        release_tag=RELEASE_TAG,
        producer_commit=PRODUCER_COMMIT,
    )


def test_package_is_canonical_bound_and_excludes_attempts(
    tmp_path: Path,
    cheap_audit_replay: list[str],
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "source")
    archive = tmp_path / "dist" / "evidence.tar"

    result = _package(evidence, archive)
    manifest_path = evidence.with_name(f"{evidence.name}.release.json")
    checksum_path = archive.with_name(f"{archive.name}.sha256")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert cheap_audit_replay.count("PARITY") == 9
    assert cheap_audit_replay.count("DEGREE_NESTING") == 6
    assert result["archive_sha256"] == packager._file_sha256(archive)
    assert checksum_path.read_text(encoding="ascii") == (
        f"{result['archive_sha256']}  {archive.name}\n"
    )
    manifest_body = {
        key: value for key, value in manifest.items() if key != "payload_sha256"
    }
    assert manifest["payload_sha256"] == packager._content_sha256(manifest_body)
    assert manifest["classification"] == "EXPLORATORY"
    assert manifest["hypothesis_status"] == "UNRESOLVED"
    assert manifest["producer_git_commit"] == PRODUCER_COMMIT
    assert manifest["archive"]["sha256"] == result["archive_sha256"]
    assert manifest["archive"]["excluded_paths"] == ["attempts/"]
    assert manifest["evidence"]["counts"] == {
        "cells": "81",
        "confirmations": "0",
        "nesting_audits": "6",
        "parity_audits": "9",
    }
    assert len(manifest["evidence"]["audits"]) == 15
    assert manifest["evidence"]["transition_summary"]["conclusion"] == (
        "NO_CERTIFIED_NEGATIVE_WITNESS_FOUND_IN_FINITE_TRANSITION_BATCH"
    )

    with tarfile.open(archive, "r:") as handle:
        members = handle.getmembers()
        names = [member.name for member in members]
        assert names[:6] == [
            "weil-transition-q7-q9-v2",
            "weil-transition-q7-q9-v2/audits",
            "weil-transition-q7-q9-v2/audits/nesting",
            "weil-transition-q7-q9-v2/audits/parity",
            "weil-transition-q7-q9-v2/cells",
            "weil-transition-q7-q9-v2/confirmations",
        ]
        file_names = [member.name for member in members if member.isfile()]
        assert file_names == sorted(file_names)
        assert not any("/attempts/" in name for name in names)
        assert "weil-transition-q7-q9-v2/confirmations" in names
        for member in members:
            assert member.mtime == 0
            assert member.uid == member.gid == 0
            assert member.uname == member.gname == ""
            assert member.mode == (0o755 if member.isdir() else 0o644)


def test_identical_trees_produce_identical_archive_sidecar_and_manifest(
    tmp_path: Path,
) -> None:
    first = _frozen_complete_tree(tmp_path / "first")
    second_parent = tmp_path / "second"
    second = second_parent / first.name
    shutil.copytree(first, second)
    first_archive = tmp_path / "dist-a" / "evidence.tar"
    second_archive = tmp_path / "dist-b" / "evidence.tar"

    first_result = _package(first, first_archive)
    second_result = _package(second, second_archive)

    assert first_archive.read_bytes() == second_archive.read_bytes()
    assert first_result["archive_sha256"] == second_result["archive_sha256"]
    assert first_result["content_set_sha256"] == second_result["content_set_sha256"]
    assert first_archive.with_name(f"{first_archive.name}.sha256").read_bytes() == (
        second_archive.with_name(f"{second_archive.name}.sha256").read_bytes()
    )
    assert first.with_name(f"{first.name}.release.json").read_bytes() == (
        second.with_name(f"{second.name}.release.json").read_bytes()
    )


@pytest.mark.parametrize(
    "mutation",
    ("incomplete", "unreferenced-cell", "unexpected-root", "unbound-confirmation"),
)
def test_incomplete_or_non_allowlisted_evidence_is_rejected(
    tmp_path: Path, mutation: str
) -> None:
    evidence = _frozen_complete_tree(tmp_path / mutation)
    index_path = evidence / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    if mutation == "incomplete":
        body = {
            key: copy.deepcopy(value)
            for key, value in index.items()
            if key != "payload_sha256"
        }
        body["cells"].pop()
        body["progress"]["completed_cells"] = "80"
        body["missing_contract_ids"] = ["missing"]
        _write_json(index_path, _hashed(body))
    elif mutation == "unreferenced-cell":
        _write_json(
            evidence / "cells" / "orphan.json",
            _hashed(
                {
                    "schema": "rh-lab/weil-search-cell/v1",
                    "classification": "EXPLORATORY",
                    "hypothesis_status": "UNRESOLVED",
                }
            ),
        )
    elif mutation == "unexpected-root":
        (evidence / "notes.txt").write_text("not allowlisted\n", encoding="utf-8")
    else:
        _write_json(
            evidence / "confirmations" / "orphan.json",
            _hashed(
                {
                    "schema": "rh-lab/weil-transition-dedicated-confirmation/v2",
                    "classification": "EXPLORATORY",
                    "hypothesis_status": "UNRESOLVED",
                }
            ),
        )

    with pytest.raises(packager.PackagingError):
        _package(evidence, tmp_path / "dist" / "evidence.tar")


def test_rehashed_fake_plan_and_promoted_conclusion_are_rejected(
    tmp_path: Path,
) -> None:
    fake_plan_tree = _frozen_complete_tree(tmp_path / "fake-plan")
    plan_path = fake_plan_tree / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_body = {
        key: copy.deepcopy(value)
        for key, value in plan.items()
        if key != "plan_sha256"
    }
    plan_body["frozen_batch"]["cell_count"] = "80"
    plan_body["plan_sha256"] = packager._content_sha256(plan_body)
    _write_json(plan_path, plan_body)
    with pytest.raises(packager.PackagingError, match="canonical|frozen"):
        _package(fake_plan_tree, tmp_path / "fake-plan.tar")

    promoted_tree = _frozen_complete_tree(tmp_path / "promoted")
    index_path = promoted_tree / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index_body = {
        key: copy.deepcopy(value)
        for key, value in index.items()
        if key != "payload_sha256"
    }
    index_body["conclusion"] = "RIEMANN_HYPOTHESIS_PROVED"
    _write_json(index_path, _hashed(index_body))
    with pytest.raises(packager.PackagingError, match="promoted or fake conclusion"):
        _package(promoted_tree, tmp_path / "promoted.tar")


@pytest.mark.parametrize("mutation", ("missing", "duplicate-target"))
def test_exact_audit_coverage_is_required(tmp_path: Path, mutation: str) -> None:
    evidence = _frozen_complete_tree(tmp_path / mutation)
    parity_paths = sorted((evidence / "audits" / "parity").glob("*.json"))
    if mutation == "missing":
        parity_paths[-1].unlink()
    else:
        first = json.loads(parity_paths[0].read_text(encoding="utf-8"))
        second = json.loads(parity_paths[1].read_text(encoding="utf-8"))
        body = {
            key: copy.deepcopy(value)
            for key, value in second.items()
            if key != "payload_sha256"
        }
        body["cell"] = first["cell"]
        _write_json(parity_paths[1], _hashed(body))

    with pytest.raises(packager.PackagingError, match="audits|cover"):
        _package(evidence, tmp_path / "dist" / "evidence.tar")


def test_live_source_replacement_after_validation_cannot_poison_archive(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "replacement")
    archive = tmp_path / "dist" / "evidence.tar"
    cell_path = sorted((evidence / "cells").glob("*.json"))[0]
    relative = cell_path.relative_to(evidence).as_posix()
    original = cell_path.read_bytes()
    original_snapshot = packager._snapshot
    replaced = False

    def replace_live_before_snapshot(
        files: Any,
    ) -> tuple[packager.FileRecord, ...]:
        nonlocal replaced
        if not replaced:
            # This is the old collect-then-read race: after semantic validation,
            # replace one live source file before the content snapshot is read.
            cell_path.write_text("not-json\n", encoding="utf-8")
            replaced = True
        return original_snapshot(files)

    monkeypatch.setattr(packager, "_snapshot", replace_live_before_snapshot)

    _package(evidence, archive)

    assert replaced
    assert cell_path.read_text(encoding="utf-8") == "not-json\n"
    with tarfile.open(archive, "r:") as handle:
        extracted = handle.extractfile(f"{packager.ARCHIVE_ROOT}/{relative}")
        assert extracted is not None
        archived = extracted.read()
    assert archived == original
    assert isinstance(json.loads(archived), dict)


def test_final_verification_rejects_archive_truncation_during_digest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "archive-tamper")
    archive = tmp_path / "dist" / "evidence.tar"
    original_digest = packager._digest_open_file
    digest_calls = 0

    def truncate_after_initial_digest(handle: Any) -> tuple[int, str]:
        nonlocal digest_calls
        result = original_digest(handle)
        digest_calls += 1
        if digest_calls == 1:
            Path(handle.name).write_bytes(b"")
        return result

    monkeypatch.setattr(packager, "_digest_open_file", truncate_after_initial_digest)

    with pytest.raises(
        packager.PackagingError,
        match="cannot reopen canonical USTAR archive",
    ):
        _package(evidence, archive)
    assert digest_calls == 1
    assert not archive.exists()
    assert not archive.with_name(f"{archive.name}.sha256").exists()
    assert not evidence.with_name(f"{evidence.name}.release.json").exists()


def test_no_late_archive_rehash_can_replace_verified_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "late-hash")
    archive = tmp_path / "dist" / "evidence.tar"
    original_hash = packager._file_sha256
    archive_hash_calls: list[Path] = []

    def truncate_if_archive_is_rehashed(path: Path) -> str:
        if path.suffix.lower() == ".tar":
            archive_hash_calls.append(path)
            path.write_bytes(b"")
        return original_hash(path)

    monkeypatch.setattr(packager, "_file_sha256", truncate_if_archive_is_rehashed)

    result = _package(evidence, archive)

    assert archive_hash_calls == []
    assert archive.stat().st_size > 0
    assert result["archive_sha256"] == original_hash(archive)
    with tarfile.open(archive, "r:") as handle:
        assert len(handle.getmembers()) == 105


def test_temporary_archive_mutation_at_old_prelink_gap_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "prelink-tamper")
    archive = tmp_path / "dist" / "evidence.tar"
    original_link = packager.os.link
    tampered = False

    def tamper_then_link(source: Any, destination: Any) -> None:
        nonlocal tampered
        data = bytearray(Path(source).read_bytes())
        assert data[-1] == 0
        data[-1] = 1
        Path(source).write_bytes(data)
        tampered = True
        original_link(source, destination)

    monkeypatch.setattr(packager.os, "link", tamper_then_link)

    with pytest.raises(packager.PackagingError, match="trailer padding"):
        _package(evidence, archive)
    assert tampered
    assert not archive.exists()
    assert not archive.with_name(f"{archive.name}.sha256").exists()
    assert not evidence.with_name(f"{evidence.name}.release.json").exists()


def test_published_archive_mutation_after_authoritative_verify_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    evidence = _frozen_complete_tree(tmp_path / "published-tamper")
    archive = tmp_path / "dist" / "evidence.tar"
    original_verify = packager._verify_written_archive
    attempted = False

    def truncate_after_verify(*args: Any, **kwargs: Any) -> tuple[int, str]:
        nonlocal attempted
        result = original_verify(*args, **kwargs)
        attempted = True
        Path(args[0]).write_bytes(b"")
        return result

    monkeypatch.setattr(
        packager,
        "_verify_written_archive",
        truncate_after_verify,
    )

    with pytest.raises(
        packager.PackagingError,
        match="unexpected byte length|exclusively lock",
    ):
        _package(evidence, archive)
    assert attempted
    assert not archive.exists()
    assert not archive.with_name(f"{archive.name}.sha256").exists()
    assert not evidence.with_name(f"{evidence.name}.release.json").exists()


@pytest.mark.parametrize("target_kind", ("root", "directory", "file"))
def test_source_reparse_points_are_portably_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    target_kind: str,
) -> None:
    evidence = _frozen_complete_tree(tmp_path / f"reparse-{target_kind}")
    target = {
        "root": evidence,
        "directory": evidence / "cells",
        "file": sorted((evidence / "cells").glob("*.json"))[0],
    }[target_kind]
    original_probe = packager._is_reparse_or_link

    def portable_reparse_probe(path: Path) -> bool:
        return path == target or original_probe(path)

    monkeypatch.setattr(packager, "_is_reparse_or_link", portable_reparse_probe)

    with pytest.raises(packager.PackagingError, match="junction|reparse point"):
        _package(evidence, tmp_path / "dist" / "evidence.tar")


@pytest.mark.parametrize("mutation", ("nonzero-trailer", "extra-zero-block"))
def test_canonical_ustar_end_padding_is_required(
    tmp_path: Path,
    mutation: str,
) -> None:
    evidence = _frozen_complete_tree(tmp_path / mutation)
    archive = tmp_path / "dist" / "evidence.tar"
    _package(evidence, archive)
    validated = packager._collect_evidence(evidence)
    records = packager._snapshot(validated.files)

    if mutation == "nonzero-trailer":
        data = bytearray(archive.read_bytes())
        assert data[-1] == 0
        data[-1] = 1
        archive.write_bytes(data)
        message = "trailer padding"
    else:
        with archive.open("ab") as handle:
            handle.write(b"\0" * tarfile.BLOCKSIZE)
        message = "unexpected byte length"

    with pytest.raises(packager.PackagingError, match=message):
        packager._verify_written_archive(archive, records)


def test_refuses_to_overwrite_any_release_output(tmp_path: Path) -> None:
    evidence = _frozen_complete_tree(tmp_path / "overwrite")
    archive = tmp_path / "dist" / "evidence.tar"
    archive.parent.mkdir(parents=True)
    archive.write_bytes(b"existing")

    with pytest.raises(packager.PackagingError, match="refusing to overwrite"):
        _package(evidence, archive)
    assert archive.read_bytes() == b"existing"
