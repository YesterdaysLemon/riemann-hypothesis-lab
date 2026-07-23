"""Build a deterministic public evidence archive for the v2 Weil batch.

The terminal cell artifacts embed their ordered numerical attempts, so the
checkpoint ``attempts/`` tree is deliberately excluded.  The attached TAR is
an uncompressed, metadata-normalized USTAR archive.  Keeping compression out
of the canonical object avoids making its digest depend on a zlib version.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import tarfile
import tempfile
from typing import Any, BinaryIO, Iterator
from urllib.parse import quote

from riemann_lab import weil_audits
from riemann_lab import weil_transition as transition
from riemann_lab import weil_transition_summary as transition_summary


ARCHIVE_ROOT = "weil-transition-q7-q9-v2"
BUILDER_ID = "canonical-weil-transition-ustar-v1"
MANIFEST_SCHEMA = "rh-lab/weil-transition-release-manifest/v1"
CLASSIFICATION = "EXPLORATORY"
HYPOTHESIS_STATUS = "UNRESOLVED"
FROZEN_PLAN_SHA256 = (
    "33185881cc619fda99b7835f5e71422b0c5cf745cb5411162638a8e08dd05336"
)
_ALLOWED_COMPLETE_CONCLUSIONS = {
    "NO_CERTIFIED_NEGATIVE_WITNESS_FOUND_IN_FINITE_TRANSITION_BATCH",
    "QUARANTINED_CONFLICTING_NEGATIVE_OBSERVATION_IN_FINITE_BATCH",
    "QUARANTINED_NEGATIVE_CANDIDATE_IN_FINITE_BATCH",
}
INDEX_SCHEMA = "rh-lab/weil-transition-index/v2"
PLAN_SCHEMA = "rh-lab/weil-transition-plan/v2"
RUN_SCHEMA = "rh-lab/weil-transition-run/v2"
CELL_SCHEMA = "rh-lab/weil-search-cell/v1"
CONFIRMATION_SCHEMA = "rh-lab/weil-transition-dedicated-confirmation/v2"

_ROOT_FILES = ("index.json", "plan.json", "run.json")
_ARCHIVE_DIRECTORIES = (
    "",
    "audits",
    "audits/nesting",
    "audits/parity",
    "cells",
    "confirmations",
)
_ALLOWED_SOURCE_DIRECTORIES = {
    ("attempts",),
    ("audits",),
    ("audits", "nesting"),
    ("audits", "parity"),
    ("cells",),
    ("confirmations",),
}
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT_RE = re.compile(r"[0-9a-fA-F]{40}\Z")
_REPOSITORY_RE = re.compile(
    r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+\Z"
)
_REPARSE_POINT_ATTRIBUTE = getattr(
    stat,
    "FILE_ATTRIBUTE_REPARSE_POINT",
    0x400,
)


class PackagingError(ValueError):
    """Raised when the evidence tree cannot be packaged without ambiguity."""


@dataclass(frozen=True)
class FileRecord:
    """A content-bound archive member, relative to the evidence root."""

    path: str
    size: int
    sha256: str

    def to_record(self) -> dict[str, str]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": str(self.size),
        }


@dataclass(frozen=True)
class EvidenceSet:
    """Validated files and semantic bindings for one complete batch."""

    source: Path
    files: tuple[tuple[str, Path], ...]
    index: dict[str, Any]
    plan: dict[str, Any]
    run: dict[str, Any]
    summary: dict[str, Any]
    audit_records: tuple[dict[str, Any], ...]
    cell_count: int
    confirmation_count: int
    parity_audit_count: int
    nesting_audit_count: int


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _content_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PackagingError(f"cannot read {label}: {path}") from exc
    if not isinstance(value, dict):
        raise PackagingError(f"{label} must be a JSON object")
    return value


def _load_referenced_object(source: Path, relative: str) -> dict[str, Any]:
    path = source / Path(*PurePosixPath(relative).parts)
    if path.is_symlink() or not path.is_file():
        raise PackagingError(f"referenced evidence file is missing: {relative}")
    return _load_object(path, relative)


def _require_exploratory(record: Mapping[str, Any], label: str) -> None:
    if record.get("classification") != CLASSIFICATION:
        raise PackagingError(f"{label} is not EXPLORATORY")
    if record.get("hypothesis_status") != HYPOTHESIS_STATUS:
        raise PackagingError(f"{label} does not keep RH UNRESOLVED")


def _require_payload_hash(record: Mapping[str, Any], label: str) -> str:
    supplied = record.get("payload_sha256")
    if not isinstance(supplied, str) or not _SHA256_RE.fullmatch(supplied):
        raise PackagingError(f"{label} lacks a canonical payload SHA-256")
    body = {key: value for key, value in record.items() if key != "payload_sha256"}
    if supplied != _content_sha256(body):
        raise PackagingError(f"{label} payload SHA-256 does not match")
    return supplied


def _require_plan_hash(plan: Mapping[str, Any]) -> str:
    supplied = plan.get("plan_sha256")
    if not isinstance(supplied, str) or not _SHA256_RE.fullmatch(supplied):
        raise PackagingError("plan lacks a canonical plan SHA-256")
    body = {key: value for key, value in plan.items() if key != "plan_sha256"}
    if supplied != _content_sha256(body):
        raise PackagingError("plan SHA-256 does not match")
    return supplied


def _canonical_count(value: Any, label: str) -> int:
    if not isinstance(value, str) or not value.isdigit():
        raise PackagingError(f"{label} must be a canonical decimal string")
    if len(value) > 1 and value.startswith("0"):
        raise PackagingError(f"{label} must be a canonical decimal string")
    return int(value)


def _reference_path(value: Any, directory: str, label: str) -> str:
    if not isinstance(value, str) or "\\" in value:
        raise PackagingError(f"{label} must be a POSIX relative JSON path")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.parts != (directory, path.name)
        or path.name in {"", ".", ".."}
        or path.suffix != ".json"
        or any(part in {".", ".."} for part in path.parts)
    ):
        raise PackagingError(f"{label} must be directly below {directory}/")
    return path.as_posix()


def _is_reparse_or_link(path: Path) -> bool:
    """Return whether a path is a symlink, junction, or Windows reparse point."""

    metadata = path.lstat()
    if stat.S_ISLNK(metadata.st_mode):
        return True
    attributes = getattr(metadata, "st_file_attributes", 0)
    if attributes & _REPARSE_POINT_ATTRIBUTE:
        return True
    is_junction = getattr(path, "is_junction", None)
    return bool(is_junction()) if callable(is_junction) else False


def _validate_source_shape(source: Path) -> tuple[set[str], set[str], set[str]]:
    """Return all non-attempt files and the two audit file sets."""

    try:
        if _is_reparse_or_link(source):
            raise PackagingError(
                "symbolic links, junctions, and reparse points are not "
                f"permitted: {source}"
            )
    except OSError as exc:
        raise PackagingError(f"cannot inspect evidence source: {source}") from exc

    actual_files: set[str] = set()
    parity: set[str] = set()
    nesting: set[str] = set()

    def visit(directory: Path) -> None:
        try:
            children = sorted(directory.iterdir(), key=lambda path: path.name)
        except OSError as exc:
            raise PackagingError(
                f"cannot enumerate evidence directory: {directory}"
            ) from exc
        for path in children:
            relative = path.relative_to(source)
            parts = relative.parts
            try:
                linklike = _is_reparse_or_link(path)
            except OSError as exc:
                raise PackagingError(
                    f"cannot inspect evidence entry: {relative}"
                ) from exc
            if linklike:
                raise PackagingError(
                    "symbolic links, junctions, and reparse points are not "
                    f"permitted: {relative}"
                )
            if path.is_dir():
                if parts == ("attempts",):
                    continue
                if parts not in _ALLOWED_SOURCE_DIRECTORIES:
                    raise PackagingError(f"unexpected directory: {relative}")
                visit(path)
                continue
            if parts[0] == "attempts":
                raise PackagingError("attempts must be a directory")
            if not path.is_file():
                raise PackagingError(f"unsupported filesystem entry: {relative}")
            relative_posix = relative.as_posix()
            actual_files.add(relative_posix)
            if len(parts) == 3 and parts[:2] == ("audits", "parity"):
                if path.suffix != ".json":
                    raise PackagingError(f"parity audit is not JSON: {relative}")
                parity.add(relative_posix)
            elif len(parts) == 3 and parts[:2] == ("audits", "nesting"):
                if path.suffix != ".json":
                    raise PackagingError(f"nesting audit is not JSON: {relative}")
                nesting.add(relative_posix)

    visit(source)
    return actual_files, parity, nesting


def _collect_evidence(source: Path) -> EvidenceSet:
    """Validate completeness, references, hashes, and the exact allowlist."""

    if not source.is_dir():
        raise PackagingError(f"evidence source is not a directory: {source}")
    for name in _ROOT_FILES:
        path = source / name
        if path.is_symlink() or not path.is_file():
            raise PackagingError(f"required evidence file is missing: {name}")

    # Reject links and off-contract paths before any semantic helper follows a
    # checkpoint reference.
    actual_files, parity_paths, nesting_paths = _validate_source_shape(source)

    index = _load_object(source / "index.json", "transition index")
    plan = _load_object(source / "plan.json", "transition plan")
    run = _load_object(source / "run.json", "run manifest")
    for record, label in (
        (index, "transition index"),
        (plan, "transition plan"),
        (run, "run manifest"),
    ):
        _require_exploratory(record, label)
    if index.get("schema") != INDEX_SCHEMA:
        raise PackagingError("unexpected transition index schema")
    if plan.get("schema") != PLAN_SCHEMA:
        raise PackagingError("unexpected transition plan schema")
    if run.get("schema") != RUN_SCHEMA:
        raise PackagingError("unexpected run manifest schema")
    _require_payload_hash(index, "transition index")
    plan_sha256 = _require_plan_hash(plan)
    _require_payload_hash(run, "run manifest")
    try:
        canonical_plan = transition.canonicalize_transition_plan(plan)
    except transition.WeilTransitionPlanError as exc:
        raise PackagingError("stored transition plan is not canonical") from exc
    if plan != canonical_plan:
        raise PackagingError("stored transition plan is not canonical")
    if (
        plan_sha256 != FROZEN_PLAN_SHA256
        or canonical_plan.get("frozen_batch") != transition.FROZEN_Q7_Q9_BATCH
        or len(canonical_plan["cells"]) != 81
    ):
        raise PackagingError("release requires the exact frozen 81-cell v2 plan")
    if run != transition._run_manifest(canonical_plan):
        raise PackagingError("run manifest is not canonical for the frozen plan")
    if index.get("conclusion") not in _ALLOWED_COMPLETE_CONCLUSIONS:
        raise PackagingError("transition index has a promoted or fake conclusion")
    try:
        summary = transition_summary.generate_weil_transition_summary(source)
    except transition_summary.WeilTransitionSummaryError as exc:
        raise PackagingError(
            "frozen transition summary cannot canonically regenerate"
        ) from exc
    _require_payload_hash(summary, "transition summary")
    if (
        summary.get("schema") != transition_summary.SUMMARY_SCHEMA
        or summary.get("classification") != CLASSIFICATION
        or summary.get("hypothesis_status") != HYPOTHESIS_STATUS
        or summary.get("frozen_batch") != transition.FROZEN_Q7_Q9_BATCH
        or summary.get("conclusion") != index.get("conclusion")
        or summary.get("counts", {}).get("cells") != "81"
        or summary.get("source", {}).get("plan", {}).get("plan_sha256")
        != FROZEN_PLAN_SHA256
        or summary.get("source", {}).get("index", {}).get("payload_sha256")
        != index["payload_sha256"]
    ):
        raise PackagingError("transition summary changed the frozen experiment")

    index_plan = index.get("plan")
    if not isinstance(index_plan, Mapping) or dict(index_plan) != {
        "path": "plan.json",
        "plan_sha256": plan_sha256,
    }:
        raise PackagingError("index does not bind the packaged plan")
    index_run = index.get("run")
    if not isinstance(index_run, Mapping) or dict(index_run) != {
        "path": "run.json"
    }:
        raise PackagingError("index does not reference the packaged run manifest")
    run_plan = run.get("plan")
    if not isinstance(run_plan, Mapping) or dict(run_plan) != {
        "path": "plan.json",
        "plan_sha256": plan_sha256,
    }:
        raise PackagingError("run manifest does not bind the packaged plan")

    progress = index.get("progress")
    if not isinstance(progress, Mapping):
        raise PackagingError("index progress is missing")
    planned = _canonical_count(progress.get("planned_cells"), "planned cells")
    completed = _canonical_count(
        progress.get("completed_cells"), "completed cells"
    )
    if planned != 81 or completed != 81:
        raise PackagingError("transition index is incomplete")
    missing = index.get("missing_contract_ids")
    if not isinstance(missing, list) or missing:
        raise PackagingError("transition index still has missing contracts")

    references = index.get("cells")
    if not isinstance(references, list) or len(references) != planned:
        raise PackagingError("index cell references do not match completed cells")
    cell_paths: set[str] = set()
    confirmation_paths: set[str] = set()
    contract_ids: set[str] = set()
    for position, reference in enumerate(references):
        if not isinstance(reference, Mapping):
            raise PackagingError(f"cell reference {position} is not an object")
        contract_id = reference.get("contract_id")
        if not isinstance(contract_id, str) or not contract_id:
            raise PackagingError(f"cell reference {position} lacks a contract ID")
        if contract_id in contract_ids:
            raise PackagingError(f"duplicate contract ID: {contract_id}")
        contract_ids.add(contract_id)
        cell_path = _reference_path(
            reference.get("path"), "cells", f"cell reference {position}"
        )
        if cell_path in cell_paths:
            raise PackagingError(f"duplicate cell path: {cell_path}")
        cell_paths.add(cell_path)
        cell = _load_referenced_object(source, cell_path)
        _require_exploratory(cell, cell_path)
        if cell.get("schema") != CELL_SCHEMA:
            raise PackagingError(f"unexpected cell schema: {cell_path}")
        cell_hash = _require_payload_hash(cell, cell_path)
        if reference.get("v1_cell_payload_sha256") != cell_hash:
            raise PackagingError(f"index hash does not bind {cell_path}")

        confirmation_reference = reference.get("dedicated_confirmation")
        if confirmation_reference is None:
            continue
        if not isinstance(confirmation_reference, Mapping):
            raise PackagingError(
                f"confirmation reference for {contract_id} is malformed"
            )
        confirmation_path = _reference_path(
            confirmation_reference.get("path"),
            "confirmations",
            f"confirmation reference for {contract_id}",
        )
        if confirmation_path in confirmation_paths:
            raise PackagingError(f"duplicate confirmation path: {confirmation_path}")
        confirmation_paths.add(confirmation_path)
        confirmation = _load_referenced_object(source, confirmation_path)
        _require_exploratory(confirmation, confirmation_path)
        if confirmation.get("schema") != CONFIRMATION_SCHEMA:
            raise PackagingError(
                f"unexpected confirmation schema: {confirmation_path}"
            )
        confirmation_hash = _require_payload_hash(
            confirmation, confirmation_path
        )
        if confirmation_reference.get("payload_sha256") != confirmation_hash:
            raise PackagingError(
                f"index hash does not bind {confirmation_path}"
            )

    status_counts = progress.get("status_counts")
    if not isinstance(status_counts, Mapping):
        raise PackagingError("index status counts are missing")
    total = sum(
        _canonical_count(value, f"status count {key}")
        for key, value in status_counts.items()
    )
    if total != completed:
        raise PackagingError("index status counts do not sum to completed cells")
    plan_cells = canonical_plan.get("cells")
    if not isinstance(plan_cells, list) or len(plan_cells) != 81:
        raise PackagingError("transition plan cell count does not match the index")

    if len(parity_paths) != 9 or len(nesting_paths) != 6:
        raise PackagingError(
            "release requires exactly 9 parity and 6 nesting audits"
        )
    audit_paths = parity_paths | nesting_paths
    expected_files = set(_ROOT_FILES) | cell_paths | confirmation_paths | audit_paths
    unexpected = sorted(actual_files - expected_files)
    missing_files = sorted(expected_files - actual_files)
    if unexpected:
        raise PackagingError(f"unexpected evidence file: {unexpected[0]}")
    if missing_files:
        raise PackagingError(f"referenced evidence file is missing: {missing_files[0]}")
    if {path for path in actual_files if path.startswith("cells/")} != cell_paths:
        raise PackagingError("cell directory does not exactly match the index")
    if {
        path for path in actual_files if path.startswith("confirmations/")
    } != confirmation_paths:
        raise PackagingError(
            "confirmation directory does not exactly match the index"
        )

    contracts = transition._contracts_from_plan(canonical_plan)
    exact_contracts = [contract for contract in contracts if contract.position == "exact"]
    if len(exact_contracts) != 9:
        raise PackagingError("frozen plan does not have nine exact transition cells")
    expected_parity: dict[str, str] = {
        _canonical_json_bytes(contract.v1_cell.to_record()).decode("utf-8"):
        contract.contract_id
        for contract in exact_contracts
    }
    exact_by_q: dict[int, list[transition.WeilTransitionCellContract]] = {}
    for contract in exact_contracts:
        exact_by_q.setdefault(contract.q, []).append(contract)
    expected_nesting: dict[str, tuple[str, str]] = {}
    for q in (7, 8, 9):
        ladder = sorted(exact_by_q.get(q, []), key=lambda contract: contract.degree)
        if len(ladder) != 3:
            raise PackagingError(f"frozen q={q} exact degree ladder changed")
        for lower, higher in zip(ladder, ladder[1:]):
            key = _canonical_json_bytes(
                {
                    "lower_cell": lower.v1_cell.to_record(),
                    "higher_cell": higher.v1_cell.to_record(),
                }
            ).decode("utf-8")
            expected_nesting[key] = (lower.contract_id, higher.contract_id)

    seen_parity: set[str] = set()
    seen_nesting: set[str] = set()
    audit_records: list[dict[str, Any]] = []
    for audit_path in sorted(audit_paths):
        audit = _load_referenced_object(source, audit_path)
        _require_exploratory(audit, audit_path)
        audit_hash = _require_payload_hash(audit, audit_path)
        is_parity = audit_path.startswith("audits/parity/")
        expected_schema = (
            weil_audits.PARITY_AUDIT_SCHEMA
            if is_parity
            else weil_audits.NESTING_AUDIT_SCHEMA
        )
        if audit.get("schema") != expected_schema:
            raise PackagingError(f"unexpected audit schema: {audit_path}")
        if is_parity:
            target_key = _canonical_json_bytes(audit.get("cell")).decode("utf-8")
            contract_id = expected_parity.get(target_key)
            if contract_id is None or target_key in seen_parity:
                raise PackagingError(
                    "parity audits do not exactly cover the nine exact cells"
                )
            seen_parity.add(target_key)
            targets = [contract_id]
            verifier = weil_audits.verify_parity_audit
            expected_kind = "PARITY"
        else:
            target_key = _canonical_json_bytes(
                {
                    "lower_cell": audit.get("lower_cell"),
                    "higher_cell": audit.get("higher_cell"),
                }
            ).decode("utf-8")
            contract_pair = expected_nesting.get(target_key)
            if contract_pair is None or target_key in seen_nesting:
                raise PackagingError(
                    "nesting audits do not exactly cover adjacent exact degrees"
                )
            seen_nesting.add(target_key)
            targets = list(contract_pair)
            verifier = weil_audits.verify_degree_nesting_audit
            expected_kind = "DEGREE_NESTING"
        try:
            verification = verifier(audit)
        except weil_audits.WeilAuditVerificationError as exc:
            raise PackagingError(f"audit does not canonically replay: {audit_path}") from exc
        if (
            not isinstance(verification, Mapping)
            or verification.get("decision") != audit.get("decision")
            or verification.get("classification")
            != "REPRODUCED_EXPLORATORY_AUDIT"
            or verification.get("hypothesis_status") != HYPOTHESIS_STATUS
            or verification.get("audit_kind") != expected_kind
            or verification.get("decision") == "AUDIT_FAILED"
        ):
            raise PackagingError(f"audit verification failed: {audit_path}")
        audit_records.append(
            {
                "path": audit_path,
                "payload_sha256": audit_hash,
                "schema": expected_schema,
                "decision": str(verification.get("decision")),
                "target_contract_ids": targets,
            }
        )
    if seen_parity != set(expected_parity) or seen_nesting != set(expected_nesting):
        raise PackagingError("audit coverage does not match the frozen exact ladders")

    files = tuple(
        (
            relative,
            source / Path(*PurePosixPath(relative).parts),
        )
        for relative in sorted(expected_files)
    )
    return EvidenceSet(
        source=source,
        files=files,
        index=index,
        plan=plan,
        run=run,
        summary=summary,
        audit_records=tuple(audit_records),
        cell_count=len(cell_paths),
        confirmation_count=len(confirmation_paths),
        parity_audit_count=len(parity_paths),
        nesting_audit_count=len(nesting_paths),
    )


def _snapshot(files: Sequence[tuple[str, Path]]) -> tuple[FileRecord, ...]:
    records: list[FileRecord] = []
    for relative, path in files:
        try:
            size = path.stat().st_size
            digest = _file_sha256(path)
        except OSError as exc:
            raise PackagingError(f"cannot hash evidence file: {relative}") from exc
        records.append(FileRecord(relative, size, digest))
    return tuple(records)


def _stage_evidence(
    source: Path,
    destination: Path,
) -> tuple[Path, tuple[FileRecord, ...]]:
    """Copy one exact non-attempt byte snapshot into a private tree.

    The live evidence tree is used only as the source of this copy.  All
    semantic validation and archive construction happens from ``destination``.
    Hashes and sizes are accumulated from the same byte stream written to the
    stage, so a later change to the live tree cannot affect the package.
    """

    if not source.is_dir():
        raise PackagingError(f"evidence source is not a directory: {source}")
    actual_files, _, _ = _validate_source_shape(source)
    try:
        destination.mkdir(mode=0o700, parents=True, exist_ok=False)
        for directory in _ARCHIVE_DIRECTORIES[1:]:
            (destination / Path(*PurePosixPath(directory).parts)).mkdir(
                mode=0o700,
                parents=True,
                exist_ok=True,
            )
    except OSError as exc:
        raise PackagingError("cannot create private evidence staging tree") from exc

    records: list[FileRecord] = []
    for relative in sorted(actual_files):
        source_path = source / Path(*PurePosixPath(relative).parts)
        staged_path = destination / Path(*PurePosixPath(relative).parts)
        try:
            linklike = _is_reparse_or_link(source_path)
        except OSError as exc:
            raise PackagingError(
                f"cannot inspect evidence file before staging: {relative}"
            ) from exc
        if linklike or not source_path.is_file():
            raise PackagingError(f"evidence file changed before staging: {relative}")
        try:
            staged_path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
            digest = hashlib.sha256()
            size = 0
            with source_path.open("rb") as source_handle, staged_path.open(
                "xb"
            ) as staged_handle:
                for block in iter(lambda: source_handle.read(1024 * 1024), b""):
                    staged_handle.write(block)
                    digest.update(block)
                    size += len(block)
        except OSError as exc:
            raise PackagingError(f"cannot stage evidence file: {relative}") from exc
        records.append(FileRecord(relative, size, digest.hexdigest()))
    return destination, tuple(records)


def _content_set_sha256(records: Sequence[FileRecord]) -> str:
    return _content_sha256([record.to_record() for record in records])


def _directory_info(name: str) -> tarfile.TarInfo:
    archive_name = ARCHIVE_ROOT if not name else f"{ARCHIVE_ROOT}/{name}"
    info = tarfile.TarInfo(f"{archive_name}/")
    info.type = tarfile.DIRTYPE
    info.mode = 0o755
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    info.size = 0
    info.pax_headers = {}
    return info


def _file_info(record: FileRecord) -> tarfile.TarInfo:
    info = tarfile.TarInfo(f"{ARCHIVE_ROOT}/{record.path}")
    info.type = tarfile.REGTYPE
    info.mode = 0o644
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    info.size = record.size
    info.linkname = ""
    info.pax_headers = {}
    return info


def _write_canonical_ustar(
    output: BinaryIO,
    evidence: EvidenceSet,
    records: Sequence[FileRecord],
) -> None:
    """Write directories first, then lexicographically ordered files."""

    paths = {relative: path for relative, path in evidence.files}
    with tarfile.open(
        fileobj=output,
        mode="w",
        format=tarfile.USTAR_FORMAT,
        encoding="utf-8",
        errors="strict",
    ) as archive:
        for directory in _ARCHIVE_DIRECTORIES:
            archive.addfile(_directory_info(directory))
        for record in records:
            with paths[record.path].open("rb") as source_handle:
                archive.addfile(_file_info(record), source_handle)


def _digest_open_file(handle: BinaryIO) -> tuple[int, str]:
    """Hash a seekable file from byte zero and return its size and digest."""

    handle.seek(0)
    digest = hashlib.sha256()
    size = 0
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
        size += len(block)
    return size, digest.hexdigest()


def _canonical_archive_sizes(
    records: Sequence[FileRecord],
) -> tuple[int, int]:
    """Return canonical member-end and total byte lengths for this USTAR."""

    member_end = len(_ARCHIVE_DIRECTORIES) * tarfile.BLOCKSIZE
    for record in records:
        padded_size = (
            (record.size + tarfile.BLOCKSIZE - 1) // tarfile.BLOCKSIZE
        ) * tarfile.BLOCKSIZE
        member_end += tarfile.BLOCKSIZE + padded_size
    minimum_total = member_end + 2 * tarfile.BLOCKSIZE
    total_size = (
        (minimum_total + tarfile.RECORDSIZE - 1) // tarfile.RECORDSIZE
    ) * tarfile.RECORDSIZE
    return member_end, total_size


def _verify_open_archive(
    raw_archive: BinaryIO,
    records: Sequence[FileRecord],
) -> tuple[int, str]:
    """Verify and fingerprint an already-open canonical archive."""

    expected_names = [
        ARCHIVE_ROOT if not name else f"{ARCHIVE_ROOT}/{name}"
        for name in _ARCHIVE_DIRECTORIES
    ] + [f"{ARCHIVE_ROOT}/{record.path}" for record in records]
    by_name = {
        f"{ARCHIVE_ROOT}/{record.path}": record for record in records
    }
    member_end, expected_size = _canonical_archive_sizes(records)
    initial_size, initial_sha256 = _digest_open_file(raw_archive)
    if initial_size != expected_size:
        raise PackagingError("canonical USTAR archive has an unexpected byte length")
    raw_archive.seek(0)
    with tarfile.open(fileobj=raw_archive, mode="r:") as archive:
        members = archive.getmembers()
        if [member.name for member in members] != expected_names:
            raise PackagingError("archive member paths or order changed")
        expected_offset = 0
        for ordinal, member in enumerate(members):
            expected_directory = ordinal < len(_ARCHIVE_DIRECTORIES)
            if (
                member.offset != expected_offset
                or member.offset_data != member.offset + tarfile.BLOCKSIZE
                or member.pax_headers
                or member.linkname != ""
                or member.devmajor != 0
                or member.devminor != 0
                or member.mtime != 0
                or member.uid != 0
                or member.gid != 0
                or member.uname != ""
                or member.gname != ""
                or member.mode != (0o755 if expected_directory else 0o644)
                or (expected_directory and not member.isdir())
                or (not expected_directory and not member.isfile())
            ):
                raise PackagingError(
                    f"archive member metadata changed: {member.name}"
                )
            raw_archive.seek(member.offset)
            header = raw_archive.read(tarfile.BLOCKSIZE)
            if (
                len(header) != tarfile.BLOCKSIZE
                or header[257:263] != b"ustar\0"
                or header[263:265] != b"00"
            ):
                raise PackagingError(
                    f"archive member is not canonical USTAR: {member.name}"
                )
            expected_offset += tarfile.BLOCKSIZE
            if expected_directory:
                if member.size != 0:
                    raise PackagingError(
                        f"archive directory has data: {member.name}"
                    )
                continue
            record = by_name[member.name]
            if member.size != record.size:
                raise PackagingError(
                    f"archive member size changed: {member.name}"
                )
            extracted = archive.extractfile(member)
            if extracted is None:
                raise PackagingError(
                    f"archive member cannot be read: {member.name}"
                )
            digest = hashlib.sha256()
            total = 0
            for block in iter(lambda: extracted.read(1024 * 1024), b""):
                total += len(block)
                digest.update(block)
            if total != record.size or digest.hexdigest() != record.sha256:
                raise PackagingError(
                    f"archive bytes differ from staged snapshot: {member.name}"
                )
            padding_size = (-record.size) % tarfile.BLOCKSIZE
            raw_archive.seek(member.offset_data + record.size)
            padding = raw_archive.read(padding_size)
            if len(padding) != padding_size or any(padding):
                raise PackagingError(
                    f"archive member padding is noncanonical: {member.name}"
                )
            expected_offset += record.size + padding_size
        if expected_offset != member_end:
            raise PackagingError("canonical USTAR member layout changed")
        raw_archive.seek(member_end)
        trailer = raw_archive.read(expected_size - member_end)
        if len(trailer) != expected_size - member_end or any(trailer):
            raise PackagingError("canonical USTAR trailer padding is not all zero")
    final_size, final_sha256 = _digest_open_file(raw_archive)
    if (final_size, final_sha256) != (initial_size, initial_sha256):
        raise PackagingError("archive changed during final verification")
    return final_size, final_sha256


def _verify_written_archive(
    archive_path: Path,
    records: Sequence[FileRecord],
    *,
    open_handle: BinaryIO | None = None,
) -> tuple[int, str]:
    """Verify/fingerprint finalized bytes, optionally through a held handle."""

    try:
        if open_handle is not None:
            return _verify_open_archive(open_handle, records)
        with archive_path.open("rb") as raw_archive:
            return _verify_open_archive(raw_archive, records)
    except PackagingError:
        raise
    except (OSError, tarfile.TarError, EOFError) as exc:
        raise PackagingError("cannot reopen canonical USTAR archive") from exc


def _require_open_path_identity(path: Path, handle: BinaryIO) -> None:
    """Reject replacement of the published path while its handle is held."""

    try:
        if _is_reparse_or_link(path):
            raise PackagingError("published archive became a link or reparse point")
        path_metadata = path.lstat()
        handle_metadata = os.fstat(handle.fileno())
    except OSError as exc:
        raise PackagingError("cannot inspect published archive identity") from exc
    if (
        not stat.S_ISREG(path_metadata.st_mode)
        or not stat.S_ISREG(handle_metadata.st_mode)
        or not os.path.samestat(path_metadata, handle_metadata)
    ):
        raise PackagingError("published archive path was replaced")


@contextmanager
def _locked_archive(path: Path) -> Iterator[BinaryIO]:
    """Hold an exclusive best-effort OS lock during finalization."""

    handle = path.open("rb")
    locked = False
    lock_size = max(os.fstat(handle.fileno()).st_size, 1)
    try:
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, lock_size)
        else:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        locked = True
    except OSError as exc:
        handle.close()
        raise PackagingError("cannot exclusively lock published archive") from exc
    try:
        yield handle
    finally:
        if locked:
            try:
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt

                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, lock_size)
                else:
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
        handle.close()


def _manifest_body(
    *,
    evidence: EvidenceSet,
    records: Sequence[FileRecord],
    archive_path: Path,
    archive_sha256: str,
    archive_size: int,
    repository: str,
    release_tag: str,
    producer_commit: str,
) -> dict[str, Any]:
    quoted_tag = quote(release_tag, safe="")
    quoted_asset = quote(archive_path.name, safe="")
    release_url = f"https://github.com/{repository}/releases/tag/{quoted_tag}"
    return {
        "schema": MANIFEST_SCHEMA,
        "classification": CLASSIFICATION,
        "hypothesis_status": HYPOTHESIS_STATUS,
        "producer_git_commit": producer_commit,
        "release": {
            "asset_url": (
                f"https://github.com/{repository}/releases/download/"
                f"{quoted_tag}/{quoted_asset}"
            ),
            "repository": repository,
            "tag": release_tag,
            "url": release_url,
        },
        "archive": {
            "builder_id": BUILDER_ID,
            "byte_size": str(archive_size),
            "checksum_sidecar": f"{archive_path.name}.sha256",
            "content_set_algorithm": (
                "sha256(canonical-json(sorted [{path,sha256,size_bytes}]))"
            ),
            "content_set_sha256": _content_set_sha256(records),
            "excluded_paths": ["attempts/"],
            "file_count": str(len(records)),
            "format": "ustar-uncompressed",
            "member_file_bytes": str(sum(record.size for record in records)),
            "metadata_profile": {
                "directory_mode": "0755",
                "file_mode": "0644",
                "gid": "0",
                "gname": "",
                "mtime": "0",
                "ordering": "fixed-directories-then-lexicographic-posix-files",
                "uid": "0",
                "uname": "",
            },
            "name": archive_path.name,
            "root": f"{ARCHIVE_ROOT}/",
            "sha256": archive_sha256,
        },
        "evidence": {
            "audits": list(evidence.audit_records),
            "counts": {
                "cells": str(evidence.cell_count),
                "confirmations": str(evidence.confirmation_count),
                "nesting_audits": str(evidence.nesting_audit_count),
                "parity_audits": str(evidence.parity_audit_count),
            },
            "index": {
                "path": "index.json",
                "payload_sha256": evidence.index["payload_sha256"],
            },
            "plan": {
                "path": "plan.json",
                "plan_sha256": evidence.plan["plan_sha256"],
            },
            "run": {
                "path": "run.json",
                "payload_sha256": evidence.run["payload_sha256"],
            },
            "transition_summary": {
                "conclusion": evidence.summary["conclusion"],
                "numerical_replay_performed": False,
                "payload_sha256": evidence.summary["payload_sha256"],
                "schema": transition_summary.SUMMARY_SCHEMA,
            },
        },
        "manifest_storage": "tracked adjacent to evidence tree; not in archive",
        "limitation": (
            "This release contains exploratory finite-dimensional evidence. "
            "It does not prove or disprove the Riemann Hypothesis."
        ),
    }


def _write_json_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    created = False
    try:
        with path.open("x", encoding="utf-8", newline="\n") as handle:
            created = True
            json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
    except BaseException:
        if created:
            try:
                path.unlink()
            except OSError:
                pass
        raise


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def package_release(
    source: Path,
    output: Path,
    *,
    repository: str,
    release_tag: str,
    producer_commit: str,
) -> dict[str, Any]:
    """Create the canonical TAR, checksum sidecar, and adjacent manifest."""

    if not _REPOSITORY_RE.fullmatch(repository):
        raise PackagingError("repository must have the form owner/name")
    if not release_tag or any(character.isspace() for character in release_tag):
        raise PackagingError("release tag must be nonempty and contain no whitespace")
    if not _COMMIT_RE.fullmatch(producer_commit):
        raise PackagingError("producer commit must be a full 40-hex Git commit")
    producer_commit = producer_commit.lower()

    source_input = Path(source)
    try:
        source_is_linklike = _is_reparse_or_link(source_input)
    except FileNotFoundError:
        source_is_linklike = False
    except OSError as exc:
        raise PackagingError("cannot inspect evidence source") from exc
    if source_is_linklike:
        raise PackagingError(
            "evidence source cannot be a symbolic link, junction, or reparse point"
        )
    source_path = source_input.resolve()
    output_input = Path(output)
    if output_input.is_symlink():
        raise PackagingError("archive output cannot be a symbolic link")
    output_path = output_input.resolve()
    if output_path.suffix.lower() != ".tar":
        raise PackagingError("canonical archive output must end in .tar")
    if _is_within(output_path, source_path):
        raise PackagingError("archive output must be outside the evidence tree")
    checksum_path = output_path.with_name(f"{output_path.name}.sha256")
    manifest_path = source_path.with_name(f"{source_path.name}.release.json")
    destinations = (output_path, checksum_path, manifest_path)
    for destination in destinations:
        if destination.exists() or destination.is_symlink():
            raise PackagingError(f"refusing to overwrite: {destination}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []
    temporary_archive: Path | None = None
    try:
        with tempfile.TemporaryDirectory(
            prefix="rh-lab-weil-transition-stage-"
        ) as staging_directory:
            staged_source, copied_records = _stage_evidence(
                source_path,
                Path(staging_directory) / "evidence",
            )
            evidence = _collect_evidence(staged_source)
            records = _snapshot(evidence.files)
            if records != copied_records:
                raise PackagingError(
                    "private evidence snapshot changed during validation"
                )

            try:
                descriptor, temporary_name = tempfile.mkstemp(
                    prefix=f".{output_path.name}.",
                    suffix=".tmp",
                    dir=output_path.parent,
                )
                temporary_archive = Path(temporary_name)
                with os.fdopen(descriptor, "wb") as output_handle:
                    _write_canonical_ustar(output_handle, evidence, records)
                try:
                    os.link(temporary_archive, output_path)
                except FileExistsError as exc:
                    raise PackagingError(
                        f"refusing to overwrite: {output_path}"
                    ) from exc
                except OSError as exc:
                    raise PackagingError(
                        "cannot atomically publish canonical USTAR archive"
                    ) from exc
                created.append(output_path)
                temporary_archive.unlink()
                temporary_archive = None
            except PackagingError:
                raise
            except (OSError, tarfile.TarError, ValueError) as exc:
                raise PackagingError("cannot write canonical USTAR archive") from exc

            with _locked_archive(output_path) as archive_handle:
                _require_open_path_identity(output_path, archive_handle)
                archive_size, archive_sha256 = _verify_written_archive(
                    output_path,
                    records,
                    open_handle=archive_handle,
                )
                manifest_body = _manifest_body(
                    evidence=evidence,
                    records=records,
                    archive_path=output_path,
                    archive_sha256=archive_sha256,
                    archive_size=archive_size,
                    repository=repository,
                    release_tag=release_tag,
                    producer_commit=producer_commit,
                )
                manifest = {
                    **manifest_body,
                    "payload_sha256": _content_sha256(manifest_body),
                }

                with checksum_path.open(
                    "x", encoding="ascii", newline="\n"
                ) as handle:
                    created.append(checksum_path)
                    handle.write(f"{archive_sha256}  {output_path.name}\n")
                _write_json_exclusive(manifest_path, manifest)
                created.append(manifest_path)

                _require_open_path_identity(output_path, archive_handle)
                final_size, final_sha256 = _verify_open_archive(
                    archive_handle,
                    records,
                )
                _require_open_path_identity(output_path, archive_handle)
                if (final_size, final_sha256) != (
                    archive_size,
                    archive_sha256,
                ):
                    raise PackagingError(
                        "published archive changed during metadata creation"
                    )
    except BaseException:
        for path in reversed(created):
            try:
                path.unlink()
            except OSError:
                pass
        raise
    finally:
        if temporary_archive is not None:
            try:
                temporary_archive.unlink()
            except OSError:
                pass

    return {
        "archive": str(output_path),
        "archive_sha256": archive_sha256,
        "checksum": str(checksum_path),
        "content_set_sha256": manifest["archive"]["content_set_sha256"],
        "manifest": str(manifest_path),
        "manifest_payload_sha256": manifest["payload_sha256"],
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Package a complete v2 Weil transition evidence tree."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repository", required=True)
    parser.add_argument("--release-tag", required=True)
    parser.add_argument("--producer-commit", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = package_release(
            args.source,
            args.output,
            repository=args.repository,
            release_tag=args.release_tag,
            producer_commit=args.producer_commit,
        )
    except (PackagingError, OSError) as exc:
        print(f"release packaging rejected: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
