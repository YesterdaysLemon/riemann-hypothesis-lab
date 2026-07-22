"""Schema and evidence-hash validation for the public claim ledger."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .artifacts import file_sha256


ALLOWED_CLASSIFICATIONS = {
    "PROVED",
    "DISPROVED",
    "CERTIFIED_FINITE",
    "CONDITIONAL",
    "EXPLORATORY",
    "FAILED",
}
RESOLUTION_STATUSES = {"PROVED", "DISPROVED"}


class ClaimValidationError(ValueError):
    """Raised when the ledger could permit an unsupported headline claim."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ClaimValidationError(message)


def _validate_file_reference(
    repository_root: Path,
    record: dict[str, Any],
    path_key: str,
    hash_key: str,
    label: str,
) -> Path:
    artifact_name = record.get(path_key)
    expected_hash = record.get(hash_key)
    _require(
        isinstance(artifact_name, str) and artifact_name,
        f"{label} path is missing",
    )
    _require(
        isinstance(expected_hash, str)
        and len(expected_hash) == 64
        and all(character in "0123456789abcdef" for character in expected_hash.lower()),
        f"{label} hash is invalid",
    )
    artifact_path = (repository_root / artifact_name).resolve()
    _require(artifact_path.is_relative_to(repository_root), f"{label} escapes repository")
    _require(artifact_path.is_file(), f"{label} is missing")
    _require(file_sha256(artifact_path) == expected_hash.lower(), f"{label} hash mismatch")
    return artifact_path


def load_and_validate_registry(path: Path) -> dict[str, Any]:
    """Load a registry, validate evidence hashes, and enforce promotion gates."""

    registry = json.loads(path.read_text(encoding="utf-8"))
    _require(registry.get("schema") == "rh-lab/claim-registry/v1", "bad schema")

    repository_root = path.parent.parent.resolve()
    hypothesis = registry.get("hypothesis", {})
    status = hypothesis.get("status")
    _require(
        status in {"UNRESOLVED", *RESOLUTION_STATUSES},
        "hypothesis status is invalid",
    )

    resolution_evidence: dict[str, Any] | None = None
    if status in RESOLUTION_STATUSES:
        evidence = hypothesis.get("resolution_evidence")
        _require(isinstance(evidence, dict), "resolved status needs evidence")
        resolution_evidence = evidence
        _require(
            isinstance(evidence.get("claim_id"), str) and evidence["claim_id"],
            "resolved status needs a linked claim id",
        )
        checker = evidence.get("formal_checker")
        _require(
            isinstance(checker, dict)
            and all(checker.get(key) for key in ("name", "version", "command")),
            "resolved status needs a pinned formal checker manifest",
        )
        _validate_file_reference(
            repository_root,
            evidence,
            "formal_proof_artifact",
            "formal_proof_sha256",
            "formal proof artifact",
        )
        _validate_file_reference(
            repository_root,
            evidence,
            "assumption_manifest",
            "assumption_manifest_sha256",
            "assumption manifest",
        )
        reproductions = evidence.get("independent_reproductions", [])
        _require(
            isinstance(reproductions, list) and len(reproductions) >= 2,
            "resolved status needs at least two independent reproductions",
        )
        reproduction_paths: set[Path] = set()
        backend_families: set[str] = set()
        for index, reproduction in enumerate(reproductions, start=1):
            _require(
                isinstance(reproduction, dict),
                f"independent reproduction {index} must be an object",
            )
            reproduction_path = _validate_file_reference(
                repository_root,
                reproduction,
                "artifact",
                "artifact_sha256",
                f"independent reproduction {index}",
            )
            backend_family = reproduction.get("backend_family")
            _require(
                isinstance(backend_family, str) and backend_family,
                f"independent reproduction {index} needs a backend family",
            )
            _require(
                bool(reproduction.get("command")),
                f"independent reproduction {index} needs a command",
            )
            _require(
                reproduction_path not in reproduction_paths,
                "independent reproductions must use distinct artifacts",
            )
            reproduction_paths.add(reproduction_path)
            backend_families.add(backend_family)
        _require(
            len(backend_families) >= 2,
            "independent reproductions must use at least two backend families",
        )
    else:
        _require(
            hypothesis.get("resolution_evidence") is None,
            "unresolved status cannot carry resolution evidence",
        )

    claims = registry.get("claims")
    _require(isinstance(claims, list), "claims must be a list")
    seen_ids: set[str] = set()
    for claim in claims:
        claim_id = claim.get("id")
        _require(isinstance(claim_id, str) and claim_id, "claim id is required")
        _require(claim_id not in seen_ids, f"duplicate claim id: {claim_id}")
        seen_ids.add(claim_id)

        classification = claim.get("classification")
        _require(
            classification in ALLOWED_CLASSIFICATIONS,
            f"invalid classification for {claim_id}",
        )
        _require(bool(claim.get("statement")), f"statement missing for {claim_id}")
        _require(bool(claim.get("limitations")), f"limitations missing for {claim_id}")

        if classification == "CERTIFIED_FINITE":
            _require(bool(claim.get("scope")), f"finite scope missing for {claim_id}")
            _require(
                claim.get("implies_rh") is False,
                f"finite claim {claim_id} must explicitly deny implying RH",
            )

        if classification in RESOLUTION_STATUSES:
            _require(
                status == classification,
                f"resolution claim {claim_id} contradicts hypothesis status",
            )
            _require(
                resolution_evidence is not None
                and resolution_evidence.get("claim_id") == claim_id,
                f"resolution claim {claim_id} is not linked from hypothesis evidence",
            )

        artifact_name = claim.get("artifact")
        expected_hash = claim.get("artifact_sha256")
        if artifact_name is not None or expected_hash is not None:
            _validate_file_reference(
                repository_root,
                claim,
                "artifact",
                "artifact_sha256",
                f"artifact for {claim_id}",
            )

        if classification in RESOLUTION_STATUSES:
            _require(
                artifact_name is not None and expected_hash is not None,
                f"resolution claim {claim_id} needs a hash-bound artifact",
            )

    if resolution_evidence is not None:
        _require(
            resolution_evidence["claim_id"] in seen_ids,
            "hypothesis resolution links to a missing claim",
        )

    return registry
