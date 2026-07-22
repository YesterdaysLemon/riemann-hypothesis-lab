from __future__ import annotations

import json
from pathlib import Path

import pytest

from riemann_lab.artifacts import content_sha256
from riemann_lab.claims import ClaimValidationError, load_and_validate_registry


def test_repository_claim_ledger_is_valid() -> None:
    registry = load_and_validate_registry(Path("claims/registry.json"))
    assert registry["hypothesis"]["status"] == "UNRESOLVED"


def test_resolution_requires_formal_and_independent_evidence(tmp_path: Path) -> None:
    claims_dir = tmp_path / "claims"
    claims_dir.mkdir()
    registry_path = claims_dir / "registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema": "rh-lab/claim-registry/v1",
                "hypothesis": {
                    "status": "PROVED",
                    "resolution_evidence": None,
                },
                "claims": [],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ClaimValidationError, match="needs evidence"):
        load_and_validate_registry(registry_path)


def test_dummy_resolution_evidence_is_rejected(tmp_path: Path) -> None:
    claims_dir = tmp_path / "claims"
    claims_dir.mkdir()
    registry_path = claims_dir / "registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema": "rh-lab/claim-registry/v1",
                "hypothesis": {
                    "status": "PROVED",
                    "resolution_evidence": {
                        "claim_id": "rh-proof",
                        "formal_checker": "x",
                        "assumption_manifest": "y",
                        "independent_reproductions": ["same", "same"],
                    },
                },
                "claims": [],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ClaimValidationError, match="pinned formal checker"):
        load_and_validate_registry(registry_path)


def test_unresolved_registry_cannot_contain_resolution_claim(tmp_path: Path) -> None:
    claims_dir = tmp_path / "claims"
    claims_dir.mkdir()
    registry_path = claims_dir / "registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema": "rh-lab/claim-registry/v1",
                "hypothesis": {
                    "status": "UNRESOLVED",
                    "resolution_evidence": None,
                },
                "claims": [
                    {
                        "id": "unsupported-proof",
                        "classification": "PROVED",
                        "statement": "RH proved",
                        "limitations": "none",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ClaimValidationError, match="contradicts hypothesis status"):
        load_and_validate_registry(registry_path)


def test_canonical_json_binding_ignores_formatting_and_line_endings(
    tmp_path: Path,
) -> None:
    claims_dir = tmp_path / "claims"
    results_dir = tmp_path / "results"
    claims_dir.mkdir()
    results_dir.mkdir()
    artifact = {"classification": "CERTIFIED_FINITE", "values": [1, 2, 3]}
    artifact_path = results_dir / "artifact.json"
    artifact_path.write_bytes(b'{\r\n  "values": [1, 2, 3],\r\n  "classification": "CERTIFIED_FINITE"\r\n}\r\n')
    registry_path = claims_dir / "registry.json"
    registry_path.write_text(
        json.dumps(
            {
                "schema": "rh-lab/claim-registry/v1",
                "hypothesis": {
                    "status": "UNRESOLVED",
                    "resolution_evidence": None,
                },
                "claims": [
                    {
                        "id": "portable-json",
                        "classification": "CERTIFIED_FINITE",
                        "statement": "finite test",
                        "limitations": "finite only",
                        "scope": {"n": "3"},
                        "implies_rh": False,
                        "artifact": "results/artifact.json",
                        "artifact_hash_mode": "canonical-json-v1",
                        "artifact_sha256": content_sha256(artifact),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    load_and_validate_registry(registry_path)
    artifact_path.write_text(json.dumps(artifact, separators=(",", ":")), encoding="utf-8")
    load_and_validate_registry(registry_path)
