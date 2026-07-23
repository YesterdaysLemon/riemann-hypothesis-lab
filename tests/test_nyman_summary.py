from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil
from typing import Any

import pytest

from riemann_lab.artifacts import content_sha256, write_json
import riemann_lab.nyman_search as search
import riemann_lab.nyman_summary as summary


REAL_CHECKPOINT = Path("results/nyman-natural-v1")


def _copy_checkpoint(tmp_path: Path) -> Path:
    target = tmp_path / "checkpoint"
    shutil.copytree(REAL_CHECKPOINT, target)
    return target


def _rehash(record: dict[str, Any], field: str = "payload_sha256") -> None:
    record[field] = content_sha256(
        {key: value for key, value in record.items() if key != field}
    )


def test_real_summary_is_compact_deterministic_and_has_no_numerical_replay(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def forbidden(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("numerical kernel or certificate replay was attempted")

    monkeypatch.setattr(search.core, "build_natural_system", forbidden)
    monkeypatch.setattr(search.core, "certify_dyadic_upper_bound", forbidden)
    monkeypatch.setattr(search.core, "certify_augmented_lower_bound", forbidden)
    monkeypatch.setattr(search, "_derive_cell", forbidden)
    monkeypatch.setattr(search, "verify_nyman_search", forbidden)

    artifact = summary.generate_nyman_summary(REAL_CHECKPOINT)
    assert artifact == summary.generate_nyman_summary(REAL_CHECKPOINT)
    assert artifact["schema"] == summary.SUMMARY_SCHEMA
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["conclusion"] == "ALL_FINITE_DISTANCE_BRACKETS_CERTIFIED"
    assert artifact["payload_sha256"] == content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )
    assert artifact["source"] == {
        "plan": {
            "path": "plan.json",
            "plan_sha256": summary.FROZEN_PLAN_SHA256,
        },
        "run": {
            "path": "run.json",
            "payload_sha256": (
                "58d51883701e8cd41f90592b3cae83c3aa71ecba791658d6a1b57a35bf9bfc15"
            ),
        },
        "index": {
            "path": "index.json",
            "payload_sha256": (
                "281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23"
            ),
        },
    }

    assert artifact["counts"] == {
        "cells": "6",
        "by_status": {
            "FINITE_DISTANCE_BRACKET_CERTIFIED": "6",
            "INCONCLUSIVE": "0",
        },
        "stored_exact_dyadic_coefficients": "504",
        "stored_ldlt_pivots": "510",
        "strict_improvement_diagnostics": "5",
        "strict_improvement_by_decision": {
            "STRICT_IMPROVEMENT_CERTIFIED": "5"
        },
    }
    assert [cell["n"] for cell in artifact["cells"]] == [
        "8",
        "16",
        "32",
        "64",
        "128",
        "256",
    ]
    assert all(
        cell["bounds"]["width"]
        == {
            "dyadic": {"numerator": "1", "denominator_exponent": "120"},
            "exact_fraction": {
                "numerator": "1",
                "denominator": "1329227995784915872903807060280344576",
            },
        }
        for cell in artifact["cells"]
    )
    assert artifact["cells"][0]["bounds"]["lower"] == {
        "dyadic": {
            "numerator": "8221705725423578681524695917750771394",
            "denominator_exponent": "128",
        },
        "exact_fraction": {
            "numerator": "4110852862711789340762347958875385697",
            "denominator": "170141183460469231731687303715884105728",
        },
    }
    assert [cell["stored_energy_display"] for cell in artifact["cells"]] == [
        "[0.02416142158589668502280894301085133744677 +/- 4.52e-42]",
        "[0.01789402347696943509905506314795876127745 +/- 4.65e-42]",
        "[0.01405194369952985983642504072212329160349 +/- 1.26e-42]",
        "[0.01137604029967365814705897002932676667883 +/- 4.18e-42]",
        "[0.009658549278111909910035392571452577088501 +/- 3.89e-43]",
        "[0.008233716277261021441598038102635434534105 +/- 3.78e-43]",
    ]
    assert artifact["generation_kernel"] == {
        "path": "kernels/generation.json",
        "precision_bits": "768",
        "max_n": "256",
        "payload_sha256": (
            "6fb54e7f63208f1bd81fb7ef2e592a8a655458dd7744081748590159b868cc26"
        ),
        "kernel_contract_sha256": (
            "bc3bc8679379a7a2c3a7c0f530829b67b7170b6ada950c893f6010a32eba06cf"
        ),
        "core_system_content_sha256": (
            "1587f9ab711c2b787f0072dcd82476cef3ee1ca0198ca167eccb0c86c95b8726"
        ),
        "prefixes": [
            {
                "n": prefix["n"],
                "core_system_content_sha256": prefix[
                    "core_system_content_sha256"
                ],
                "prefix_kernel_sha256": prefix["prefix_kernel_sha256"],
            }
            for prefix in json.loads(
                (REAL_CHECKPOINT / "kernels" / "generation.json").read_text(
                    encoding="utf-8"
                )
            )["prefixes"]
        ],
    }
    assert artifact["strict_improvement"]["role"] == (
        "separate exact diagnostic; never a cell-validity gate"
    )
    assert all(
        item["decision"] == "STRICT_IMPROVEMENT_CERTIFIED"
        and item["cell_validity_gate"] is False
        for item in artifact["strict_improvement"]["diagnostics"]
    )
    assert artifact["summary_method"] == {
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
        "stored_numerical_values_treated_as_hash_bound_records": True,
        "content_hashes_checked": [
            "plan.plan_sha256",
            "run.payload_sha256",
            "index.payload_sha256",
            "generation_kernel.payload_sha256",
            "cell.payload_sha256",
            "candidate.candidate_sha256",
            "upper_certificate.payload_sha256",
            "lower_certificate.payload_sha256",
        ],
    }
    assert "does not prove" in artifact["limitation"]


def test_verifier_accepts_mapping_and_path_and_exactly_regenerates(
    tmp_path: Path,
) -> None:
    artifact = summary.generate_nyman_summary(REAL_CHECKPOINT)
    summary_path = tmp_path / "summary.json"
    write_json(summary_path, artifact)

    expected = {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": artifact["payload_sha256"],
        "plan_sha256": summary.FROZEN_PLAN_SHA256,
        "index_payload_sha256": artifact["source"]["index"]["payload_sha256"],
        "verified_cells": "6",
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }
    assert summary.verify_nyman_summary(artifact, REAL_CHECKPOINT) == expected
    assert summary.verify_nyman_summary(summary_path, REAL_CHECKPOINT) == expected


def test_summary_verifier_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    canonical = json.dumps(
        summary.generate_nyman_summary(REAL_CHECKPOINT),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    path = tmp_path / "duplicate-summary-key.json"
    path.write_text('{"schema":"ignored",' + canonical[1:], encoding="utf-8")

    with pytest.raises(
        summary.NymanSummaryVerificationError,
        match="duplicate JSON object key: schema",
    ):
        summary.verify_nyman_summary(path, REAL_CHECKPOINT)


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_summary_verifier_rejects_nonstandard_json_constants(
    constant: str,
    tmp_path: Path,
) -> None:
    canonical = json.dumps(
        summary.generate_nyman_summary(REAL_CHECKPOINT),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    path = tmp_path / "nonstandard-summary-constant.json"
    path.write_text('{"extra":' + constant + "," + canonical[1:], encoding="utf-8")

    with pytest.raises(
        summary.NymanSummaryVerificationError,
        match=f"nonstandard JSON constant: {constant}",
    ):
        summary.verify_nyman_summary(path, REAL_CHECKPOINT)


def test_summary_generator_rejects_duplicate_checkpoint_keys(
    tmp_path: Path,
) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    path = checkpoint / "index.json"
    canonical = path.read_text(encoding="utf-8")
    path.write_text('{"schema":"ignored",' + canonical.lstrip()[1:], encoding="utf-8")

    with pytest.raises(
        summary.NymanSummaryError,
        match="duplicate JSON object key: schema",
    ):
        summary.generate_nyman_summary(checkpoint)


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_summary_generator_rejects_nonstandard_checkpoint_constants(
    constant: str,
    tmp_path: Path,
) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    path = checkpoint / "index.json"
    canonical = path.read_text(encoding="utf-8")
    path.write_text(
        '{"extra":' + constant + "," + canonical.lstrip()[1:],
        encoding="utf-8",
    )

    with pytest.raises(
        summary.NymanSummaryError,
        match=f"nonstandard JSON constant: {constant}",
    ):
        summary.generate_nyman_summary(checkpoint)


def test_verifier_rejects_rehashed_summary_tampering() -> None:
    artifact = summary.generate_nyman_summary(REAL_CHECKPOINT)
    tampered = copy.deepcopy(artifact)
    tampered["counts"]["cells"] = "5"
    _rehash(tampered)

    with pytest.raises(
        summary.NymanSummaryVerificationError,
        match="does not canonically regenerate",
    ):
        summary.verify_nyman_summary(tampered, REAL_CHECKPOINT)


def test_generator_rejects_a_cell_payload_mutation(tmp_path: Path) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    path = checkpoint / "cells" / "n-0008.json"
    artifact = json.loads(path.read_text(encoding="utf-8"))
    artifact["generation"]["upper_certificate"]["energy"]["display"] = "changed"
    write_json(path, artifact)

    with pytest.raises(summary.NymanSummaryError, match="cell payload hash mismatch"):
        summary.generate_nyman_summary(checkpoint)


def test_generator_rejects_a_rehashed_unknown_certificate_field(
    tmp_path: Path,
) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    cell_path = checkpoint / "cells" / "n-0008.json"
    artifact = json.loads(cell_path.read_text(encoding="utf-8"))
    upper = artifact["generation"]["upper_certificate"]
    upper["unknown_field"] = "must fail closed"
    _rehash(upper)
    _rehash(artifact)
    write_json(cell_path, artifact)

    index_path = checkpoint / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index["cells"][0]["payload_sha256"] = artifact["payload_sha256"]
    _rehash(index)
    write_json(index_path, index)

    with pytest.raises(
        summary.NymanSummaryError, match="upper certificate has noncanonical fields"
    ):
        summary.generate_nyman_summary(checkpoint)


def test_generator_rejects_rehashed_kernel_unknown_fields(tmp_path: Path) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    kernel_path = checkpoint / "kernels" / "generation.json"
    kernel = json.loads(kernel_path.read_text(encoding="utf-8"))
    kernel["unknown_field"] = "must fail closed"
    _rehash(kernel)
    write_json(kernel_path, kernel)

    index_path = checkpoint / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index["generation_kernel"]["payload_sha256"] = kernel["payload_sha256"]
    _rehash(index)
    write_json(index_path, index)

    with pytest.raises(
        summary.NymanSummaryError, match="generation kernel has noncanonical fields"
    ):
        summary.generate_nyman_summary(checkpoint)


@pytest.mark.parametrize(
    "path",
    ("../outside.json", "cells/../../outside.json", "C:/outside.json"),
)
def test_generator_rejects_rehashed_cell_path_escapes(
    path: str, tmp_path: Path
) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    index_path = checkpoint / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index["cells"][0]["path"] = path
    _rehash(index)
    write_json(index_path, index)

    with pytest.raises(summary.NymanSummaryError, match="cell path changed"):
        summary.generate_nyman_summary(checkpoint)


def test_generator_rejects_rehashed_index_semantic_mutation(tmp_path: Path) -> None:
    checkpoint = _copy_checkpoint(tmp_path)
    index_path = checkpoint / "index.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    index["strict_improvement_diagnostics"][0]["decision"] = (
        "STRICT_IMPROVEMENT_NOT_CERTIFIED"
    )
    _rehash(index)
    write_json(index_path, index)

    with pytest.raises(summary.NymanSummaryError, match="index is not canonical"):
        summary.generate_nyman_summary(checkpoint)
