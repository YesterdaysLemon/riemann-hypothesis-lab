from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from riemann_lab.cli import _integer_witness, main


def test_integer_witness_parser_is_exact_and_rejects_empty_or_zero() -> None:
    assert _integer_witness("1, -2, 3") == (1, -2, 3)
    for malformed in ("", "1,,2", "1,x,2", "0,0"):
        with pytest.raises(argparse.ArgumentTypeError):
            _integer_witness(malformed)


def test_transition_cli_writes_and_replays_zero_cell_checkpoint(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    checkpoint = tmp_path / "transition"
    exit_code = main(
        [
            "weil-transition-search",
            "--plan",
            "plans/weil-transition-q7-q9-v2.json",
            "--checkpoint-dir",
            str(checkpoint),
            "--max-cells",
            "0",
        ]
    )
    generated = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert generated["classification"] == "EXPLORATORY"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["progress"]["planned_cells"] == "81"
    assert generated["progress"]["completed_cells"] == "0"

    exit_code = main(
        [
            "verify-weil-transition-search",
            "--index",
            str(checkpoint / "index.json"),
            "--checkpoint-dir",
            str(checkpoint),
        ]
    )
    replay = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert replay["classification"] == (
        "REPRODUCED_EXPLORATORY_TRANSITION_SEARCH"
    )
    assert replay["verified_cells"] == "0"
    assert replay["hypothesis_status"] == "UNRESOLVED"


def test_nyman_cli_writes_and_replays_zero_cell_checkpoint(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    checkpoint = tmp_path / "nyman"
    exit_code = main(
        [
            "nyman-search",
            "--plan",
            "plans/nyman-natural-v1.json",
            "--checkpoint-dir",
            str(checkpoint),
            "--max-cells",
            "0",
        ]
    )
    generated = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert generated["classification"] == "EXPLORATORY"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["progress"]["planned_cells"] == "6"
    assert generated["progress"]["completed_cells"] == "0"

    exit_code = main(
        [
            "verify-nyman-search",
            "--index",
            str(checkpoint / "index.json"),
            "--checkpoint-dir",
            str(checkpoint),
            "--bits",
            "1536",
        ]
    )
    replay = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert replay["classification"] == "REPRODUCED_EXPLORATORY_NYMAN_SEARCH"
    assert replay["verified_cells"] == "0"
    assert replay["hypothesis_status"] == "UNRESOLVED"


def test_nyman_normalization_cli_generates_and_exactly_regenerates(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    artifact = tmp_path / "normalization.json"
    assert main(
        ["nyman-normalization-audit", "--output", str(artifact)]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "EXPLORATORY"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["audit_outcome"] == (
        "NORMALIZATION_AUDIT_PASSED_EXPLORATORY"
    )

    assert main(
        [
            "verify-nyman-normalization-audit",
            "--artifact",
            str(artifact),
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed["classification"] == (
        "REPRODUCED_EXPLORATORY_NYMAN_NORMALIZATION_AUDIT"
    )
    assert replayed["verified_ratio_comparisons"] == "14"
    assert replayed["harmonic_comparison_reproduced"] is True
    assert replayed["harmonic_comparison_passed"] is True


def test_nyman_summary_cli_generates_and_structurally_regenerates(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    summary = tmp_path / "summary.json"
    checkpoint = Path("results/nyman-natural-v1")
    assert main(
        [
            "summarize-nyman",
            "--checkpoint-dir",
            str(checkpoint),
            "--output",
            str(summary),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "EXPLORATORY"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["counts"]["cells"] == "6"

    assert main(
        [
            "verify-nyman-summary",
            "--summary",
            str(summary),
            "--checkpoint-dir",
            str(checkpoint),
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed["classification"] == (
        "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY"
    )
    assert replayed["verified_cells"] == "6"
    assert replayed["numerical_replay_performed"] is False


def test_parity_and_nesting_cli_artifacts_replay(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    parity = tmp_path / "parity.json"
    assert main(
        [
            "weil-parity-audit",
            "--cutoff-numerator",
            "3",
            "--cutoff-denominator",
            "2",
            "--degree",
            "0",
            "--bits",
            "96",
            "--output",
            str(parity),
        ]
    ) == 0
    generated_parity = json.loads(capsys.readouterr().out)
    assert generated_parity["decision"] == "AUDIT_PASSED_EXPLORATORY"
    assert main(
        ["verify-weil-parity-audit", "--artifact", str(parity)]
    ) == 0
    replayed_parity = json.loads(capsys.readouterr().out)
    assert replayed_parity["audit_kind"] == "PARITY"

    nesting = tmp_path / "nesting.json"
    assert main(
        [
            "weil-nesting-audit",
            "--cutoff-numerator",
            "3",
            "--cutoff-denominator",
            "2",
            "--lower-degree",
            "0",
            "--higher-degree",
            "1",
            "--bits",
            "96",
            "--replay-bits",
            "128",
            "--output",
            str(nesting),
        ]
    ) == 0
    generated_nesting = json.loads(capsys.readouterr().out)
    assert generated_nesting["decision"] == "AUDIT_PASSED_EXPLORATORY"
    assert main(
        ["verify-weil-nesting-audit", "--artifact", str(nesting)]
    ) == 0
    replayed_nesting = json.loads(capsys.readouterr().out)
    assert replayed_nesting["audit_kind"] == "DEGREE_NESTING"
