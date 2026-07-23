from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

import riemann_lab.cli as cli
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


def test_nyman_beta2_cli_routes_candidate_audit_and_verifier(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    candidate_path = tmp_path / "candidate.json"
    audit_path = tmp_path / "audit.json"
    fake_candidate = {
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "payload_sha256": "a" * 64,
        "solver_role": "approximate-untrusted-candidate-generator-only",
    }
    fake_audit = {
        "classification": "CERTIFIED_FINITE",
        "audit_outcome": "N512_STRONG_FINITE_CONTRACTION_CERTIFIED",
        "hypothesis_status": "UNRESOLVED",
        "payload_sha256": "b" * 64,
    }
    calls: list[tuple[object, ...]] = []

    monkeypatch.setattr(
        cli, "propose_nyman_beta2_candidate", lambda: fake_candidate
    )

    def fake_generate(*args: object) -> dict[str, str]:
        calls.append(args)
        return fake_audit

    def fake_verify(*args: object, **kwargs: object) -> dict[str, object]:
        calls.append((*args, kwargs))
        return {
            "classification": "REPRODUCED_CERTIFIED_FINITE_NYMAN_CONTRACTION",
            "hypothesis_status": "UNRESOLVED",
        }

    monkeypatch.setattr(cli, "generate_nyman_beta2_audit", fake_generate)
    monkeypatch.setattr(cli, "verify_nyman_beta2_audit", fake_verify)

    assert main(
        ["propose-nyman-beta2-n512", "--output", str(candidate_path)]
    ) == 0
    proposed = json.loads(capsys.readouterr().out)
    assert proposed["classification"] == "EXPLORATORY"
    assert json.loads(candidate_path.read_text(encoding="utf-8")) == fake_candidate

    assert main(
        [
            "nyman-beta2-n512-audit",
            "--candidate",
            str(candidate_path),
            "--summary",
            "summary.json",
            "--checkpoint-dir",
            "checkpoint",
            "--output",
            str(audit_path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "CERTIFIED_FINITE"
    assert json.loads(audit_path.read_text(encoding="utf-8")) == fake_audit

    assert main(
        [
            "verify-nyman-beta2-n512-audit",
            "--artifact",
            str(audit_path),
            "--candidate",
            str(candidate_path),
            "--summary",
            "summary.json",
            "--checkpoint-dir",
            "checkpoint",
            "--bits",
            "1536",
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed["classification"] == (
        "REPRODUCED_CERTIFIED_FINITE_NYMAN_CONTRACTION"
    )
    assert calls[-1][-1] == {"replay_precision_bits": 1536}


def test_nyman_trial_cli_routes_audit_and_verifier(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    audit_path = tmp_path / "trial-audit.json"
    fake_audit = {
        "classification": "CERTIFIED_FINITE",
        "audit_outcome": (
            "CANONICAL_EIGHT_COLUMN_TRIAL_SUBSPACE_REJECTED_AT_N512"
        ),
        "hypothesis_status": "UNRESOLVED",
        "payload_sha256": "c" * 64,
    }
    calls: list[tuple[object, ...]] = []

    def fake_generate(*args: object) -> dict[str, str]:
        calls.append(args)
        return fake_audit

    def fake_verify(*args: object, **kwargs: object) -> dict[str, object]:
        calls.append((*args, kwargs))
        return {
            "classification": (
                "REPRODUCED_CERTIFIED_FINITE_NYMAN_TRIAL_SUBSPACE_REJECTION"
            ),
            "hypothesis_status": "UNRESOLVED",
        }

    monkeypatch.setattr(
        cli, "generate_nyman_trial_subspace_audit", fake_generate
    )
    monkeypatch.setattr(
        cli, "verify_nyman_trial_subspace_audit", fake_verify
    )

    assert main(
        [
            "nyman-trial-subspace-audit",
            "--summary",
            "summary.json",
            "--checkpoint-dir",
            "checkpoint",
            "--output",
            str(audit_path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "CERTIFIED_FINITE"
    assert json.loads(audit_path.read_text(encoding="utf-8")) == fake_audit

    assert main(
        [
            "verify-nyman-trial-subspace-audit",
            "--artifact",
            str(audit_path),
            "--summary",
            "summary.json",
            "--checkpoint-dir",
            "checkpoint",
            "--bits",
            "1536",
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed["classification"] == (
        "REPRODUCED_CERTIFIED_FINITE_NYMAN_TRIAL_SUBSPACE_REJECTION"
    )
    assert calls[-1][-1] == {"replay_precision_bits": 1536}


def test_vasyunin_greedy_cli_generates_and_exactly_replays(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    artifact_path = tmp_path / "vasyunin-audit.json"

    assert main(
        [
            "nyman-vasyunin-greedy-audit",
            "--limit",
            "32",
            "--output",
            str(artifact_path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "CERTIFIED_FINITE"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["limit"] == "32"

    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert artifact["payload_sha256"] == generated["payload_sha256"]
    assert artifact["finite_scope"]["integer_intervals"] == "[1,33)"

    assert main(
        [
            "verify-nyman-vasyunin-greedy-audit",
            "--artifact",
            str(artifact_path),
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed == {
        "classification": (
            "REPRODUCED_CERTIFIED_FINITE_VASYUNIN_GREEDY_PREFIX"
        ),
        "hypothesis_status": "UNRESOLVED",
        "infinite_divergence_machine_reproved": False,
        "payload_sha256": artifact["payload_sha256"],
        "verified_limit": "32",
    }

    duplicate_path = tmp_path / "duplicate-vasyunin-audit.json"
    duplicate_path.write_text(
        '{"schema":"first","schema":"second"}',
        encoding="utf-8",
    )
    assert main(
        [
            "verify-nyman-vasyunin-greedy-audit",
            "--artifact",
            str(duplicate_path),
        ]
    ) == 2
    assert "duplicate JSON object key: schema" in capsys.readouterr().out

    deeply_nested_path = tmp_path / "deeply-nested-vasyunin-audit.json"
    deeply_nested_path.write_text(
        "[" * 10_000 + "0" + "]" * 10_000,
        encoding="utf-8",
    )
    assert main(
        [
            "verify-nyman-vasyunin-greedy-audit",
            "--artifact",
            str(deeply_nested_path),
        ]
    ) == 2
    assert "cannot read or parse Vasyunin audit" in capsys.readouterr().out


def test_alias_sharp_truncation_cli_generates_and_exactly_replays(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    artifact_path = tmp_path / "alias-audit.json"

    assert main(
        [
            "nyman-alias-sharp-truncation-audit",
            "--limit",
            "64",
            "--output",
            str(artifact_path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "CERTIFIED_FINITE"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["limit"] == "64"

    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert artifact["payload_sha256"] == generated["payload_sha256"]
    assert artifact["frozen_shell"]["target_step"] == (
        "f_y(t)=-1 on [9,16) and 0 elsewhere"
    )

    assert main(
        [
            "verify-nyman-alias-sharp-truncation-audit",
            "--artifact",
            str(artifact_path),
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed == {
        "classification": "REPRODUCED_CERTIFIED_FINITE_ALIAS_IDENTITIES",
        "hypothesis_status": "UNRESOLVED",
        "infinite_sharp_truncation_obstruction_machine_reproved": False,
        "payload_sha256": artifact["payload_sha256"],
        "verified_limit": "64",
    }

    duplicate_path = tmp_path / "duplicate-alias-audit.json"
    duplicate_path.write_text(
        '{"schema":"first","schema":"second"}',
        encoding="utf-8",
    )
    assert main(
        [
            "verify-nyman-alias-sharp-truncation-audit",
            "--artifact",
            str(duplicate_path),
        ]
    ) == 2
    assert "duplicate JSON object key: schema" in capsys.readouterr().out


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
