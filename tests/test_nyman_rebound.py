from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from flint import ctx
import pytest

from riemann_lab.artifacts import content_sha256, write_json
from riemann_lab.balls import arb_from_dyadic
from riemann_lab.cli import main
import riemann_lab.nyman_rebound as rebound


REAL_SUMMARY = Path("results/nyman-natural-v1-summary.json")
REAL_CHECKPOINT = Path("results/nyman-natural-v1")
REAL_REBOUND = Path("results/nyman-forced-rebound-v1.json")
FROZEN_REBOUND_PAYLOAD_SHA256 = (
    "4f5ecef5798cd273774aa4e2d05e7ac3d6b31cf241cadfe632742ea879bff825"
)


def _artifact() -> dict[str, Any]:
    return rebound.generate_nyman_rebound_audit(
        REAL_SUMMARY,
        REAL_CHECKPOINT,
    )


def _rehash(record: dict[str, Any]) -> None:
    record["payload_sha256"] = content_sha256(
        {key: value for key, value in record.items() if key != "payload_sha256"}
    )


def _fraction(record: dict[str, str]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def test_real_rebound_audit_is_deterministic_and_binds_frozen_sources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls = 0
    real_verify = rebound.verify_nyman_summary

    def tracked(*args: Any, **kwargs: Any) -> dict[str, Any]:
        nonlocal calls
        calls += 1
        return real_verify(*args, **kwargs)

    monkeypatch.setattr(rebound, "verify_nyman_summary", tracked)
    previous_precision = ctx.prec
    try:
        ctx.prec = 333
        artifact = _artifact()
        assert ctx.prec == 333
        assert artifact == _artifact()
        assert ctx.prec == 333
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()

    assert calls == 2
    assert artifact["schema"] == rebound.REBOUND_SCHEMA
    assert artifact["audit_id"] == rebound.FROZEN_AUDIT_ID
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["audit_outcome"] == (
        "FINITE_FORCED_REBOUND_PRECONDITIONS_CERTIFIED"
    )
    assert artifact["payload_sha256"] == content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )
    assert artifact["payload_sha256"] == FROZEN_REBOUND_PAYLOAD_SHA256
    assert artifact["source_verification"] == {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": rebound.FROZEN_SUMMARY_PAYLOAD_SHA256,
        "plan_sha256": rebound.FROZEN_PLAN_SHA256,
        "index_payload_sha256": rebound.FROZEN_INDEX_PAYLOAD_SHA256,
        "verified_cells": "6",
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }
    assert artifact["frozen_source_bounds"] == {
        "n128_lower": {
            "numerator": "3286634009378442248435358862600663947",
            "denominator": "340282366920938463463374607431768211456",
        },
        "n256_lower": {
            "numerator": "2801788463381838392697853210438663515",
            "denominator": "340282366920938463463374607431768211456",
        },
        "n256_upper": {
            "numerator": "2801788463381838392697853210438663771",
            "denominator": "340282366920938463463374607431768211456",
        },
    }
    assert "does not locate a rebound" in artifact["limitation"]
    assert "UNRESOLVED" in artifact["limitation"]


def test_all_five_scaled_declines_are_exact_positive_rational_implications() -> None:
    artifact = _artifact()
    bundle = artifact["five_exact_scaled_declines"]
    assert bundle["count"] == "5"
    assert bundle["all_certified"] is True
    assert [record["k"] for record in bundle["records"]] == [
        "3",
        "4",
        "5",
        "6",
        "7",
    ]

    for record in bundle["records"]:
        k = int(record["k"])
        n = 1 << k
        lower = _fraction(record["source_lower_bound"])
        next_upper = _fraction(record["source_next_upper_bound"])
        factor = _fraction(record["contraction_factor"])
        margin = _fraction(record["exact_positive_margin"])
        assert record["n"] == str(n)
        assert record["next_n"] == str(2 * n)
        assert factor == Fraction(k, k + 1)
        assert margin == factor * lower - next_upper
        assert margin > 0
        assert record["decision"] == "SCALED_DECLINE_CERTIFIED"
        assert record["derived_strict_inequality"] == (
            f"{k + 1}*d_{2 * n}^2 < {k}*d_{n}^2"
        )


def test_n256_arb_chain_has_two_strict_positive_margins() -> None:
    gap = _artifact()["n256_scaled_gap"]
    assert gap["precision_bits"] == "256"
    assert _fraction(gap["rational_separator"]) == Fraction(23, 500)
    assert gap["checks"] == {
        "d256_squared_at_most_exact_source_upper": True,
        "scaled_upper_strictly_below_23_over_500": True,
        "23_over_500_strictly_below_c0": True,
        "c0_strictly_below_9_l256_log_2": True,
    }
    assert arb_from_dyadic(gap["finite_positive_margin"]["dyadic"]) > 0
    assert arb_from_dyadic(gap["asymptotic_positive_margin"]["dyadic"]) > 0
    assert arb_from_dyadic(gap["beta_two_anchor_floor"]["dyadic"]) > 0
    assert arb_from_dyadic(gap["beta_two_compatibility_margin"]["dyadic"]) > 0
    assert gap["certified_chain"] == "d_256^2*log(256) < 23/500 < C0"


def test_dependency_manifest_is_explicit_about_published_theorems_and_backend() -> None:
    manifest = _artifact()["dependency_manifest"]
    assert manifest["finite_inputs"] == {
        "summary_schema": "rh-lab/nyman-natural-summary/v1",
        "summary_payload_sha256": rebound.FROZEN_SUMMARY_PAYLOAD_SHA256,
        "index_payload_sha256": rebound.FROZEN_INDEX_PAYLOAD_SHA256,
        "plan_sha256": rebound.FROZEN_PLAN_SHA256,
        "required_certified_cells": ["8", "16", "32", "64", "128", "256"],
    }
    assert manifest["numeric_backend"]["precision_bits"] == "256"
    assert manifest["numeric_backend"]["python_flint"]
    assert manifest["numeric_backend"]["flint"]
    citations = manifest["published_analytic_dependencies"]
    assert [item["url"] for item in citations] == [
        "https://arxiv.org/abs/math/0202141",
        "https://arxiv.org/abs/math/0103058",
        "https://arxiv.org/abs/1211.5191",
    ]
    assert all(item["reproved_by_this_artifact"] is False for item in citations)
    inference = manifest["analytic_inference"]
    assert inference["statement"] == (
        "The published-theorem dichotomy gives "
        "liminf_N d_N^2*log(N) >= C0; together with the finite chain, "
        "a later dyadic scaled value must exceed the N=256 value."
    )
    assert len(inference["derivation_steps"]) == 4
    assert "Hadamard-product zero-sum identity" in inference["derivation_steps"][2]
    assert inference["machine_checked_by_this_artifact"] is False
    assert inference["locates_rebound"] is False
    assert inference["resolves_rh"] is False


def test_committed_rebound_artifact_exactly_regenerates() -> None:
    supplied = json.loads(REAL_REBOUND.read_text(encoding="utf-8"))
    generated = _artifact()
    assert supplied == generated
    assert supplied["payload_sha256"] == FROZEN_REBOUND_PAYLOAD_SHA256
    assert rebound.verify_nyman_rebound_audit(
        REAL_REBOUND,
        REAL_SUMMARY,
        REAL_CHECKPOINT,
    )["artifact_payload_sha256"] == generated["payload_sha256"]


def test_rebound_verifier_accepts_mapping_and_path_and_exactly_regenerates(
    tmp_path: Path,
) -> None:
    artifact = _artifact()
    path = tmp_path / "rebound.json"
    write_json(path, artifact)
    expected = {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_REBOUND_AUDIT",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "FINITE_FORCED_REBOUND_PRECONDITIONS_CERTIFIED",
        "artifact_payload_sha256": artifact["payload_sha256"],
        "summary_payload_sha256": rebound.FROZEN_SUMMARY_PAYLOAD_SHA256,
        "index_payload_sha256": rebound.FROZEN_INDEX_PAYLOAD_SHA256,
        "plan_sha256": rebound.FROZEN_PLAN_SHA256,
        "verified_scaled_declines": "5",
        "precision_bits": "256",
    }
    assert rebound.verify_nyman_rebound_audit(
        artifact,
        REAL_SUMMARY,
        REAL_CHECKPOINT,
    ) == expected
    assert rebound.verify_nyman_rebound_audit(
        path,
        REAL_SUMMARY,
        REAL_CHECKPOINT,
    ) == expected


def test_rebound_verifier_rejects_duplicate_json_keys(tmp_path: Path) -> None:
    canonical = json.dumps(
        _artifact(),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"ignored",' + canonical[1:], encoding="utf-8")
    with pytest.raises(
        rebound.NymanReboundVerificationError,
        match="duplicate JSON object key: schema",
    ):
        rebound.verify_nyman_rebound_audit(
            path,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_rebound_verifier_rejects_nonstandard_json_constants(
    constant: str,
    tmp_path: Path,
) -> None:
    canonical = json.dumps(
        _artifact(),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    path = tmp_path / "nonstandard.json"
    path.write_text(
        '{"extra":' + constant + "," + canonical[1:],
        encoding="utf-8",
    )
    with pytest.raises(
        rebound.NymanReboundVerificationError,
        match=f"nonstandard JSON constant: {constant}",
    ):
        rebound.verify_nyman_rebound_audit(
            path,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )


def test_rebound_verifier_rejects_non_finite_mapping_value() -> None:
    artifact = _artifact()
    artifact["extra"] = float("nan")
    with pytest.raises(
        rebound.NymanReboundVerificationError,
        match="non-finite JSON number",
    ):
        rebound.verify_nyman_rebound_audit(
            artifact,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )


@pytest.mark.parametrize(
    ("mutator", "expected"),
    [
        (
            lambda artifact: artifact.__setitem__("unknown", "field"),
            "does not canonically regenerate",
        ),
        (
            lambda artifact: artifact["n256_scaled_gap"]["checks"].__setitem__(
                "23_over_500_strictly_below_c0",
                False,
            ),
            "does not canonically regenerate",
        ),
        (
            lambda artifact: artifact["five_exact_scaled_declines"]["records"][
                0
            ]["exact_positive_margin"].__setitem__("numerator", "1"),
            "does not canonically regenerate",
        ),
    ],
)
def test_rebound_verifier_rejects_rehashed_semantic_or_field_mutations(
    mutator: Any,
    expected: str,
) -> None:
    artifact = copy.deepcopy(_artifact())
    mutator(artifact)
    _rehash(artifact)
    with pytest.raises(
        rebound.NymanReboundVerificationError,
        match=expected,
    ):
        rebound.verify_nyman_rebound_audit(
            artifact,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )


def test_generator_rejects_duplicate_source_summary_keys(tmp_path: Path) -> None:
    canonical = REAL_SUMMARY.read_text(encoding="utf-8")
    path = tmp_path / "duplicate-summary.json"
    path.write_text(
        '{"schema":"ignored",' + canonical.lstrip()[1:],
        encoding="utf-8",
    )
    with pytest.raises(
        rebound.NymanReboundError,
        match="duplicate JSON object key: schema",
    ):
        rebound.generate_nyman_rebound_audit(path, REAL_CHECKPOINT)


def test_generator_uses_the_same_deep_snapshot_that_it_verifies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = json.loads(REAL_SUMMARY.read_text(encoding="utf-8"))
    real_verify = rebound.verify_nyman_summary

    def verify_then_mutate_caller(
        snapshot: dict[str, Any],
        checkpoint_dir: Path,
    ) -> dict[str, Any]:
        result = real_verify(snapshot, checkpoint_dir)
        source["cells"][-1]["bounds"]["upper"]["exact_fraction"][
            "numerator"
        ] = "0"
        return result

    monkeypatch.setattr(rebound, "verify_nyman_summary", verify_then_mutate_caller)
    artifact = rebound.generate_nyman_rebound_audit(source, REAL_CHECKPOINT)
    assert artifact["frozen_source_bounds"]["n256_upper"]["numerator"] == (
        "2801788463381838392697853210438663771"
    )


def test_generator_normalizes_hostile_nested_mapping_before_verification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class StickyMapping(dict[str, Any]):
        def __deepcopy__(self, memo: dict[int, Any]) -> StickyMapping:
            return self

    source = json.loads(REAL_SUMMARY.read_text(encoding="utf-8"))
    exact_upper = source["cells"][-1]["bounds"]["upper"]["exact_fraction"]
    hostile = StickyMapping(exact_upper)
    source["cells"][-1]["bounds"]["upper"]["exact_fraction"] = hostile
    real_verify = rebound.verify_nyman_summary

    def verify_then_mutate_hostile_source(
        snapshot: dict[str, Any],
        checkpoint_dir: Path,
    ) -> dict[str, Any]:
        result = real_verify(snapshot, checkpoint_dir)
        hostile["numerator"] = "1"
        return result

    monkeypatch.setattr(
        rebound,
        "verify_nyman_summary",
        verify_then_mutate_hostile_source,
    )
    artifact = rebound.generate_nyman_rebound_audit(source, REAL_CHECKPOINT)
    assert artifact["frozen_source_bounds"]["n256_upper"]["numerator"] == (
        "2801788463381838392697853210438663771"
    )


def test_rebound_cli_generates_and_exactly_regenerates(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    artifact_path = tmp_path / "rebound.json"
    assert main(
        [
            "nyman-rebound-audit",
            "--summary",
            str(REAL_SUMMARY),
            "--checkpoint-dir",
            str(REAL_CHECKPOINT),
            "--output",
            str(artifact_path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["classification"] == "EXPLORATORY"
    assert generated["hypothesis_status"] == "UNRESOLVED"
    assert generated["scaled_declines"] == "5"

    assert main(
        [
            "verify-nyman-rebound-audit",
            "--artifact",
            str(artifact_path),
            "--summary",
            str(REAL_SUMMARY),
            "--checkpoint-dir",
            str(REAL_CHECKPOINT),
        ]
    ) == 0
    replayed = json.loads(capsys.readouterr().out)
    assert replayed["classification"] == (
        "REPRODUCED_EXPLORATORY_NYMAN_REBOUND_AUDIT"
    )
    assert replayed["verified_scaled_declines"] == "5"
