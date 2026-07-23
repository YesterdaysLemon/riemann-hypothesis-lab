from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

from flint import arb, ctx
import pytest

from riemann_lab.artifacts import content_sha256, write_json
from riemann_lab.balls import arb_from_dyadic
from riemann_lab.cli import main
import riemann_lab.nyman_mobius as mobius


REAL_SUMMARY = Path("results/nyman-natural-v1-summary.json")
REAL_CHECKPOINT = Path("results/nyman-natural-v1")
REAL_ARTIFACT = Path("results/nyman-mobius-core-tail-v1.json")
FROZEN_PAYLOAD_SHA256 = (
    "b6579bf5774e5413ee5f05b2d3f125d1afcd6f8e519a3942b7352924d01acb8e"
)
FROZEN_REPLAY_KERNEL_SHA256 = (
    "fdb64334dbb6b234d92b38f2b2956d19bf1afe1077c37562171be9aebf57aae5"
)


@pytest.fixture(scope="module")
def artifact() -> dict[str, Any]:
    return mobius.generate_nyman_mobius_audit(REAL_SUMMARY, REAL_CHECKPOINT)


def _fraction(record: dict[str, str]) -> Fraction:
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def _rehash(record: dict[str, Any]) -> None:
    record["payload_sha256"] = content_sha256(
        {key: value for key, value in record.items() if key != "payload_sha256"}
    )


def test_artifact_is_deterministic_frozen_and_restores_precision(
    artifact: dict[str, Any],
) -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 333
        regenerated = mobius.generate_nyman_mobius_audit(
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )
        assert ctx.prec == 333
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()

    assert regenerated == artifact
    assert artifact["schema"] == mobius.MOBIUS_SCHEMA
    assert artifact["audit_id"] == mobius.FROZEN_AUDIT_ID
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["audit_outcome"] == "FINITE_EXPLICIT_FAMILY_ENERGY_CERTIFIED"
    assert artifact["payload_sha256"] == FROZEN_PAYLOAD_SHA256
    assert artifact["payload_sha256"] == content_sha256(
        {key: value for key, value in artifact.items() if key != "payload_sha256"}
    )


def test_exact_mobius_divisor_manifest_covers_every_k(
    artifact: dict[str, Any],
) -> None:
    exact = artifact["generation"]["exact_arithmetic"]
    assert exact["max_n"] == "256"
    assert len(exact["mobius_values_n_1_through_256"]) == 256
    assert len(exact["prime_power_base_n_1_through_256"]) == 256
    assert exact["checks"] == {
        "checked_k": "256",
        "prime_power_count": "70",
        "sum_mu_over_divisors_is_delta_k1": True,
        "product_d_power_mu_d_is_inverse_prime_exactly_at_prime_powers": True,
        "derived_log_identity": "sum_(d|k) mu(d)*log(d) = -Lambda(k)",
        "derived_bcf_divisor_identity_through_cutoff": (
            "h_N(k)=1_(k=1)+Lambda(k)/log(N) for k<=N"
        ),
    }
    assert exact["mobius_values_n_1_through_256"][:10] == [
        "1",
        "-1",
        "-1",
        "0",
        "-1",
        "1",
        "-1",
        "0",
        "0",
        "1",
    ]
    assert exact["prime_power_base_n_1_through_256"][1:10] == [
        "2",
        "3",
        "2",
        "5",
        None,
        "7",
        "2",
        "3",
        None,
    ]


def test_two_core_paths_agree_and_core_tail_split_is_strict(
    artifact: dict[str, Any],
) -> None:
    records = artifact["generation"]["records"]
    assert [record["n"] for record in records] == [
        "8",
        "16",
        "32",
        "64",
        "128",
        "256",
    ]
    for record in records:
        checks = record["checks"]
        assert all(checks.values())
        total = arb_from_dyadic(record["candidate_total_error_squared"]["dyadic"])
        general = arb_from_dyadic(
            record["general_divisor_core_error_squared"]["dyadic"]
        )
        prime = arb_from_dyadic(
            record["prime_counting_core_error_squared"]["dyadic"]
        )
        core = arb_from_dyadic(record["certified_core_intersection"]["dyadic"])
        tail = arb_from_dyadic(
            record["tail_error_squared_by_total_minus_core"]["dyadic"]
        )
        assert general.overlaps(prime)
        # Reconstructing an Arb dyadic radius rounds outward once more, so
        # overlap is the stable serialized check for the stored intersection.
        assert general.overlaps(core)
        assert prime.overlaps(core)
        assert total > core > 0
        assert tail > 0
        assert (total - core).overlaps(tail)


def test_frozen_n256_values_and_core_share(
    artifact: dict[str, Any],
) -> None:
    record = artifact["generation"]["records"][-1]
    assert record["candidate_total_error_squared"]["display"].startswith(
        "[0.1242877077163881918488587124202360151830"
    )
    assert record["certified_core_intersection"]["display"].startswith(
        "[0.1227498440627963540023256270925295630070"
    )
    assert record["tail_error_squared_by_total_minus_core"]["display"].startswith(
        "[0.001537863653591837846533085327706452176070"
    )
    share = arb_from_dyadic(record["core_fraction_of_candidate_total"]["dyadic"])
    assert share * 100 > 98


def test_all_finite_grid_trends_are_certified_but_not_extrapolated(
    artifact: dict[str, Any],
) -> None:
    trends = artifact["generation"]["finite_grid_trends"]
    assert len(trends) == 7
    assert all(trends.values())
    assert "finite dimensions" in artifact["limitation"]
    assert "Finite decreases or fitted rates cannot prove" in artifact["limitation"]
    assert artifact["dependency_manifest"]["decisive_all_scale_target"][
        "proved_by_this_artifact"
    ] is False


def test_candidate_is_strictly_worse_than_the_certified_optimum(
    artifact: dict[str, Any],
) -> None:
    expected_factors = [34, 26, 21, 19, 16, 15]
    comparisons = artifact["candidate_vs_certified_optimum"]
    assert len(comparisons) == 6
    for comparison, factor in zip(comparisons, expected_factors, strict=True):
        lower = _fraction(comparison["source_optimal_lower_bound"])
        upper = _fraction(comparison["source_optimal_upper_bound"])
        assert lower < upper
        assert comparison["certified_strict_integer_factor"] == str(factor)
        assert arb_from_dyadic(comparison["strict_factor_margin"]["dyadic"]) > 0
        assert comparison["derived_inequality"] == (
            f"E_BCF({comparison['n']}) > {factor}*d_{comparison['n']}^2"
        )


def test_dependency_manifest_names_the_exact_unproved_bridge(
    artifact: dict[str, Any],
) -> None:
    manifest = artifact["dependency_manifest"]
    assert [source["url"] for source in manifest["published_sources"]] == [
        "https://arxiv.org/abs/math/0202141",
        "https://arxiv.org/abs/1211.5191",
        "https://arxiv.org/abs/math/0306251",
    ]
    assert all(
        source["reproved_by_this_artifact"] is False
        for source in manifest["published_sources"]
    )
    identity = manifest["manual_exact_identity"]
    assert "for 1<=t<N" in identity["bcf_core_residual"]
    assert "fails" in identity["tail_boundary"]
    target = manifest["decisive_all_scale_target"]
    assert "for every N>=N0" in target["statement"]
    assert "proves RH" in target["consequence"]
    assert target["proved_by_this_artifact"] is False


def test_direct_reciprocal_integral_controls_sign_jacobian_and_cutoff(
    artifact: dict[str, Any],
) -> None:
    """Cross-check the Gram total by directly integrating unit intervals."""

    n = 8
    truncation = 512
    with mobius._clean_precision(mobius.FROZEN_GENERATION_BITS):
        _, mu, _ = mobius._arithmetic_manifest()
        coefficients = mobius._bcf_coefficients(mu, n)
        _, boundary_primes = mobius._mobius_sieve(truncation)
        boundary_bases = mobius._prime_power_bases(truncation, boundary_primes)
        slope = sum(
            (coefficients[index - 1] / index for index in range(1, n + 1)),
            arb(0),
        )
        direct_partial = slope * slope
        cumulative_h = arb(0)
        psi = arb(0)
        q_at_n = None
        q_at_ten = None
        psi_q_at_n = None
        psi_q_at_ten = None
        for interval in range(1, truncation):
            h_value = sum(
                (
                    coefficients[divisor - 1]
                    for divisor in range(1, n + 1)
                    if interval % divisor == 0
                ),
                arb(0),
            )
            cumulative_h += h_value
            q_value = 1 - cumulative_h
            # This psi path is only a boundary control, not the tail path.
            if boundary_bases[interval]:
                psi += arb(boundary_bases[interval]).log()
            psi_q = -psi / arb(n).log()
            if interval == n:
                q_at_n = q_value
                psi_q_at_n = psi_q
            if interval == 10:
                q_at_ten = q_value
                psi_q_at_ten = psi_q
            interval_log = (arb(interval + 1) / interval).log()
            direct_partial += (
                slope * slope
                + 2 * slope * q_value * interval_log
                + q_value * q_value / (interval * (interval + 1))
            )

        assert q_at_n is not None and psi_q_at_n is not None
        assert (q_at_n - psi_q_at_n).contains(0)
        assert q_at_ten is not None and psi_q_at_ten is not None
        assert not q_at_ten.overlaps(psi_q_at_ten)

        absolute_bound = 1 + sum((abs(value) for value in coefficients), arb(0))
        remainder_bound = absolute_bound * absolute_bound / truncation
        total = arb_from_dyadic(
            artifact["generation"]["records"][0][
                "candidate_total_error_squared"
            ]["dyadic"]
        )
        assert total > direct_partial
        assert direct_partial + remainder_bound > total


def test_committed_artifact_exactly_regenerates(
    artifact: dict[str, Any],
) -> None:
    supplied = json.loads(REAL_ARTIFACT.read_text(encoding="utf-8"))
    assert supplied == artifact


def test_verifier_performs_512_bit_same_backend_replay(
    artifact: dict[str, Any],
) -> None:
    result = mobius.verify_nyman_mobius_audit(
        artifact,
        REAL_SUMMARY,
        REAL_CHECKPOINT,
    )
    assert result == {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_MOBIUS_AUDIT",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "FINITE_EXPLICIT_FAMILY_ENERGY_CERTIFIED",
        "artifact_payload_sha256": FROZEN_PAYLOAD_SHA256,
        "summary_payload_sha256": mobius.FROZEN_SUMMARY_PAYLOAD_SHA256,
        "verified_cells": "6",
        "generation_precision_bits": "256",
        "replay_precision_bits": "512",
        "replay_kernel_system_content_sha256": FROZEN_REPLAY_KERNEL_SHA256,
        "all_generation_enclosures_contain_replay": True,
        "all_comparison_enclosures_contain_replay": True,
        "same_backend_replay_only": True,
    }


def test_verifier_rejects_duplicate_keys_and_rehashed_mutations(
    artifact: dict[str, Any],
    tmp_path: Path,
) -> None:
    canonical = json.dumps(
        artifact,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text('{"schema":"ignored",' + canonical[1:], encoding="utf-8")
    with pytest.raises(
        mobius.NymanMobiusVerificationError,
        match="duplicate JSON object key: schema",
    ):
        mobius.verify_nyman_mobius_audit(
            duplicate,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )

    changed = copy.deepcopy(artifact)
    changed["generation"]["finite_grid_trends"][
        "tail_strictly_decreases_on_grid"
    ] = False
    _rehash(changed)
    with pytest.raises(
        mobius.NymanMobiusVerificationError,
        match="does not canonically regenerate",
    ):
        mobius.verify_nyman_mobius_audit(
            changed,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )


def test_generator_uses_one_detached_source_snapshot(
    artifact: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = json.loads(REAL_SUMMARY.read_text(encoding="utf-8"))
    real_verify = mobius.verify_nyman_summary

    def verify_then_mutate_caller(
        snapshot: dict[str, Any],
        checkpoint_dir: Path,
    ) -> dict[str, Any]:
        result = real_verify(snapshot, checkpoint_dir)
        source["cells"][-1]["bounds"]["upper"]["exact_fraction"][
            "numerator"
        ] = "0"
        return result

    monkeypatch.setattr(mobius, "verify_nyman_summary", verify_then_mutate_caller)
    regenerated = mobius.generate_nyman_mobius_audit(source, REAL_CHECKPOINT)
    assert regenerated == artifact


def test_cli_generates_and_verifies(
    artifact: dict[str, Any],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    path = tmp_path / "mobius.json"
    assert main(
        [
            "nyman-mobius-core-tail-audit",
            "--summary",
            str(REAL_SUMMARY),
            "--checkpoint-dir",
            str(REAL_CHECKPOINT),
            "--output",
            str(path),
        ]
    ) == 0
    generated = json.loads(capsys.readouterr().out)
    assert generated["payload_sha256"] == FROZEN_PAYLOAD_SHA256
    assert generated["certified_cells"] == "6"
    assert json.loads(path.read_text(encoding="utf-8")) == artifact

    assert main(
        [
            "verify-nyman-mobius-core-tail-audit",
            "--artifact",
            str(path),
            "--summary",
            str(REAL_SUMMARY),
            "--checkpoint-dir",
            str(REAL_CHECKPOINT),
        ]
    ) == 0
    verified = json.loads(capsys.readouterr().out)
    assert verified["artifact_payload_sha256"] == FROZEN_PAYLOAD_SHA256
    assert verified["all_generation_enclosures_contain_replay"] is True
    assert verified["all_comparison_enclosures_contain_replay"] is True


def test_write_json_round_trip_is_canonical(
    artifact: dict[str, Any],
    tmp_path: Path,
) -> None:
    path = tmp_path / "round-trip.json"
    write_json(path, artifact)
    assert json.loads(path.read_text(encoding="utf-8")) == artifact
