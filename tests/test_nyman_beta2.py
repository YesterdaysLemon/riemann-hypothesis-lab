from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path

import pytest
from flint import ctx

from riemann_lab.artifacts import content_sha256
import riemann_lab.nyman_beta2 as beta2


REAL_CANDIDATE = Path("results/nyman-beta2-n512-candidate-v1.json")
REAL_ARTIFACT = Path("results/nyman-beta2-n512-v1.json")
REAL_SUMMARY = Path("results/nyman-natural-v1-summary.json")
REAL_CHECKPOINT = Path("results/nyman-natural-v1")


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rehash(record: dict[str, object]) -> dict[str, object]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _fraction(record: object) -> Fraction:
    assert isinstance(record, dict)
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def test_frozen_candidate_is_exact_and_hash_bound() -> None:
    candidate = beta2._load_object(
        REAL_CANDIDATE,
        beta2.NymanBeta2VerificationError,
        label="candidate",
    )
    numerators = beta2._validate_candidate(
        candidate,
        beta2.NymanBeta2VerificationError,
        require_frozen=True,
    )

    assert len(numerators) == 512
    assert candidate["payload_sha256"] == beta2.FROZEN_CANDIDATE_PAYLOAD_SHA256
    assert (
        candidate["coefficient_vector_sha256"]
        == beta2.FROZEN_COEFFICIENT_VECTOR_SHA256
    )
    assert candidate["kernel"] == {
        "core_system_content_sha256": beta2.FROZEN_GENERATION_KERNEL_SHA256,
        "old_prefix_matches_frozen_v1_numeric_kernel": True,
        "prefix_256_numeric_sha256": beta2.FROZEN_N256_PREFIX_KERNEL_SHA256,
        "prefix_512_numeric_sha256": beta2.FROZEN_GENERATION_PREFIX_512_SHA256,
    }
    assert candidate["solver_role"] == (
        "approximate-untrusted-candidate-generator-only"
    )


def test_frozen_audit_has_only_the_claimed_finite_scope() -> None:
    candidate = _load(REAL_CANDIDATE)
    artifact = _load(REAL_ARTIFACT)
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}

    assert artifact["payload_sha256"] == content_sha256(body)
    assert artifact["payload_sha256"] == beta2.FROZEN_AUDIT_PAYLOAD_SHA256
    assert artifact["schema"] == beta2.AUDIT_SCHEMA
    assert artifact["audit_id"] == beta2.FROZEN_AUDIT_ID
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["candidate"] == candidate
    assert artifact["audit_outcome"] == (
        "N512_STRONG_FINITE_CONTRACTION_CERTIFIED"
    )
    assert artifact["n512_lower_bound_attempted"] is False
    assert artifact["theorem_bridge"]["resolves_rh"] is False
    assert "CERTIFIED_FINITE" in artifact["limitation"]
    assert "no all-scale" in artifact["limitation"]
    assert "no proof or disproof" in artifact["limitation"]


def test_exact_endpoint_chain_certifies_stronger_than_beta_two() -> None:
    comparison = _load(REAL_ARTIFACT)["generation"]["exact_comparisons"]
    u512 = _fraction(comparison["stored_upper_bound_u512"]["exact_fraction"])
    l256 = _fraction(comparison["source_lower_bound_l256"])
    u256 = _fraction(comparison["source_upper_bound_u256"])
    strong_threshold = _fraction(
        comparison["strong_threshold_449_over_500_times_l256"]
    )
    beta_two_threshold = _fraction(
        comparison["beta_two_threshold_9_over_10_times_l256"]
    )
    schur_safe_threshold = _fraction(
        comparison["schur_safe_threshold_l256_minus_u256_over_10"]
    )

    assert u512 == beta2.FROZEN_U512
    assert l256 == beta2.FROZEN_L256
    assert u256 == beta2.FROZEN_U256
    assert strong_threshold == Fraction(449, 500) * l256
    assert beta_two_threshold == Fraction(9, 10) * l256
    assert schur_safe_threshold == l256 - u256 / 10
    assert u512 < strong_threshold < beta_two_threshold
    assert u512 < schur_safe_threshold
    assert _fraction(comparison["strong_exact_positive_margin"]) == (
        strong_threshold - u512
    )
    assert _fraction(comparison["beta_two_exact_positive_margin"]) == (
        beta_two_threshold - u512
    )
    assert _fraction(comparison["schur_safe_exact_positive_margin"]) == (
        schur_safe_threshold - u512
    )
    assert comparison["checks"] == {
        "449_over_500_strictly_below_9_over_10": True,
        "all_comparisons_exact_rational": True,
        "u512_strictly_below_449_over_500_l256": True,
        "u512_strictly_below_9_over_10_l256": True,
        "u512_strictly_below_l256_minus_u256_over_10": True,
    }


def test_stored_terminal_certificates_are_strict_and_hash_bound() -> None:
    artifact = _load(REAL_ARTIFACT)
    generation = artifact["generation"]
    lower = generation["source_n256_lower_certificate_replay"]
    upper = generation["n512_direct_upper_certificate"]

    for certificate in (lower, upper):
        body = {
            key: value
            for key, value in certificate.items()
            if key != "payload_sha256"
        }
        assert certificate["payload_sha256"] == content_sha256(body)
        assert certificate["hypothesis_status"] == "UNRESOLVED"
    assert lower["decision"] == "LOWER_BOUND_CERTIFIED"
    assert lower["checks"]["fixed_order_interval_ldlt_positive"] is True
    assert len(lower["augmented_ldlt"]["pivots"]) == 257
    assert upper["decision"] == "UPPER_BOUND_CERTIFIED"
    assert upper["checks"]["energy_at_most_claimed_upper_bound"] is True
    assert upper["checks"]["energy_strictly_below_claimed_upper_bound"] is True
    assert upper["checks"]["approximate_solve_used_as_evidence"] is False
    assert len(upper["coefficients"]["numerators"]) == 512
    assert _fraction(upper["claimed_upper_bound"]) == beta2.FROZEN_U512


def test_frozen_source_summary_regenerates_structurally() -> None:
    summary = beta2._load_object(
        REAL_SUMMARY,
        beta2.NymanBeta2Error,
        label="source summary",
    )
    verification, cell = beta2._validate_source_summary(summary, REAL_CHECKPOINT)

    assert verification["summary_payload_sha256"] == (
        beta2.FROZEN_SUMMARY_PAYLOAD_SHA256
    )
    assert verification["numerical_replay_performed"] is False
    assert cell["hashes"]["cell_payload_sha256"] == (
        beta2.FROZEN_N256_CELL_PAYLOAD_SHA256
    )


def test_candidate_rejects_coefficient_mutation_even_when_outer_hash_is_repaired() -> None:
    candidate = _load(REAL_CANDIDATE)
    candidate["coefficients"]["numerators"][0] = str(
        int(candidate["coefficients"]["numerators"][0]) + 1
    )
    candidate = _rehash(candidate)

    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="coefficient vector hash mismatch",
    ):
        beta2._validate_candidate(
            candidate,
            beta2.NymanBeta2VerificationError,
            require_frozen=True,
        )


def test_candidate_rejects_fully_rehashed_alternate_vector() -> None:
    candidate = _load(REAL_CANDIDATE)
    candidate["coefficients"]["numerators"][0] = str(
        int(candidate["coefficients"]["numerators"][0]) + 1
    )
    candidate["coefficient_vector_sha256"] = content_sha256(
        candidate["coefficients"]
    )
    candidate = _rehash(candidate)

    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="frozen candidate payload hash changed",
    ):
        beta2._validate_candidate(
            candidate,
            beta2.NymanBeta2VerificationError,
            require_frozen=True,
        )


def test_strict_loader_rejects_duplicate_keys(tmp_path: Path) -> None:
    raw = REAL_CANDIDATE.read_text(encoding="utf-8")
    path = tmp_path / "duplicate.json"
    path.write_text('{"schema":"forged",' + raw.lstrip()[1:], encoding="utf-8")

    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="duplicate JSON object key: schema",
    ):
        beta2._load_object(
            path,
            beta2.NymanBeta2VerificationError,
            label="candidate",
        )


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_strict_loader_rejects_nonstandard_constants(
    constant: str, tmp_path: Path
) -> None:
    raw = REAL_CANDIDATE.read_text(encoding="utf-8")
    path = tmp_path / "constant.json"
    path.write_text('{"extra":' + constant + "," + raw.lstrip()[1:], encoding="utf-8")

    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match=f"nonstandard JSON constant: {constant}",
    ):
        beta2._load_object(
            path,
            beta2.NymanBeta2VerificationError,
            label="candidate",
        )


def test_mapping_loader_rejects_float() -> None:
    candidate = _load(REAL_CANDIDATE)
    candidate["extra"] = 0.898

    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="floating-point JSON number is noncanonical",
    ):
        beta2._load_object(
            candidate,
            beta2.NymanBeta2VerificationError,
            label="candidate",
        )


@pytest.mark.parametrize("bits", [True, 0, 1535, 1537])
def test_verifier_rejects_weak_replay_precision_before_io(bits: object) -> None:
    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="replay precision must equal the frozen 1536 bits",
    ):
        beta2.verify_nyman_beta2_audit(
            Path("missing-audit.json"),
            Path("missing-candidate.json"),
            Path("missing-summary.json"),
            Path("missing-checkpoint"),
            replay_precision_bits=bits,  # type: ignore[arg-type]
        )


def test_comparison_builder_fails_closed_if_endpoint_is_weakened(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        beta2,
        "FROZEN_U512",
        beta2.STRONG_CONTRACTION * beta2.FROZEN_L256,
    )

    with pytest.raises(
        beta2.NymanBeta2Error,
        match="misses a contraction threshold",
    ):
        beta2._comparison_record()


def test_rehashed_audit_mutation_is_rejected_before_expensive_regeneration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _load(REAL_ARTIFACT)
    artifact["terminal_statement"] = "forged"
    artifact = _rehash(artifact)
    called = False

    def forbidden_derive(*args: object, **kwargs: object) -> object:
        nonlocal called
        called = True
        raise AssertionError("expensive regeneration must not run")

    monkeypatch.setattr(beta2, "_derive_audit", forbidden_derive)
    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="frozen audit payload hash changed",
    ):
        beta2.verify_nyman_beta2_audit(
            artifact,
            REAL_CANDIDATE,
            REAL_SUMMARY,
            REAL_CHECKPOINT,
        )
    assert called is False


def test_clean_precision_restores_context_after_exception() -> None:
    original = ctx.prec
    with pytest.raises(RuntimeError, match="control"):
        with beta2._clean_precision(333):
            assert ctx.prec == 333
            raise RuntimeError("control")
    assert ctx.prec == original


def _arb_record(
    mid_mantissa: int,
    mid_exponent: int,
    radius_mantissa: int,
    radius_exponent: int,
) -> dict[str, object]:
    return {
        "display": "test-only",
        "dyadic": {
            "mid_mantissa": str(mid_mantissa),
            "mid_exponent": str(mid_exponent),
            "radius_mantissa": str(radius_mantissa),
            "radius_exponent": str(radius_exponent),
        },
        "is_exact": radius_mantissa == 0,
    }


def test_replay_containment_accepts_subset_and_rejects_escape() -> None:
    generation = _arb_record(0, 0, 2, 0)
    contained = _arb_record(0, 0, 1, 0)
    escaped = _arb_record(2, 0, 1, 0)

    beta2._assert_replay_contained(generation, contained, "test")
    with pytest.raises(
        beta2.NymanBeta2VerificationError,
        match="escaped the generation enclosure",
    ):
        beta2._assert_replay_contained(generation, escaped, "test")


def test_candidate_generation_wraps_backend_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_build(*args: object, **kwargs: object) -> object:
        raise ArithmeticError("backend control")

    monkeypatch.setattr(beta2.core, "build_natural_system", fail_build)
    with pytest.raises(
        beta2.NymanBeta2Error,
        match="canonical N=512 candidate generation failed",
    ):
        beta2.propose_nyman_beta2_candidate()
