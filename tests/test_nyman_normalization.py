from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Callable

from flint import arb, ctx
import pytest

from riemann_lab.artifacts import content_sha256
import riemann_lab.nyman_normalization as normalization
import riemann_lab.nyman_oracle as oracle


FROZEN_BUNDLE_SHA256 = (
    "021a060fb3eb0121c325a0af13cc39095c0c6f26bc5a4b235269611a66a9442c"
)


@pytest.fixture(scope="module")
def frozen_bundle() -> dict[str, Any]:
    return normalization.generate_nyman_normalization_bundle()


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def test_frozen_bundle_passes_both_normalization_checks_and_is_hash_bound(
    frozen_bundle: dict[str, Any],
) -> None:
    bundle = frozen_bundle
    assert bundle["schema"] == normalization.BUNDLE_SCHEMA
    assert bundle["classification"] == "EXPLORATORY"
    assert bundle["hypothesis_status"] == "UNRESOLVED"
    assert bundle["bundle_id"] == normalization.FROZEN_BUNDLE_ID
    assert bundle["audit_outcome"] == (
        "NORMALIZATION_AUDIT_PASSED_EXPLORATORY"
    )
    assert bundle["payload_sha256"] == FROZEN_BUNDLE_SHA256
    assert bundle["payload_sha256"] == content_sha256(
        {key: value for key, value in bundle.items() if key != "payload_sha256"}
    )

    scope = bundle["frozen_scope"]
    assert scope == normalization._frozen_scope()
    assert scope["ratio_count"] == "14"
    assert scope["truncation"] == "4096"
    assert scope["harmonic_last_index"] == "4096"
    assert scope["precision_bits"] == "128"
    assert scope["generic_tail_upper_bound"] == {
        "numerator": "1",
        "denominator": "4096",
    }

    formula = bundle["formula_vs_truncated_integral"]
    assert formula["audit_outcome"] == (
        "ALL_14_CONSISTENT_WITH_TRUNCATION_ORACLE"
    )
    assert formula["evaluator_label"] == normalization.FORMULA_EVALUATOR_LABEL
    assert len(formula["comparisons"]) == 14
    assert {item["outcome"] for item in formula["comparisons"]} == {
        "CONSISTENT_WITH_TRUNCATION_ORACLE"
    }
    assert formula["comparisons"][5]["ratio"] == {
        "supplied": {"numerator": "4", "denominator": "6"},
        "reduced": {"numerator": "2", "denominator": "3"},
        "supplied_was_reduced": False,
        "canonicalization": (
            "fractions.Fraction gcd reduction with positive denominator"
        ),
    }
    assert formula["payload_sha256"] == content_sha256(
        {key: value for key, value in formula.items() if key != "payload_sha256"}
    )

    harmonic = bundle["one_minus_euler_gamma"]
    assert harmonic["schema"] == normalization.HARMONIC_AUDIT_SCHEMA
    assert harmonic["audit_outcome"] == (
        "CORE_ONE_MINUS_EULER_GAMMA_OVERLAPS_HARMONIC_ORACLE"
    )
    assert harmonic["comparison"] == "rigorous Arb enclosure overlap"
    assert harmonic["independent_harmonic_oracle"]["last_harmonic_index"] == (
        "4096"
    )
    assert harmonic["payload_sha256"] == content_sha256(
        {key: value for key, value in harmonic.items() if key != "payload_sha256"}
    )
    harmonic_oracle = harmonic["independent_harmonic_oracle"]
    assert harmonic_oracle["payload_sha256"] == content_sha256(
        {
            key: value
            for key, value in harmonic_oracle.items()
            if key != "payload_sha256"
        }
    )


def test_bundle_states_resolution_and_independence_limits(
    frozen_bundle: dict[str, Any],
) -> None:
    limitation = frozen_bundle["limitation"]
    assert "about 12 bits" in limitation
    assert "generic tail is 1/4096" in limitation
    assert "same Arb/FLINT backend" in limitation
    assert "not a clean-room" in limitation
    assert "UNRESOLVED" in limitation
    independence = frozen_bundle["independence"]
    assert independence["clean_room_or_alternate_backend"] is False
    assert independence["shared_arithmetic_backend"] == "python-flint/Arb"
    assert "no nyman formula import" in independence[
        "truncated_integral_implementation"
    ]


def test_generation_is_deterministic_across_precision_and_cache_perturbation(
    frozen_bundle: dict[str, Any],
) -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 512
        for index in range(1, 12):
            _ = (arb.pi() * index).sin() + arb(index + 1).log()
        regenerated = normalization.generate_nyman_normalization_bundle()
        assert ctx.prec == 512
        assert regenerated == frozen_bundle
    finally:
        ctx.prec = previous_precision


def test_exact_regeneration_verifier_accepts_mapping_and_path(
    frozen_bundle: dict[str, Any],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reproduced = normalization.verify_nyman_normalization_bundle(frozen_bundle)
    assert reproduced == {
        "classification": (
            "REPRODUCED_EXPLORATORY_NYMAN_NORMALIZATION_AUDIT"
        ),
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "NORMALIZATION_AUDIT_PASSED_EXPLORATORY",
        "verified_ratio_comparisons": "14",
        "harmonic_comparison_reproduced": True,
        "harmonic_comparison_passed": True,
        "payload_sha256": FROZEN_BUNDLE_SHA256,
        "exact_regeneration": True,
        "same_arb_backend_not_clean_room": True,
    }

    path = tmp_path / "normalization.json"
    path.write_text(
        json.dumps(frozen_bundle, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        normalization,
        "generate_nyman_normalization_bundle",
        lambda: copy.deepcopy(frozen_bundle),
    )
    assert normalization.verify_nyman_normalization_bundle(path) == reproduced


def _mutate_formula_result(bundle: dict[str, Any]) -> None:
    formula = bundle["formula_vs_truncated_integral"]
    encoded = formula["comparisons"][0]["candidate_enclosure"]["dyadic"]
    encoded["mid_mantissa"] = str(int(encoded["mid_mantissa"]) + 1)
    bundle["formula_vs_truncated_integral"] = _rehash(formula)


def _mutate_harmonic_result(bundle: dict[str, Any]) -> None:
    harmonic = bundle["one_minus_euler_gamma"]
    encoded = harmonic["core_enclosure"]["dyadic"]
    encoded["mid_mantissa"] = str(int(encoded["mid_mantissa"]) + 1)
    bundle["one_minus_euler_gamma"] = _rehash(harmonic)


def _mutate_classification(bundle: dict[str, Any]) -> None:
    bundle["classification"] = "PROVED"


@pytest.mark.parametrize(
    ("label", "mutate", "message"),
    [
        (
            "formula result",
            _mutate_formula_result,
            "does not exactly regenerate",
        ),
        (
            "harmonic result",
            _mutate_harmonic_result,
            "does not exactly regenerate",
        ),
        (
            "classification",
            _mutate_classification,
            "improperly promoted",
        ),
    ],
)
def test_rehashed_semantic_mutations_fail_closed(
    label: str,
    mutate: Callable[[dict[str, Any]], None],
    message: str,
    frozen_bundle: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutated = copy.deepcopy(frozen_bundle)
    mutate(mutated)
    mutated = _rehash(mutated)
    monkeypatch.setattr(
        normalization,
        "generate_nyman_normalization_bundle",
        lambda: copy.deepcopy(frozen_bundle),
    )
    with pytest.raises(
        normalization.NymanNormalizationVerificationError,
        match=message,
    ):
        normalization.verify_nyman_normalization_bundle(mutated)


def test_payload_hash_mutation_fails_before_regeneration(
    frozen_bundle: dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mutated = copy.deepcopy(frozen_bundle)
    mutated["payload_sha256"] = "0" * 64
    monkeypatch.setattr(
        normalization,
        "generate_nyman_normalization_bundle",
        lambda: (_ for _ in ()).throw(AssertionError("regenerated")),
    )
    with pytest.raises(
        normalization.NymanNormalizationVerificationError,
        match="payload hash mismatch",
    ):
        normalization.verify_nyman_normalization_bundle(mutated)


@pytest.mark.parametrize("duplicate_key", ["schema", "payload_sha256"])
def test_duplicate_json_keys_are_rejected_at_parse_boundary(
    duplicate_key: str,
    frozen_bundle: dict[str, Any],
    tmp_path: Path,
) -> None:
    canonical = json.dumps(
        frozen_bundle, ensure_ascii=False, separators=(",", ":"), sort_keys=True
    )
    path = tmp_path / f"duplicate-{duplicate_key}.json"
    path.write_text(
        "{" + json.dumps(duplicate_key) + ":\"bogus\"," + canonical[1:],
        encoding="utf-8",
    )
    with pytest.raises(
        normalization.NymanNormalizationVerificationError,
        match=f"duplicate JSON object key: {duplicate_key}",
    ):
        normalization.verify_nyman_normalization_bundle(path)


def test_nonstandard_json_constants_are_rejected(
    frozen_bundle: dict[str, Any],
    tmp_path: Path,
) -> None:
    canonical = json.dumps(
        frozen_bundle, ensure_ascii=False, separators=(",", ":"), sort_keys=True
    )
    path = tmp_path / "nan.json"
    path.write_text("{\"extra\":NaN," + canonical[1:], encoding="utf-8")
    with pytest.raises(
        normalization.NymanNormalizationVerificationError,
        match="nonstandard JSON constant: NaN",
    ):
        normalization.verify_nyman_normalization_bundle(path)


def test_clean_precision_restores_state_when_generation_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    previous_precision = ctx.prec
    try:
        ctx.prec = 333
        monkeypatch.setattr(
            normalization.oracle,
            "audit_frozen_ratio_normalization",
            lambda *args, **kwargs: (_ for _ in ()).throw(
                RuntimeError("injected generation failure")
            ),
        )
        with pytest.raises(RuntimeError, match="injected generation failure"):
            normalization.generate_nyman_normalization_bundle()
        assert ctx.prec == 333
    finally:
        ctx.prec = previous_precision


def test_invalid_bundle_inputs_fail_closed(tmp_path: Path) -> None:
    invalid = tmp_path / "invalid.json"
    invalid.write_text("not-json", encoding="utf-8")
    with pytest.raises(
        normalization.NymanNormalizationVerificationError,
        match="cannot read",
    ):
        normalization.verify_nyman_normalization_bundle(invalid)
    with pytest.raises(TypeError, match="mapping or Path"):
        normalization.verify_nyman_normalization_bundle([])  # type: ignore[arg-type]
