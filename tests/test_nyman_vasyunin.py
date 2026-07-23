from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import pytest

import riemann_lab.nyman_vasyunin as vasyunin


def _rehash(artifact: dict[str, Any]) -> None:
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    encoded = json.dumps(
        body,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    artifact["payload_sha256"] = hashlib.sha256(encoded).hexdigest()


def test_seed_is_the_exact_floor_difference_and_binary_parity() -> None:
    expected_h3 = [0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0]
    assert [
        vasyunin.vasyunin_seed(3, m) for m in range(1, 13)
    ] == expected_h3

    for n in range(1, 20):
        for m in range(1, 80):
            value = vasyunin.vasyunin_seed(n, m)
            assert value == m // n - 2 * (m // (2 * n))
            assert value == (m // n) % 2
            assert value in {0, 1}
        assert all(
            vasyunin.vasyunin_seed(n, prior_m) == 0
            for prior_m in range(1, n)
        )
        assert vasyunin.vasyunin_seed(n, n) == 1


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (1, 1),
        (2, -1),
        (3, -1),
        (4, 0),
        (5, -1),
        (6, 1),
        (10, 1),
        (15, 1),
        (30, -1),
        (49, 0),
        (210, 1),
    ],
)
def test_exact_mobius_values(n: int, expected: int) -> None:
    assert vasyunin.mobius(n) == expected


def test_greedy_recurrence_matches_two_adic_mobius_formula() -> None:
    coefficients = vasyunin.vasyunin_greedy_coefficients(256)
    assert coefficients[:16] == (
        1,
        1,
        -1,
        2,
        -1,
        -1,
        -1,
        4,
        0,
        -1,
        -1,
        -2,
        -1,
        -1,
        1,
        8,
    )
    for n, coefficient in enumerate(coefficients, start=1):
        recurrence = 1 - sum(
            coefficients[k - 1] * vasyunin.vasyunin_seed(k, n)
            for k in range(1, n)
        )
        assert coefficient == recurrence
        assert coefficient == vasyunin.vasyunin_formula_coefficient(n)


def test_every_prefix_interpolates_one_on_all_of_its_integer_intervals() -> None:
    coefficients = vasyunin.vasyunin_greedy_coefficients(96)
    for n in range(1, len(coefficients) + 1):
        prefix = coefficients[:n]
        assert all(
            vasyunin.vasyunin_phi(prefix, m) == 1
            for m in range(1, n + 1)
        )


def test_finite_audit_has_strict_schema_and_exact_power_of_two_factors() -> None:
    artifact = vasyunin.audit_vasyunin_greedy(32)
    assert set(artifact) == {
        "schema",
        "classification",
        "hypothesis_status",
        "audit_outcome",
        "finite_scope",
        "definitions",
        "coefficient_records",
        "power_of_two_increment_records",
        "checks",
        "external_theorem_bridge",
        "limitation",
        "payload_sha256",
    }
    assert artifact["schema"] == vasyunin.VASYUNIN_AUDIT_SCHEMA
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["finite_scope"]["limit"] == "32"
    assert artifact["finite_scope"]["integer_intervals"] == "[1,33)"
    assert len(artifact["coefficient_records"]) == 32
    assert all(artifact["checks"].values())

    powers = artifact["power_of_two_increment_records"]
    assert [record["n"] for record in powers] == ["2", "4", "8", "16", "32"]
    for record in powers:
        n = int(record["n"])
        assert int(record["coefficient_c_n"]) == n // 2
        assert record["l1_increment"] == {
            "expression": "(1/2)*log(2)",
            "symbolic_rational_multiplier": "1/2",
            "exact_rational_multiplier": {
                "numerator": "1",
                "denominator": "2",
            },
        }
        l2 = record["l2_squared_increment"]
        assert l2["expression"] == f"({n}/4)*log(2)"
        assert l2["symbolic_rational_multiplier"] == "n/4"
        assert l2["substituted_rational_multiplier"] == f"{n}/4"
        numerator = int(l2["exact_reduced_rational_multiplier"]["numerator"])
        denominator = int(
            l2["exact_reduced_rational_multiplier"]["denominator"]
        )
        assert numerator * 4 == n * denominator


def test_general_increment_records_use_reduced_exact_rationals() -> None:
    records = vasyunin.audit_vasyunin_greedy(12)["coefficient_records"]
    by_n = {int(record["n"]): record for record in records}
    assert by_n[3]["recurrence_coefficient"] == "-1"
    assert by_n[3]["l1_increment_log2_multiplier"] == {
        "numerator": "1",
        "denominator": "3",
    }
    assert by_n[3]["l2_squared_increment_log2_multiplier"] == {
        "numerator": "1",
        "denominator": "3",
    }
    assert by_n[4]["l1_increment_log2_multiplier"] == {
        "numerator": "1",
        "denominator": "2",
    }
    assert by_n[4]["l2_squared_increment_log2_multiplier"] == {
        "numerator": "1",
        "denominator": "1",
    }
    assert by_n[9]["recurrence_coefficient"] == "0"
    assert by_n[9]["l1_increment_log2_multiplier"] == {
        "numerator": "0",
        "denominator": "1",
    }


def test_published_bridge_does_not_promote_finite_enumeration() -> None:
    artifact = vasyunin.audit_vasyunin_greedy(16)
    bridge = artifact["external_theorem_bridge"]
    assert bridge["citation"]["primary_url"] == (
        "https://arxiv.org/abs/math/0506318"
    )
    scope = bridge["machine_scope"]
    assert scope == {
        "finite_prefix_recurrence_and_closed_form_checked": True,
        "finite_prefix_interval_interpolation_checked": True,
        "finite_prefix_rational_increment_multipliers_checked": True,
        "log2_integral_identity_reproved_by_finite_enumeration": False,
        "infinite_closed_form_reproved_by_finite_enumeration": False,
        "infinite_l1_divergence_reproved_by_finite_enumeration": False,
        "all_vasyunin_corrections_covered": False,
        "resolves_rh": False,
    }
    assert "finite prefix only" in artifact["limitation"]
    assert "UNRESOLVED" in artifact["limitation"]


def test_audit_exactly_regenerates_and_returns_a_fresh_tree() -> None:
    artifact = vasyunin.audit_vasyunin_greedy(24)
    assert vasyunin.verify_vasyunin_greedy_audit(artifact) == {
        "classification": (
            "REPRODUCED_CERTIFIED_FINITE_VASYUNIN_GREEDY_PREFIX"
        ),
        "hypothesis_status": "UNRESOLVED",
        "verified_limit": "24",
        "payload_sha256": artifact["payload_sha256"],
        "infinite_divergence_machine_reproved": False,
    }

    artifact["coefficient_records"][0]["recurrence_coefficient"] = "999"
    regenerated = vasyunin.audit_vasyunin_greedy(24)
    assert regenerated["coefficient_records"][0]["recurrence_coefficient"] == "1"


def _mutate_coefficient(artifact: dict[str, Any]) -> None:
    artifact["coefficient_records"][7]["recurrence_coefficient"] = "8"


def _mutate_infinite_claim(artifact: dict[str, Any]) -> None:
    artifact["external_theorem_bridge"]["machine_scope"][
        "infinite_l1_divergence_reproved_by_finite_enumeration"
    ] = True


def _mutate_power_factor(artifact: dict[str, Any]) -> None:
    artifact["power_of_two_increment_records"][2]["l1_increment"][
        "exact_rational_multiplier"
    ]["numerator"] = "2"


@pytest.mark.parametrize(
    "mutator",
    [_mutate_coefficient, _mutate_infinite_claim, _mutate_power_factor],
)
def test_rehashed_semantic_mutations_fail_exact_regeneration(
    mutator: Callable[[dict[str, Any]], None],
) -> None:
    artifact = copy.deepcopy(vasyunin.audit_vasyunin_greedy(16))
    mutator(artifact)
    _rehash(artifact)
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="does not exactly regenerate",
    ):
        vasyunin.verify_vasyunin_greedy_audit(artifact)


def test_unhashed_mutation_fails_payload_binding() -> None:
    artifact = vasyunin.audit_vasyunin_greedy(8)
    artifact["coefficient_records"][0]["recurrence_coefficient"] = "0"
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="payload hash mismatch",
    ):
        vasyunin.verify_vasyunin_greedy_audit(artifact)


def test_audit_detects_a_closed_form_implementation_mutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_formula = vasyunin.vasyunin_formula_coefficient

    def mutated_formula(n: int) -> int:
        value = real_formula(n)
        return value + 1 if n == 8 else value

    monkeypatch.setattr(
        vasyunin,
        "vasyunin_formula_coefficient",
        mutated_formula,
    )
    with pytest.raises(
        vasyunin.NymanVasyuninError,
        match="recurrence and closed form disagree at n=8",
    ):
        vasyunin.audit_vasyunin_greedy(8)


@pytest.mark.parametrize("bad", [True, False, 0, -1, 1.0, "8", None])
def test_positive_integer_inputs_are_strictly_validated(bad: Any) -> None:
    expected = (
        TypeError
        if isinstance(bad, (bool, float, str)) or bad is None
        else ValueError
    )
    with pytest.raises(expected):
        vasyunin.mobius(bad)
    with pytest.raises(expected):
        vasyunin.vasyunin_seed(bad, 1)
    with pytest.raises(expected):
        vasyunin.vasyunin_seed(1, bad)
    with pytest.raises(expected):
        vasyunin.vasyunin_formula_coefficient(bad)
    with pytest.raises(expected):
        vasyunin.vasyunin_greedy_coefficients(bad)
    with pytest.raises(expected):
        vasyunin.audit_vasyunin_greedy(bad)


def test_audit_limit_and_phi_coefficient_validation() -> None:
    with pytest.raises(ValueError, match="at most"):
        vasyunin.audit_vasyunin_greedy(vasyunin.MAX_AUDIT_LIMIT + 1)
    with pytest.raises(TypeError, match="nonempty integer sequence"):
        vasyunin.vasyunin_phi("1,1,-1", 3)
    with pytest.raises(ValueError, match="nonempty"):
        vasyunin.vasyunin_phi([], 1)
    with pytest.raises(TypeError, match="coefficient 2"):
        vasyunin.vasyunin_phi([1, True], 2)
    with pytest.raises(TypeError, match="coefficient 2"):
        vasyunin.vasyunin_phi([1, 1.0], 2)


def test_verifier_rejects_noncanonical_or_malformed_inputs() -> None:
    artifact = vasyunin.audit_vasyunin_greedy(8)

    with pytest.raises(TypeError, match="mapping"):
        vasyunin.verify_vasyunin_greedy_audit([])  # type: ignore[arg-type]

    malformed = copy.deepcopy(artifact)
    malformed["finite_scope"]["limit"] = "08"
    _rehash(malformed)
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="not canonical",
    ):
        vasyunin.verify_vasyunin_greedy_audit(malformed)

    non_json = copy.deepcopy(artifact)
    non_json["extra"] = object()
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="non-JSON value",
    ):
        vasyunin.verify_vasyunin_greedy_audit(non_json)

    floating = copy.deepcopy(artifact)
    floating["finite_scope"]["limit"] = 8.0
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="floating-point value is noncanonical",
    ):
        vasyunin.verify_vasyunin_greedy_audit(floating)


def test_path_verifier_rejects_duplicate_keys_and_nonstandard_constants(
    tmp_path: Path,
) -> None:
    artifact = vasyunin.audit_vasyunin_greedy(8)
    canonical_path = tmp_path / "canonical.json"
    canonical_path.write_text(
        json.dumps(artifact, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    assert vasyunin.verify_vasyunin_greedy_audit(canonical_path)[
        "verified_limit"
    ] == "8"

    duplicate_path = tmp_path / "duplicate.json"
    duplicate_path.write_text(
        '{"schema":"first","schema":"second"}',
        encoding="utf-8",
    )
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="duplicate JSON object key: schema",
    ):
        vasyunin.verify_vasyunin_greedy_audit(duplicate_path)

    nan_path = tmp_path / "nan.json"
    nan_path.write_text('{"value":NaN}', encoding="utf-8")
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="nonstandard JSON constant: NaN",
    ):
        vasyunin.verify_vasyunin_greedy_audit(nan_path)

    deeply_nested_path = tmp_path / "deeply-nested.json"
    deeply_nested_path.write_text(
        "[" * 10_000 + "0" + "]" * 10_000,
        encoding="utf-8",
    )
    with pytest.raises(
        vasyunin.NymanVasyuninVerificationError,
        match="cannot read or parse Vasyunin audit",
    ):
        vasyunin.verify_vasyunin_greedy_audit(deeply_nested_path)
