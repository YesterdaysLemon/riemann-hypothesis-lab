from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import pytest

import riemann_lab.nyman_alias as alias


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


def test_linear_mobius_sieve_matches_known_prefix() -> None:
    assert alias.mobius_sieve(16) == (
        0,
        1,
        -1,
        -1,
        0,
        -1,
        1,
        -1,
        0,
        0,
        1,
        -1,
        0,
        -1,
        1,
        1,
        0,
    )


def test_alias_convolution_and_divisor_inversion_are_exact() -> None:
    shell = {9: 1, 16: -1}
    coefficients = alias.alias_coefficients(shell, 64)
    assert coefficients[9] == 1
    assert coefficients[16] == -1
    assert coefficients[18] == -1
    assert coefficients[27] == -1
    assert coefficients[32] == 1
    assert coefficients[36] == 0

    inverted = alias.divisor_sum_prefix(coefficients)
    assert [
        (index, value)
        for index, value in enumerate(inverted)
        if value
    ] == [(9, 1), (16, -1)]


def test_finite_audit_has_exact_slope_records_and_honest_bridge() -> None:
    artifact = alias.audit_alias_sharp_truncation(64)
    assert artifact["schema"] == alias.ALIAS_AUDIT_SCHEMA
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["finite_scope"]["limit"] == "64"
    assert all(artifact["checks"].values())

    records = artifact["checkpoint_records"]
    assert [record["X"] for record in records] == ["16", "32", "64"]
    assert records[0]["S_X"] == {"numerator": "7", "denominator": "144"}
    assert records[0]["X_times_S_X_squared"] == {
        "numerator": "49",
        "denominator": "1296",
    }

    bridge = artifact["external_theorem_bridge"]
    assert bridge["primary_sources"][0]["primary_url"] == (
        "https://arxiv.org/abs/math/0504402"
    )
    assert bridge["primary_sources"][1]["primary_url"] == (
        "https://doi.org/10.1515/crll.1989.399.1"
    )
    scope = bridge["machine_scope"]
    assert scope["finite_mobius_convolution_checked"] is True
    assert scope["infinite_non_little_o_conclusion_machine_proved"] is False
    assert scope["scale_dependent_regularizations_covered"] is False
    assert scope["resolves_rh"] is False
    assert "UNRESOLVED" in artifact["limitation"]


def test_audit_exactly_regenerates_and_returns_a_fresh_tree() -> None:
    artifact = alias.audit_alias_sharp_truncation(128)
    assert alias.verify_alias_sharp_truncation_audit(artifact) == {
        "classification": "REPRODUCED_CERTIFIED_FINITE_ALIAS_IDENTITIES",
        "hypothesis_status": "UNRESOLVED",
        "verified_limit": "128",
        "payload_sha256": artifact["payload_sha256"],
        "infinite_sharp_truncation_obstruction_machine_reproved": False,
    }

    artifact["checkpoint_records"][0]["S_X"]["numerator"] = "999"
    regenerated = alias.audit_alias_sharp_truncation(128)
    assert regenerated["checkpoint_records"][0]["S_X"]["numerator"] == "7"


def _mutate_slope(artifact: dict[str, Any]) -> None:
    artifact["checkpoint_records"][1]["S_X"]["numerator"] = "0"


def _mutate_machine_scope(artifact: dict[str, Any]) -> None:
    artifact["external_theorem_bridge"]["machine_scope"][
        "infinite_non_little_o_conclusion_machine_proved"
    ] = True


def _mutate_limitation(artifact: dict[str, Any]) -> None:
    artifact["limitation"] = "RH proved"


@pytest.mark.parametrize(
    "mutator",
    [_mutate_slope, _mutate_machine_scope, _mutate_limitation],
)
def test_rehashed_semantic_mutations_fail_exact_regeneration(
    mutator: Callable[[dict[str, Any]], None],
) -> None:
    artifact = copy.deepcopy(alias.audit_alias_sharp_truncation(64))
    mutator(artifact)
    _rehash(artifact)
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="does not exactly regenerate",
    ):
        alias.verify_alias_sharp_truncation_audit(artifact)


def test_unhashed_mutation_fails_payload_binding() -> None:
    artifact = alias.audit_alias_sharp_truncation(64)
    artifact["checks"]["mobius_sieve_completed"] = False
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="payload hash mismatch",
    ):
        alias.verify_alias_sharp_truncation_audit(artifact)


@pytest.mark.parametrize("bad", [True, False, 0, -1, 1.0, "64", None])
def test_strict_integer_inputs_are_rejected(bad: Any) -> None:
    expected = (
        TypeError
        if isinstance(bad, (bool, float, str)) or bad is None
        else ValueError
    )
    with pytest.raises(expected):
        alias.mobius_sieve(bad)
    with pytest.raises(expected):
        alias.audit_alias_sharp_truncation(bad)


def test_shell_and_coefficient_validation() -> None:
    with pytest.raises(TypeError, match="mapping"):
        alias.alias_coefficients([(9, 1)], 16)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="nonzero"):
        alias.alias_coefficients({9: 0}, 16)
    with pytest.raises(TypeError, match="coefficients"):
        alias.alias_coefficients({9: True}, 16)
    with pytest.raises(ValueError, match="complete frozen shell"):
        alias.audit_alias_sharp_truncation(15)
    with pytest.raises(ValueError, match="at most"):
        alias.audit_alias_sharp_truncation(alias.MAX_AUDIT_LIMIT + 1)
    with pytest.raises(TypeError, match="integer sequence"):
        alias.divisor_sum_prefix("0,1,-1")
    with pytest.raises(ValueError, match="index zero"):
        alias.divisor_sum_prefix([1, 0])


def test_path_verifier_rejects_duplicate_keys_and_nonstandard_values(
    tmp_path: Path,
) -> None:
    artifact = alias.audit_alias_sharp_truncation(64)
    canonical_path = tmp_path / "canonical.json"
    canonical_path.write_text(
        json.dumps(artifact, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    assert alias.verify_alias_sharp_truncation_audit(canonical_path)[
        "verified_limit"
    ] == "64"

    duplicate_path = tmp_path / "duplicate.json"
    duplicate_path.write_text(
        '{"schema":"first","schema":"second"}', encoding="utf-8"
    )
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="duplicate JSON object key: schema",
    ):
        alias.verify_alias_sharp_truncation_audit(duplicate_path)

    nan_path = tmp_path / "nan.json"
    nan_path.write_text('{"value":NaN}', encoding="utf-8")
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="nonstandard JSON constant: NaN",
    ):
        alias.verify_alias_sharp_truncation_audit(nan_path)

    malformed = copy.deepcopy(artifact)
    malformed["finite_scope"]["limit"] = "064"
    _rehash(malformed)
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="finite limit is not canonical",
    ):
        alias.verify_alias_sharp_truncation_audit(malformed)

    deeply_nested_path = tmp_path / "deeply-nested.json"
    deeply_nested_path.write_text(
        '{"value":' + "[" * 1_000 + "0" + "]" * 1_000 + "}",
        encoding="utf-8",
    )
    with pytest.raises(
        alias.NymanAliasVerificationError,
        match="nesting is too deep|cannot read or parse alias audit",
    ):
        alias.verify_alias_sharp_truncation_audit(deeply_nested_path)
