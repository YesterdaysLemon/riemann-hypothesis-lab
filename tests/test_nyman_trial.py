from __future__ import annotations

import copy
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb
import pytest

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
from riemann_lab import nyman as core
import riemann_lab.nyman_trial as trial


def _rehash(record: dict[str, Any]) -> dict[str, Any]:
    body = {
        key: copy.deepcopy(value)
        for key, value in record.items()
        if key != "payload_sha256"
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _exact(value: Any, integer: int) -> bool:
    return bool(value.is_exact() and value == integer)


def test_canonical_basis_order_and_zero_policies() -> None:
    with trial.beta2._clean_precision(192):
        basis = trial._build_trial_basis()
        manifest = trial._basis_manifest(
            basis,
            old_n=trial.FROZEN_OLD_N,
            new_n=trial.FROZEN_NEW_N,
        )

    assert trial.TRIAL_COLUMN_FORMULAS == (
        "mu(j)*log(512/j)",
        "mu(j)",
        "mu(j)*log(512/j)^2",
        "abs(mu(j))*log(512/j)",
        "abs(mu(j))",
        "1-abs(mu(j))",
        "1_{2|j}",
        "1_{3|j}",
    )
    assert len(basis) == 256
    assert all(len(row) == 8 for row in basis)

    # 257 is prime, so mu(257)=-1.
    row_257 = basis[0]
    assert row_257[0] < 0
    assert _exact(row_257[1], -1)
    assert row_257[2] < 0
    assert row_257[3] > 0
    assert _exact(row_257[4], 1)
    assert _exact(row_257[5], 0)
    assert _exact(row_257[6], 0)
    assert _exact(row_257[7], 0)

    # 260 contains 2^2.  Its Mobius value is zero, and the log-bearing
    # columns are formed as exact zero without evaluating 0*log.
    row_260 = basis[260 - 257]
    assert all(_exact(row_260[index], 0) for index in (0, 1, 2, 3, 4))
    assert _exact(row_260[5], 1)
    assert _exact(row_260[6], 1)
    assert _exact(row_260[7], 0)

    row_512 = basis[-1]
    assert all(_exact(row_512[index], 0) for index in (0, 1, 2, 3, 4))
    assert _exact(row_512[5], 1)
    assert _exact(row_512[6], 1)
    assert _exact(row_512[7], 0)
    assert manifest["checks"] == {
        "column_order_frozen": True,
        "endpoint_j_equals_512_log_columns_are_exact_zero": True,
        "zero_mobius_rows_skip_log_multiplication": True,
        "binary_float_used": False,
    }
    assert len(manifest["coefficient_matrix_content_sha256"]) == 64


def test_aggregate_construction_and_augmented_ldlt_on_identity_fixture() -> None:
    identity = tuple(
        tuple(arb(1 if row == column else 0) for column in range(4))
        for row in range(4)
    )
    system = core.NaturalSystem(
        (1, 2, 3, 4),
        identity,
        (arb(0), arb(0), arb(0), arb(0)),
        precision_bits=192,
    )
    basis = (
        (arb(1), arb(0)),
        (arb(0), arb(1)),
    )
    aggregate = trial._build_aggregate_system(system, basis, old_n=2)

    assert aggregate.dimension == 4
    assert aggregate.old_dimension == 2
    assert aggregate.trial_dimension == 2
    assert all(
        _exact(aggregate.gram[row][column], 1 if row == column else 0)
        for row in range(4)
        for column in range(4)
    )
    assert all(_exact(value, 0) for value in aggregate.target)

    certificate = trial._certify_aggregate_lower_bound(
        aggregate,
        Fraction(1, 2),
    )
    assert certificate["decision"] == "LOWER_BOUND_CERTIFIED"
    assert certificate["matrix_order"] == (
        "2 old aggregate columns; 2 trial columns in supplied order; "
        "target Schur coordinate"
    )
    assert certificate["strict_statement"] == (
        "restricted squared distance > 1/2"
    )
    assert certificate["checks"] == {
        "fixed_order_interval_ldlt_positive": True,
        "aggregate_gram_positive_as_leading_principal_submatrix": True,
        "trial_columns_independent_modulo_old_span": True,
        "approximate_solve_used_as_evidence": False,
    }
    assert len(certificate["augmented_ldlt"]["pivots"]) == 5


def test_small_end_to_end_trial_evidence_builds_augmented_ldlt(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Exercise the real natural-kernel -> aggregate -> augmented-LDL path on
    # ten dilates.  The eight test columns are an exact identity on j=3..10,
    # so this fixture tests construction rather than the canonical ansatz.
    monkeypatch.setattr(trial, "FROZEN_OLD_N", 2)
    monkeypatch.setattr(trial, "FROZEN_NEW_N", 10)
    monkeypatch.setattr(
        trial,
        "FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND",
        Fraction(1, 1_000_000),
    )

    def identity_new_block(
        old_n: int = 2,
        new_n: int = 10,
    ) -> tuple[tuple[Any, ...], ...]:
        assert (old_n, new_n) == (2, 10)
        return tuple(
            tuple(arb(1 if row == column else 0) for column in range(8))
            for row in range(8)
        )

    monkeypatch.setattr(trial, "_build_trial_basis", identity_new_block)
    evidence = trial._build_trial_evidence(192)

    assert evidence["aggregate_system"] == {
        "content_sha256": evidence["aggregate_system"]["content_sha256"],
        "dimension": "10",
        "old_dimension": "2",
        "trial_dimension": "8",
        "ordering": "rho_1,...,rho_256 followed by the eight manifest columns",
    }
    certificate = evidence["restricted_distance_lower_certificate"]
    assert certificate["decision"] == "LOWER_BOUND_CERTIFIED"
    assert certificate["augmented_ldlt"]["classification"] == "POSITIVE_DEFINITE"
    assert len(certificate["augmented_ldlt"]["pivots"]) == 11


def test_exact_threshold_implies_strict_gain_ceiling() -> None:
    thresholds = trial._exact_threshold_record()
    lower = trial._source_fraction(
        thresholds["restricted_distance_lower_bound"],
        "restricted distance threshold",
    )
    source_upper = trial._source_fraction(
        thresholds["source_upper_bound_u256"],
        "source U256",
    )
    assert lower == Fraction(9, 10) * source_upper
    assert lower == trial.FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND
    assert trial.GAIN_CEILING_FACTOR == Fraction(1, 10)

    # For every admissible d^2<=U and every certified F>L, subtraction gives
    # d^2-F < d^2-(9/10)U <= d^2/10.
    for d_squared in (source_upper, source_upper / 2, source_upper / 10):
        assert d_squared - lower <= d_squared / 10
    bridge = trial._theorem_bridge()
    assert bridge["all_scale_trial_subspace_lemma_proved"] is False
    assert bridge["all_trial_subspaces_rejected"] is False
    assert bridge["resolves_rh"] is False


def _fake_lower_certificate(
    aggregate_sha256: str,
) -> dict[str, Any]:
    dimension = trial.FROZEN_OLD_N + len(trial.TRIAL_COLUMN_IDS)
    body = {
        "schema": trial.LOWER_CERTIFICATE_SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(trial.FROZEN_GENERATION_BITS),
        "normalization": core.NYMAN_NORMALIZATION,
        "aggregate_system_content_sha256": aggregate_sha256,
        "aggregate_dimension": str(dimension),
        "old_dimension": str(trial.FROZEN_OLD_N),
        "trial_dimension": str(len(trial.TRIAL_COLUMN_IDS)),
        "claimed_lower_bound": trial._fraction_record(
            trial.FROZEN_RESTRICTED_DISTANCE_LOWER_BOUND
        ),
        "matrix_order": (
            "rho_1,...,rho_256; eight frozen trial columns in manifest order; "
            "target Schur coordinate"
        ),
        "strict_statement": "F_V > (9/10)*U_256",
        "decision": "LOWER_BOUND_CERTIFIED",
        "augmented_ldlt": {
            "classification": "POSITIVE_DEFINITE",
            "fixed_order": [str(index) for index in range(dimension + 1)],
            "failed_pivot_index": None,
            "pivots": [
                {"index": str(index), "value": arb_record(arb(1))}
                for index in range(dimension + 1)
            ],
        },
        "checks": {
            "fixed_order_interval_ldlt_positive": True,
            "aggregate_gram_positive_as_leading_principal_submatrix": True,
            "trial_columns_independent_modulo_old_span": True,
            "approximate_solve_used_as_evidence": False,
        },
        "limitation": (
            "A strict finite lower bound for one aggregate subspace is not an "
            "all-scale Nyman result and does not resolve RH."
        ),
    }
    return trial._with_payload_hash(body)


def _fake_evidence() -> dict[str, Any]:
    basis_sha256 = "a" * 64
    aggregate_sha256 = "b" * 64
    mu = trial._mobius_sieve(trial.FROZEN_NEW_N)
    return {
        "precision_bits": str(trial.FROZEN_GENERATION_BITS),
        "backend": trial._backend_record(),
        "max_512_kernel": {
            "core_system_content_sha256": (
                trial.beta2.FROZEN_GENERATION_KERNEL_SHA256
            ),
            "prefix_256_numeric_sha256": (
                trial.beta2.FROZEN_N256_PREFIX_KERNEL_SHA256
            ),
            "prefix_512_numeric_sha256": (
                trial.beta2.FROZEN_GENERATION_PREFIX_512_SHA256
            ),
            "prefix_nesting_check": True,
        },
        "trial_basis": {
            "new_block": {
                "first_j": str(trial.FROZEN_OLD_N + 1),
                "last_j": str(trial.FROZEN_NEW_N),
                "row_count": str(trial.FROZEN_NEW_N - trial.FROZEN_OLD_N),
            },
            "ordered_columns": [
                {
                    "index": str(index),
                    "column_id": column_id,
                    "formula": formula,
                }
                for index, (column_id, formula) in enumerate(
                    zip(
                        trial.TRIAL_COLUMN_IDS,
                        trial.TRIAL_COLUMN_FORMULAS,
                        strict=True,
                    )
                )
            ],
            "coefficient_encoding": (
                "Arb exact-dyadic enclosure per coefficient"
            ),
            "coefficient_matrix_content_sha256": basis_sha256,
            "mobius_values_content_sha256": content_sha256(
                {
                    "first_j": str(trial.FROZEN_OLD_N + 1),
                    "last_j": str(trial.FROZEN_NEW_N),
                    "values": [
                        str(mu[j])
                        for j in range(
                            trial.FROZEN_OLD_N + 1,
                            trial.FROZEN_NEW_N + 1,
                        )
                    ],
                }
            ),
            "checks": {
                "column_order_frozen": True,
                "endpoint_j_equals_512_log_columns_are_exact_zero": True,
                "zero_mobius_rows_skip_log_multiplication": True,
                "binary_float_used": False,
            },
        },
        "aggregate_system": {
            "content_sha256": aggregate_sha256,
            "dimension": str(trial.FROZEN_OLD_N + len(trial.TRIAL_COLUMN_IDS)),
            "old_dimension": str(trial.FROZEN_OLD_N),
            "trial_dimension": str(len(trial.TRIAL_COLUMN_IDS)),
            "ordering": (
                "rho_1,...,rho_256 followed by the eight manifest columns"
            ),
        },
        "restricted_distance_lower_certificate": _fake_lower_certificate(
            aggregate_sha256
        ),
    }


def _fake_artifact(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    monkeypatch.setattr(trial, "FROZEN_AUDIT_PAYLOAD_SHA256", "")
    monkeypatch.setattr(trial, "FROZEN_BASIS_COEFFICIENTS_SHA256", "")
    monkeypatch.setattr(trial, "FROZEN_AGGREGATE_SYSTEM_SHA256", "")
    monkeypatch.setattr(
        trial,
        "FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256",
        "",
    )
    verification = {
        "classification": "test-only-source-verification",
    }
    source_cell = {
        "terminal_statement": "test-only frozen source statement",
    }
    monkeypatch.setattr(
        trial,
        "_validate_source_summary",
        lambda summary, checkpoint_dir: (verification, source_cell),
    )
    monkeypatch.setattr(
        trial,
        "_build_trial_evidence",
        lambda precision_bits: _fake_evidence(),
    )
    return trial.generate_nyman_trial_subspace_audit(
        {"test": "source"},
        Path("test-checkpoint"),
    )


def test_generation_contract_with_monkeypatched_small_evidence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _fake_artifact(monkeypatch)

    assert artifact["schema"] == trial.AUDIT_SCHEMA
    assert artifact["audit_id"] == trial.FROZEN_AUDIT_ID
    assert artifact["classification"] == "CERTIFIED_FINITE"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert artifact["audit_outcome"] == trial.AUDIT_OUTCOME
    assert artifact["terminal_statement"] == trial.TERMINAL_STATEMENT
    assert artifact["payload_sha256"] == content_sha256(
        {
            key: value
            for key, value in artifact.items()
            if key != "payload_sha256"
        }
    )
    trial._preflight_audit(artifact)


def test_preflight_rejects_rehashed_mutation_before_regeneration(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _fake_artifact(monkeypatch)
    artifact["terminal_statement"] = "forged"
    artifact = _rehash(artifact)
    called = False

    def forbidden_derive(*args: object, **kwargs: object) -> object:
        nonlocal called
        called = True
        raise AssertionError("expensive regeneration must not run")

    monkeypatch.setattr(trial, "_derive_audit", forbidden_derive)
    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="audit terminal_statement changed",
    ):
        trial.verify_nyman_trial_subspace_audit(
            artifact,
            {"test": "source"},
            Path("test-checkpoint"),
        )
    assert called is False


def test_preflight_rejects_nonpositive_rehashed_pivot(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    artifact = _fake_artifact(monkeypatch)
    certificate = artifact["generation"]["restricted_distance_lower_certificate"]
    certificate["augmented_ldlt"]["pivots"][-1]["value"] = arb_record(arb(-1))
    artifact["generation"]["restricted_distance_lower_certificate"] = _rehash(
        certificate
    )
    artifact = _rehash(artifact)

    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="stores a nonpositive pivot",
    ):
        trial._preflight_audit(artifact)


def test_strict_loader_rejects_duplicate_keys_and_float(tmp_path: Path) -> None:
    duplicate = tmp_path / "duplicate.json"
    duplicate.write_text(
        '{"schema":"first","schema":"second"}',
        encoding="utf-8",
    )
    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="duplicate JSON object key: schema",
    ):
        trial._load_object(
            duplicate,
            trial.NymanTrialVerificationError,
            label="trial artifact",
        )


def test_loader_and_hash_check_wrap_pathological_mapping_errors() -> None:
    cyclic: dict[str, Any] = {}
    cyclic["self"] = cyclic
    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="is not a finite canonical JSON tree",
    ):
        trial._load_object(
            cyclic,
            trial.NymanTrialVerificationError,
            label="trial artifact",
        )

    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="is not canonical JSON",
    ):
        trial._require_payload_hash(
            {
                "oversized": 10**20_000,
                "payload_sha256": "0" * 64,
            },
            trial.NymanTrialVerificationError,
            label="trial artifact",
        )

    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="floating-point JSON number is noncanonical",
    ):
        trial._load_object(
            {"schema": trial.AUDIT_SCHEMA, "forged": 0.9},
            trial.NymanTrialVerificationError,
            label="trial artifact",
        )


@pytest.mark.parametrize("constant", ["NaN", "Infinity", "-Infinity"])
def test_strict_loader_rejects_nonstandard_constants(
    constant: str,
    tmp_path: Path,
) -> None:
    path = tmp_path / "constant.json"
    path.write_text('{"forged":' + constant + "}", encoding="utf-8")
    with pytest.raises(
        trial.NymanTrialVerificationError,
        match=f"nonstandard JSON constant: {constant}",
    ):
        trial._load_object(
            path,
            trial.NymanTrialVerificationError,
            label="trial artifact",
        )


@pytest.mark.parametrize("bits", [True, 0, 1535, 1537])
def test_verifier_requires_exact_replay_precision_before_io(bits: object) -> None:
    with pytest.raises(
        trial.NymanTrialVerificationError,
        match="replay precision must equal the frozen 1536 bits",
    ):
        trial.verify_nyman_trial_subspace_audit(
            Path("missing-artifact.json"),
            Path("missing-summary.json"),
            Path("missing-checkpoint"),
            replay_precision_bits=bits,  # type: ignore[arg-type]
        )


def test_public_api_and_canonical_hashes_are_frozen() -> None:
    assert callable(trial.generate_nyman_trial_subspace_audit)
    assert callable(trial.verify_nyman_trial_subspace_audit)
    assert "generate_nyman_trial_subspace_audit" in trial.__all__
    assert "verify_nyman_trial_subspace_audit" in trial.__all__
    assert trial.FROZEN_AUDIT_PAYLOAD_SHA256 == (
        "467d6d819700a87f656e917bcb3e63008d7243fa5f412f75b6eec4ade5a67956"
    )
    assert trial.FROZEN_BASIS_COEFFICIENTS_SHA256 == (
        "ec9ea0b1a3cb7439f2e30c552c38e2e5a4620bf5a5ec1285d68901a2a29a7f99"
    )
    assert trial.FROZEN_AGGREGATE_SYSTEM_SHA256 == (
        "dd9a571997cc10a32fe332ad2127d58913c021054f058aee90253fdb85422fec"
    )
    assert trial.FROZEN_LOWER_CERTIFICATE_PAYLOAD_SHA256 == (
        "a58186360f701a26fede8872e463e28d57f4e2fd6f06a001137495216b29eea0"
    )
