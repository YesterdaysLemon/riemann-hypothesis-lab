from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from tools import generate_nyman_balanced_scout
from tools.scout_nyman_balanced import (
    _candidate_data,
    _require_finite,
    divisor_cumulative_columns,
    ideal_shell,
    old_residual_interval_masses,
    scout,
    solve_balanced,
    target_cumulative,
)


ROOT = Path(__file__).resolve().parents[1]


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def test_ideal_shell_reconstructs_documented_n8_values() -> None:
    path = ROOT / "results/nyman-natural-v1/cells/n-0008.json"
    candidate, old_energy = _candidate_data(path, 8)
    indices, shell, local_norm = ideal_shell(candidate)

    assert indices.tolist() == list(range(9, 17))
    assert float(np.sum(shell)) == pytest.approx(0.0, abs=1e-15)
    assert float(np.sum(np.abs(shell))) == pytest.approx(
        4.733961161269381, rel=2e-15
    )
    assert float(np.sum(shell / indices)) == pytest.approx(
        -0.002466416276690815, rel=2e-14
    )
    assert local_norm == pytest.approx(0.00828945337577692, rel=2e-15)
    assert old_energy == pytest.approx(0.024161421585896685, rel=2e-15)


def test_divisor_columns_match_direct_floor_definition() -> None:
    indices = np.asarray([2, 3], dtype=np.int64)
    shell = np.asarray([1.25, -1.25], dtype=np.float64)
    interval_limit = 64
    multiplier_limit = 7
    columns = divisor_cumulative_columns(
        indices, shell, multiplier_limit, interval_limit
    )

    for m in range(1, interval_limit + 1):
        for k in range(1, multiplier_limit + 1):
            direct = sum(
                value * (m // (int(index) * k))
                for index, value in zip(indices, shell, strict=True)
            )
            assert columns[m - 1, k - 1] == pytest.approx(direct, abs=1e-14)


def test_scout_rejects_unknown_objective() -> None:
    with pytest.raises(ValueError, match="objective"):
        scout(8, 4, 64, ROOT, "unknown")


def test_scout_rejects_nonfinite_output() -> None:
    with pytest.raises(ArithmeticError, match="cell.value"):
        _require_finite({"value": float("nan")})


@pytest.mark.parametrize("direct_gain", [False, True])
def test_balanced_solve_satisfies_constraints_and_interval_identity(
    direct_gain: bool,
) -> None:
    path = ROOT / "results/nyman-natural-v1/cells/n-0008.json"
    candidate, _ = _candidate_data(path, 8)
    indices, shell, _ = ideal_shell(candidate)
    interval_limit = 4096
    columns = divisor_cumulative_columns(indices, shell, 8, interval_limit)
    original_columns = columns.copy()
    target = target_cumulative(indices, shell, interval_limit)
    masses = (
        old_residual_interval_masses(candidate, interval_limit)
        if direct_gain
        else None
    )

    coefficients, residual, condition = solve_balanced(
        columns, target, masses
    )

    assert coefficients[0] == 1.0
    assert float(
        np.sum(coefficients / np.arange(1, 9, dtype=np.float64))
    ) == pytest.approx(0.0, abs=3e-15)
    assert residual == pytest.approx(
        target - original_columns @ coefficients,
        abs=2e-13,
    )
    assert condition > 1.0


def test_frozen_exploratory_artifact_is_hash_bound_and_honest() -> None:
    path = ROOT / "results/nyman-balanced-multiplier-scout-v1.json"
    artifact = json.loads(path.read_text(encoding="utf-8"))
    payload_hash = artifact.pop("payload_sha256")

    assert payload_hash == (
        "f31c330b373b3e8616fb1837193d8fbac6e31ce56dd4986364342dc4369f3144"
    )
    assert _canonical_sha256(artifact) == payload_hash
    assert artifact["classification"] == "EXPLORATORY"
    assert artifact["hypothesis_status"] == "UNRESOLVED"
    assert generate_nyman_balanced_scout.SCHEMA == artifact["schema"]
    for binding in artifact["source_bindings"]:
        source = ROOT / binding["path"]
        source_artifact = json.loads(source.read_text(encoding="utf-8"))
        assert hashlib.sha256(source.read_bytes()).hexdigest() == binding[
            "raw_sha256"
        ]
        assert source_artifact["candidate"]["candidate_sha256"] == binding[
            "candidate_sha256"
        ]
    assert artifact["checks"] == {
        "all_cells_keep_rh_unresolved": True,
        "all_cells_remain_exploratory": True,
        "all_harmonic_balance_residuals_below_2e_15": True,
        "exact_dyadic_large_scale_scout_vectors_stored": False,
        "fixed_k64_truncated_direct_gain_positive_on_grid": True,
        "full_infinite_alias_tail_certified": False,
        "full_signed_direct_gain_tail_certified": False,
        "resolves_rh": False,
        "uniform_all_scale_estimate_proved": False,
    }
    assert "neither a lower nor an upper bound" in artifact["limitation"]

    fixed = artifact["fixed_width_grid"]["cells"]
    all_cells = [
        *fixed,
        *artifact["n256_width_sweep"]["cells"],
        artifact["n256_k64_convergence"],
    ]
    root_bindings = {
        int(binding["n"]): binding for binding in artifact["source_bindings"]
    }
    for cell in all_cells:
        cell_body = {
            key: value for key, value in cell.items() if key != "payload_sha256"
        }
        assert cell["schema"] == (
            "rh-lab/nyman-balanced-multiplier-scout-cell/v1"
        )
        assert _canonical_sha256(cell_body) == cell["payload_sha256"]
        source_binding = cell["source_binding"]
        root_binding = root_bindings[cell["n"]]
        assert source_binding["path"] == root_binding["path"]
        assert source_binding["raw_sha256"] == root_binding["raw_sha256"]
        assert source_binding["candidate_sha256"] == root_binding[
            "candidate_sha256"
        ]
    assert [cell["n"] for cell in fixed] == [8, 16, 32, 64, 128, 256]
    assert all(
        cell["partial_norms"][-1][
            "direct_gain_fraction_of_old_energy_bracket_midpoint"
        ]
        > 0.0
        for cell in fixed
    )
    assert artifact["n256_k64_convergence"]["partial_norms"][-1][
        "direct_gain_fraction_of_old_energy_bracket_midpoint"
    ] == pytest.approx(0.07075747744826753, rel=2e-15)
