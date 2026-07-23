"""Exploratory truncated-norm scout for balanced Nyman multipliers.

This script is deliberately not a certificate.  It reads a stored exact-dyadic
Nyman candidate, constructs the associated ideal first-shell step target, and
solves the finite constrained least-squares problem

    c_1 = 1,                 sum_(k <= K) c_k / k = 0

against the exact integer-interval error formula, truncated at ``M``.  The
truncated alias norm is monotone below the full Hilbert norm.  The truncated
direct gain also has a signed cross-term tail, so it has no one-sided meaning.
The floating-point solve and stored candidate make every reported value
EXPLORATORY.

NumPy is required only for this research scout; it is not a package runtime
dependency.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "rh-lab/nyman-balanced-multiplier-scout-cell/v1"
ENGINE = "binary64-divisor-interval-normal-equations/v1"


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _require_finite(value: Any, path: str = "cell") -> None:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ArithmeticError(f"{path} is not finite")
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _require_finite(item, f"{path}[{index}]")
        return
    if isinstance(value, dict):
        for key, item in value.items():
            _require_finite(item, f"{path}.{key}")


def _dyadic_float(record: dict[str, Any]) -> float:
    return math.ldexp(
        float(record["numerator"]),
        -int(record["denominator_exponent"]),
    )


def _candidate_data(path: Path, n: int) -> tuple[np.ndarray, float]:
    artifact: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    candidate = artifact["candidate"]
    record = candidate["coefficients"]
    exponent = int(record["denominator_exponent"])
    numerators = record["numerators"]
    if len(numerators) != n:
        raise ValueError(f"stored candidate dimension is not N={n}")
    scale = math.ldexp(1.0, -exponent)
    coefficients = np.asarray([float(value) * scale for value in numerators])
    lower = _dyadic_float(candidate["lower_bound"])
    upper = _dyadic_float(candidate["upper_bound"])
    return coefficients, (lower + upper) / 2.0


def ideal_shell(coefficients: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return shell indices, coefficients, and approximate local norm."""

    n = int(coefficients.size)
    source_indices = np.arange(1, n + 1, dtype=np.float64)
    slope = -float(np.sum(coefficients / source_indices))

    cumulative: dict[int, float] = {}
    local_norm = 0.0
    for m in range(n + 1, 2 * n):
        floors = np.floor(float(m) / source_indices)
        intercept = 1.0 + float(np.dot(coefficients, floors))
        logarithmic_mass = math.log1p(1.0 / m)
        interval_weight = 1.0 / (m * (m + 1.0))
        averaged_residual = (
            intercept + slope * logarithmic_mass / interval_weight
        )
        cumulative[m] = -averaged_residual
        local_norm += averaged_residual * averaged_residual * interval_weight

    indices = np.arange(n + 1, 2 * n + 1, dtype=np.int64)
    shell = np.empty(n, dtype=np.float64)
    shell[0] = cumulative[n + 1]
    for offset, m in enumerate(range(n + 2, 2 * n), start=1):
        shell[offset] = cumulative[m] - cumulative[m - 1]
    shell[-1] = -cumulative[2 * n - 1]
    return indices, shell, local_norm


def divisor_cumulative_columns(
    indices: np.ndarray,
    shell: np.ndarray,
    multiplier_limit: int,
    interval_limit: int,
) -> np.ndarray:
    """Build F_k(M)=sum_j y_j floor(M/(j*k)) for every M and k."""

    columns = np.empty((interval_limit, multiplier_limit), dtype=np.float64)
    for k in range(1, multiplier_limit + 1):
        jumps = np.zeros(interval_limit + 1, dtype=np.float64)
        for j, value in zip(indices, shell, strict=True):
            period = int(j) * k
            if period <= interval_limit:
                jumps[period::period] += value
        columns[:, k - 1] = np.cumsum(jumps)[1:]
    return columns


def target_cumulative(
    indices: np.ndarray,
    shell: np.ndarray,
    interval_limit: int,
) -> np.ndarray:
    jumps = np.zeros(interval_limit + 1, dtype=np.float64)
    jumps[indices] = shell
    return np.cumsum(jumps)[1:]


def solve_balanced(
    columns: np.ndarray,
    target: np.ndarray,
    residual_interval_masses: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve with c1=1 after eliminating c2 by harmonic balance.

    With no interval masses, minimize the truncated alias norm.  Otherwise
    maximize the signed direct gain against the supplied old residual.
    """

    interval_limit, multiplier_limit = columns.shape
    if multiplier_limit < 2:
        raise ValueError("a balanced multiplier needs K >= 2")
    m = np.arange(1, interval_limit + 1, dtype=np.float64)
    weight = 1.0 / (m * (m + 1.0))
    base = target - columns[:, 0] + 2.0 * columns[:, 1]

    coefficients = np.zeros(multiplier_limit, dtype=np.float64)
    coefficients[0] = 1.0
    coefficients[1] = -2.0
    if (
        residual_interval_masses is not None
        and residual_interval_masses.shape != target.shape
    ):
        raise ValueError("old-residual interval masses changed dimension")
    if multiplier_limit == 2:
        return coefficients, base, 1.0

    variable_indices = np.arange(3, multiplier_limit + 1, dtype=np.float64)
    # Reuse the column allocation and accumulate normal equations in blocks.
    # This keeps memory linear in M*K and makes K proportional to N practical
    # for the scout.  It is intentionally a floating-point diagnostic, not a
    # certificate path.
    effective = columns[:, 2:]
    effective -= (2.0 / variable_indices)[None, :] * columns[:, 1, None]
    variable_count = multiplier_limit - 2
    normal = np.zeros((variable_count, variable_count), dtype=np.float64)
    right_hand_side = np.zeros(variable_count, dtype=np.float64)
    block_size = 100_000
    for first in range(0, interval_limit, block_size):
        last = min(first + block_size, interval_limit)
        block = effective[first:last]
        block_weight = weight[first:last]
        normal += block.T @ (block * block_weight[:, None])
        right = base[first:last] * block_weight
        if residual_interval_masses is not None:
            right = right - residual_interval_masses[first:last]
        right_hand_side += block.T @ right
    eigenvalues = np.linalg.eigvalsh(normal)
    if eigenvalues[0] <= 0.0:
        raise ArithmeticError("truncated normal matrix is not positive definite")
    solution = np.linalg.solve(normal, right_hand_side)
    coefficients[2:] = solution
    coefficients[1] -= float(np.sum(2.0 * solution / variable_indices))
    residual = base - effective @ solution
    condition = float(math.sqrt(eigenvalues[-1] / eigenvalues[0]))
    return coefficients, residual, condition


def old_residual_interval_masses(
    coefficients: np.ndarray,
    interval_limit: int,
) -> np.ndarray:
    """Return integral of the old residual on each weighted unit interval."""

    jumps = np.zeros(interval_limit + 1, dtype=np.float64)
    for index, value in enumerate(coefficients, start=1):
        jumps[index::index] += value
    intercept = 1.0 + np.cumsum(jumps)[1:]
    slope = -float(
        np.sum(
            coefficients
            / np.arange(1, coefficients.size + 1, dtype=np.float64)
        )
    )
    m = np.arange(1, interval_limit + 1, dtype=np.float64)
    return slope * np.log1p(1.0 / m) + intercept / (m * (m + 1.0))


def scout(
    n: int,
    k: int,
    interval_limit: int,
    root: Path,
    objective: str,
) -> dict[str, Any]:
    if n < 1:
        raise ValueError("N must be positive")
    if k < 2:
        raise ValueError("K must be at least 2")
    if interval_limit < 2 * n:
        raise ValueError("the interval limit must include the complete shell")
    if objective not in {"alias-norm", "direct-gain"}:
        raise ValueError("objective must be 'alias-norm' or 'direct-gain'")
    candidate_path = root / "results" / "nyman-natural-v1" / "cells" / f"n-{n:04d}.json"
    candidate_artifact: dict[str, Any] = json.loads(
        candidate_path.read_text(encoding="utf-8")
    )
    candidate, old_energy_midpoint = _candidate_data(candidate_path, n)
    indices, shell, local_norm = ideal_shell(candidate)
    columns = divisor_cumulative_columns(indices, shell, k, interval_limit)
    target = target_cumulative(indices, shell, interval_limit)
    interval_masses = (
        old_residual_interval_masses(candidate, interval_limit)
        if objective == "direct-gain"
        else None
    )
    multiplier, residual, condition = solve_balanced(
        columns,
        target,
        interval_masses,
    )

    m = np.arange(1, interval_limit + 1, dtype=np.float64)
    summands = residual * residual / (m * (m + 1.0))
    cumulative = np.cumsum(summands)
    direct_gain = (
        local_norm
        + 2.0 * np.cumsum(residual * interval_masses)
        - cumulative
        if interval_masses is not None
        else None
    )
    checkpoints = sorted(
        {
            max(2 * n, interval_limit // 8),
            max(2 * n, interval_limit // 4),
            max(2 * n, interval_limit // 2),
            interval_limit,
        }
    )
    harmonic_balance = float(
        np.sum(multiplier / np.arange(1, k + 1, dtype=np.float64))
    )
    body = {
        "schema": SCHEMA,
        "engine": ENGINE,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "source_binding": {
            "path": candidate_path.relative_to(root).as_posix(),
            "raw_sha256": hashlib.sha256(candidate_path.read_bytes()).hexdigest(),
            "candidate_sha256": candidate_artifact["candidate"][
                "candidate_sha256"
            ],
            "coefficient_denominator_exponent": candidate_artifact["candidate"][
                "coefficients"
            ]["denominator_exponent"],
        },
        "n": n,
        "multiplier_limit": k,
        "interval_limit": interval_limit,
        "objective": objective,
        "old_candidate_energy_bracket_midpoint": old_energy_midpoint,
        "local_gain_fraction_of_old_energy_bracket_midpoint": (
            local_norm / old_energy_midpoint
        ),
        "local_norm": local_norm,
        "shell_sum": float(np.sum(shell)),
        "shell_harmonic_sum": float(np.sum(shell / indices)),
        "shell_critical_line_zero_height_value": float(
            np.sum(shell / np.sqrt(indices))
        ),
        "shell_l1": float(np.sum(np.abs(shell))),
        "harmonic_balance_residual": harmonic_balance,
        "least_squares_condition": condition,
        "coefficient_l1": float(np.sum(np.abs(multiplier))),
        "coefficient_linf": float(np.max(np.abs(multiplier))),
        "coefficients": multiplier.tolist(),
        "partial_norms": [
            {
                "last_interval": checkpoint,
                "squared_norm": float(cumulative[checkpoint - 1]),
                "ratio_to_local_norm": float(
                    cumulative[checkpoint - 1] / local_norm
                ),
                **(
                    {
                        "direct_gain": float(direct_gain[checkpoint - 1]),
                        "direct_gain_fraction_of_old_energy_bracket_midpoint": float(
                            direct_gain[checkpoint - 1] / old_energy_midpoint
                        ),
                        "direct_gain_fraction_of_local_gain": float(
                            direct_gain[checkpoint - 1] / local_norm
                        ),
                    }
                    if direct_gain is not None
                    else {}
                ),
            }
            for checkpoint in checkpoints
        ],
        "limitation": (
            "The ideal shell comes from one stored exact-dyadic candidate, "
            "energy fractions use the midpoint of its certified bracket, "
            "and the solve uses binary64 NumPy. The alias norm omits a "
            "positive tail beyond the declared interval limit; the displayed "
            "direct gain also omits a signed cross-term tail and is therefore "
            "neither a lower nor an upper bound. This is a convergence "
            "diagnostic, not a certificate and not evidence resolving RH."
        ),
    }
    _require_finite(body)
    return {**body, "payload_sha256": _canonical_sha256(body)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--k", type=int, required=True)
    parser.add_argument("--interval-limit", type=int, required=True)
    parser.add_argument(
        "--objective",
        choices=("alias-norm", "direct-gain"),
        default="alias-norm",
    )
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    result = scout(
        args.n,
        args.k,
        args.interval_limit,
        args.root,
        args.objective,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
