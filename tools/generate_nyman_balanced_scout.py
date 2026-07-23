"""Generate the frozen exploratory balanced-multiplier scout grid."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np

if __package__:
    from tools.scout_nyman_balanced import scout
else:
    from scout_nyman_balanced import scout


SCHEMA = "rh-lab/nyman-balanced-multiplier-scout/v1"
GRID_N = (8, 16, 32, 64, 128, 256)
WIDTHS_N256 = (2, 4, 8, 16, 32, 64, 128, 256)
FIXED_WIDTH = 64
GRID_INTERVAL_LIMIT = 1_000_000
CONVERGENCE_INTERVAL_LIMIT = 4_000_000


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _raw_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate_binding(root: Path, n: int) -> dict[str, str]:
    path = root / "results" / "nyman-natural-v1" / "cells" / f"n-{n:04d}.json"
    artifact = json.loads(path.read_text(encoding="utf-8"))
    return {
        "n": str(n),
        "path": path.relative_to(root).as_posix(),
        "raw_sha256": _raw_sha256(path),
        "candidate_sha256": artifact["candidate"]["candidate_sha256"],
    }


def generate(root: Path) -> dict[str, Any]:
    fixed_width = [
        scout(
            n,
            FIXED_WIDTH,
            GRID_INTERVAL_LIMIT,
            root,
            "direct-gain",
        )
        for n in GRID_N
    ]
    by_width: list[dict[str, Any]] = []
    for width in WIDTHS_N256:
        if width == FIXED_WIDTH:
            by_width.append(fixed_width[-1])
        else:
            by_width.append(
                scout(
                    256,
                    width,
                    GRID_INTERVAL_LIMIT,
                    root,
                    "direct-gain",
                )
            )
    convergence = scout(
        256,
        FIXED_WIDTH,
        CONVERGENCE_INTERVAL_LIMIT,
        root,
        "direct-gain",
    )

    body = {
        "schema": SCHEMA,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "as_of": "2026-07-23",
        "reporting_plan": {
            "status": "RETROSPECTIVE_EXPLORATORY_GRID_NOT_PREREGISTERED",
            "fixed_width": {
                "N": list(GRID_N),
                "K": FIXED_WIDTH,
                "interval_limit": GRID_INTERVAL_LIMIT,
                "objective": "direct-gain",
            },
            "width_sweep": {
                "N": 256,
                "K": list(WIDTHS_N256),
                "interval_limit": GRID_INTERVAL_LIMIT,
                "objective": "direct-gain",
            },
            "convergence": {
                "N": 256,
                "K": FIXED_WIDTH,
                "interval_limit": CONVERGENCE_INTERVAL_LIMIT,
                "objective": "direct-gain",
                "interpretation": (
                    "one vector optimized at the full cutoff and evaluated "
                    "on four fixed prefixes"
                ),
            },
        },
        "method": {
            "shell": (
                "Ideal first-shell projection reconstructed from each stored "
                "exact-dyadic finite Nyman candidate in binary64"
            ),
            "constraints": "c_1=1 and sum_(k<=K)c_k/k=0",
            "interval_error": (
                "D_c(M)=Y(M)-sum_k c_k sum_j y_j floor(M/(j*k))"
            ),
            "alias_norm": "sum_M D_c(M)^2/(M*(M+1))",
            "direct_gain": "G+2*sum_M D_c(M)*R_N(M)-alias_norm",
            "optimization": (
                "binary64 blocked normal equations after exact affine "
                "elimination of both multiplier constraints"
            ),
            "support_bound": "max support <=2*N*K",
            "energy_fraction_denominator": (
                "midpoint of the stored candidate's certified energy bracket"
            ),
        },
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "platform": platform.platform(),
            "byteorder": sys.byteorder,
        },
        "source_bindings": [_candidate_binding(root, n) for n in GRID_N],
        "fixed_width_grid": {
            "K": str(FIXED_WIDTH),
            "interval_limit": str(GRID_INTERVAL_LIMIT),
            "cells": fixed_width,
        },
        "n256_width_sweep": {
            "N": "256",
            "interval_limit": str(GRID_INTERVAL_LIMIT),
            "cells": by_width,
        },
        "n256_k64_convergence": convergence,
        "checks": {
            "all_cells_remain_exploratory": all(
                cell["classification"] == "EXPLORATORY"
                for cell in [*fixed_width, *by_width, convergence]
            ),
            "all_cells_keep_rh_unresolved": all(
                cell["hypothesis_status"] == "UNRESOLVED"
                for cell in [*fixed_width, *by_width, convergence]
            ),
            "fixed_k64_truncated_direct_gain_positive_on_grid": all(
                cell["partial_norms"][-1][
                    "direct_gain_fraction_of_old_energy_bracket_midpoint"
                ]
                > 0.0
                for cell in fixed_width
            ),
            "all_harmonic_balance_residuals_below_2e_15": all(
                abs(cell["harmonic_balance_residual"]) < 2e-15
                for cell in [*fixed_width, *by_width, convergence]
            ),
            "full_infinite_alias_tail_certified": False,
            "full_signed_direct_gain_tail_certified": False,
            "uniform_all_scale_estimate_proved": False,
            "exact_dyadic_large_scale_scout_vectors_stored": False,
            "resolves_rh": False,
        },
        "limitation": (
            "This artifact freezes binary64 finite-cutoff diagnostics. The "
            "alias norm omits a positive tail; the direct gain also omits a "
            "signed cross-term tail and is neither a lower nor an upper "
            "bound. Platform linear algebra may change final floating-point "
            "digits. No uniform estimate is proved, and RH remains unresolved."
        ),
    }
    return {**body, "payload_sha256": _canonical_sha256(body)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/nyman-balanced-multiplier-scout-v1.json"),
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output
    if not output.is_absolute():
        output = root / output
    artifact = generate(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(output)
    print(
        json.dumps(
            {
                "classification": artifact["classification"],
                "hypothesis_status": artifact["hypothesis_status"],
                "output": str(output),
                "payload_sha256": artifact["payload_sha256"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
