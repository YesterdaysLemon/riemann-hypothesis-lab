"""Certify a six-scale fixed-width balanced-correction ladder."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
from typing import Any, Mapping

import flint
import numpy as np

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record

if __package__:
    from tools.certify_nyman_balanced_prefix import (
        PrefixBound,
        _arb_from_fraction,
        certified_prefix_gain,
    )
    from tools.certify_nyman_balanced_tail import (
        ExactTailInputs,
        build_exact_vectors,
        coefficient_sum,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        strict_json_loads,
    )
    from tools.generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as V1_TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record,
        _strict_lower_proved,
    )
    from tools.generate_nyman_large_sieve_tail_certificate import (
        REFERENCE,
        _tail_record,
        _verify_tail_record,
    )
else:
    from certify_nyman_balanced_prefix import (
        PrefixBound,
        _arb_from_fraction,
        certified_prefix_gain,
    )
    from certify_nyman_balanced_tail import (
        ExactTailInputs,
        build_exact_vectors,
        coefficient_sum,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        strict_json_loads,
    )
    from generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as V1_TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record,
        _strict_lower_proved,
    )
    from generate_nyman_large_sieve_tail_certificate import (
        REFERENCE,
        _tail_record,
        _verify_tail_record,
    )


SCHEMA = "rh-lab/nyman-large-sieve-scaling-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-large-sieve-scaling-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "e974817e59e33efd84b18aa6e81d0592e96c7bb5fd5584dfb83f7763de8c0320"
)
N_VALUES = (8, 16, 32, 64, 128, 256)
CLAIMED_GAIN_LOWER_BOUNDS = {
    8: Fraction(1, 150),
    16: Fraction(1, 300),
    32: Fraction(1, 625),
    64: Fraction(1, 750),
    128: Fraction(1, 1_600),
    256: Fraction(1, 5_000),
}
CONFIG = {
    "multiplier_limit": 64,
    "shell_bits": 9,
    "z_bits": 9,
    "old_bits": 16,
    "cutoff_per_n": 1 << 15,
    "generation_block_size": 1 << 18,
    "generation_precision_bits": 256,
    "replay_block_size": 1 << 17,
    "replay_precision_bits": 384,
    "scout_section": "fixed_width_grid",
}
AS_OF = "2026-07-23"
STATEMENT = (
    "For each N in {8,16,32,64,128,256}, one explicit exact-dyadic "
    "balanced K=64 correction has strictly positive complete infinite "
    "weighted-L2 direct gain. Every prefix ends at the common-rule cutoff "
    "T=32768*N, and every omitted tail uses the Fourier/Farey large-sieve "
    "bound."
)
LIMITATION = (
    "This certifies six finite contractions selected from one frozen numerical "
    "scout. It does not prove that the construction continues, that the lower "
    "bounds have a uniform asymptotic size, or that iterating these "
    "corrections drives a Nyman distance to zero. RH remains unresolved."
)
METHOD = {
    "selection": "frozen direct-gain fixed_width_grid cells with K=64",
    "exactization": (
        "old coefficients on 2^-16; shells and multiplier harmonic variables "
        "on 2^-9; exact balance and Dirichlet convolution"
    ),
    "uniform_cutoff_rule": "T=32768*N",
    "periodic_tail": (
        "C_LS=Q(Q-1)(rho+tau), tau=rho+mu-2nu, followed by Abel summation"
    ),
    "prefix": (
        "exact int64 divisor recurrences, certified binary64 reductions, and "
        "Arb logarithmic compression"
    ),
    "interpretation": (
        "the six signs are certified finite; any scaling trend or continuation "
        "past N=256 is exploratory"
    ),
}
TRUSTED_COMPUTING_BASE = dict(V1_TRUSTED_COMPUTING_BASE)
EXPECTED_CHECKS = {
    "all_six_cells_present": True,
    "all_uniform_cutoffs": True,
    "all_corrections_balanced": True,
    "all_complete_gains_strictly_above_claims": True,
    "resolves_rh": False,
}
TOP_LEVEL_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "as_of",
    "statement",
    "configuration",
    "cells",
    "trend_diagnostics",
    "method",
    "reference",
    "trusted_computing_base",
    "checks",
    "limitation",
    "payload_sha256",
}


def _cutoff(n: int) -> int:
    return CONFIG["cutoff_per_n"] * n


def _construct_vectors(root: Path, n: int) -> ExactTailInputs:
    return build_exact_vectors(
        root,
        n=n,
        multiplier_limit=CONFIG["multiplier_limit"],
        shell_bits=CONFIG["shell_bits"],
        z_bits=CONFIG["z_bits"],
        old_bits=CONFIG["old_bits"],
        scout_section=CONFIG["scout_section"],
    )


def _vector_commitment(values: Mapping[int, Fraction]) -> dict[str, str]:
    full_record = _dyadic_vector_record(values)
    return {
        "denominator_exponent": str(full_record["denominator_exponent"]),
        "support_count": str(full_record["support_count"]),
        "maximum_index": str(full_record["maximum_index"]),
        "payload_sha256": content_sha256(full_record),
    }


def _vector_commitments(vectors: ExactTailInputs) -> dict[str, object]:
    return {
        "old": _vector_commitment(vectors.old_coefficients),
        "shell": _vector_commitment(vectors.shell),
        "multiplier": _vector_commitment(vectors.multiplier),
        "added": _vector_commitment(vectors.added_coefficients),
    }


def _full_gain_record(full_gain: Any) -> dict[str, object]:
    """Store only the authoritative Arb enclosure, without helper endpoints."""

    return {"enclosure": arb_record(full_gain)}


def _cell_checks(
    n: int,
    vectors: ExactTailInputs,
    full_gain: Any,
) -> dict[str, bool]:
    threshold = _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUNDS[n])
    return {
        "old_support_at_most_n": max(vectors.old_coefficients) <= n,
        "added_support_at_most_2nk": (
            max(vectors.added_coefficients)
            <= 2 * n * CONFIG["multiplier_limit"]
        ),
        "added_coefficient_sum_zero": (
            coefficient_sum(vectors.added_coefficients) == 0
        ),
        "added_harmonic_sum_zero": (
            harmonic_sum(vectors.added_coefficients) == 0
        ),
        "uniform_cutoff_rule": _cutoff(n) == CONFIG["cutoff_per_n"] * n,
        "full_gain_strictly_above_claim": _strict_lower_proved(
            full_gain,
            threshold,
        ),
        "resolves_rh": False,
    }


def _generate_cell(root: Path, n: int) -> tuple[dict[str, object], Any]:
    vectors = _construct_vectors(root, n)
    cutoff = _cutoff(n)
    prefix = certified_prefix_gain(
        vectors.old_coefficients,
        vectors.added_coefficients,
        cutoff,
        block_size=CONFIG["generation_block_size"],
        precision_bits=CONFIG["generation_precision_bits"],
    )
    tail = large_sieve_tail_center_radius(
        vectors.old_coefficients,
        vectors.added_coefficients,
        cutoff,
    )
    full_gain = _full_gain(prefix, tail, CONFIG["generation_precision_bits"])
    checks = _cell_checks(n, vectors, full_gain)
    if not all(value is True for key, value in checks.items() if key != "resolves_rh"):
        raise ArithmeticError(f"N={n} scaling cell did not close")
    if checks["resolves_rh"] is not False:
        raise ArithmeticError("scaling cell incorrectly resolves RH")
    record: dict[str, object] = {
        "n": str(n),
        "cutoff": str(cutoff),
        "claimed_gain_lower_bound": _fraction_record(
            CLAIMED_GAIN_LOWER_BOUNDS[n]
        ),
        "source": {
            "scout_section": CONFIG["scout_section"],
            **_source_record(vectors),
        },
        "vector_commitments": _vector_commitments(vectors),
        "prefix": _prefix_record(prefix),
        "tail": _tail_record(tail),
        "full_gain": _full_gain_record(full_gain),
        "checks": checks,
    }
    return record, tail


def _trend_diagnostics(
    cells_and_tails: list[tuple[dict[str, object], Any]],
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for cell, tail in cells_and_tails:
        rows.append(
            {
                "n": cell["n"],
                "support_limit": str(tail.support_limit),
                "cutoff_to_support_ratio": _fraction_record(
                    Fraction(int(cell["cutoff"]), tail.support_limit)
                ),
                "old_plus_new_periodic_mean_square": _fraction_record(
                    tail.old_periodic_mean_square
                    + tail.new_periodic_mean_square
                ),
                "periodic_tail_radius": _fraction_record(
                    tail.periodic_radius
                ),
            }
        )
    return {
        "classification": "EXPLORATORY",
        "rows": rows,
        "observation": (
            "All six complete gains are positive at an approximately constant "
            "cutoff/support ratio, but the periodic energy grows and the "
            "certified gain fraction deteriorates across this finite range."
        ),
        "limitation": (
            "Six points cannot establish an asymptotic energy or gain law."
        ),
    }


def _verify_prefix_record(
    record: object,
    replay: PrefixBound,
    *,
    cutoff: int,
) -> None:
    expected_fields = set(_prefix_record(replay))
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError("scaling prefix record fields changed")
    expected_blocks = (
        cutoff + CONFIG["generation_block_size"] - 1
    ) // CONFIG["generation_block_size"]
    integer_expectations = {
        "cutoff": cutoff,
        "block_size": CONFIG["generation_block_size"],
        "block_count": expected_blocks,
        "added_denominator": replay.added_denominator,
        "old_denominator": replay.old_denominator,
        "last_added_step_numerator": replay.last_added_step_numerator,
        "last_old_intercept_numerator": replay.last_old_intercept_numerator,
        "maximum_added_step_numerator": replay.maximum_added_step_numerator,
        "maximum_old_intercept_numerator": replay.maximum_old_intercept_numerator,
        "maximum_interval_numerator": replay.maximum_interval_numerator,
    }
    for name, expected in integer_expectations.items():
        supplied = record[name]
        if not isinstance(supplied, str) or str(int(supplied)) != supplied:
            raise ValueError(f"scaling prefix {name} is not canonical")
        if int(supplied) != expected:
            raise ValueError(f"scaling prefix {name} changed")
    midpoint = _parse_fraction(record["rational_midpoint"], "prefix midpoint")
    error = _parse_fraction(record["rational_error"], "prefix error")
    lower = _parse_fraction(record["rational_lower"], "prefix lower")
    upper = _parse_fraction(record["rational_upper"], "prefix upper")
    if error < 0 or lower != midpoint - error or upper != midpoint + error:
        raise ValueError("scaling rational prefix identities failed")
    if lower > replay.rational_lower or upper < replay.rational_upper:
        raise ValueError("scaling rational prefix does not contain replay")
    stored_absolute = _parse_fraction(
        record["computed_absolute_term_sum_upper"],
        "prefix absolute term sum",
    )
    if stored_absolute < replay.computed_absolute_term_sum_upper:
        raise ValueError("scaling absolute term sum does not contain replay")
    stored_prefix = _ball_from_record(record["prefix_gain"], "stored prefix")
    if not stored_prefix.contains(replay.prefix_gain):
        raise ValueError("scaling prefix ball does not contain replay")
    stored_log = _ball_from_record(record["logarithmic_gain"], "stored log")
    if not stored_log.contains(replay.logarithmic_gain):
        raise ValueError("scaling logarithmic ball does not contain replay")


def _verify_full_gain_record(
    record: object,
    replay: Any,
    threshold: Fraction,
) -> None:
    expected_fields = {"enclosure"}
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError("scaling full gain record is missing")
    stored = _ball_from_record(record["enclosure"], "stored full gain")
    threshold_ball = _arb_from_fraction(threshold)
    if not _strict_lower_proved(stored, threshold_ball):
        raise ValueError("stored scaling gain does not prove its claim")
    if not _strict_lower_proved(replay, threshold_ball):
        raise ValueError("replayed scaling gain does not prove its claim")
    if not stored.contains(replay):
        raise ValueError("stored scaling gain does not contain replay")


def generate(root: Path) -> dict[str, object]:
    _require_environment()
    cells_and_tails = [_generate_cell(root, n) for n in N_VALUES]
    cells = [cell for cell, _ in cells_and_tails]
    checks = {
        "all_six_cells_present": len(cells) == 6,
        "all_uniform_cutoffs": all(
            int(cell["cutoff"]) == _cutoff(int(cell["n"]))
            for cell in cells
        ),
        "all_corrections_balanced": all(
            cell["checks"]["added_coefficient_sum_zero"]
            and cell["checks"]["added_harmonic_sum_zero"]
            for cell in cells
        ),
        "all_complete_gains_strictly_above_claims": all(
            cell["checks"]["full_gain_strictly_above_claim"]
            for cell in cells
        ),
        "resolves_rh": False,
    }
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("scaling certificate checks did not close")
    body: dict[str, object] = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": {
            **{key: str(value) for key, value in CONFIG.items()},
            "n_values": [str(n) for n in N_VALUES],
        },
        "cells": cells,
        "trend_diagnostics": _trend_diagnostics(cells_and_tails),
        "method": METHOD,
        "reference": REFERENCE,
        "trusted_computing_base": TRUSTED_COMPUTING_BASE,
        "checks": checks,
        "limitation": LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    supplied = strict_json_loads(artifact_path.read_text(encoding="utf-8"))
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("scaling certificate top-level fields changed")
    body = {key: value for key, value in supplied.items() if key != "payload_sha256"}
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("scaling certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("scaling certificate is not the frozen v1 payload")
    expected_scalars = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "limitation": LIMITATION,
    }
    for name, expected in expected_scalars.items():
        if supplied.get(name) != expected:
            raise ValueError(f"scaling certificate {name} changed")
    expected_configuration = {
        **{key: str(value) for key, value in CONFIG.items()},
        "n_values": [str(n) for n in N_VALUES],
    }
    if supplied.get("configuration") != expected_configuration:
        raise ValueError("scaling certificate configuration changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("scaling certificate method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("scaling certificate trusted-computing-base changed")
    stored_cells = supplied.get("cells")
    if not isinstance(stored_cells, list) or len(stored_cells) != len(N_VALUES):
        raise ValueError("scaling certificate cells changed")

    replay_rows: list[dict[str, object]] = []
    cells_and_tails: list[tuple[dict[str, object], Any]] = []
    for expected_n, cell in zip(N_VALUES, stored_cells):
        cell_fields = {
            "n",
            "cutoff",
            "claimed_gain_lower_bound",
            "source",
            "vector_commitments",
            "prefix",
            "tail",
            "full_gain",
            "checks",
        }
        if not isinstance(cell, Mapping) or set(cell) != cell_fields:
            raise ValueError("scaling cell fields changed")
        if str(expected_n) != cell.get("n"):
            raise ValueError("scaling cell order or N changed")
        cutoff = _cutoff(expected_n)
        if cell.get("cutoff") != str(cutoff):
            raise ValueError(f"N={expected_n} cutoff changed")
        threshold = _parse_fraction(
            cell.get("claimed_gain_lower_bound"),
            f"N={expected_n} claimed gain",
        )
        if threshold != CLAIMED_GAIN_LOWER_BOUNDS[expected_n]:
            raise ValueError(f"N={expected_n} claimed gain changed")
        vectors = _construct_vectors(root, expected_n)
        expected_source = {
            "scout_section": CONFIG["scout_section"],
            **_source_record(vectors),
        }
        if cell.get("source") != expected_source:
            raise ValueError(f"N={expected_n} source binding changed")
        if cell.get("vector_commitments") != _vector_commitments(vectors):
            raise ValueError(f"N={expected_n} vector commitments changed")
        require_balanced(vectors.added_coefficients, f"N={expected_n} added")
        replay_prefix = certified_prefix_gain(
            vectors.old_coefficients,
            vectors.added_coefficients,
            cutoff,
            block_size=CONFIG["replay_block_size"],
            precision_bits=CONFIG["replay_precision_bits"],
        )
        replay_tail = large_sieve_tail_center_radius(
            vectors.old_coefficients,
            vectors.added_coefficients,
            cutoff,
        )
        _verify_prefix_record(cell.get("prefix"), replay_prefix, cutoff=cutoff)
        _verify_tail_record(cell.get("tail"), replay_tail)
        replay_full = _full_gain(
            replay_prefix,
            replay_tail,
            CONFIG["replay_precision_bits"],
        )
        _verify_full_gain_record(cell.get("full_gain"), replay_full, threshold)
        expected_cell_checks = _cell_checks(expected_n, vectors, replay_full)
        if cell.get("checks") != expected_cell_checks:
            raise ValueError(f"N={expected_n} checks changed")
        cells_and_tails.append((dict(cell), replay_tail))
        replay_rows.append(
            {
                "n": str(expected_n),
                "lower": replay_full.lower().str(18),
                "upper": replay_full.upper().str(18),
            }
        )
    if supplied.get("trend_diagnostics") != _trend_diagnostics(cells_and_tails):
        raise ValueError("scaling trend diagnostics changed")
    if supplied.get("checks") != EXPECTED_CHECKS:
        raise ValueError("scaling certificate check ledger changed")
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "replay_cells": replay_rows,
        "replay_environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "python_flint": flint.__version__,
            "numpy": np.__version__,
        },
        "verified": True,
    }


def _raw_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    artifact_path = args.artifact
    if not artifact_path.is_absolute():
        artifact_path = root / artifact_path
    if args.verify:
        print(json.dumps(verify(root, artifact_path), sort_keys=True))
        return
    artifact = generate(root)
    artifact_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = artifact_path.with_name(f".{artifact_path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    temporary.replace(artifact_path)
    print(
        json.dumps(
            {
                "artifact": str(artifact_path),
                "payload_sha256": artifact["payload_sha256"],
                "raw_sha256": _raw_sha256(artifact_path),
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
