"""Generate and independently replay a genuine two-step Nyman chain.

The first update is the frozen N=8, K=64 cell from the large-sieve scaling
certificate.  The second update is constructed from the *updated* first
vector, rather than from another independently selected old vector.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
from typing import Any, Mapping

import flint
from flint import ctx
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
        convolve_coefficients,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        rounded_ideal_shell,
        strict_json_loads,
    )
    from tools.generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record as _exact_source_record,
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
        convolve_coefficients,
        harmonic_sum,
        large_sieve_tail_center_radius,
        require_balanced,
        rounded_ideal_shell,
        strict_json_loads,
    )
    from generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _full_gain,
        _parse_fraction,
        _prefix_record,
        _require_environment,
        _source_record as _exact_source_record,
        _strict_lower_proved,
    )
    from generate_nyman_large_sieve_tail_certificate import (
        REFERENCE,
        _tail_record,
        _verify_tail_record,
    )


SCHEMA = "rh-lab/nyman-nested-chain-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-nested-chain-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "1171b416a8aa949ce380d5687cf8865d5550891ac24ab0d781cfa8261627bcad"
)

FROZEN_SCALING_PATH = Path("results/nyman-large-sieve-scaling-v1.json")
FROZEN_SCALING_SCHEMA = "rh-lab/nyman-large-sieve-scaling-certificate/v1"
FROZEN_SCALING_RAW_SHA256 = (
    "02886cec84010235b3a6e564b7aba8d42e57bbcf2f4e28783dc91894d01f9cce"
)
FROZEN_SCALING_PAYLOAD_SHA256 = (
    "e974817e59e33efd84b18aa6e81d0592e96c7bb5fd5584dfb83f7763de8c0320"
)
FROZEN_SCALING_N8_CELL_SHA256 = (
    "0f9b0bb2c83319204464cb3a6f7c1d7f72312d4cc7d24f5c01bbeb4fac72b8a6"
)

STEP1_CLAIM = Fraction(1, 150)
STEP2_CLAIM = Fraction(1, 9_000)
TOTAL_CLAIM = STEP1_CLAIM + STEP2_CLAIM
STEP2_ALPHA = Fraction(1, 4)
STEP2_MULTIPLIER_VALUES = (
    Fraction(1),
    -Fraction(317, 256),
    -Fraction(591, 512),
    Fraction(5, 128),
    -Fraction(255, 512),
    Fraction(141, 128),
    -Fraction(133, 256),
    -Fraction(1, 8),
)

CONFIG = {
    "step1_n": 8,
    "step1_multiplier_limit": 64,
    "step1_shell_bits": 9,
    "step1_z_bits": 9,
    "step1_old_bits": 16,
    "step1_cutoff": 1 << 18,
    "step2_shell_cutoff": 1 << 10,
    "step2_shell_bits": 9,
    "step2_multiplier_limit": 8,
    "step2_cutoff": 1 << 25,
    "generation_block_size": 1 << 18,
    "generation_precision_bits": 256,
    "replay_block_size": 196_613,
    "replay_precision_bits": 448,
    "scout_section": "fixed_width_grid",
}
AS_OF = "2026-07-23"
STATEMENT = (
    "For explicit exact-dyadic Nyman vectors p0, p1=p0+a1, and p2=p1+a2, "
    "the complete infinite weighted-L2 direct gain from p0 to p1 is greater "
    "than 1/150 and the complete gain from that same p1 to p2 is greater "
    "than 1/9000. Thus the two-step telescoped gain is greater than 61/9000."
)
LIMITATION = (
    "This is one genuine nested chain of length two. It proves neither a "
    "uniform recurrence nor convergence of Nyman distances, and it does not "
    "prove or disprove the Riemann Hypothesis. RH remains unresolved."
)
METHOD = {
    "step_1": (
        "reconstruct the frozen N=8, K=64 fixed-width scaling cell exactly"
    ),
    "nesting": "set p1=p0+a1 before constructing or certifying the second step",
    "step_2_shell": (
        "construct the sparse p1 ideal shell with explicit support cutoff "
        "1024 on indices 1025..2048, and round to exact multiples of 2^-9"
    ),
    "step_2_multiplier": (
        "use the prescribed exact harmonic-balanced K=8 multiplier and outer "
        "scale alpha=1/4"
    ),
    "prefix": (
        "exact int64 divisor recurrences, pinned binary64 reductions with an "
        "exact forward-error radius, and Arb logarithmic compression"
    ),
    "tail": (
        "Fourier/Farey large-sieve periodic discrepancy followed by Abel "
        "summation"
    ),
    "selection_boundary": (
        "the step-2 multiplier was selected in an exploratory numerical pilot; "
        "this certificate freezes it and recomputes both steps from exact data"
    ),
}

EXPECTED_CHECKS = {
    "has_exactly_two_steps": True,
    "step_1_is_frozen_scaling_n8_cell": True,
    "step_1_update_identity": True,
    "step_2_old_equals_step_1_updated": True,
    "step_2_update_identity": True,
    "both_corrections_balanced": True,
    "both_complete_gains_strictly_above_claims": True,
    "telescoped_gain_strictly_above_sum_of_claims": True,
    "resolves_rh": False,
}
TOP_LEVEL_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "as_of",
    "statement",
    "configuration",
    "claims",
    "source",
    "chain",
    "telescoped_full_gain",
    "method",
    "reference",
    "trusted_computing_base",
    "checks",
    "limitation",
    "payload_sha256",
}
STEP_FIELDS = {
    "name",
    "relation",
    "cutoff",
    "claimed_gain_lower_bound",
    "vector_commitments",
    "prefix",
    "tail",
    "full_gain",
    "checks",
}
VECTOR_NAMES = {"old", "shell", "multiplier", "added", "updated"}


@dataclass(frozen=True)
class NestedChain:
    """Exact vectors and source bindings for the two linked updates."""

    step1_inputs: ExactTailInputs
    p0: Mapping[int, Fraction]
    a1: Mapping[int, Fraction]
    p1: Mapping[int, Fraction]
    shell2: Mapping[int, Fraction]
    multiplier2: Mapping[int, Fraction]
    a2: Mapping[int, Fraction]
    p2: Mapping[int, Fraction]
    source: Mapping[str, object]


def _clean(values: Mapping[int, Fraction]) -> dict[int, Fraction]:
    return {
        int(index): Fraction(value)
        for index, value in values.items()
        if Fraction(value)
    }


def _add_vectors(
    left: Mapping[int, Fraction],
    right: Mapping[int, Fraction],
) -> dict[int, Fraction]:
    result = dict(_clean(left))
    for index, value in right.items():
        result[int(index)] = result.get(int(index), Fraction(0)) + Fraction(value)
    return _clean(result)


def _configuration_record() -> dict[str, object]:
    return {
        **{name: str(value) for name, value in CONFIG.items()},
        "step2_alpha": _fraction_record(STEP2_ALPHA),
        "step2_multiplier": [
            _fraction_record(value) for value in STEP2_MULTIPLIER_VALUES
        ],
    }


def _claims_record() -> dict[str, object]:
    return {
        "step_1_gain_strictly_greater_than": _fraction_record(STEP1_CLAIM),
        "step_2_gain_strictly_greater_than": _fraction_record(STEP2_CLAIM),
        "telescoped_gain_strictly_greater_than": _fraction_record(TOTAL_CLAIM),
        "hypothesis_status": "UNRESOLVED",
    }


def _vector_commitment(values: Mapping[int, Fraction]) -> dict[str, str]:
    record = _dyadic_vector_record(values)
    return {
        "denominator_exponent": str(record["denominator_exponent"]),
        "support_count": str(record["support_count"]),
        "maximum_index": str(record["maximum_index"]),
        "payload_sha256": content_sha256(record),
    }


def _vector_commitments(
    old: Mapping[int, Fraction],
    shell: Mapping[int, Fraction],
    multiplier: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    updated: Mapping[int, Fraction],
) -> dict[str, object]:
    return {
        "old": _vector_commitment(old),
        "shell": _vector_commitment(shell),
        "multiplier": _vector_commitment(multiplier),
        "added": _vector_commitment(added),
        "updated": _vector_commitment(updated),
    }


def _load_scaling_source(
    root: Path,
    step1: ExactTailInputs,
) -> dict[str, object]:
    path = root / FROZEN_SCALING_PATH
    raw = path.read_bytes()
    raw_sha = hashlib.sha256(raw).hexdigest()
    if raw_sha != FROZEN_SCALING_RAW_SHA256:
        raise ValueError("frozen scaling certificate raw SHA-256 mismatch")
    artifact = strict_json_loads(raw.decode("utf-8"))
    if not isinstance(artifact, Mapping):
        raise ValueError("frozen scaling certificate must be an object")
    if artifact.get("schema") != FROZEN_SCALING_SCHEMA:
        raise ValueError("frozen scaling certificate schema changed")
    body = {
        key: value
        for key, value in artifact.items()
        if key != "payload_sha256"
    }
    if artifact.get("payload_sha256") != content_sha256(body):
        raise ValueError("frozen scaling certificate payload hash mismatch")
    if artifact.get("payload_sha256") != FROZEN_SCALING_PAYLOAD_SHA256:
        raise ValueError("unexpected frozen scaling certificate payload")
    cells = artifact.get("cells")
    if not isinstance(cells, list):
        raise ValueError("frozen scaling certificate has no cells")
    matches = [
        cell
        for cell in cells
        if isinstance(cell, Mapping) and cell.get("n") == "8"
    ]
    if len(matches) != 1:
        raise ValueError("frozen scaling certificate has no unique N=8 cell")
    cell = matches[0]
    if content_sha256(cell) != FROZEN_SCALING_N8_CELL_SHA256:
        raise ValueError("frozen N=8 scaling cell changed")
    expected_step1_commitments = {
        name: _vector_commitment(values)
        for name, values in {
            "old": step1.old_coefficients,
            "shell": step1.shell,
            "multiplier": step1.multiplier,
            "added": step1.added_coefficients,
        }.items()
    }
    if cell.get("vector_commitments") != expected_step1_commitments:
        raise ValueError("frozen N=8 scaling cell vectors do not reconstruct")
    return {
        "scaling_certificate": {
            "path": FROZEN_SCALING_PATH.as_posix(),
            "raw_sha256": FROZEN_SCALING_RAW_SHA256,
            "payload_sha256": FROZEN_SCALING_PAYLOAD_SHA256,
            "n8_cell_payload_sha256": FROZEN_SCALING_N8_CELL_SHA256,
        },
        "step_1_exact_sources": _exact_source_record(step1),
    }


def _construct_chain(root: Path) -> NestedChain:
    step1 = build_exact_vectors(
        root,
        n=CONFIG["step1_n"],
        multiplier_limit=CONFIG["step1_multiplier_limit"],
        shell_bits=CONFIG["step1_shell_bits"],
        z_bits=CONFIG["step1_z_bits"],
        old_bits=CONFIG["step1_old_bits"],
        scout_section=CONFIG["scout_section"],
    )
    p0 = _clean(step1.old_coefficients)
    a1 = _clean(step1.added_coefficients)
    require_balanced(a1, "step-1 correction")
    p1 = _add_vectors(p0, a1)
    shell_cutoff = CONFIG["step2_shell_cutoff"]
    if max(p1, default=0) > shell_cutoff:
        raise ValueError("step-1 updated vector exceeds step-2 shell cutoff")
    shell2 = rounded_ideal_shell(
        p1,
        CONFIG["step2_shell_bits"],
        support_cutoff=shell_cutoff,
    )
    multiplier2 = {
        index: value
        for index, value in enumerate(STEP2_MULTIPLIER_VALUES, start=1)
        if value
    }
    if harmonic_sum(multiplier2):
        raise ArithmeticError("prescribed step-2 multiplier is not balanced")
    a2 = _clean(
        {
            index: STEP2_ALPHA * value
            for index, value in convolve_coefficients(
                shell2,
                multiplier2,
            ).items()
        }
    )
    require_balanced(a2, "step-2 correction")
    p2 = _add_vectors(p1, a2)
    source = _load_scaling_source(root, step1)
    return NestedChain(
        step1_inputs=step1,
        p0=p0,
        a1=a1,
        p1=p1,
        shell2=shell2,
        multiplier2=multiplier2,
        a2=a2,
        p2=p2,
        source=source,
    )


def _full_gain_record(ball: Any) -> dict[str, object]:
    return {"enclosure": arb_record(ball)}


def _add_balls(left: Any, right: Any, precision_bits: int) -> Any:
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        return left + right
    finally:
        ctx.prec = previous_precision


def _step1_checks(chain: NestedChain, full_gain: Any) -> dict[str, bool]:
    return {
        "old_support_at_most_8": max(chain.p0) <= 8,
        "added_support_at_most_1024": max(chain.a1) <= 1024,
        "added_coefficient_sum_zero": coefficient_sum(chain.a1) == 0,
        "added_harmonic_sum_zero": harmonic_sum(chain.a1) == 0,
        "updated_equals_old_plus_added": chain.p1 == _add_vectors(
            chain.p0,
            chain.a1,
        ),
        "full_gain_strictly_above_claim": _strict_lower_proved(
            full_gain,
            _arb_from_fraction(STEP1_CLAIM),
        ),
        "resolves_rh": False,
    }


def _step2_checks(chain: NestedChain, full_gain: Any) -> dict[str, bool]:
    shell_min = min(chain.shell2, default=0)
    shell_max = max(chain.shell2, default=0)
    return {
        "old_is_step_1_updated": chain.p1 == _add_vectors(chain.p0, chain.a1),
        "old_support_at_most_1024": max(chain.p1) <= 1024,
        "shell_support_within_1025_through_2048": (
            shell_min >= 1025 and shell_max <= 2048
        ),
        "shell_coefficient_sum_zero": coefficient_sum(chain.shell2) == 0,
        "multiplier_is_prescribed": chain.multiplier2
        == {
            index: value
            for index, value in enumerate(STEP2_MULTIPLIER_VALUES, start=1)
            if value
        },
        "multiplier_harmonic_sum_zero": harmonic_sum(chain.multiplier2) == 0,
        "added_support_at_most_16384": max(chain.a2) <= 16384,
        "added_coefficient_sum_zero": coefficient_sum(chain.a2) == 0,
        "added_harmonic_sum_zero": harmonic_sum(chain.a2) == 0,
        "updated_equals_old_plus_added": chain.p2 == _add_vectors(
            chain.p1,
            chain.a2,
        ),
        "full_gain_strictly_above_claim": _strict_lower_proved(
            full_gain,
            _arb_from_fraction(STEP2_CLAIM),
        ),
        "resolves_rh": False,
    }


def _step_record(
    *,
    name: str,
    relation: str,
    cutoff: int,
    claim: Fraction,
    old: Mapping[int, Fraction],
    shell: Mapping[int, Fraction],
    multiplier: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    updated: Mapping[int, Fraction],
    checks: Mapping[str, bool],
    prefix: PrefixBound,
    tail: Any,
    full_gain: Any,
) -> dict[str, object]:
    return {
        "name": name,
        "relation": relation,
        "cutoff": str(cutoff),
        "claimed_gain_lower_bound": _fraction_record(claim),
        "vector_commitments": _vector_commitments(
            old,
            shell,
            multiplier,
            added,
            updated,
        ),
        "prefix": _prefix_record(prefix),
        "tail": _tail_record(tail),
        "full_gain": _full_gain_record(full_gain),
        "checks": dict(checks),
    }


def _certify_step(
    old: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    cutoff: int,
) -> tuple[PrefixBound, Any, Any]:
    prefix = certified_prefix_gain(
        old,
        added,
        cutoff,
        block_size=CONFIG["generation_block_size"],
        precision_bits=CONFIG["generation_precision_bits"],
    )
    tail = large_sieve_tail_center_radius(old, added, cutoff)
    full_gain = _full_gain(
        prefix,
        tail,
        CONFIG["generation_precision_bits"],
    )
    return prefix, tail, full_gain


def generate(root: Path) -> dict[str, object]:
    _require_environment()
    chain = _construct_chain(root)
    prefix1, tail1, gain1 = _certify_step(
        chain.p0,
        chain.a1,
        CONFIG["step1_cutoff"],
    )
    prefix2, tail2, gain2 = _certify_step(
        chain.p1,
        chain.a2,
        CONFIG["step2_cutoff"],
    )
    checks1 = _step1_checks(chain, gain1)
    checks2 = _step2_checks(chain, gain2)
    if not all(value for key, value in checks1.items() if key != "resolves_rh"):
        raise ArithmeticError("step-1 checks did not close")
    if not all(value for key, value in checks2.items() if key != "resolves_rh"):
        raise ArithmeticError("step-2 checks did not close")
    total_gain = _add_balls(
        gain1,
        gain2,
        CONFIG["generation_precision_bits"],
    )
    checks = {
        "has_exactly_two_steps": True,
        "step_1_is_frozen_scaling_n8_cell": True,
        "step_1_update_identity": chain.p1
        == _add_vectors(chain.p0, chain.a1),
        "step_2_old_equals_step_1_updated": True,
        "step_2_update_identity": chain.p2
        == _add_vectors(chain.p1, chain.a2),
        "both_corrections_balanced": (
            coefficient_sum(chain.a1) == 0
            and harmonic_sum(chain.a1) == 0
            and coefficient_sum(chain.a2) == 0
            and harmonic_sum(chain.a2) == 0
        ),
        "both_complete_gains_strictly_above_claims": (
            checks1["full_gain_strictly_above_claim"]
            and checks2["full_gain_strictly_above_claim"]
        ),
        "telescoped_gain_strictly_above_sum_of_claims": _strict_lower_proved(
            total_gain,
            _arb_from_fraction(TOTAL_CLAIM),
        ),
        "resolves_rh": False,
    }
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("nested-chain check ledger did not close")
    steps = [
        _step_record(
            name="step_1",
            relation="p1=p0+a1",
            cutoff=CONFIG["step1_cutoff"],
            claim=STEP1_CLAIM,
            old=chain.p0,
            shell=chain.step1_inputs.shell,
            multiplier=chain.step1_inputs.multiplier,
            added=chain.a1,
            updated=chain.p1,
            checks=checks1,
            prefix=prefix1,
            tail=tail1,
            full_gain=gain1,
        ),
        _step_record(
            name="step_2",
            relation="p2=p1+a2",
            cutoff=CONFIG["step2_cutoff"],
            claim=STEP2_CLAIM,
            old=chain.p1,
            shell=chain.shell2,
            multiplier=chain.multiplier2,
            added=chain.a2,
            updated=chain.p2,
            checks=checks2,
            prefix=prefix2,
            tail=tail2,
            full_gain=gain2,
        ),
    ]
    body: dict[str, object] = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": _configuration_record(),
        "claims": _claims_record(),
        "source": chain.source,
        "chain": steps,
        "telescoped_full_gain": _full_gain_record(total_gain),
        "method": METHOD,
        "reference": REFERENCE,
        "trusted_computing_base": TRUSTED_COMPUTING_BASE,
        "checks": checks,
        "limitation": LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _verify_prefix_record(
    record: object,
    replay: PrefixBound,
    *,
    cutoff: int,
) -> None:
    expected_fields = set(_prefix_record(replay))
    if not isinstance(record, Mapping) or set(record) != expected_fields:
        raise ValueError("nested-chain prefix record fields changed")
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
            raise ValueError(f"nested-chain prefix {name} is not canonical")
        if int(supplied) != expected:
            raise ValueError(f"nested-chain prefix {name} changed")
    midpoint = _parse_fraction(record["rational_midpoint"], "prefix midpoint")
    error = _parse_fraction(record["rational_error"], "prefix error")
    lower = _parse_fraction(record["rational_lower"], "prefix lower")
    upper = _parse_fraction(record["rational_upper"], "prefix upper")
    if error < 0 or lower != midpoint - error or upper != midpoint + error:
        raise ValueError("nested-chain rational prefix identities failed")
    if lower > replay.rational_lower or upper < replay.rational_upper:
        raise ValueError("stored nested-chain prefix does not contain replay")
    absolute_sum = _parse_fraction(
        record["computed_absolute_term_sum_upper"],
        "prefix absolute term sum",
    )
    if absolute_sum < replay.computed_absolute_term_sum_upper:
        raise ValueError("stored prefix absolute sum does not contain replay")
    stored_prefix = _ball_from_record(record["prefix_gain"], "stored prefix")
    if not stored_prefix.contains(replay.prefix_gain):
        raise ValueError("stored prefix Arb ball does not contain replay")
    stored_log = _ball_from_record(record["logarithmic_gain"], "stored log")
    if not stored_log.contains(replay.logarithmic_gain):
        raise ValueError("stored logarithmic Arb ball does not contain replay")


def _verify_full_gain_record(
    record: object,
    replay: Any,
    threshold: Fraction,
    name: str,
) -> Any:
    if not isinstance(record, Mapping) or set(record) != {"enclosure"}:
        raise ValueError(f"{name} full-gain record fields changed")
    stored = _ball_from_record(record["enclosure"], f"{name} full gain")
    threshold_ball = _arb_from_fraction(threshold)
    if not _strict_lower_proved(stored, threshold_ball):
        raise ValueError(f"stored {name} gain does not prove its claim")
    if not _strict_lower_proved(replay, threshold_ball):
        raise ValueError(f"replayed {name} gain does not prove its claim")
    if not stored.contains(replay):
        raise ValueError(f"stored {name} gain does not contain replay")
    return stored


def _step_specs(chain: NestedChain) -> tuple[dict[str, object], ...]:
    return (
        {
            "name": "step_1",
            "relation": "p1=p0+a1",
            "cutoff": CONFIG["step1_cutoff"],
            "claim": STEP1_CLAIM,
            "old": chain.p0,
            "shell": chain.step1_inputs.shell,
            "multiplier": chain.step1_inputs.multiplier,
            "added": chain.a1,
            "updated": chain.p1,
            "checks": _step1_checks,
        },
        {
            "name": "step_2",
            "relation": "p2=p1+a2",
            "cutoff": CONFIG["step2_cutoff"],
            "claim": STEP2_CLAIM,
            "old": chain.p1,
            "shell": chain.shell2,
            "multiplier": chain.multiplier2,
            "added": chain.a2,
            "updated": chain.p2,
            "checks": _step2_checks,
        },
    )


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    supplied = strict_json_loads(artifact_path.read_text(encoding="utf-8"))
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("nested-chain certificate top-level fields changed")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("nested-chain certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("nested-chain certificate is not the frozen v1 payload")
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
            raise ValueError(f"nested-chain certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("nested-chain configuration changed")
    if supplied.get("claims") != _claims_record():
        raise ValueError("nested-chain claims changed")
    for name, expected in (
        ("step_1_gain_strictly_greater_than", STEP1_CLAIM),
        ("step_2_gain_strictly_greater_than", STEP2_CLAIM),
        ("telescoped_gain_strictly_greater_than", TOTAL_CLAIM),
    ):
        if _parse_fraction(supplied["claims"][name], name) != expected:
            raise ValueError(f"nested-chain {name} changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("nested-chain method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("nested-chain trusted-computing-base changed")

    chain = _construct_chain(root)
    if supplied.get("source") != chain.source:
        raise ValueError("nested-chain source binding changed")
    stored_steps = supplied.get("chain")
    if not isinstance(stored_steps, list) or len(stored_steps) != 2:
        raise ValueError("nested-chain step list changed")

    replay_gains: list[Any] = []
    replay_rows: list[dict[str, str]] = []
    for stored, spec in zip(stored_steps, _step_specs(chain), strict=True):
        if not isinstance(stored, Mapping) or set(stored) != STEP_FIELDS:
            raise ValueError("nested-chain step fields changed")
        name = str(spec["name"])
        cutoff = int(spec["cutoff"])
        claim = Fraction(spec["claim"])
        if stored.get("name") != name or stored.get("relation") != spec["relation"]:
            raise ValueError(f"{name} identity changed")
        if stored.get("cutoff") != str(cutoff):
            raise ValueError(f"{name} cutoff changed")
        if _parse_fraction(
            stored.get("claimed_gain_lower_bound"),
            f"{name} claimed gain",
        ) != claim:
            raise ValueError(f"{name} claimed gain changed")
        expected_commitments = _vector_commitments(
            spec["old"],
            spec["shell"],
            spec["multiplier"],
            spec["added"],
            spec["updated"],
        )
        commitments = stored.get("vector_commitments")
        if not isinstance(commitments, Mapping) or set(commitments) != VECTOR_NAMES:
            raise ValueError(f"{name} vector commitment fields changed")
        if commitments != expected_commitments:
            raise ValueError(f"{name} vector commitments changed")
        replay_prefix = certified_prefix_gain(
            spec["old"],
            spec["added"],
            cutoff,
            block_size=CONFIG["replay_block_size"],
            precision_bits=CONFIG["replay_precision_bits"],
        )
        replay_tail = large_sieve_tail_center_radius(
            spec["old"],
            spec["added"],
            cutoff,
        )
        _verify_prefix_record(stored.get("prefix"), replay_prefix, cutoff=cutoff)
        _verify_tail_record(stored.get("tail"), replay_tail)
        replay_gain = _full_gain(
            replay_prefix,
            replay_tail,
            CONFIG["replay_precision_bits"],
        )
        _verify_full_gain_record(
            stored.get("full_gain"),
            replay_gain,
            claim,
            name,
        )
        expected_step_checks = spec["checks"](chain, replay_gain)
        if stored.get("checks") != expected_step_checks:
            raise ValueError(f"{name} check ledger changed")
        replay_gains.append(replay_gain)
        replay_rows.append(
            {
                "name": name,
                "lower": replay_gain.lower().str(20),
                "upper": replay_gain.upper().str(20),
            }
        )

    replay_total = _add_balls(
        replay_gains[0],
        replay_gains[1],
        CONFIG["replay_precision_bits"],
    )
    _verify_full_gain_record(
        supplied.get("telescoped_full_gain"),
        replay_total,
        TOTAL_CLAIM,
        "telescoped",
    )
    replay_checks = {
        "has_exactly_two_steps": len(stored_steps) == 2,
        "step_1_is_frozen_scaling_n8_cell": True,
        "step_1_update_identity": chain.p1
        == _add_vectors(chain.p0, chain.a1),
        "step_2_old_equals_step_1_updated": (
            _vector_commitment(chain.p1)
            == stored_steps[1]["vector_commitments"]["old"]
        ),
        "step_2_update_identity": chain.p2
        == _add_vectors(chain.p1, chain.a2),
        "both_corrections_balanced": (
            coefficient_sum(chain.a1) == 0
            and harmonic_sum(chain.a1) == 0
            and coefficient_sum(chain.a2) == 0
            and harmonic_sum(chain.a2) == 0
        ),
        "both_complete_gains_strictly_above_claims": (
            _strict_lower_proved(
                replay_gains[0],
                _arb_from_fraction(STEP1_CLAIM),
            )
            and _strict_lower_proved(
                replay_gains[1],
                _arb_from_fraction(STEP2_CLAIM),
            )
        ),
        "telescoped_gain_strictly_above_sum_of_claims": _strict_lower_proved(
            replay_total,
            _arb_from_fraction(TOTAL_CLAIM),
        ),
        "resolves_rh": False,
    }
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("replayed nested-chain checks did not close")
    if supplied.get("checks") != replay_checks:
        raise ValueError("nested-chain top-level check ledger changed")
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "replay_steps": replay_rows,
        "replay_telescoped_gain": {
            "lower": replay_total.lower().str(20),
            "upper": replay_total.upper().str(20),
        },
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
