"""Certify a fixed-modulus refinement of the third rebased Nyman step.

This certificate keeps the exact published ``p2`` and ``p3`` vectors and the
same finite prefix.  Its only change is a sharper exact-rational upper bound
for the denominator-spacing tail constant, obtained by grouping active
Fourier denominators through ``gcd(d, 840)``.  The resulting finite claim is
neither a proof nor a disproof of RH; the hypothesis remains unresolved.
"""

from __future__ import annotations

import argparse
from dataclasses import fields
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any, Mapping

import flint
from flint import arb, ctx
import numpy as np


if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
from tools import generate_nyman_rebased_third_step_certificate as third
from tools.certify_nyman_fixed_modulus_tail import (
    FixedModulusTail,
    fixed_modulus_tail_upper,
)


SCHEMA = "rh-lab/nyman-fixed-modulus-refinement-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-fixed-modulus-refinement-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "de69623dee4550796471850d62bf78ec3162fdec34620526f21d0ac578f35f92"
)

PARENT_ARTIFACT_PATH = Path("results/nyman-rebased-third-step-v1.json")
PARENT_ARTIFACT_RAW_LF_SHA256 = (
    "e5bd17bf94f2ce95611ad09bf75acaecc156db003ce4481b4a6e41c3504340bd"
)
PARENT_ARTIFACT_CANONICAL_SHA256 = (
    "811ebaf03089288c99d8155171a07e4b5e5731d53678265354ed2738157fa1ca"
)
PARENT_ARTIFACT_PAYLOAD_SHA256 = (
    "bbafc2c92acb954c9c9e550f781612965bc5c76fe59ebd4a1adc5b596c3bee02"
)
PARENT_P3_PAYLOAD_SHA256 = third.P3_PAYLOAD_SHA256

CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 17_500)
MODULUS = 840
CUTOFF = third.CUTOFF
GENERATION_BLOCK_SIZE = third.GENERATION_BLOCK_SIZE
GENERATION_PRECISION_BITS = third.GENERATION_PRECISION_BITS
REPLAY_BLOCK_SIZE = third.REPLAY_BLOCK_SIZE
REPLAY_PRECISION_BITS = third.REPLAY_PRECISION_BITS

AS_OF = "2026-07-23"
STATEMENT = (
    "For the same explicit exact-dyadic rebased Nyman vectors p2 and p3 "
    "published in nyman-rebased-third-step-v1, the complete weighted-L2 "
    "energy satisfies E(p2)-E(p3)>1/17500."
)
LIMITATION = (
    "This sharpens one finite p2-to-p3 contraction by a grouped local-spacing "
    "tail bound. It supplies no fourth step, all-scale recurrence, convergence "
    "theorem, proof, or disproof. The Riemann Hypothesis remains unresolved."
)
METHOD = {
    "parent_binding": (
        "pin the published third-step artifact by raw LF, canonical, payload, "
        "schema, semantic status, and exact p3-vector commitments"
    ),
    "active_denominators": (
        "compute P_d=sum_(d|n)p3_n/n exactly and retain precisely the "
        "denominators with P_d nonzero"
    ),
    "farey_gap": (
        "for distinct reduced a/d and b/e use circular gap at least "
        "gcd(d,e)/(d*e)"
    ),
    "fixed_modulus_grouping": (
        "put m=gcd(d,840) and R_m=max_active_e e/gcd(m,e), so the reciprocal "
        "local gap at denominator d is at most d*R_m"
    ),
    "weighted_cosecant": (
        "apply Montgomery-Vaughan's 3/2 weighted periodic cosecant theorem "
        "and aggregate exact squared Fourier mass by denominator"
    ),
    "decision": (
        "subtract the fixed-modulus p3 prefix-plus-tail upper bound from the "
        "nonnegative-tail p2 prefix lower bound at the same cutoff"
    ),
    "selection_boundary": (
        "the published optimizer proposal is inherited unchanged and is not "
        "evidence; exact reconstruction and interval arithmetic decide the sign"
    ),
}
REFERENCE = third.REFERENCE
TRUSTED_COMPUTING_BASE = {
    **third.TRUSTED_COMPUTING_BASE,
    "parent_artifact": (
        "published nyman-rebased-third-step-v1 artifact pinned by raw LF and "
        "canonical content before exact reconstruction"
    ),
    "fixed_modulus_tail": (
        "exact divisor sums, active-denominator filtering, gcd classes modulo "
        "840, Jordan J2, and Montgomery-Vaughan's published 3/2 theorem"
    ),
}

TOP_LEVEL_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "as_of",
    "statement",
    "configuration",
    "claimed_gain_lower_bound",
    "source",
    "proof",
    "method",
    "reference",
    "trusted_computing_base",
    "checks",
    "limitation",
    "payload_sha256",
}
EXPECTED_CHECKS = {
    "parent_artifact_commitments_match": True,
    "parent_semantics_remain_certified_finite_and_unresolved": True,
    "parent_p3_commitment_matches_exact_reconstruction": True,
    "active_denominator_count_is_17157": True,
    "modulus_840_has_32_active_gcd_classes": True,
    "fixed_modulus_constant_strictly_improves_active_q_constant": True,
    "old_complete_energy_lower_uses_only_p2_prefix": True,
    "new_complete_energy_upper_uses_same_cutoff": True,
    "optimizer_used_as_evidence": False,
    "complete_gain_strictly_above_claim": True,
    "resolves_rh": False,
}


def _load_parent_artifact(root: Path) -> Mapping[str, object]:
    path = root / PARENT_ARTIFACT_PATH
    raw = path.read_bytes()
    if b"\r" in raw or not raw.endswith(b"\n"):
        raise ValueError("parent third-step artifact is not LF-only")
    if hashlib.sha256(raw).hexdigest() != PARENT_ARTIFACT_RAW_LF_SHA256:
        raise ValueError("parent third-step raw LF SHA-256 mismatch")
    supplied = third.parent.strict_json_loads(raw.decode("utf-8"))
    if not isinstance(supplied, Mapping) or set(supplied) != third.TOP_LEVEL_FIELDS:
        raise ValueError("parent third-step artifact fields changed")
    if content_sha256(supplied) != PARENT_ARTIFACT_CANONICAL_SHA256:
        raise ValueError("parent third-step canonical SHA-256 mismatch")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("parent third-step body SHA-256 mismatch")
    if supplied.get("payload_sha256") != PARENT_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("parent third-step payload commitment changed")
    if (
        supplied.get("schema") != third.SCHEMA
        or supplied.get("classification") != "CERTIFIED_FINITE"
        or supplied.get("hypothesis_status") != "UNRESOLVED"
        or supplied.get("claimed_gain_lower_bound")
        != third.parent._fraction_record(third.CLAIMED_GAIN_LOWER_BOUND)
    ):
        raise ValueError("parent third-step semantics changed")
    vectors = supplied.get("vectors")
    p3_record = vectors.get("p3") if isinstance(vectors, Mapping) else None
    if (
        not isinstance(vectors, Mapping)
        or not isinstance(p3_record, Mapping)
        or p3_record.get("payload_sha256")
        != PARENT_P3_PAYLOAD_SHA256
    ):
        raise ValueError("parent third-step p3 commitment changed")
    return supplied


def _source_record(parent_artifact: Mapping[str, object]) -> dict[str, object]:
    return {
        "parent_artifact": {
            "path": PARENT_ARTIFACT_PATH.as_posix(),
            "raw_lf_sha256": PARENT_ARTIFACT_RAW_LF_SHA256,
            "canonical_sha256": PARENT_ARTIFACT_CANONICAL_SHA256,
            "payload_sha256": PARENT_ARTIFACT_PAYLOAD_SHA256,
            "schema": third.SCHEMA,
            "p3_payload_sha256": PARENT_P3_PAYLOAD_SHA256,
        },
        "inherited_source": parent_artifact["source"],
    }


def _configuration_record() -> dict[str, str]:
    return {
        "cutoff": str(CUTOFF),
        "fixed_modulus": str(MODULUS),
        "generation_block_size": str(GENERATION_BLOCK_SIZE),
        "generation_precision_bits": str(GENERATION_PRECISION_BITS),
        "replay_block_size": str(REPLAY_BLOCK_SIZE),
        "replay_precision_bits": str(REPLAY_PRECISION_BITS),
    }


def _tail_record(tail: FixedModulusTail) -> dict[str, object]:
    integer_fields = {
        "cutoff",
        "modulus",
        "support_limit",
        "active_denominator_count",
        "active_denominator_maximum",
    }
    result: dict[str, object] = {}
    for field in fields(FixedModulusTail):
        value = getattr(tail, field.name)
        if field.name in integer_fields:
            result[field.name] = str(value)
        elif field.name == "group_reciprocal_bounds":
            result[field.name] = [
                {
                    "gcd_class": str(group),
                    "reciprocal_factor": str(reciprocal),
                }
                for group, reciprocal in value
            ]
        else:
            result[field.name] = third.parent._fraction_record(value)
    return result


def _combine_proof(
    old_prefix: third.parent.AbsoluteEnergyPrefix,
    new_prefix: third.parent.AbsoluteEnergyPrefix,
    tail: FixedModulusTail,
    precision_bits: int,
) -> tuple[Any, Any]:
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        new_upper = (
            new_prefix.prefix_energy
            + third.parent._arb_from_fraction(tail.upper_bound)
        )
        gain = old_prefix.prefix_energy - new_upper
    finally:
        ctx.prec = previous_precision
    return new_upper, gain


def _generate_components(
    root: Path,
    *,
    block_size: int,
    precision_bits: int,
) -> tuple[
    third.ThirdStepWitness,
    third.parent.AbsoluteEnergyPrefix,
    third.parent.AbsoluteEnergyPrefix,
    FixedModulusTail,
    Any,
    Any,
]:
    witness = third._construct_witness(root)
    old_prefix = third.parent.absolute_energy_prefix(
        witness.parent_witness.p2,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    new_prefix = third.parent.absolute_energy_prefix(
        witness.p3,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    tail = fixed_modulus_tail_upper(witness.p3, CUTOFF, MODULUS)
    new_upper, gain = _combine_proof(
        old_prefix,
        new_prefix,
        tail,
        precision_bits,
    )
    return witness, old_prefix, new_prefix, tail, new_upper, gain


def _proof_record(
    old_prefix: third.parent.AbsoluteEnergyPrefix,
    new_prefix: third.parent.AbsoluteEnergyPrefix,
    tail: FixedModulusTail,
    new_upper: Any,
    gain: Any,
) -> dict[str, object]:
    return {
        "old_complete_energy_lower_from_prefix": {
            "relation": "E(p2)>=E_prefix(p2;T) by nonnegative omitted square integral",
            "prefix": third.parent._prefix_record(old_prefix),
            "prefix_lower_endpoint": arb_record(old_prefix.prefix_energy.lower()),
        },
        "new_prefix": third.parent._prefix_record(new_prefix),
        "new_fixed_modulus_tail_upper": _tail_record(tail),
        "new_complete_energy_upper_computation": {
            "enclosure": arb_record(new_upper),
        },
        "gain_lower_bound_computation": {
            "enclosure": arb_record(gain),
        },
    }


def _checks(
    witness: third.ThirdStepWitness,
    old_prefix: third.parent.AbsoluteEnergyPrefix,
    new_prefix: third.parent.AbsoluteEnergyPrefix,
    tail: FixedModulusTail,
    gain: Any,
) -> dict[str, bool]:
    return {
        "parent_artifact_commitments_match": True,
        "parent_semantics_remain_certified_finite_and_unresolved": True,
        "parent_p3_commitment_matches_exact_reconstruction": (
            third.parent._vector_commitment(witness.p3)["payload_sha256"]
            == PARENT_P3_PAYLOAD_SHA256
        ),
        "active_denominator_count_is_17157": (
            tail.active_denominator_count == 17_157
        ),
        "modulus_840_has_32_active_gcd_classes": (
            tail.modulus == MODULUS
            and len(tail.group_reciprocal_bounds) == 32
        ),
        "fixed_modulus_constant_strictly_improves_active_q_constant": (
            tail.fixed_modulus_spacing_constant
            < tail.active_q_spacing_constant
        ),
        "old_complete_energy_lower_uses_only_p2_prefix": (
            old_prefix.cutoff == CUTOFF
        ),
        "new_complete_energy_upper_uses_same_cutoff": (
            new_prefix.cutoff == CUTOFF and tail.cutoff == CUTOFF
        ),
        "optimizer_used_as_evidence": False,
        "complete_gain_strictly_above_claim": (
            gain.lower()
            > third.parent._arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND)
        ),
        "resolves_rh": False,
    }


def generate(root: Path) -> dict[str, object]:
    third.parent._require_environment()
    parent_artifact = _load_parent_artifact(root)
    components = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    checks = _checks(
        components[0],
        components[1],
        components[2],
        components[3],
        components[5],
    )
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("fixed-modulus refinement checks did not close")
    body = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": _configuration_record(),
        "claimed_gain_lower_bound": third.parent._fraction_record(
            CLAIMED_GAIN_LOWER_BOUND
        ),
        "source": _source_record(parent_artifact),
        "proof": _proof_record(*components[1:]),
        "method": METHOD,
        "reference": REFERENCE,
        "trusted_computing_base": TRUSTED_COMPUTING_BASE,
        "checks": checks,
        "limitation": LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _verify_ledger_shape(stored: object) -> None:
    if not isinstance(stored, Mapping) or set(stored) != set(EXPECTED_CHECKS):
        raise ValueError("fixed-modulus check ledger fields changed")
    if any(type(value) is not bool for value in stored.values()):
        raise ValueError("fixed-modulus check ledger values are not booleans")


def _require_overlap(
    generation: tuple[Any, ...],
    replay: tuple[Any, ...],
) -> None:
    pairs = (
        ("old prefix", generation[1].prefix_energy, replay[1].prefix_energy),
        ("new prefix", generation[2].prefix_energy, replay[2].prefix_energy),
        ("new complete upper", generation[4], replay[4]),
        ("gain", generation[5], replay[5]),
    )
    for label, first, second in pairs:
        if not first.overlaps(second):
            raise ArithmeticError(
                f"generation and independent replay disagree on {label}"
            )


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    third.parent._require_environment()
    supplied = third.parent.strict_json_loads(
        artifact_path.read_text(encoding="utf-8")
    )
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("fixed-modulus certificate top-level fields changed")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("fixed-modulus certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("fixed-modulus certificate is not the frozen v1 payload")
    scalars = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "limitation": LIMITATION,
    }
    for name, expected in scalars.items():
        if supplied.get(name) != expected:
            raise ValueError(f"fixed-modulus certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("fixed-modulus certificate configuration changed")
    if supplied.get("claimed_gain_lower_bound") != third.parent._fraction_record(
        CLAIMED_GAIN_LOWER_BOUND
    ):
        raise ValueError("fixed-modulus claimed gain changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("fixed-modulus method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("fixed-modulus trusted-computing-base changed")

    parent_artifact = _load_parent_artifact(root)
    if supplied.get("source") != _source_record(parent_artifact):
        raise ValueError("fixed-modulus parent or source binding changed")
    _verify_ledger_shape(supplied.get("checks"))

    generation = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    generated_proof = _proof_record(*generation[1:])
    third._verify_generated_proof(supplied.get("proof"), generated_proof)
    generated_checks = _checks(
        generation[0],
        generation[1],
        generation[2],
        generation[3],
        generation[5],
    )
    if generated_checks != EXPECTED_CHECKS:
        raise ArithmeticError("fixed-modulus generation checks did not close")
    if supplied.get("checks") != generated_checks:
        raise ValueError("fixed-modulus check ledger changed")

    replay = _generate_components(
        root,
        block_size=REPLAY_BLOCK_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    replay_checks = _checks(
        replay[0],
        replay[1],
        replay[2],
        replay[3],
        replay[5],
    )
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("independent fixed-modulus replay did not close")
    _require_overlap(generation, replay)
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "fixed_modulus_spacing_constant": third.parent._fraction_record(
            replay[3].fixed_modulus_spacing_constant
        ),
        "replay_gain_lower_bound_computation": {
            "lower": replay[5].lower().str(30),
            "upper": replay[5].upper().str(30),
        },
        "replay_environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "python_flint": flint.__version__,
            "numpy": np.__version__,
            "block_size": str(REPLAY_BLOCK_SIZE),
            "precision_bits": str(REPLAY_PRECISION_BITS),
        },
        "verified": True,
    }


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
        result = verify(root, artifact_path)
    else:
        result = generate(root)
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        artifact_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
