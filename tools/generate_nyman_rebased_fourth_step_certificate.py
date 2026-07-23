"""Certify one finite fourth step of the rebased Nyman construction.

The fourth shell is reconstructed from the exact published ``p3``.  Only its
first multiplier is retained, so the new basis has 33 coordinates and support
through 65,536 rather than the eight-multiplier exploratory support through
524,288.  The optimizer proposes one exact-dyadic vector; independent prefix
and complete-tail enclosure decides the sign.  This is one finite contraction
only, and RH remains unresolved.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
from functools import wraps
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any, Mapping

import flint
from flint import ctx
import numpy as np


if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
from tools import generate_nyman_rebased_third_step_certificate as third


SCHEMA = "rh-lab/nyman-rebased-fourth-step-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-rebased-fourth-step-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "a807f63928869c6451fa397edcc294ff6b3553b69a54099e8053722bd2b594b0"
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

CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 100_000)
SHELL_BITS = third.SHELL_BITS
WEIGHT_BITS = third.WEIGHT_BITS
SHELL_4_MULTIPLIER_LIMIT = 1
CUTOFF = 1 << 26
GENERATION_BLOCK_SIZE = 1 << 20
GENERATION_PRECISION_BITS = 512
REPLAY_BLOCK_SIZE = 250_003
REPLAY_PRECISION_BITS = 640
MAX_INTEGER_DECIMAL_DIGITS = 100_000

WEIGHT_4_NUMERATORS = (
    -62_704,
    60_496,
    60_567,
    7_832,
    49_590,
    -32_697,
    40_718,
    3_700,
    52_259,
    -30_480,
    -31_835,
    931,
    -17_330,
    19_824,
    -15_002,
    847,
    41_140,
    -28_229,
    -33_848,
    -3_095,
    -14_364,
    30_892,
    -10_955,
    5_776,
    35_212,
    -28_199,
    -28_837,
    949,
    -11_048,
    28_333,
    -8_579,
    -710,
    2_390,
)
WEIGHTS_4 = tuple(
    Fraction(numerator, 1 << WEIGHT_BITS)
    for numerator in WEIGHT_4_NUMERATORS
)

SHELL_4_PAYLOAD_SHA256 = (
    "05d1aa0b6a1f7d2f8b8e303f6a699725f632b4211462331dee6ab3c846b88fad"
)
WEIGHTS_4_PAYLOAD_SHA256 = (
    "eca0b2afabc8ef71888292d21f0d704728957bd740cf532445bebbaf8f7c7436"
)
P4_PAYLOAD_SHA256 = (
    "bcc33e26af682975e7c711b74cfd954f246f0fb8756e3d93de6e02e065f801f1"
)
TAIL_RECORD_SHA256 = (
    "8f3cc3d85fd790886c8b5bd8af7ed43f3d4711ab095f5c433e9c1496b1b70cd1"
)

AS_OF = "2026-07-23"
STATEMENT = (
    "For the explicit exact-dyadic rebased Nyman vectors p3 and p4, where "
    "the fourth shell is reconstructed from exact published p3 and only its "
    "first multiplier is retained, the complete weighted-L2 energy satisfies "
    "E(p3)-E(p4)>1/100000."
)
LIMITATION = (
    "This certifies one finite fourth rebased contraction. The frozen fourth "
    "weight list was selected by an exploratory optimizer. No fifth step, "
    "all-scale recurrence, convergence theorem, proof, or disproof follows. "
    "The Riemann Hypothesis remains unresolved."
)
METHOD = {
    "selection_boundary": (
        "the optimizer proposed one 33-coordinate dyadic vector; exact "
        "reconstruction and independent interval proofs decide the sign"
    ),
    "parent_binding": (
        "pin the published third-step artifact by raw LF, canonical, payload, "
        "source-lineage, and p3-vector commitments"
    ),
    "fourth_shell": (
        "derive the ideal shell on indices 32769..65536 from exact p3, prove "
        "unique nearest 2^-9 bins, and impose final exact balance"
    ),
    "reduced_basis": (
        "reoptimize the 32 inherited direct/shell coordinates and retain only "
        "multiplier one of shell four, avoiding the exploratory K=8 support"
    ),
    "old_lower": (
        "use nonnegativity of the omitted integral to bound complete E(p3) "
        "below by its arbitrary-slope prefix through T=2^26"
    ),
    "new_upper": (
        "bound complete E(p4) above by its prefix at the identical cutoff plus "
        "the exact denominator-weighted local-spacing tail"
    ),
    "decision": (
        "subtract the p4 complete-energy upper bound from the p3 prefix lower "
        "bound and test the exact rational threshold"
    ),
}
REFERENCE = third.REFERENCE
TRUSTED_COMPUTING_BASE = {
    **third.TRUSTED_COMPUTING_BASE,
    "parent_artifact": (
        "published nyman-rebased-third-step-v1 artifact pinned byte-for-byte "
        "and by canonical content before exact p3 reconstruction"
    ),
    "fourth_basis": (
        "exact 33-coordinate basis with one fourth-shell multiplier; no "
        "optimizer objective is accepted as evidence"
    ),
    "old_energy_lower": (
        "nonnegativity of the weighted-L2 integrand outside the certified p3 "
        "prefix at the same cutoff used for p4"
    ),
    "integer_serialization": (
        "CPython decimal integer conversion limit raised to 100000 digits for "
        "the exact local-spacing rational; inputs remain frozen artifacts"
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
    "vectors",
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
    "parent_source_lineage_matches_exact_reconstruction": True,
    "parent_p3_commitment_matches": True,
    "fourth_shell_is_exact_p3_derived_2^-9_rule": True,
    "fourth_shell_is_exactly_balanced": True,
    "fourth_weights_are_exact_2^-16_list": True,
    "p4_has_47345_nonzero_coefficients_through_65536": True,
    "old_p3_has_13537_changed_coefficients": True,
    "update_escapes_frozen_prefix_no_go_hypotheses": True,
    "old_complete_energy_lower_uses_only_p3_prefix": True,
    "new_tail_uses_exact_local_spacing_constant": True,
    "local_spacing_constant_improves_global_spacing_constant": True,
    "optimizer_used_as_evidence": False,
    "complete_gain_strictly_above_claim": True,
    "resolves_rh": False,
}


@dataclass(frozen=True)
class FourthStepWitness:
    parent_witness: third.ThirdStepWitness
    shell4: Mapping[int, Fraction]
    weights4: Mapping[int, Fraction]
    p4: Mapping[int, Fraction]
    source: Mapping[str, object]


def _require_environment() -> None:
    third.parent._require_environment()


@contextmanager
def _integer_serialization_context() -> Any:
    previous_limit = sys.get_int_max_str_digits()
    needs_raise = (
        previous_limit != 0
        and previous_limit < MAX_INTEGER_DECIMAL_DIGITS
    )
    if needs_raise:
        sys.set_int_max_str_digits(MAX_INTEGER_DECIMAL_DIGITS)
    try:
        yield
    finally:
        if needs_raise:
            sys.set_int_max_str_digits(previous_limit)


def _with_integer_serialization_context(function: Any) -> Any:
    @wraps(function)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        with _integer_serialization_context():
            return function(*args, **kwargs)

    return wrapped


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
    ):
        raise ValueError("parent third-step semantics changed")
    vectors = supplied.get("vectors")
    p3_record = vectors.get("p3") if isinstance(vectors, Mapping) else None
    if (
        not isinstance(p3_record, Mapping)
        or p3_record.get("payload_sha256") != PARENT_P3_PAYLOAD_SHA256
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


def _aggregate_reduced_basis(
    weights: tuple[Fraction, ...],
    shells: tuple[Mapping[int, Fraction], ...],
) -> dict[int, Fraction]:
    if len(weights) != 33 or len(shells) != 4:
        raise ValueError("fourth-step reduced basis contract changed")
    result = third.parent._aggregate_basis(weights[:32], shells[:3])
    shell_weight = weights[32]
    for index, value in shells[3].items():
        third.parent._add_coefficient(
            result,
            int(index),
            shell_weight * Fraction(value),
        )
    return third.parent._clean(result)


def _construct_witness(root: Path) -> FourthStepWitness:
    parent_artifact = _load_parent_artifact(root)
    parent_witness = third._construct_witness(root)
    if third.parent._vector_commitment(parent_witness.p3)["payload_sha256"] != (
        PARENT_P3_PAYLOAD_SHA256
    ):
        raise ArithmeticError("exact parent p3 commitment changed")
    if parent_artifact.get("source") != parent_witness.source:
        raise ArithmeticError("parent source lineage differs from reconstruction")
    if parent_artifact.get("vectors") != third._vectors_record(parent_witness):
        raise ArithmeticError("parent vectors differ from reconstruction")
    shell4 = third.parent.rounded_ideal_shell(
        parent_witness.p3,
        SHELL_BITS,
        support_cutoff=32_768,
    )
    weights4 = {
        index: value for index, value in enumerate(WEIGHTS_4, start=1) if value
    }
    p4 = _aggregate_reduced_basis(
        WEIGHTS_4,
        (
            parent_witness.parent_witness.shell1,
            parent_witness.parent_witness.shell2,
            parent_witness.shell3,
            shell4,
        ),
    )
    commitments = {
        "shell4": third.parent._vector_commitment(shell4)["payload_sha256"],
        "weights4": third.parent._vector_commitment(weights4)["payload_sha256"],
        "p4": third.parent._vector_commitment(p4)["payload_sha256"],
    }
    expected = {
        "shell4": SHELL_4_PAYLOAD_SHA256,
        "weights4": WEIGHTS_4_PAYLOAD_SHA256,
        "p4": P4_PAYLOAD_SHA256,
    }
    if commitments != expected:
        raise ArithmeticError("fourth-step vector commitment changed")
    if (
        len(shell4) != 32_685
        or min(shell4) != 32_769
        or max(shell4) != 65_536
        or third.parent._coefficient_sum(shell4)
    ):
        raise ArithmeticError("fourth shell support or balance changed")
    if len(p4) != 47_345 or min(p4) != 1 or max(p4) != 65_536:
        raise ArithmeticError("rebased p4 support contract changed")
    return FourthStepWitness(
        parent_witness=parent_witness,
        shell4=shell4,
        weights4=weights4,
        p4=p4,
        source=_source_record(parent_artifact),
    )


def _vectors_record(witness: FourthStepWitness) -> dict[str, object]:
    p3_commitment = third.parent._vector_commitment(witness.parent_witness.p3)
    return {
        **third._vectors_record(witness.parent_witness),
        "shell_4": third.parent._vector_commitment(witness.shell4),
        "weights_4": third.parent._vector_commitment(witness.weights4),
        "p4": third.parent._vector_commitment(witness.p4),
        "shell_4_parent_p3_payload_sha256": p3_commitment["payload_sha256"],
    }


def _configuration_record() -> dict[str, object]:
    return {
        "cutoff": str(CUTOFF),
        "shell_bits": str(SHELL_BITS),
        "weight_bits": str(WEIGHT_BITS),
        "shell_4_multiplier_limit": str(SHELL_4_MULTIPLIER_LIMIT),
        "generation_block_size": str(GENERATION_BLOCK_SIZE),
        "generation_precision_bits": str(GENERATION_PRECISION_BITS),
        "replay_block_size": str(REPLAY_BLOCK_SIZE),
        "replay_precision_bits": str(REPLAY_PRECISION_BITS),
        "weights_4": [
            third.parent._fraction_record(value) for value in WEIGHTS_4
        ],
    }


def _combine_proof(
    old_prefix: third.parent.AbsoluteEnergyPrefix,
    new_prefix: third.parent.AbsoluteEnergyPrefix,
    tail: third.LocalSpacingTail,
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
    FourthStepWitness,
    third.parent.AbsoluteEnergyPrefix,
    third.parent.AbsoluteEnergyPrefix,
    third.LocalSpacingTail,
    Any,
    Any,
]:
    witness = _construct_witness(root)
    old_prefix = third.parent.absolute_energy_prefix(
        witness.parent_witness.p3,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    new_prefix = third.parent.absolute_energy_prefix(
        witness.p4,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    tail = third.local_spacing_tail_upper(witness.p4, CUTOFF)
    if content_sha256(third._tail_record(tail)) != TAIL_RECORD_SHA256:
        raise ArithmeticError("fourth-step exact tail record changed")
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
    tail: third.LocalSpacingTail,
    new_upper: Any,
    gain: Any,
) -> dict[str, object]:
    return {
        "old_complete_energy_lower_from_prefix": {
            "relation": "E(p3)>=E_prefix(p3;T) by nonnegative omitted square integral",
            "prefix": third.parent._prefix_record(old_prefix),
            "prefix_lower_endpoint": arb_record(old_prefix.prefix_energy.lower()),
        },
        "new_prefix": third.parent._prefix_record(new_prefix),
        "new_tail_upper": third._tail_record(tail),
        "new_complete_energy_upper_computation": {
            "enclosure": arb_record(new_upper),
        },
        "gain_lower_bound_computation": {
            "enclosure": arb_record(gain),
        },
    }


def _checks(
    witness: FourthStepWitness,
    old_prefix: third.parent.AbsoluteEnergyPrefix,
    new_prefix: third.parent.AbsoluteEnergyPrefix,
    tail: third.LocalSpacingTail,
    gain: Any,
) -> dict[str, bool]:
    delta = {
        index: witness.p4.get(index, Fraction())
        - witness.parent_witness.p3.get(index, Fraction())
        for index in set(witness.p4) | set(witness.parent_witness.p3)
    }
    delta = third.parent._clean(delta)
    return {
        "parent_artifact_commitments_match": True,
        "parent_source_lineage_matches_exact_reconstruction": True,
        "parent_p3_commitment_matches": (
            third.parent._vector_commitment(witness.parent_witness.p3)[
                "payload_sha256"
            ]
            == PARENT_P3_PAYLOAD_SHA256
        ),
        "fourth_shell_is_exact_p3_derived_2^-9_rule": (
            witness.shell4
            == third.parent.rounded_ideal_shell(
                witness.parent_witness.p3,
                SHELL_BITS,
                support_cutoff=32_768,
            )
        ),
        "fourth_shell_is_exactly_balanced": (
            third.parent._coefficient_sum(witness.shell4) == 0
        ),
        "fourth_weights_are_exact_2^-16_list": (
            tuple(witness.weights4[index] for index in range(1, 34))
            == WEIGHTS_4
        ),
        "p4_has_47345_nonzero_coefficients_through_65536": (
            len(witness.p4) == 47_345 and max(witness.p4) == 65_536
        ),
        "old_p3_has_13537_changed_coefficients": (
            sum(
                witness.p4.get(index, Fraction()) != value
                for index, value in witness.parent_witness.p3.items()
            )
            == 13_537
        ),
        "update_escapes_frozen_prefix_no_go_hypotheses": (
            third.parent._coefficient_sum(delta) != 0
            and third.parent._harmonic_sum(delta) != 0
        ),
        "old_complete_energy_lower_uses_only_p3_prefix": (
            old_prefix.cutoff == CUTOFF
        ),
        "new_tail_uses_exact_local_spacing_constant": (
            tail.local_spacing_constant
            == Fraction(3, 2) * tail.support_limit * tail.sigma
            and new_prefix.cutoff == tail.cutoff == CUTOFF
        ),
        "local_spacing_constant_improves_global_spacing_constant": (
            tail.local_spacing_constant < tail.global_spacing_constant
        ),
        "optimizer_used_as_evidence": False,
        "complete_gain_strictly_above_claim": (
            gain.lower()
            > third.parent._arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND)
        ),
        "resolves_rh": False,
    }


@_with_integer_serialization_context
def generate(root: Path) -> dict[str, object]:
    _require_environment()
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
        raise ArithmeticError("fourth-step certificate checks did not close")
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
        "source": components[0].source,
        "vectors": _vectors_record(components[0]),
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
        raise ValueError("fourth-step check ledger fields changed")
    if any(type(value) is not bool for value in stored.values()):
        raise ValueError("fourth-step check ledger values are not booleans")


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


@_with_integer_serialization_context
def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    supplied = third.parent.strict_json_loads(
        artifact_path.read_text(encoding="utf-8")
    )
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("fourth-step certificate top-level fields changed")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("fourth-step certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("fourth-step certificate is not the frozen v1 payload")
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
            raise ValueError(f"fourth-step certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("fourth-step certificate configuration changed")
    if supplied.get("claimed_gain_lower_bound") != third.parent._fraction_record(
        CLAIMED_GAIN_LOWER_BOUND
    ):
        raise ValueError("fourth-step claimed gain changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("fourth-step method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("fourth-step trusted-computing-base changed")

    parent_artifact = _load_parent_artifact(root)
    if supplied.get("source") != _source_record(parent_artifact):
        raise ValueError("fourth-step parent or source binding changed")
    vectors = supplied.get("vectors")
    if not isinstance(vectors, Mapping):
        raise ValueError("fourth-step vector commitments are not an object")
    _verify_ledger_shape(supplied.get("checks"))

    generation = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    if vectors != _vectors_record(generation[0]):
        raise ValueError("fourth-step vector commitments changed")
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
        raise ArithmeticError("fourth-step generation checks did not close")
    if supplied.get("checks") != generated_checks:
        raise ValueError("fourth-step check ledger changed")

    replay = _generate_components(
        root,
        block_size=REPLAY_BLOCK_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    if _vectors_record(replay[0]) != _vectors_record(generation[0]):
        raise ArithmeticError("independent replay reconstructed new vectors")
    replay_checks = _checks(
        replay[0],
        replay[1],
        replay[2],
        replay[3],
        replay[5],
    )
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("independent fourth-step replay did not close")
    _require_overlap(generation, replay)
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "replay_old_prefix_energy": {
            "lower": replay[1].prefix_energy.lower().str(30),
            "upper": replay[1].prefix_energy.upper().str(30),
        },
        "replay_new_complete_energy_upper_computation": {
            "lower": replay[4].lower().str(30),
            "upper": replay[4].upper().str(30),
        },
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
