"""Certify one finite third step of the rebased Nyman construction.

The frozen optimizer output is proposal-only.  This module binds to the
published two-step certificate, reconstructs its exact p2 vector, derives a
third rounded shell, and proves a complete-energy decrease without reusing an
optimizer objective value.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
from typing import Any, Mapping, Sequence

import flint
from flint import ctx
import numpy as np

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record

if __package__:
    from tools import generate_nyman_rebased_schur_certificate as parent
else:
    import generate_nyman_rebased_schur_certificate as parent


SCHEMA = "rh-lab/nyman-rebased-third-step-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-rebased-third-step-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "bbafc2c92acb954c9c9e550f781612965bc5c76fe59ebd4a1adc5b596c3bee02"
)

PARENT_ARTIFACT_PATH = Path("results/nyman-rebased-schur-v1.json")
PARENT_ARTIFACT_RAW_LF_SHA256 = (
    "a2e68b81dc37a252b35070f73579cea73705a063060c4dab517d42e2738ef0aa"
)
PARENT_ARTIFACT_CANONICAL_SHA256 = (
    "f26192ff4f2646cd74dc3f2e55582580618641722e2576ad0edd4a634891671c"
)
PARENT_ARTIFACT_PAYLOAD_SHA256 = (
    "a6e1119fb457e18b521e70203483101a77be7509d52d4f41889951399ce95edd"
)
PARENT_P2_PAYLOAD_SHA256 = (
    "b392a0c8d9ea2bf09b5fdf98575906bb8c17b75044c4a01586bdc12fa299520a"
)

CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 20_000)
DIRECT_DIMENSION = parent.DIRECT_DIMENSION
MULTIPLIER_LIMIT = parent.MULTIPLIER_LIMIT
SHELL_BITS = parent.SHELL_BITS
WEIGHT_BITS = parent.WEIGHT_BITS
CUTOFF = 1 << 23
GENERATION_BLOCK_SIZE = 65_536
GENERATION_PRECISION_BITS = 320
REPLAY_BLOCK_SIZE = 250_003
REPLAY_PRECISION_BITS = 512

WEIGHT_3_NUMERATORS = (
    -62_698,
    60_485,
    60_554,
    7_849,
    49_569,
    -32_652,
    40_703,
    3_693,
    52_082,
    -30_277,
    -31_529,
    928,
    -17_135,
    19_550,
    -14_768,
    831,
    39_574,
    -26_774,
    -31_574,
    -2_862,
    -13_563,
    28_601,
    -9_984,
    4_575,
    23_635,
    -16_402,
    -14_389,
    949,
    -4_871,
    11_964,
    -3_181,
    19,
)
WEIGHTS_3 = tuple(
    Fraction(numerator, 1 << WEIGHT_BITS)
    for numerator in WEIGHT_3_NUMERATORS
)

SHELL_3_PAYLOAD_SHA256 = (
    "290a8f9dcb6d9d5731a24d97b2e6848ca4a0778bc16b9890e9505ae83f205335"
)
WEIGHTS_3_PAYLOAD_SHA256 = (
    "8df8d30ce08214f66aa746204dbb1dadaae95871878d03058f5d632f7bce418f"
)
P3_PAYLOAD_SHA256 = (
    "4797b563591fd11321c2ff2ef7adbe6cbe6ff35f7ccf7d639eb780ef50ca9c97"
)

AS_OF = "2026-07-23"
STATEMENT = (
    "For the explicit exact-dyadic rebased Nyman vectors p2 and p3, where "
    "the third shell is reconstructed from the exact published p2, the "
    "complete weighted-L2 energy satisfies E(p2)-E(p3)>1/20000."
)
LIMITATION = (
    "This certifies one finite third rebased contraction. The frozen third "
    "weight list was selected by an exploratory optimizer, and no all-scale "
    "recurrence or convergence theorem is proved. The Riemann Hypothesis "
    "remains unresolved."
)
METHOD = {
    "selection_boundary": (
        "the optimizer only proposed one dyadic weight list; the verifier "
        "reconstructs every vector and proves the inequality afresh"
    ),
    "parent_binding": (
        "strictly pin the published rebased-v1 artifact by raw LF, canonical, "
        "payload, source-lineage, and p2-vector commitments"
    ),
    "third_shell": (
        "derive the ideal shell on indices 2049..4096 from the exact p2, "
        "prove unique nearest 2^-9 bins, and impose final exact balance"
    ),
    "p3": (
        "reoptimize the eight direct coordinates and eight dilates of each "
        "of shells one, two, and three using the frozen 2^-16 weight list"
    ),
    "old_lower": (
        "use nonnegativity of the omitted integral to bound complete E(p2) "
        "below by its arbitrary-slope prefix over (0,1) and M=1..2^23"
    ),
    "new_upper": (
        "bound complete E(p3) above by the same finite prefix plus the exact "
        "denominator-weighted local-spacing tail starting at T+1"
    ),
    "local_spacing": (
        "use sigma=c0^2+(1/12)sum_d d*J2(d)*P_d^2 and the periodic "
        "Montgomery-Vaughan inequality to obtain C_loc=(3/2)Q*sigma"
    ),
    "decision": (
        "subtract the certified p3 complete-energy upper bound from the p2 "
        "prefix enclosure and test the resulting lower endpoint"
    ),
}
REFERENCE = {
    "authors": "H. L. Montgomery and R. C. Vaughan",
    "title": "Hilbert's inequality",
    "journal": "Journal of the London Mathematical Society (2) 8",
    "year": "1974",
    "pages": "73-82",
    "doi": "10.1112/jlms/s2-8.1.73",
    "url": "https://doi.org/10.1112/jlms/s2-8.1.73",
    "used_result": (
        "Theorem 1 in its periodic weighted cosecant form, applied after "
        "grouping Farey frequencies by denominator."
    ),
}
TRUSTED_COMPUTING_BASE = {
    **parent.TRUSTED_COMPUTING_BASE,
    "parent_artifact": (
        "published nyman-rebased-schur-v1 artifact, pinned byte-for-byte and "
        "by canonical content before exact p2 reconstruction"
    ),
    "old_energy_lower": (
        "nonnegativity of the weighted-L2 integrand outside the certified "
        "p2 prefix"
    ),
    "local_spacing_tail": (
        "exact denominator grouping with "
        "sigma=c0^2+(1/12)sum d*J2(d)*P_d^2 and C_loc=(3/2)Q*sigma"
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
VECTOR_NAMES = parent.VECTOR_NAMES | {
    "shell_3",
    "weights_3",
    "p3",
    "shell_3_parent_p2_payload_sha256",
}
EXPECTED_CHECKS = {
    "parent_artifact_commitments_match": True,
    "parent_source_lineage_matches_exact_reconstruction": True,
    "parent_p2_commitment_matches": True,
    "third_shell_is_exact_p2_derived_2^-9_rule": True,
    "third_shell_is_exactly_balanced": True,
    "third_weights_are_exact_2^-16_list": True,
    "p3_has_14660_nonzero_coefficients_through_32768": True,
    "old_complete_energy_lower_uses_only_p2_prefix": True,
    "new_tail_uses_exact_local_spacing_constant": True,
    "local_spacing_constant_improves_global_spacing_constant": True,
    "optimizer_used_as_evidence": False,
    "complete_gain_strictly_above_claim": True,
    "resolves_rh": False,
}


@dataclass(frozen=True)
class ThirdStepWitness:
    parent_witness: parent.RebasedWitness
    shell3: Mapping[int, Fraction]
    weights3: Mapping[int, Fraction]
    p3: Mapping[int, Fraction]
    source: Mapping[str, object]


@dataclass(frozen=True)
class LocalSpacingTail:
    cutoff: int
    support_limit: int
    coefficient_sum: Fraction
    harmonic_sum: Fraction
    c0: Fraction
    rho: Fraction
    sigma: Fraction
    local_spacing_constant: Fraction
    global_spacing_constant: Fraction
    global_to_local_improvement: Fraction
    absolute_residual_bound: Fraction
    mean_term: Fraction
    discrepancy_term: Fraction
    cross_term: Fraction
    slope_square_term: Fraction
    upper_bound: Fraction


def _load_parent_artifact(root: Path) -> Mapping[str, object]:
    path = root / PARENT_ARTIFACT_PATH
    raw = path.read_bytes()
    if b"\r" in raw:
        raise ValueError("parent rebased artifact is not LF-only")
    if hashlib.sha256(raw).hexdigest() != PARENT_ARTIFACT_RAW_LF_SHA256:
        raise ValueError("parent rebased artifact raw LF SHA-256 mismatch")
    supplied = parent.strict_json_loads(raw.decode("utf-8"))
    if (
        not isinstance(supplied, Mapping)
        or set(supplied) != parent.TOP_LEVEL_FIELDS
    ):
        raise ValueError("parent rebased artifact fields changed")
    if content_sha256(supplied) != PARENT_ARTIFACT_CANONICAL_SHA256:
        raise ValueError("parent rebased artifact canonical SHA-256 mismatch")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if content_sha256(body) != PARENT_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("parent rebased artifact body SHA-256 mismatch")
    if supplied.get("payload_sha256") != PARENT_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("parent rebased artifact payload commitment changed")
    if (
        supplied.get("schema") != parent.SCHEMA
        or supplied.get("classification") != "CERTIFIED_FINITE"
        or supplied.get("hypothesis_status") != "UNRESOLVED"
    ):
        raise ValueError("parent rebased artifact semantics changed")
    vectors = supplied.get("vectors")
    if (
        not isinstance(vectors, Mapping)
        or set(vectors) != parent.VECTOR_NAMES
    ):
        raise ValueError("parent rebased artifact vector fields changed")
    p2 = vectors.get("p2")
    if (
        not isinstance(p2, Mapping)
        or p2.get("payload_sha256") != PARENT_P2_PAYLOAD_SHA256
    ):
        raise ValueError("parent rebased artifact p2 commitment changed")
    return supplied


def _source_record(
    parent_artifact: Mapping[str, object],
) -> dict[str, object]:
    return {
        "parent_artifact": {
            "path": PARENT_ARTIFACT_PATH.as_posix(),
            "raw_lf_sha256": PARENT_ARTIFACT_RAW_LF_SHA256,
            "canonical_sha256": PARENT_ARTIFACT_CANONICAL_SHA256,
            "payload_sha256": PARENT_ARTIFACT_PAYLOAD_SHA256,
            "schema": parent.SCHEMA,
            "p2_payload_sha256": PARENT_P2_PAYLOAD_SHA256,
        },
        "inherited_natural_source": parent_artifact["source"],
    }


def _construct_witness(root: Path) -> ThirdStepWitness:
    parent_artifact = _load_parent_artifact(root)
    parent_witness = parent._construct_witness(root)
    if parent_artifact.get("source") != parent_witness.source:
        raise ValueError("parent source lineage differs from exact reconstruction")
    if parent_artifact.get("vectors") != parent._vectors_record(parent_witness):
        raise ValueError("parent vectors differ from exact reconstruction")
    if (
        parent._vector_commitment(parent_witness.p2)["payload_sha256"]
        != PARENT_P2_PAYLOAD_SHA256
    ):
        raise ArithmeticError("reconstructed parent p2 commitment changed")

    shell3 = parent.rounded_ideal_shell(
        parent_witness.p2,
        SHELL_BITS,
        support_cutoff=2048,
    )
    weights3 = {
        index: value
        for index, value in enumerate(WEIGHTS_3, start=1)
        if value
    }
    p3 = parent._aggregate_basis(
        WEIGHTS_3,
        (parent_witness.shell1, parent_witness.shell2, shell3),
    )
    if (
        len(shell3) != 2048
        or min(shell3) != 2049
        or max(shell3) != 4096
        or parent._coefficient_sum(shell3)
    ):
        raise ArithmeticError("third-shell support or balance contract changed")
    if len(p3) != 14_660 or min(p3) != 1 or max(p3) != 32_768:
        raise ArithmeticError("rebased p3 support contract changed")
    commitments = {
        "shell3": parent._vector_commitment(shell3)["payload_sha256"],
        "weights3": parent._vector_commitment(weights3)["payload_sha256"],
        "p3": parent._vector_commitment(p3)["payload_sha256"],
    }
    if commitments != {
        "shell3": SHELL_3_PAYLOAD_SHA256,
        "weights3": WEIGHTS_3_PAYLOAD_SHA256,
        "p3": P3_PAYLOAD_SHA256,
    }:
        raise ArithmeticError("third-step vector commitment changed")
    return ThirdStepWitness(
        parent_witness=parent_witness,
        shell3=shell3,
        weights3=weights3,
        p3=p3,
        source=_source_record(parent_artifact),
    )


def _vectors_record(witness: ThirdStepWitness) -> dict[str, object]:
    result = dict(parent._vectors_record(witness.parent_witness))
    result.update(
        {
            "shell_3": parent._vector_commitment(witness.shell3),
            "weights_3": parent._vector_commitment(witness.weights3),
            "p3": parent._vector_commitment(witness.p3),
            "shell_3_parent_p2_payload_sha256": PARENT_P2_PAYLOAD_SHA256,
        }
    )
    return result


def _configuration_record() -> dict[str, object]:
    return {
        "direct_dimension": str(DIRECT_DIMENSION),
        "multiplier_limit": str(MULTIPLIER_LIMIT),
        "shell_bits": str(SHELL_BITS),
        "weight_bits": str(WEIGHT_BITS),
        "cutoff": str(CUTOFF),
        "generation_block_size": str(GENERATION_BLOCK_SIZE),
        "generation_precision_bits": str(GENERATION_PRECISION_BITS),
        "replay_block_size": str(REPLAY_BLOCK_SIZE),
        "replay_precision_bits": str(REPLAY_PRECISION_BITS),
        "weights_3": [
            parent._fraction_record(value) for value in WEIGHTS_3
        ],
    }


def local_spacing_tail_upper(
    coefficients: Mapping[int, Fraction],
    cutoff: int,
) -> LocalSpacingTail:
    """Return the audited denominator-weighted absolute tail bound."""

    if cutoff < 1:
        raise ValueError("tail cutoff must be positive")
    cleaned = parent._clean(coefficients)
    support_limit = max(cleaned, default=0)
    coefficient_sum = parent._coefficient_sum(cleaned)
    harmonic_sum = parent._harmonic_sum(cleaned)
    c0 = Fraction(1) - coefficient_sum / 2
    if support_limit:
        divisor_sums = parent._divisor_harmonic_sums(
            cleaned,
            support_limit,
        )
        jordan = parent._jordan_j2_sieve(support_limit)
        rho_sum = sum(
            (
                jordan[divisor]
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in range(2, support_limit + 1)
            ),
            start=Fraction(),
        )
        sigma_sum = sum(
            (
                divisor
                * jordan[divisor]
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in range(2, support_limit + 1)
            ),
            start=Fraction(),
        )
        rho = c0 * c0 + rho_sum / 12
        sigma = c0 * c0 + sigma_sum / 12
    else:
        rho = c0 * c0
        sigma = c0 * c0
    if rho < 0 or sigma < 0:
        raise ArithmeticError("periodic residual moment became negative")
    local_spacing_constant = (
        Fraction(3, 2) * support_limit * sigma
    )
    global_spacing_constant = (
        support_limit * (support_limit - 1) * rho
        if support_limit >= 2
        else Fraction()
    )
    global_to_local_improvement = (
        global_spacing_constant / local_spacing_constant
        if local_spacing_constant
        else Fraction()
    )
    absolute_residual_bound = abs(c0) + sum(
        (abs(Fraction(value)) for value in cleaned.values()),
        start=Fraction(),
    ) / 2
    mean_term = rho / (cutoff + 1)
    discrepancy_term = local_spacing_constant / (
        (cutoff + 1) * (cutoff + 2)
    )
    cross_term = (
        abs(harmonic_sum) * absolute_residual_bound / (cutoff + 1)
    )
    slope_square_term = (
        harmonic_sum * harmonic_sum / (4 * (cutoff + 1))
    )
    upper_bound = (
        mean_term
        + discrepancy_term
        + cross_term
        + slope_square_term
    )
    return LocalSpacingTail(
        cutoff=cutoff,
        support_limit=support_limit,
        coefficient_sum=coefficient_sum,
        harmonic_sum=harmonic_sum,
        c0=c0,
        rho=rho,
        sigma=sigma,
        local_spacing_constant=local_spacing_constant,
        global_spacing_constant=global_spacing_constant,
        global_to_local_improvement=global_to_local_improvement,
        absolute_residual_bound=absolute_residual_bound,
        mean_term=mean_term,
        discrepancy_term=discrepancy_term,
        cross_term=cross_term,
        slope_square_term=slope_square_term,
        upper_bound=upper_bound,
    )


def _tail_record(tail: LocalSpacingTail) -> dict[str, object]:
    return {
        field.name: (
            str(getattr(tail, field.name))
            if field.name in {"cutoff", "support_limit"}
            else parent._fraction_record(getattr(tail, field.name))
        )
        for field in fields(LocalSpacingTail)
    }


def _combine_proof(
    old_prefix: parent.AbsoluteEnergyPrefix,
    new_prefix: parent.AbsoluteEnergyPrefix,
    new_tail: LocalSpacingTail,
    precision_bits: int,
) -> tuple[Any, Any]:
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        new_complete_upper = (
            new_prefix.prefix_energy
            + parent._arb_from_fraction(new_tail.upper_bound)
        )
        gain_lower_computation = (
            old_prefix.prefix_energy - new_complete_upper
        )
    finally:
        ctx.prec = previous_precision
    return new_complete_upper, gain_lower_computation


def _proof_record(
    old_prefix: parent.AbsoluteEnergyPrefix,
    new_prefix: parent.AbsoluteEnergyPrefix,
    new_tail: LocalSpacingTail,
    new_complete_upper: Any,
    gain_lower_computation: Any,
) -> dict[str, object]:
    return {
        "old_complete_energy_lower_from_prefix": {
            "relation": "E(p2)>=E_prefix(p2;T)",
            "prefix": parent._prefix_record(old_prefix),
            "prefix_lower_endpoint": arb_record(
                old_prefix.prefix_energy.lower()
            ),
        },
        "new_prefix": parent._prefix_record(new_prefix),
        "new_tail_upper": _tail_record(new_tail),
        "new_complete_energy_upper_computation": {
            "enclosure": arb_record(new_complete_upper)
        },
        "gain_lower_bound_computation": {
            "enclosure": arb_record(gain_lower_computation)
        },
    }


def _checks(
    witness: ThirdStepWitness,
    old_prefix: parent.AbsoluteEnergyPrefix,
    new_tail: LocalSpacingTail,
    gain_lower_computation: Any,
) -> dict[str, bool]:
    return {
        "parent_artifact_commitments_match": True,
        "parent_source_lineage_matches_exact_reconstruction": True,
        "parent_p2_commitment_matches": (
            parent._vector_commitment(witness.parent_witness.p2)[
                "payload_sha256"
            ]
            == PARENT_P2_PAYLOAD_SHA256
        ),
        "third_shell_is_exact_p2_derived_2^-9_rule": (
            witness.shell3
            == parent.rounded_ideal_shell(
                witness.parent_witness.p2,
                SHELL_BITS,
                support_cutoff=2048,
            )
        ),
        "third_shell_is_exactly_balanced": (
            parent._coefficient_sum(witness.shell3) == 0
        ),
        "third_weights_are_exact_2^-16_list": tuple(
            witness.weights3[index] for index in range(1, 33)
        )
        == WEIGHTS_3,
        "p3_has_14660_nonzero_coefficients_through_32768": (
            len(witness.p3) == 14_660 and max(witness.p3) == 32_768
        ),
        "old_complete_energy_lower_uses_only_p2_prefix": (
            old_prefix.cutoff == CUTOFF
        ),
        "new_tail_uses_exact_local_spacing_constant": (
            new_tail.cutoff == CUTOFF
            and new_tail.local_spacing_constant
            == Fraction(3, 2) * new_tail.support_limit * new_tail.sigma
        ),
        "local_spacing_constant_improves_global_spacing_constant": (
            new_tail.local_spacing_constant
            < new_tail.global_spacing_constant
        ),
        "optimizer_used_as_evidence": False,
        "complete_gain_strictly_above_claim": parent._strict_lower_proved(
            gain_lower_computation,
            parent._arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND),
        ),
        "resolves_rh": False,
    }


def _generate_components(
    root: Path,
    *,
    block_size: int,
    precision_bits: int,
) -> tuple[
    ThirdStepWitness,
    parent.AbsoluteEnergyPrefix,
    parent.AbsoluteEnergyPrefix,
    LocalSpacingTail,
    Any,
    Any,
]:
    witness = _construct_witness(root)
    old_prefix = parent.absolute_energy_prefix(
        witness.parent_witness.p2,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    new_prefix = parent.absolute_energy_prefix(
        witness.p3,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    new_tail = local_spacing_tail_upper(witness.p3, CUTOFF)
    new_complete_upper, gain_lower_computation = _combine_proof(
        old_prefix,
        new_prefix,
        new_tail,
        precision_bits,
    )
    return (
        witness,
        old_prefix,
        new_prefix,
        new_tail,
        new_complete_upper,
        gain_lower_computation,
    )


def generate(root: Path) -> dict[str, object]:
    parent._require_environment()
    components = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    witness, old_prefix = components[:2]
    checks = _checks(witness, old_prefix, components[3], components[-1])
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("third-step certificate checks did not close")
    body = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": _configuration_record(),
        "claimed_gain_lower_bound": parent._fraction_record(
            CLAIMED_GAIN_LOWER_BOUND
        ),
        "source": witness.source,
        "vectors": _vectors_record(witness),
        "proof": _proof_record(*components[1:]),
        "method": METHOD,
        "reference": REFERENCE,
        "trusted_computing_base": TRUSTED_COMPUTING_BASE,
        "checks": checks,
        "limitation": LIMITATION,
    }
    return {**body, "payload_sha256": content_sha256(body)}


def _verify_generated_proof(
    stored: object,
    generated: Mapping[str, object],
) -> None:
    if stored != generated:
        raise ValueError("stored third-step generation proof record changed")
    if not isinstance(stored, Mapping):
        raise ValueError("stored third-step proof is not an object")
    old = stored["old_complete_energy_lower_from_prefix"]
    if not isinstance(old, Mapping) or set(old) != {
        "relation",
        "prefix",
        "prefix_lower_endpoint",
    }:
        raise ValueError("stored old-prefix lower proof fields changed")
    parent._ball_from_record(
        old["prefix_lower_endpoint"],
        "old prefix lower endpoint",
    )
    old_prefix = old["prefix"]
    if not isinstance(old_prefix, Mapping):
        raise ValueError("stored old prefix is not an object")
    new_prefix = stored["new_prefix"]
    if not isinstance(new_prefix, Mapping):
        raise ValueError("stored new prefix is not an object")
    for label, record in (
        ("old", old_prefix),
        ("new", new_prefix),
    ):
        for name in (
            "logarithmic_sum",
            "logarithmic_cross_term",
            "prefix_energy",
        ):
            parent._ball_from_record(
                record[name],
                f"{label} prefix {name}",
            )
    for name in (
        "new_complete_energy_upper_computation",
        "gain_lower_bound_computation",
    ):
        record = stored[name]
        if not isinstance(record, Mapping) or set(record) != {"enclosure"}:
            raise ValueError(f"{name} fields changed")
        parent._ball_from_record(record["enclosure"], name)


def _verify_check_ledger_shape(stored: object) -> None:
    if not isinstance(stored, Mapping) or set(stored) != set(EXPECTED_CHECKS):
        raise ValueError("third-step check ledger fields changed")
    if any(type(value) is not bool for value in stored.values()):
        raise ValueError("third-step check ledger values are not booleans")


def _verify_check_ledger(
    stored: object,
    generated: Mapping[str, bool],
) -> None:
    _verify_check_ledger_shape(stored)
    if stored != generated:
        raise ValueError("third-step check ledger changed")


def _require_replay_overlap(
    generation: tuple[Any, ...],
    replay: tuple[Any, ...],
) -> None:
    pairs = (
        (
            "old prefix energy",
            generation[1].prefix_energy,
            replay[1].prefix_energy,
        ),
        (
            "new prefix energy",
            generation[2].prefix_energy,
            replay[2].prefix_energy,
        ),
        (
            "new complete energy upper computation",
            generation[4],
            replay[4],
        ),
        (
            "gain lower-bound computation",
            generation[5],
            replay[5],
        ),
    )
    for label, generation_ball, replay_ball in pairs:
        if not generation_ball.overlaps(replay_ball):
            raise ArithmeticError(
                f"generation and independent replay disagree on {label}"
            )


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    parent._require_environment()
    supplied = parent.strict_json_loads(
        artifact_path.read_text(encoding="utf-8")
    )
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("third-step certificate top-level fields changed")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("third-step certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("third-step certificate is not the frozen v1 payload")
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
            raise ValueError(f"third-step certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("third-step certificate configuration changed")
    claim = parent._parse_fraction(
        supplied.get("claimed_gain_lower_bound"),
        "third-step claimed gain",
    )
    if claim != CLAIMED_GAIN_LOWER_BOUND:
        raise ValueError("third-step claimed gain changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("third-step method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("third-step trusted-computing-base changed")

    parent_artifact = _load_parent_artifact(root)
    if supplied.get("source") != _source_record(parent_artifact):
        raise ValueError("third-step parent or source binding changed")
    vectors = supplied.get("vectors")
    if not isinstance(vectors, Mapping) or set(vectors) != VECTOR_NAMES:
        raise ValueError("third-step vector commitment fields changed")
    _verify_check_ledger_shape(supplied.get("checks"))

    generation = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    generation_witness = generation[0]
    if vectors != _vectors_record(generation_witness):
        raise ValueError("third-step vector commitments changed")
    generated_proof = _proof_record(*generation[1:])
    _verify_generated_proof(supplied.get("proof"), generated_proof)
    generated_checks = _checks(
        generation_witness,
        generation[1],
        generation[3],
        generation[-1],
    )
    if generated_checks != EXPECTED_CHECKS:
        raise ArithmeticError("third-step generation replay checks did not close")
    _verify_check_ledger(supplied.get("checks"), generated_checks)

    replay = _generate_components(
        root,
        block_size=REPLAY_BLOCK_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    replay_witness = replay[0]
    if _vectors_record(replay_witness) != _vectors_record(generation_witness):
        raise ArithmeticError("independent replay reconstructed new vectors")
    replay_checks = _checks(
        replay_witness,
        replay[1],
        replay[3],
        replay[-1],
    )
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("independent third-step replay did not close")
    _require_replay_overlap(generation, replay)
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "replay_old_prefix_energy": {
            "lower": replay[1].prefix_energy.lower().str(30),
            "upper": replay[1].prefix_energy.upper().str(30),
        },
        "replay_new_complete_energy_upper_computation": {
            "lower": replay[-2].lower().str(30),
            "upper": replay[-2].upper().str(30),
        },
        "replay_gain_lower_bound_computation": {
            "lower": replay[-1].lower().str(30),
            "upper": replay[-1].upper().str(30),
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
