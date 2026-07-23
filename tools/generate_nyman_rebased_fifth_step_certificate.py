"""Certify one finite fifth step of the rebased Nyman construction.

The fifth shell is reconstructed from the exact published ``p4``.  As in the
published fourth step, only the first multiplier of the new shell is retained.
The resulting exact-dyadic candidate has 34 coordinates and support through
131,072.

The prefix cutoff is deliberately large: ``T = 2^30``.  Prefixes are enclosed
by the versioned large-cutoff backend in
``tools.certify_nyman_absolute_prefix``.  That backend streams the exact
integer recurrence and attaches an explicit binary64 forward-error radius to
its sequential term evaluation and pairwise reduction tree.  The omitted p5
tail is bounded by the exact denominator-weighted local-spacing formula.

This program certifies one finite contraction only.  It does not establish an
all-scale recurrence or resolve the Riemann Hypothesis.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass, fields
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
from tools import certify_nyman_absolute_prefix as absolute_prefix
from tools import generate_nyman_rebased_fourth_step_certificate as fourth


SCHEMA = "rh-lab/nyman-rebased-fifth-step-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-rebased-fifth-step-v1.json")

# This is intentionally unset until the expensive generation result has been
# separately audited.  ``verify`` refuses to operate until the exact
# payload commitment is frozen here.
FROZEN_ARTIFACT_PAYLOAD_SHA256: str | None = None

PARENT_ARTIFACT_PATH = Path("results/nyman-rebased-fourth-step-v1.json")
PARENT_ARTIFACT_RAW_LF_SHA256 = (
    "16bd8e0f0fea73599a8a2fcd0cf45bed60163b5fab0c52314c43cefbb1beb560"
)
PARENT_ARTIFACT_CANONICAL_SHA256 = (
    "48b142aa8e21bcc7ac0ec980765f30ca13da831fad54c99a6d9e655db9258f3d"
)
PARENT_ARTIFACT_PAYLOAD_SHA256 = (
    "a807f63928869c6451fa397edcc294ff6b3553b69a54099e8053722bd2b594b0"
)
PARENT_P4_PAYLOAD_SHA256 = fourth.P4_PAYLOAD_SHA256

CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 300_000)
SHELL_BITS = fourth.SHELL_BITS
WEIGHT_BITS = fourth.WEIGHT_BITS
SHELL_5_MULTIPLIER_LIMIT = 1
CUTOFF = 1 << 30

GENERATION_CHUNK_SIZE = 1 << 24
GENERATION_REDUCTION_LEAF_SIZE = 1 << 20
GENERATION_PRECISION_BITS = 768

# The distinct-layout replay changes both partitions and precision.  Chunk
# boundaries affect the exact recurrence schedule; leaf boundaries define a
# genuinely different floating-point reduction tree.  It is a diversity check
# within the same implementation and trusted computing base, not an
# independent implementation.
REPLAY_CHUNK_SIZE = 16_000_003
REPLAY_REDUCTION_LEAF_SIZE = 1_000_003
REPLAY_PRECISION_BITS = 896

MAX_INTEGER_DECIMAL_DIGITS = 500_000

WEIGHT_5_NUMERATORS = (
    -62_705,
    60_498,
    60_570,
    7_832,
    49_592,
    -32_705,
    40_717,
    3_694,
    52_281,
    -30_495,
    -31_874,
    919,
    -17_366,
    19_839,
    -15_054,
    849,
    41_394,
    -28_295,
    -34_191,
    -3_174,
    -14_657,
    31_119,
    -11_369,
    5_734,
    35_933,
    -28_568,
    -29_901,
    830,
    -11_013,
    31_034,
    -8_915,
    -198,
    8_913,
    9_783,
)
WEIGHTS_5 = tuple(
    Fraction(numerator, 1 << WEIGHT_BITS)
    for numerator in WEIGHT_5_NUMERATORS
)

SHELL_5_PAYLOAD_SHA256 = (
    "0cb080416c91cf8dc07b5b5ad154e5ec1974aa4d9bfa1bd42621fca8fe4fc759"
)
WEIGHTS_5_PAYLOAD_SHA256 = (
    "0560615ddbffb655c3b5f016f4fbf97c8797484e2cbaeb48dee4cf06aef68d73"
)
P5_PAYLOAD_SHA256 = (
    "a3ebe954f5a8f72f761fdb220a67e2d34c08ddda4d3e6ad362e69150ae3c92f5"
)

SHELL_5_SUPPORT_COUNT = 65_489
SHELL_5_MINIMUM_INDEX = 65_537
SHELL_5_MAXIMUM_INDEX = 131_072
SHELL_5_L1_NORM = Fraction(186_133, 8)

P5_SUPPORT_COUNT = 112_834
P5_MAXIMUM_INDEX = 131_072
P5_DENOMINATOR_EXPONENT = 25
P5_COEFFICIENT_SUM = Fraction(127_493, 65_536)
P5_C0 = Fraction(3_579, 131_072)
P5_L1_NORM = Fraction(118_651_826_999, 16_777_216)

# Exact safety scalars from a separate full-cutoff scout.  The certified
# prefix run must reproduce all of them before its sign is accepted.
P5_RECURRENCE_ABSOLUTE_BOUND = 118_831_420_454_094_570
P5_LAST_Q_NUMERATOR = 207_489_199_598_882
P5_MAXIMUM_Q_NUMERATOR = 207_489_384_506_602
P5_MAXIMUM_SQUARE_NUMERATOR = 43_051_844_682_928_530_082_981_586_404
P5_INT64_HEADROOM = 9_104_540_616_400_681_237

PARENT_SUPPORT_COUNT = 47_345
PARENT_CHANGED_COUNT = 47_344
PARENT_UNCHANGED_COUNT = 1
DELTA_SUPPORT_COUNT = 112_833
DELTA_COEFFICIENT_SUM = Fraction(-9, 65_536)

GENERATION_CHUNK_COUNT = 64
GENERATION_REDUCTION_LEAF_COUNT = 1_024
GENERATION_MAXIMUM_LEAF_LENGTH = 1 << 20
GENERATION_MAXIMUM_LEAF_DEPTH = 20
GENERATION_ROOT_REDUCTION_DEPTH = 10
GENERATION_REDUCTION_DEPTH_BOUND = 30

REPLAY_CHUNK_COUNT = 68
REPLAY_REDUCTION_LEAF_COUNT = 1_074
REPLAY_MAXIMUM_LEAF_LENGTH = 1_000_003
REPLAY_MAXIMUM_LEAF_DEPTH = 20
REPLAY_ROOT_REDUCTION_DEPTH = 11
REPLAY_REDUCTION_DEPTH_BOUND = 31

# Frozen by exact p5 reconstruction and local-tail hashing.  Unlike the
# prefix, this computation does not scale with T and is covered by
# ``--self-check``.
TAIL_RECORD_SHA256 = (
    "2594584fbf58095038da310b2e9120836016b41163f8054f5206a19cc979a4a5"
)

AS_OF = "2026-07-23"
STATEMENT = (
    "For the explicit exact-dyadic rebased Nyman vectors p4 and p5, where "
    "the fifth shell is reconstructed from exact published p4 and only its "
    "first multiplier is retained, the complete weighted-L2 energy satisfies "
    "E(p4)-E(p5)>1/300000."
)
LIMITATION = (
    "This certifies one finite fifth rebased contraction. The frozen fifth "
    "weight list was selected by an exploratory optimizer. No sixth step, "
    "all-scale recurrence, convergence theorem, proof, or disproof follows. "
    "The Riemann Hypothesis remains unresolved."
)
METHOD = {
    "selection_boundary": (
        "the optimizer proposed one 34-coordinate dyadic vector; exact "
        "reconstruction and audited interval enclosures decide the sign"
    ),
    "parent_binding": (
        "pin the published fourth-step artifact by raw LF, canonical, "
        "payload, source-lineage, and p4-vector commitments"
    ),
    "fifth_shell": (
        "derive the ideal shell on indices 65537..131072 from exact p4, "
        "prove unique nearest 2^-9 bins, and impose final exact balance"
    ),
    "reduced_basis": (
        "reoptimize the 33 inherited direct/shell coordinates and retain only "
        "multiplier one of shell five"
    ),
    "old_lower": (
        "use nonnegativity of the omitted integral to bound complete E(p4) "
        "below by its arbitrary-slope prefix through T=2^30"
    ),
    "new_upper": (
        "bound complete E(p5) above by its prefix at the identical cutoff plus "
        "the exact denominator-weighted local-spacing tail"
    ),
    "large_cutoff_prefix": (
        "stream the exact signed-int64 Q_M recurrence; evaluate nonnegative "
        "terms by exact dyadic scaling, squaring, and sequential divisions; "
        "then cover every binary64 operation with an exact forward-error "
        "radius for a declared pairwise reduction tree"
    ),
    "decision": (
        "subtract the p5 complete-energy upper bound from the p4 prefix lower "
        "bound and test the exact rational threshold"
    ),
}
REFERENCE = fourth.REFERENCE
TRUSTED_COMPUTING_BASE = {
    **fourth.TRUSTED_COMPUTING_BASE,
    "parent_artifact": (
        "published nyman-rebased-fourth-step-v1 artifact pinned byte-for-byte "
        "and by canonical content before exact p4 reconstruction"
    ),
    "fifth_basis": (
        "exact 34-coordinate reduced basis with one fifth-shell multiplier; "
        "no optimizer objective is accepted as evidence"
    ),
    "large_cutoff_prefix_backend": (
        "tools/certify_nyman_absolute_prefix.py; exact int64 recurrence, "
        "sequential binary64 terms, explicit pairwise tree, exact rational "
        "forward-error radius, and Arb logarithmic compression"
    ),
    "old_energy_lower": (
        "nonnegativity of the weighted-L2 integrand outside the certified p4 "
        "prefix at the same cutoff used for p5"
    ),
    "distinct_layout_replay": (
        "different chunk size, different reduction-leaf size, and higher Arb "
        "precision in the same backend must close the same claim, overlap "
        "generation enclosures, and expose different layout/depth metadata"
    ),
    "integer_serialization": (
        "CPython decimal integer conversion limit is raised only inside a "
        "scoped context and restored on exit"
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
    "parent_p4_commitment_matches": True,
    "fifth_shell_is_exact_p4_derived_2^-9_rule": True,
    "fifth_shell_is_exactly_balanced": True,
    "fifth_shell_support_and_l1_match": True,
    "fifth_weights_are_exact_2^-16_list": True,
    "p5_vector_commitment_and_support_match": True,
    "p5_exact_sum_c0_and_l1_match": True,
    "parent_support_is_preserved_through_65536": True,
    "parent_coefficients_have_47344_changed_and_1_unchanged": True,
    "delta_has_112833_nonzero_coefficients_and_exact_sum": True,
    "update_escapes_frozen_prefix_no_go_hypotheses": True,
    "prefixes_use_versioned_large_cutoff_backend": True,
    "prefixes_use_declared_schedule": True,
    "prefixes_use_exact_schedule_and_depth_metadata": True,
    "p5_recurrence_safety_scalars_match": True,
    "p5_exact_max_square_and_int64_headroom_match": True,
    "old_complete_energy_lower_uses_only_p4_prefix": True,
    "new_tail_uses_exact_local_spacing_constant": True,
    "local_spacing_constant_improves_global_spacing_constant": True,
    "optimizer_used_as_evidence": False,
    "complete_gain_strictly_above_claim": True,
    "resolves_rh": False,
}

PREFIX_TEXT_FIELDS = {
    "algorithm_version",
    "term_evaluation_algorithm",
    "reduction_algorithm",
}
PREFIX_INTEGER_FIELDS = {
    "cutoff",
    "chunk_size",
    "chunk_count",
    "reduction_leaf_size",
    "reduction_leaf_count",
    "maximum_leaf_length",
    "maximum_leaf_reduction_depth",
    "root_reduction_depth",
    "reduction_depth_bound",
    "coefficient_denominator",
    "coefficient_denominator_exponent",
    "recurrence_absolute_bound",
    "int64_headroom",
    "term_rounding_factor_count",
    "last_q_numerator",
    "maximum_q_numerator",
    "maximum_square_numerator",
}
PREFIX_FRACTION_FIELDS = {
    "term_relative_error",
    "reduction_gamma",
    "rational_midpoint",
    "rational_error",
    "rational_lower",
    "rational_upper",
    "computed_absolute_term_sum_upper",
    "term_evaluation_error",
    "reduction_error",
    "slope_square_term",
}
PREFIX_BALL_FIELDS = {
    "logarithmic_sum",
    "logarithmic_cross_term",
    "prefix_energy",
}
PROOF_FIELDS = {
    "old_complete_energy_lower_from_prefix",
    "new_prefix",
    "new_tail_upper",
    "new_complete_energy_upper_computation",
    "gain_lower_bound_computation",
}
OLD_PREFIX_PROOF_FIELDS = {
    "relation",
    "prefix",
    "prefix_lower_endpoint",
}
EXPECTED_OLD_PREFIX_RELATION = (
    "E(p4)>=E_prefix(p4;T) by nonnegative omitted square integral"
)
LAYOUT_INDEPENDENT_PREFIX_FIELDS = {
    "algorithm_version",
    "term_evaluation_algorithm",
    "reduction_algorithm",
    "cutoff",
    "coefficient_denominator",
    "coefficient_denominator_exponent",
    "recurrence_absolute_bound",
    "int64_headroom",
    "term_rounding_factor_count",
    "term_relative_error",
    "slope_square_term",
    "last_q_numerator",
    "maximum_q_numerator",
    "maximum_square_numerator",
}


@dataclass(frozen=True)
class FifthStepWitness:
    parent_witness: fourth.FourthStepWitness
    shell5: Mapping[int, Fraction]
    weights5: Mapping[int, Fraction]
    p5: Mapping[int, Fraction]
    source: Mapping[str, object]


def _require_environment() -> None:
    fourth._require_environment()
    absolute_prefix.require_binary64_tcb()


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
        raise ValueError("parent fourth-step artifact is not LF-only")
    if hashlib.sha256(raw).hexdigest() != PARENT_ARTIFACT_RAW_LF_SHA256:
        raise ValueError("parent fourth-step raw LF SHA-256 mismatch")
    supplied = fourth.third.parent.strict_json_loads(raw.decode("utf-8"))
    if not isinstance(supplied, Mapping) or set(supplied) != fourth.TOP_LEVEL_FIELDS:
        raise ValueError("parent fourth-step artifact fields changed")
    if content_sha256(supplied) != PARENT_ARTIFACT_CANONICAL_SHA256:
        raise ValueError("parent fourth-step canonical SHA-256 mismatch")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("parent fourth-step body SHA-256 mismatch")
    if supplied.get("payload_sha256") != PARENT_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("parent fourth-step payload commitment changed")
    if (
        supplied.get("schema") != fourth.SCHEMA
        or supplied.get("classification") != "CERTIFIED_FINITE"
        or supplied.get("hypothesis_status") != "UNRESOLVED"
    ):
        raise ValueError("parent fourth-step semantics changed")
    vectors = supplied.get("vectors")
    p4_record = vectors.get("p4") if isinstance(vectors, Mapping) else None
    if (
        not isinstance(p4_record, Mapping)
        or p4_record.get("payload_sha256") != PARENT_P4_PAYLOAD_SHA256
    ):
        raise ValueError("parent fourth-step p4 commitment changed")
    return supplied


def _source_record(parent_artifact: Mapping[str, object]) -> dict[str, object]:
    return {
        "parent_artifact": {
            "path": PARENT_ARTIFACT_PATH.as_posix(),
            "raw_lf_sha256": PARENT_ARTIFACT_RAW_LF_SHA256,
            "canonical_sha256": PARENT_ARTIFACT_CANONICAL_SHA256,
            "payload_sha256": PARENT_ARTIFACT_PAYLOAD_SHA256,
            "schema": fourth.SCHEMA,
            "p4_payload_sha256": PARENT_P4_PAYLOAD_SHA256,
        },
        "inherited_source": parent_artifact["source"],
    }


def _l1_norm(values: Mapping[int, Fraction]) -> Fraction:
    return sum(
        (abs(Fraction(value)) for value in values.values()),
        start=Fraction(),
    )


def _aggregate_reduced_basis(
    weights: tuple[Fraction, ...],
    shells: tuple[Mapping[int, Fraction], ...],
) -> dict[int, Fraction]:
    if len(weights) != 34 or len(shells) != 5:
        raise ValueError("fifth-step reduced basis contract changed")
    result = fourth._aggregate_reduced_basis(weights[:33], shells[:4])
    shell_weight = weights[33]
    for index, value in shells[4].items():
        fourth.third.parent._add_coefficient(
            result,
            int(index),
            shell_weight * Fraction(value),
        )
    return fourth.third.parent._clean(result)


def _construct_witness(root: Path) -> FifthStepWitness:
    parent_artifact = _load_parent_artifact(root)
    parent_witness = fourth._construct_witness(root)
    if fourth.third.parent._vector_commitment(parent_witness.p4)[
        "payload_sha256"
    ] != PARENT_P4_PAYLOAD_SHA256:
        raise ArithmeticError("exact parent p4 commitment changed")
    if parent_artifact.get("source") != parent_witness.source:
        raise ArithmeticError("parent source lineage differs from reconstruction")
    if parent_artifact.get("vectors") != fourth._vectors_record(parent_witness):
        raise ArithmeticError("parent vectors differ from reconstruction")

    shell5 = fourth.third.parent.rounded_ideal_shell(
        parent_witness.p4,
        SHELL_BITS,
        support_cutoff=65_536,
    )
    weights5 = {
        index: value
        for index, value in enumerate(WEIGHTS_5, start=1)
        if value
    }
    inherited = parent_witness.parent_witness
    base = inherited.parent_witness
    p5 = _aggregate_reduced_basis(
        WEIGHTS_5,
        (
            base.shell1,
            base.shell2,
            inherited.shell3,
            parent_witness.shell4,
            shell5,
        ),
    )

    commitments = {
        "shell5": fourth.third.parent._vector_commitment(shell5),
        "weights5": fourth.third.parent._vector_commitment(weights5),
        "p5": fourth.third.parent._vector_commitment(p5),
    }
    expected = {
        "shell5": {
            "denominator_exponent": str(SHELL_BITS),
            "support_count": str(SHELL_5_SUPPORT_COUNT),
            "maximum_index": str(SHELL_5_MAXIMUM_INDEX),
            "payload_sha256": SHELL_5_PAYLOAD_SHA256,
        },
        "weights5": {
            "denominator_exponent": str(WEIGHT_BITS),
            "support_count": str(len(WEIGHTS_5)),
            "maximum_index": str(len(WEIGHTS_5)),
            "payload_sha256": WEIGHTS_5_PAYLOAD_SHA256,
        },
        "p5": {
            "denominator_exponent": str(P5_DENOMINATOR_EXPONENT),
            "support_count": str(P5_SUPPORT_COUNT),
            "maximum_index": str(P5_MAXIMUM_INDEX),
            "payload_sha256": P5_PAYLOAD_SHA256,
        },
    }
    if commitments != expected:
        raise ArithmeticError("fifth-step vector commitment changed")
    if (
        len(shell5) != SHELL_5_SUPPORT_COUNT
        or min(shell5) != SHELL_5_MINIMUM_INDEX
        or max(shell5) != SHELL_5_MAXIMUM_INDEX
        or fourth.third.parent._coefficient_sum(shell5)
        or _l1_norm(shell5) != SHELL_5_L1_NORM
    ):
        raise ArithmeticError("fifth-shell support, balance, or l1 changed")
    if (
        len(p5) != P5_SUPPORT_COUNT
        or min(p5) != 1
        or max(p5) != P5_MAXIMUM_INDEX
        or fourth.third.parent._coefficient_sum(p5) != P5_COEFFICIENT_SUM
        or _l1_norm(p5) != P5_L1_NORM
    ):
        raise ArithmeticError("rebased p5 exact vector invariants changed")
    parent_support = set(parent_witness.p4)
    if (
        len(parent_support) != PARENT_SUPPORT_COUNT
        or {index for index in p5 if index <= 65_536} != parent_support
    ):
        raise ArithmeticError("fifth-step parent support preservation changed")
    changed_parent_count = sum(
        p5[index] != value
        for index, value in parent_witness.p4.items()
    )
    if (
        changed_parent_count != PARENT_CHANGED_COUNT
        or PARENT_SUPPORT_COUNT - changed_parent_count
        != PARENT_UNCHANGED_COUNT
    ):
        raise ArithmeticError("fifth-step parent change count changed")
    delta = fourth.third.parent._clean(
        {
            index: p5.get(index, Fraction())
            - parent_witness.p4.get(index, Fraction())
            for index in set(p5) | parent_support
        }
    )
    if (
        len(delta) != DELTA_SUPPORT_COUNT
        or fourth.third.parent._coefficient_sum(delta)
        != DELTA_COEFFICIENT_SUM
    ):
        raise ArithmeticError("fifth-step exact delta invariants changed")
    return FifthStepWitness(
        parent_witness=parent_witness,
        shell5=shell5,
        weights5=weights5,
        p5=p5,
        source=_source_record(parent_artifact),
    )


def _vectors_record(witness: FifthStepWitness) -> dict[str, object]:
    p4_commitment = fourth.third.parent._vector_commitment(
        witness.parent_witness.p4
    )
    return {
        **fourth._vectors_record(witness.parent_witness),
        "shell_5": fourth.third.parent._vector_commitment(witness.shell5),
        "weights_5": fourth.third.parent._vector_commitment(witness.weights5),
        "p5": fourth.third.parent._vector_commitment(witness.p5),
        "shell_5_parent_p4_payload_sha256": p4_commitment["payload_sha256"],
    }


def _configuration_record() -> dict[str, object]:
    return {
        "cutoff": str(CUTOFF),
        "shell_bits": str(SHELL_BITS),
        "weight_bits": str(WEIGHT_BITS),
        "shell_5_multiplier_limit": str(SHELL_5_MULTIPLIER_LIMIT),
        "prefix_algorithm_version": absolute_prefix.PREFIX_ALGORITHM_VERSION,
        "term_evaluation_algorithm": (
            absolute_prefix.TERM_EVALUATION_ALGORITHM
        ),
        "reduction_algorithm": absolute_prefix.REDUCTION_ALGORITHM,
        "generation_chunk_size": str(GENERATION_CHUNK_SIZE),
        "generation_chunk_count": str(GENERATION_CHUNK_COUNT),
        "generation_reduction_leaf_size": str(
            GENERATION_REDUCTION_LEAF_SIZE
        ),
        "generation_reduction_leaf_count": str(
            GENERATION_REDUCTION_LEAF_COUNT
        ),
        "generation_maximum_leaf_length": str(
            GENERATION_MAXIMUM_LEAF_LENGTH
        ),
        "generation_maximum_leaf_depth": str(
            GENERATION_MAXIMUM_LEAF_DEPTH
        ),
        "generation_root_reduction_depth": str(
            GENERATION_ROOT_REDUCTION_DEPTH
        ),
        "generation_reduction_depth_bound": str(
            GENERATION_REDUCTION_DEPTH_BOUND
        ),
        "generation_precision_bits": str(GENERATION_PRECISION_BITS),
        "replay_chunk_size": str(REPLAY_CHUNK_SIZE),
        "replay_chunk_count": str(REPLAY_CHUNK_COUNT),
        "replay_reduction_leaf_size": str(REPLAY_REDUCTION_LEAF_SIZE),
        "replay_reduction_leaf_count": str(REPLAY_REDUCTION_LEAF_COUNT),
        "replay_maximum_leaf_length": str(REPLAY_MAXIMUM_LEAF_LENGTH),
        "replay_maximum_leaf_depth": str(REPLAY_MAXIMUM_LEAF_DEPTH),
        "replay_root_reduction_depth": str(
            REPLAY_ROOT_REDUCTION_DEPTH
        ),
        "replay_reduction_depth_bound": str(
            REPLAY_REDUCTION_DEPTH_BOUND
        ),
        "replay_precision_bits": str(REPLAY_PRECISION_BITS),
        "weights_5": [
            fourth.third.parent._fraction_record(value)
            for value in WEIGHTS_5
        ],
    }


def _prefix_record(
    prefix: absolute_prefix.AbsoluteEnergyPrefix,
) -> dict[str, object]:
    field_names = {field.name for field in fields(prefix)}
    declared = (
        PREFIX_TEXT_FIELDS
        | PREFIX_INTEGER_FIELDS
        | PREFIX_FRACTION_FIELDS
        | PREFIX_BALL_FIELDS
    )
    if field_names != declared:
        raise ArithmeticError("large-cutoff prefix record schema changed")
    result: dict[str, object] = {}
    for field in fields(prefix):
        name = field.name
        value = getattr(prefix, name)
        if name in PREFIX_TEXT_FIELDS:
            if not isinstance(value, str):
                raise TypeError(f"prefix {name} is not text")
            result[name] = value
        elif name in PREFIX_INTEGER_FIELDS:
            if type(value) is not int:
                raise TypeError(f"prefix {name} is not an exact integer")
            result[name] = str(value)
        elif name in PREFIX_FRACTION_FIELDS:
            result[name] = fourth.third.parent._fraction_record(
                Fraction(value)
            )
        else:
            result[name] = arb_record(value)
    return result


def _combine_proof(
    old_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    new_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    tail: fourth.third.LocalSpacingTail,
    precision_bits: int,
) -> tuple[Any, Any]:
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        new_upper = (
            new_prefix.prefix_energy
            + absolute_prefix._arb_from_fraction(tail.upper_bound)
        )
        gain = old_prefix.prefix_energy - new_upper
    finally:
        ctx.prec = previous_precision
    return new_upper, gain


def _prepare_exact_tail(
    witness: FifthStepWitness,
) -> fourth.third.LocalSpacingTail:
    """Compute and bind the cutoff-independent exact p5 tail record."""

    tail = fourth.third.local_spacing_tail_upper(witness.p5, CUTOFF)
    tail_record = fourth.third._tail_record(tail)
    if content_sha256(tail_record) != TAIL_RECORD_SHA256:
        raise ArithmeticError("fifth-step exact tail record changed")
    return tail


def _prepare_exact_data(
    root: Path,
) -> tuple[FifthStepWitness, fourth.third.LocalSpacingTail]:
    """Reconstruct and bind every cutoff-independent exact object."""

    witness = _construct_witness(root)
    tail = _prepare_exact_tail(witness)
    return witness, tail


def _generate_components(
    witness: FifthStepWitness,
    tail: fourth.third.LocalSpacingTail,
    *,
    chunk_size: int,
    reduction_leaf_size: int,
    precision_bits: int,
) -> tuple[
    FifthStepWitness,
    absolute_prefix.AbsoluteEnergyPrefix,
    absolute_prefix.AbsoluteEnergyPrefix,
    fourth.third.LocalSpacingTail,
    Any,
    Any,
]:
    old_prefix = absolute_prefix.certified_absolute_energy_prefix(
        witness.parent_witness.p4,
        CUTOFF,
        chunk_size=chunk_size,
        reduction_leaf_size=reduction_leaf_size,
        precision_bits=precision_bits,
    )
    new_prefix = absolute_prefix.certified_absolute_energy_prefix(
        witness.p5,
        CUTOFF,
        chunk_size=chunk_size,
        reduction_leaf_size=reduction_leaf_size,
        precision_bits=precision_bits,
    )
    new_upper, gain = _combine_proof(
        old_prefix,
        new_prefix,
        tail,
        precision_bits,
    )
    return witness, old_prefix, new_prefix, tail, new_upper, gain


def _proof_record(
    old_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    new_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    tail: fourth.third.LocalSpacingTail,
    new_upper: Any,
    gain: Any,
) -> dict[str, object]:
    return {
        "old_complete_energy_lower_from_prefix": {
            "relation": EXPECTED_OLD_PREFIX_RELATION,
            "prefix": _prefix_record(old_prefix),
            "prefix_lower_endpoint": arb_record(
                old_prefix.prefix_energy.lower()
            ),
        },
        "new_prefix": _prefix_record(new_prefix),
        "new_tail_upper": fourth.third._tail_record(tail),
        "new_complete_energy_upper_computation": {
            "enclosure": arb_record(new_upper),
        },
        "gain_lower_bound_computation": {
            "enclosure": arb_record(gain),
        },
    }


def _prefix_layout_signature(
    prefix: absolute_prefix.AbsoluteEnergyPrefix,
) -> tuple[int, int, int, int, int]:
    return (
        prefix.chunk_size,
        prefix.chunk_count,
        prefix.reduction_leaf_size,
        prefix.reduction_leaf_count,
        prefix.maximum_leaf_length,
    )


def _prefix_depth_signature(
    prefix: absolute_prefix.AbsoluteEnergyPrefix,
) -> tuple[int, int, int]:
    return (
        prefix.maximum_leaf_reduction_depth,
        prefix.root_reduction_depth,
        prefix.reduction_depth_bound,
    )


def _expected_schedule_signatures(
    chunk_size: int,
    reduction_leaf_size: int,
    precision_bits: int,
) -> tuple[tuple[int, int, int, int, int], tuple[int, int, int]] | None:
    if (
        chunk_size,
        reduction_leaf_size,
        precision_bits,
    ) == (
        GENERATION_CHUNK_SIZE,
        GENERATION_REDUCTION_LEAF_SIZE,
        GENERATION_PRECISION_BITS,
    ):
        return (
            (
                GENERATION_CHUNK_SIZE,
                GENERATION_CHUNK_COUNT,
                GENERATION_REDUCTION_LEAF_SIZE,
                GENERATION_REDUCTION_LEAF_COUNT,
                GENERATION_MAXIMUM_LEAF_LENGTH,
            ),
            (
                GENERATION_MAXIMUM_LEAF_DEPTH,
                GENERATION_ROOT_REDUCTION_DEPTH,
                GENERATION_REDUCTION_DEPTH_BOUND,
            ),
        )
    if (
        chunk_size,
        reduction_leaf_size,
        precision_bits,
    ) == (
        REPLAY_CHUNK_SIZE,
        REPLAY_REDUCTION_LEAF_SIZE,
        REPLAY_PRECISION_BITS,
    ):
        return (
            (
                REPLAY_CHUNK_SIZE,
                REPLAY_CHUNK_COUNT,
                REPLAY_REDUCTION_LEAF_SIZE,
                REPLAY_REDUCTION_LEAF_COUNT,
                REPLAY_MAXIMUM_LEAF_LENGTH,
            ),
            (
                REPLAY_MAXIMUM_LEAF_DEPTH,
                REPLAY_ROOT_REDUCTION_DEPTH,
                REPLAY_REDUCTION_DEPTH_BOUND,
            ),
        )
    return None


def _checks(
    witness: FifthStepWitness,
    old_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    new_prefix: absolute_prefix.AbsoluteEnergyPrefix,
    tail: fourth.third.LocalSpacingTail,
    gain: Any,
    *,
    chunk_size: int,
    reduction_leaf_size: int,
    precision_bits: int,
) -> dict[str, bool]:
    delta = {
        index: witness.p5.get(index, Fraction())
        - witness.parent_witness.p4.get(index, Fraction())
        for index in set(witness.p5) | set(witness.parent_witness.p4)
    }
    delta = fourth.third.parent._clean(delta)
    parent_support = set(witness.parent_witness.p4)
    p5_parent_range_support = {
        index for index in witness.p5 if index <= 65_536
    }
    changed_parent_count = sum(
        witness.p5[index] != value
        for index, value in witness.parent_witness.p4.items()
    )
    unchanged_parent_count = PARENT_SUPPORT_COUNT - changed_parent_count
    p5_commitment = fourth.third.parent._vector_commitment(witness.p5)
    p5_sum = fourth.third.parent._coefficient_sum(witness.p5)
    expected_schedule = _expected_schedule_signatures(
        chunk_size,
        reduction_leaf_size,
        precision_bits,
    )
    return {
        "parent_artifact_commitments_match": True,
        "parent_source_lineage_matches_exact_reconstruction": True,
        "parent_p4_commitment_matches": (
            fourth.third.parent._vector_commitment(
                witness.parent_witness.p4
            )["payload_sha256"]
            == PARENT_P4_PAYLOAD_SHA256
        ),
        "fifth_shell_is_exact_p4_derived_2^-9_rule": (
            witness.shell5
            == fourth.third.parent.rounded_ideal_shell(
                witness.parent_witness.p4,
                SHELL_BITS,
                support_cutoff=65_536,
            )
        ),
        "fifth_shell_is_exactly_balanced": (
            fourth.third.parent._coefficient_sum(witness.shell5) == 0
        ),
        "fifth_shell_support_and_l1_match": (
            len(witness.shell5) == SHELL_5_SUPPORT_COUNT
            and min(witness.shell5) == SHELL_5_MINIMUM_INDEX
            and max(witness.shell5) == SHELL_5_MAXIMUM_INDEX
            and _l1_norm(witness.shell5) == SHELL_5_L1_NORM
        ),
        "fifth_weights_are_exact_2^-16_list": (
            tuple(witness.weights5[index] for index in range(1, 35))
            == WEIGHTS_5
        ),
        "p5_vector_commitment_and_support_match": (
            p5_commitment["payload_sha256"] == P5_PAYLOAD_SHA256
            and p5_commitment["denominator_exponent"]
            == str(P5_DENOMINATOR_EXPONENT)
            and len(witness.p5) == P5_SUPPORT_COUNT
            and max(witness.p5) == P5_MAXIMUM_INDEX
        ),
        "p5_exact_sum_c0_and_l1_match": (
            p5_sum == P5_COEFFICIENT_SUM
            and 1 - p5_sum / 2 == P5_C0
            and _l1_norm(witness.p5) == P5_L1_NORM
        ),
        "parent_support_is_preserved_through_65536": (
            len(parent_support) == PARENT_SUPPORT_COUNT
            and p5_parent_range_support == parent_support
        ),
        "parent_coefficients_have_47344_changed_and_1_unchanged": (
            changed_parent_count == PARENT_CHANGED_COUNT
            and unchanged_parent_count == PARENT_UNCHANGED_COUNT
        ),
        "delta_has_112833_nonzero_coefficients_and_exact_sum": (
            len(delta) == DELTA_SUPPORT_COUNT
            and fourth.third.parent._coefficient_sum(delta)
            == DELTA_COEFFICIENT_SUM
        ),
        "update_escapes_frozen_prefix_no_go_hypotheses": (
            fourth.third.parent._coefficient_sum(delta) != 0
            and fourth.third.parent._harmonic_sum(delta) != 0
        ),
        "prefixes_use_versioned_large_cutoff_backend": (
            old_prefix.algorithm_version
            == new_prefix.algorithm_version
            == absolute_prefix.PREFIX_ALGORITHM_VERSION
            and old_prefix.term_evaluation_algorithm
            == new_prefix.term_evaluation_algorithm
            == absolute_prefix.TERM_EVALUATION_ALGORITHM
            and old_prefix.reduction_algorithm
            == new_prefix.reduction_algorithm
            == absolute_prefix.REDUCTION_ALGORITHM
        ),
        "prefixes_use_declared_schedule": (
            old_prefix.cutoff == new_prefix.cutoff == CUTOFF
            and old_prefix.chunk_size
            == new_prefix.chunk_size
            == chunk_size
            and old_prefix.reduction_leaf_size
            == new_prefix.reduction_leaf_size
            == reduction_leaf_size
            and precision_bits >= 64
        ),
        "prefixes_use_exact_schedule_and_depth_metadata": (
            expected_schedule is not None
            and _prefix_layout_signature(old_prefix)
            == _prefix_layout_signature(new_prefix)
            == expected_schedule[0]
            and _prefix_depth_signature(old_prefix)
            == _prefix_depth_signature(new_prefix)
            == expected_schedule[1]
        ),
        "p5_recurrence_safety_scalars_match": (
            new_prefix.coefficient_denominator_exponent
            == P5_DENOMINATOR_EXPONENT
            and new_prefix.recurrence_absolute_bound
            == P5_RECURRENCE_ABSOLUTE_BOUND
            and new_prefix.last_q_numerator == P5_LAST_Q_NUMERATOR
            and new_prefix.maximum_q_numerator
            == P5_MAXIMUM_Q_NUMERATOR
            and new_prefix.maximum_square_numerator
            == P5_MAXIMUM_SQUARE_NUMERATOR
            and new_prefix.int64_headroom == P5_INT64_HEADROOM
        ),
        "p5_exact_max_square_and_int64_headroom_match": (
            P5_MAXIMUM_SQUARE_NUMERATOR
            == P5_MAXIMUM_Q_NUMERATOR * P5_MAXIMUM_Q_NUMERATOR
            and P5_INT64_HEADROOM
            == absolute_prefix.INT64_MAX - P5_RECURRENCE_ABSOLUTE_BOUND
            and new_prefix.maximum_square_numerator
            == P5_MAXIMUM_SQUARE_NUMERATOR
            and new_prefix.int64_headroom == P5_INT64_HEADROOM
        ),
        "old_complete_energy_lower_uses_only_p4_prefix": (
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
            > absolute_prefix._arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND)
        ),
        "resolves_rh": False,
    }


@_with_integer_serialization_context
def self_check(root: Path) -> dict[str, object]:
    """Reconstruct frozen finite data without evaluating the huge prefix."""

    _require_environment()
    witness, tail = _prepare_exact_data(root)
    tail_hash = content_sha256(fourth.third._tail_record(tail))
    p5_sum = fourth.third.parent._coefficient_sum(witness.p5)
    if 1 - p5_sum / 2 != P5_C0:
        raise ArithmeticError("fifth-step c0 identity changed")
    return {
        "mode": "EXACT_DATA_PREFLIGHT",
        "parent_payload_sha256": PARENT_ARTIFACT_PAYLOAD_SHA256,
        "shell_5_payload_sha256": SHELL_5_PAYLOAD_SHA256,
        "weights_5_payload_sha256": WEIGHTS_5_PAYLOAD_SHA256,
        "p5_payload_sha256": P5_PAYLOAD_SHA256,
        "tail_record_sha256": tail_hash,
        "p5_support_count": str(len(witness.p5)),
        "p5_maximum_index": str(max(witness.p5)),
        "p5_denominator_exponent": str(P5_DENOMINATOR_EXPONENT),
        "cutoff": str(CUTOFF),
        "prefix_evaluated": False,
        "hypothesis_status": "UNRESOLVED",
        "self_check_passed": True,
        "certificate_verified": False,
    }


@_with_integer_serialization_context
def generate(root: Path) -> dict[str, object]:
    _require_environment()
    exact_data = _prepare_exact_data(root)
    components = _generate_components(
        *exact_data,
        chunk_size=GENERATION_CHUNK_SIZE,
        reduction_leaf_size=GENERATION_REDUCTION_LEAF_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    checks = _checks(
        components[0],
        components[1],
        components[2],
        components[3],
        components[5],
        chunk_size=GENERATION_CHUNK_SIZE,
        reduction_leaf_size=GENERATION_REDUCTION_LEAF_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("fifth-step certificate checks did not close")
    body = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": _configuration_record(),
        "claimed_gain_lower_bound": fourth.third.parent._fraction_record(
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
    result = {**body, "payload_sha256": content_sha256(body)}
    if (
        FROZEN_ARTIFACT_PAYLOAD_SHA256 is not None
        and result["payload_sha256"] != FROZEN_ARTIFACT_PAYLOAD_SHA256
    ):
        raise ArithmeticError(
            "generated fifth-step payload differs from the frozen commitment"
        )
    return result


def _verify_ledger_shape(stored: object) -> None:
    if not isinstance(stored, Mapping) or set(stored) != set(EXPECTED_CHECKS):
        raise ValueError("fifth-step check ledger fields changed")
    if any(type(value) is not bool for value in stored.values()):
        raise ValueError("fifth-step check ledger values are not booleans")


def _parse_canonical_integer(value: object, name: str) -> int:
    return fourth.third.parent._canonical_int_string(value, name)


def _verify_prefix_record_shape(
    record: object,
    label: str,
    *,
    is_new: bool,
) -> None:
    if not isinstance(record, Mapping):
        raise ValueError(f"stored {label} prefix is not an object")
    expected_prefix_fields = {
        field.name for field in fields(absolute_prefix.AbsoluteEnergyPrefix)
    }
    if set(record) != expected_prefix_fields:
        raise ValueError(f"stored {label} prefix fields changed")

    expected_text = {
        "algorithm_version": absolute_prefix.PREFIX_ALGORITHM_VERSION,
        "term_evaluation_algorithm": (
            absolute_prefix.TERM_EVALUATION_ALGORITHM
        ),
        "reduction_algorithm": absolute_prefix.REDUCTION_ALGORITHM,
    }
    for name, expected in expected_text.items():
        if record[name] != expected:
            raise ValueError(f"stored {label} prefix {name} changed")

    integers = {
        name: _parse_canonical_integer(
            record[name],
            f"stored {label} prefix {name}",
        )
        for name in PREFIX_INTEGER_FIELDS
    }
    fractions = {
        name: fourth.third.parent._parse_fraction(
            record[name],
            f"stored {label} prefix {name}",
        )
        for name in PREFIX_FRACTION_FIELDS
    }
    for name in PREFIX_BALL_FIELDS:
        fourth.third.parent._ball_from_record(
            record[name],
            f"stored {label} prefix {name}",
        )

    if (
        integers["cutoff"] != CUTOFF
        or integers["chunk_size"] != GENERATION_CHUNK_SIZE
        or integers["reduction_leaf_size"]
        != GENERATION_REDUCTION_LEAF_SIZE
    ):
        raise ValueError(f"stored {label} prefix generation layout changed")
    expected_chunk_count = (
        CUTOFF + GENERATION_CHUNK_SIZE - 1
    ) // GENERATION_CHUNK_SIZE
    if (
        expected_chunk_count != GENERATION_CHUNK_COUNT
        or integers["chunk_count"] != GENERATION_CHUNK_COUNT
        or integers["reduction_leaf_count"]
        != GENERATION_REDUCTION_LEAF_COUNT
        or integers["maximum_leaf_length"]
        != GENERATION_MAXIMUM_LEAF_LENGTH
        or integers["maximum_leaf_reduction_depth"]
        != GENERATION_MAXIMUM_LEAF_DEPTH
        or integers["root_reduction_depth"]
        != GENERATION_ROOT_REDUCTION_DEPTH
        or integers["reduction_depth_bound"]
        != GENERATION_REDUCTION_DEPTH_BOUND
    ):
        raise ValueError(f"stored {label} prefix schedule metadata changed")
    if not (
        1
        <= integers["maximum_leaf_length"]
        <= integers["reduction_leaf_size"]
        <= integers["chunk_size"]
    ):
        raise ValueError(f"stored {label} prefix leaf layout is invalid")
    if integers["reduction_leaf_count"] < 1:
        raise ValueError(f"stored {label} prefix has no reduction leaves")
    if integers["maximum_leaf_reduction_depth"] != (
        integers["maximum_leaf_length"] - 1
    ).bit_length():
        raise ValueError(f"stored {label} prefix leaf depth changed")
    if integers["root_reduction_depth"] != (
        integers["reduction_leaf_count"] - 1
    ).bit_length():
        raise ValueError(f"stored {label} prefix root depth changed")
    if integers["reduction_depth_bound"] != (
        integers["maximum_leaf_reduction_depth"]
        + integers["root_reduction_depth"]
    ):
        raise ValueError(f"stored {label} prefix depth bound changed")
    denominator = integers["coefficient_denominator"]
    exponent = integers["coefficient_denominator_exponent"]
    if exponent < 0 or denominator != 1 << exponent:
        raise ValueError(f"stored {label} prefix denominator is not dyadic")
    if integers["recurrence_absolute_bound"] < 0:
        raise ValueError(f"stored {label} prefix recurrence bound is negative")
    if integers["int64_headroom"] != (
        absolute_prefix.INT64_MAX
        - integers["recurrence_absolute_bound"]
    ):
        raise ValueError(f"stored {label} prefix int64 headroom changed")
    if integers["term_rounding_factor_count"] != (
        absolute_prefix.BINARY64_TERM_ROUNDING_FACTOR_COUNT
    ):
        raise ValueError(f"stored {label} prefix term model changed")
    maximum_q = integers["maximum_q_numerator"]
    if maximum_q < 0 or integers["maximum_square_numerator"] != (
        maximum_q * maximum_q
    ):
        raise ValueError(f"stored {label} prefix maximum square changed")
    if fractions["rational_error"] < 0:
        raise ValueError(f"stored {label} prefix rational error is negative")
    if fractions["rational_lower"] != (
        fractions["rational_midpoint"] - fractions["rational_error"]
    ):
        raise ValueError(f"stored {label} prefix rational lower changed")
    if fractions["rational_upper"] != (
        fractions["rational_midpoint"] + fractions["rational_error"]
    ):
        raise ValueError(f"stored {label} prefix rational upper changed")
    if fractions["rational_error"] != (
        fractions["term_evaluation_error"]
        + fractions["reduction_error"]
    ):
        raise ValueError(f"stored {label} prefix error decomposition changed")
    expected_term_relative_error = (
        (1 + absolute_prefix.BINARY64_UNIT_ROUNDOFF)
        ** integers["term_rounding_factor_count"]
        - 1
    )
    if fractions["term_relative_error"] != expected_term_relative_error:
        raise ValueError(f"stored {label} prefix term error model changed")
    expected_reduction_gamma = absolute_prefix._gamma(
        integers["reduction_depth_bound"]
    )
    if fractions["reduction_gamma"] != expected_reduction_gamma:
        raise ValueError(f"stored {label} prefix reduction gamma changed")
    expected_absolute_upper = fractions["rational_midpoint"] / (
        1 - expected_reduction_gamma
    )
    if (
        fractions["rational_midpoint"] < 0
        or fractions["computed_absolute_term_sum_upper"]
        != expected_absolute_upper
    ):
        raise ValueError(f"stored {label} prefix midpoint bound changed")
    if fractions["reduction_error"] != (
        expected_reduction_gamma * expected_absolute_upper
    ):
        raise ValueError(f"stored {label} prefix reduction error changed")
    if fractions["term_evaluation_error"] != (
        expected_term_relative_error
        / (1 - expected_term_relative_error)
        * expected_absolute_upper
    ):
        raise ValueError(f"stored {label} prefix term error changed")
    if (
        fractions["term_relative_error"] < 0
        or fractions["reduction_gamma"] < 0
        or fractions["computed_absolute_term_sum_upper"] < 0
    ):
        raise ValueError(f"stored {label} prefix error model is invalid")
    if is_new and (
        exponent != P5_DENOMINATOR_EXPONENT
        or integers["recurrence_absolute_bound"]
        != P5_RECURRENCE_ABSOLUTE_BOUND
        or integers["last_q_numerator"] != P5_LAST_Q_NUMERATOR
        or maximum_q != P5_MAXIMUM_Q_NUMERATOR
        or integers["maximum_square_numerator"]
        != P5_MAXIMUM_SQUARE_NUMERATOR
        or integers["int64_headroom"] != P5_INT64_HEADROOM
    ):
        raise ValueError("stored new prefix p5 safety scalars changed")


def _verify_tail_record_shape(record: object) -> None:
    if not isinstance(record, Mapping):
        raise ValueError("stored fifth-step tail is not an object")
    expected_fields = {
        field.name for field in fields(fourth.third.LocalSpacingTail)
    }
    if set(record) != expected_fields:
        raise ValueError("stored fifth-step tail fields changed")
    if content_sha256(record) != TAIL_RECORD_SHA256:
        raise ValueError("stored fifth-step exact tail hash changed")

    cutoff = _parse_canonical_integer(
        record["cutoff"],
        "stored fifth-step tail cutoff",
    )
    support_limit = _parse_canonical_integer(
        record["support_limit"],
        "stored fifth-step tail support limit",
    )
    exact = {
        name: fourth.third.parent._parse_fraction(
            record[name],
            f"stored fifth-step tail {name}",
        )
        for name in expected_fields - {"cutoff", "support_limit"}
    }
    if cutoff != CUTOFF or support_limit != P5_MAXIMUM_INDEX:
        raise ValueError("stored fifth-step tail cutoff or support changed")
    if (
        exact["coefficient_sum"] != P5_COEFFICIENT_SUM
        or exact["c0"] != P5_C0
        or exact["c0"] != 1 - exact["coefficient_sum"] / 2
    ):
        raise ValueError("stored fifth-step tail c0 identity changed")
    if exact["rho"] < 0 or exact["sigma"] < 0:
        raise ValueError("stored fifth-step tail moment is negative")
    expected_local = Fraction(3, 2) * support_limit * exact["sigma"]
    expected_global = (
        support_limit * (support_limit - 1) * exact["rho"]
        if support_limit >= 2
        else Fraction()
    )
    if exact["local_spacing_constant"] != expected_local:
        raise ValueError("stored fifth-step local spacing constant changed")
    if exact["global_spacing_constant"] != expected_global:
        raise ValueError("stored fifth-step global spacing constant changed")
    expected_improvement = (
        expected_global / expected_local
        if expected_local
        else Fraction()
    )
    if exact["global_to_local_improvement"] != expected_improvement:
        raise ValueError("stored fifth-step spacing ratio changed")
    if exact["mean_term"] != exact["rho"] / (cutoff + 1):
        raise ValueError("stored fifth-step tail mean term changed")
    if exact["discrepancy_term"] != expected_local / (
        (cutoff + 1) * (cutoff + 2)
    ):
        raise ValueError("stored fifth-step tail discrepancy term changed")
    if exact["cross_term"] != (
        abs(exact["harmonic_sum"])
        * exact["absolute_residual_bound"]
        / (cutoff + 1)
    ):
        raise ValueError("stored fifth-step tail cross term changed")
    if exact["slope_square_term"] != (
        exact["harmonic_sum"]
        * exact["harmonic_sum"]
        / (4 * (cutoff + 1))
    ):
        raise ValueError("stored fifth-step tail slope term changed")
    if exact["upper_bound"] != (
        exact["mean_term"]
        + exact["discrepancy_term"]
        + exact["cross_term"]
        + exact["slope_square_term"]
    ):
        raise ValueError("stored fifth-step tail upper bound changed")


def _verify_stored_proof_shape(stored: object) -> None:
    if not isinstance(stored, Mapping) or set(stored) != PROOF_FIELDS:
        raise ValueError("stored fifth-step proof fields changed")
    old = stored["old_complete_energy_lower_from_prefix"]
    if not isinstance(old, Mapping) or set(old) != OLD_PREFIX_PROOF_FIELDS:
        raise ValueError("stored old-prefix lower proof fields changed")
    if old["relation"] != EXPECTED_OLD_PREFIX_RELATION:
        raise ValueError("stored old-prefix lower relation changed")
    fourth.third.parent._ball_from_record(
        old["prefix_lower_endpoint"],
        "old prefix lower endpoint",
    )
    _verify_prefix_record_shape(old["prefix"], "old", is_new=False)
    _verify_prefix_record_shape(stored["new_prefix"], "new", is_new=True)
    _verify_tail_record_shape(stored["new_tail_upper"])
    for name in (
        "new_complete_energy_upper_computation",
        "gain_lower_bound_computation",
    ):
        record = stored[name]
        if not isinstance(record, Mapping) or set(record) != {"enclosure"}:
            raise ValueError(f"{name} fields changed")
        fourth.third.parent._ball_from_record(record["enclosure"], name)


def _verify_generated_proof(
    stored: object,
    generated: Mapping[str, object],
) -> None:
    if stored != generated:
        raise ValueError("stored fifth-step generation proof record changed")


def _require_distinct_layout_consistency(
    generation: tuple[Any, ...],
    replay: tuple[Any, ...],
) -> None:
    expected_generation = _expected_schedule_signatures(
        GENERATION_CHUNK_SIZE,
        GENERATION_REDUCTION_LEAF_SIZE,
        GENERATION_PRECISION_BITS,
    )
    expected_replay = _expected_schedule_signatures(
        REPLAY_CHUNK_SIZE,
        REPLAY_REDUCTION_LEAF_SIZE,
        REPLAY_PRECISION_BITS,
    )
    if (
        expected_generation is None
        or expected_replay is None
        or expected_generation == expected_replay
    ):
        raise ArithmeticError("declared distinct-layout schedules are invalid")
    if _vectors_record(generation[0]) != _vectors_record(replay[0]):
        raise ArithmeticError("distinct-layout replay reconstructed new vectors")
    if fourth.third._tail_record(generation[3]) != fourth.third._tail_record(
        replay[3]
    ):
        raise ArithmeticError("distinct-layout replay changed the exact tail")

    for label, generated_prefix, replayed_prefix in (
        ("old", generation[1], replay[1]),
        ("new", generation[2], replay[2]),
    ):
        for name in LAYOUT_INDEPENDENT_PREFIX_FIELDS:
            if getattr(generated_prefix, name) != getattr(
                replayed_prefix,
                name,
            ):
                raise ArithmeticError(
                    "distinct-layout replay changed "
                    f"{label} prefix recurrence metadata {name}"
                )
        if max(
            generated_prefix.rational_lower,
            replayed_prefix.rational_lower,
        ) > min(
            generated_prefix.rational_upper,
            replayed_prefix.rational_upper,
        ):
            raise ArithmeticError(
                f"distinct-layout {label} rational intervals do not overlap"
            )
        for name in (
            "logarithmic_sum",
            "logarithmic_cross_term",
            "prefix_energy",
        ):
            if not getattr(generated_prefix, name).overlaps(
                getattr(replayed_prefix, name)
            ):
                raise ArithmeticError(
                    f"distinct-layout {label} prefix {name} does not overlap"
                )

        generated_layout = _prefix_layout_signature(generated_prefix)
        replayed_layout = _prefix_layout_signature(replayed_prefix)
        if generated_layout != expected_generation[0]:
            raise ArithmeticError(
                f"generation {label} prefix layout metadata changed"
            )
        if replayed_layout != expected_replay[0]:
            raise ArithmeticError(
                f"distinct-layout {label} replay metadata changed"
            )
        if generated_layout == replayed_layout:
            raise ArithmeticError(
                f"distinct-layout {label} replay reused the generation layout"
            )
        generated_depths = _prefix_depth_signature(generated_prefix)
        replayed_depths = _prefix_depth_signature(replayed_prefix)
        if generated_depths != expected_generation[1]:
            raise ArithmeticError(
                f"generation {label} prefix depth metadata changed"
            )
        if replayed_depths != expected_replay[1]:
            raise ArithmeticError(
                f"distinct-layout {label} replay depth metadata changed"
            )
        if generated_depths == replayed_depths:
            raise ArithmeticError(
                f"distinct-layout {label} replay reused depth metadata"
            )

    for label, first, second in (
        ("new complete upper", generation[4], replay[4]),
        ("gain", generation[5], replay[5]),
    ):
        if not first.overlaps(second):
            raise ArithmeticError(
                f"generation and distinct-layout replay disagree on {label}"
            )


@_with_integer_serialization_context
def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    if FROZEN_ARTIFACT_PAYLOAD_SHA256 is None:
        raise RuntimeError(
            "fifth-step artifact payload is not frozen; a separate audit "
            "must set FROZEN_ARTIFACT_PAYLOAD_SHA256 first"
        )
    supplied = fourth.third.parent.strict_json_loads(
        artifact_path.read_text(encoding="utf-8")
    )
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("fifth-step certificate top-level fields changed")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("fifth-step certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("fifth-step certificate is not the frozen v1 payload")
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
            raise ValueError(f"fifth-step certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("fifth-step certificate configuration changed")
    if supplied.get(
        "claimed_gain_lower_bound"
    ) != fourth.third.parent._fraction_record(CLAIMED_GAIN_LOWER_BOUND):
        raise ValueError("fifth-step claimed gain changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("fifth-step method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("fifth-step trusted-computing-base changed")

    vectors = supplied.get("vectors")
    if not isinstance(vectors, Mapping):
        raise ValueError("fifth-step vector commitments are not an object")
    _verify_ledger_shape(supplied.get("checks"))
    stored_proof = supplied.get("proof")
    _verify_stored_proof_shape(stored_proof)

    parent_artifact = _load_parent_artifact(root)
    if supplied.get("source") != _source_record(parent_artifact):
        raise ValueError("fifth-step parent or source binding changed")

    witness = _construct_witness(root)
    if vectors != _vectors_record(witness):
        raise ValueError("fifth-step vector commitments changed")
    exact_tail = _prepare_exact_tail(witness)
    if not isinstance(stored_proof, Mapping):
        raise ValueError("stored fifth-step proof is not an object")
    exact_tail_record = fourth.third._tail_record(exact_tail)
    if stored_proof["new_tail_upper"] != exact_tail_record:
        raise ValueError("stored fifth-step exact tail differs from reconstruction")
    exact_data = (witness, exact_tail)

    generation = _generate_components(
        *exact_data,
        chunk_size=GENERATION_CHUNK_SIZE,
        reduction_leaf_size=GENERATION_REDUCTION_LEAF_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    generated_proof = _proof_record(*generation[1:])
    _verify_generated_proof(stored_proof, generated_proof)
    generated_checks = _checks(
        generation[0],
        generation[1],
        generation[2],
        generation[3],
        generation[5],
        chunk_size=GENERATION_CHUNK_SIZE,
        reduction_leaf_size=GENERATION_REDUCTION_LEAF_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    if generated_checks != EXPECTED_CHECKS:
        raise ArithmeticError("fifth-step generation checks did not close")
    if supplied.get("checks") != generated_checks:
        raise ValueError("fifth-step check ledger changed")

    replay = _generate_components(
        *exact_data,
        chunk_size=REPLAY_CHUNK_SIZE,
        reduction_leaf_size=REPLAY_REDUCTION_LEAF_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    replay_checks = _checks(
        replay[0],
        replay[1],
        replay[2],
        replay[3],
        replay[5],
        chunk_size=REPLAY_CHUNK_SIZE,
        reduction_leaf_size=REPLAY_REDUCTION_LEAF_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("distinct-layout fifth-step replay did not close")
    _require_distinct_layout_consistency(generation, replay)
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
            "chunk_size": str(REPLAY_CHUNK_SIZE),
            "reduction_leaf_size": str(REPLAY_REDUCTION_LEAF_SIZE),
            "precision_bits": str(REPLAY_PRECISION_BITS),
        },
        "verified": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUTPUT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    artifact_path = args.artifact
    if not artifact_path.is_absolute():
        artifact_path = root / artifact_path
    if args.self_check:
        result = self_check(root)
        display = result
    elif args.verify:
        result = verify(root, artifact_path)
        display = result
    else:
        result = generate(root)
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        artifact_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        display = {
            "artifact": str(artifact_path),
            "artifact_payload_sha256": result["payload_sha256"],
            "classification": result["classification"],
            "hypothesis_status": result["hypothesis_status"],
            "claimed_gain_lower_bound": result["claimed_gain_lower_bound"],
            "written": True,
        }
    print(json.dumps(display, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
