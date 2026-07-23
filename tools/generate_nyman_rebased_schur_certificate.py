"""Generate and independently replay a finite rebased Nyman certificate.

The numerical optimizer is proposal-only.  This module reconstructs two
explicit exact-dyadic vectors from a frozen N=8 source, proves the first
vector's complete energy by a direct Arb Gram evaluation, and bounds the
second vector's complete energy by an arbitrary-slope interval prefix plus
an absolute large-sieve tail.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
from typing import Any, Mapping, Sequence

import flint
from flint import arb, ctx
import numpy as np

from riemann_lab.artifacts import content_sha256
from riemann_lab.balls import arb_record
from riemann_lab.nyman import (
    build_natural_system,
    evaluate_natural_distance,
    natural_system_content_sha256,
)

if __package__:
    from tools.certify_nyman_balanced_prefix import (
        BINARY64_UNIT_ROUNDOFF,
        INT64_MAX,
        _arb_from_fraction,
        _radius_ball,
        common_dyadic_numerators,
        require_binary64_tcb,
    )
    from tools.certify_nyman_balanced_tail import (
        rounded_ideal_shell,
        strict_json_loads,
    )
    from tools.generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as BASE_TCB,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _parse_fraction,
        _require_environment,
        _strict_lower_proved,
    )
else:
    from certify_nyman_balanced_prefix import (
        BINARY64_UNIT_ROUNDOFF,
        INT64_MAX,
        _arb_from_fraction,
        _radius_ball,
        common_dyadic_numerators,
        require_binary64_tcb,
    )
    from certify_nyman_balanced_tail import (
        rounded_ideal_shell,
        strict_json_loads,
    )
    from generate_nyman_balanced_tail_certificate import (
        TRUSTED_COMPUTING_BASE as BASE_TCB,
        _ball_from_record,
        _dyadic_vector_record,
        _fraction_record,
        _parse_fraction,
        _require_environment,
        _strict_lower_proved,
    )


SCHEMA = "rh-lab/nyman-rebased-schur-certificate/v1"
DEFAULT_OUTPUT = Path("results/nyman-rebased-schur-v1.json")
FROZEN_ARTIFACT_PAYLOAD_SHA256 = (
    "a6e1119fb457e18b521e70203483101a77be7509d52d4f41889951399ce95edd"
)

SOURCE_PATH = Path("results/nyman-natural-v1/cells/n-0008.json")
SOURCE_SCHEMA = "rh-lab/nyman-natural-cell/v1"
SOURCE_RAW_SHA256 = (
    "a4f24c707c3dd7b1a1ca39313abad588647908f302f6b29cfdfaaba09c918af7"
)
SOURCE_PAYLOAD_SHA256 = (
    "07cabf14353ae6ce31b973b2184e00eec087f232aeb7ca2fe9d5352319fa14fa"
)
SOURCE_CANDIDATE_SHA256 = (
    "e5da7950782bd701b3c1b272713e3642afba24d8e1f563cee95deef211f890b0"
)
SOURCE_CELL_CONTRACT_SHA256 = (
    "9f6feeedd9c5832ba0bf3d162e66c3df4c8c5f878b871da2746fdd9da41ae1b9"
)

KERNEL_PATH = Path("results/nyman-natural-v1/kernels/generation.json")
KERNEL_SCHEMA = "rh-lab/nyman-natural-shared-kernel/v1"
KERNEL_RAW_SHA256 = (
    "fe7720a7ed139d57d8f433c53b0353c50b926f0c4c5b0a4ef19f323727bd4d6f"
)
KERNEL_PAYLOAD_SHA256 = (
    "6fb54e7f63208f1bd81fb7ef2e592a8a655458dd7744081748590159b868cc26"
)
KERNEL_CONTRACT_SHA256 = (
    "bc3bc8679379a7a2c3a7c0f530829b67b7170b6ada950c893f6010a32eba06cf"
)
KERNEL_N8_PREFIX_SHA256 = (
    "eb4533a1ffc11d2cc831e241a4e70d22edffe4fd0a483dcfa4f866003cb77aca"
)

CLAIMED_GAIN_LOWER_BOUND = Fraction(1, 5_000)
DIRECT_DIMENSION = 8
MULTIPLIER_LIMIT = 8
SHELL_BITS = 9
WEIGHT_BITS = 16
CUTOFF = 1 << 17
GENERATION_BLOCK_SIZE = 1 << 16
GENERATION_PRECISION_BITS = 256
REPLAY_BLOCK_SIZE = 32_749
REPLAY_PRECISION_BITS = 448

WEIGHTS_1 = tuple(
    Fraction(value, 1 << WEIGHT_BITS)
    for value in (
        -62593,
        60323,
        60359,
        8121,
        49226,
        -32009,
        40467,
        3590,
        49069,
        -27314,
        -27116,
        1877,
        -13291,
        18210,
        -9392,
        -1947,
    )
)
WEIGHTS_2 = tuple(
    Fraction(value, 1 << WEIGHT_BITS)
    for value in (
        -62689,
        60472,
        60537,
        7871,
        49541,
        -32597,
        40683,
        3686,
        51831,
        -30054,
        -31131,
        998,
        -16782,
        19346,
        -14278,
        815,
        37397,
        -25285,
        -29420,
        -2386,
        -12530,
        28264,
        -8007,
        1248,
    )
)

AS_OF = "2026-07-23"
STATEMENT = (
    "For two explicit exact-dyadic rebased Nyman vectors p1 and p2, where "
    "the second shell is reconstructed from the exact p1, the complete "
    "weighted-L2 energy satisfies E(p1)-E(p2)>1/5000."
)
LIMITATION = (
    "This certifies one finite rebased contraction. The two frozen weight "
    "lists were selected by an exploratory optimizer, and no all-scale "
    "recurrence or convergence theorem is proved. The Riemann Hypothesis "
    "remains unresolved."
)
METHOD = {
    "selection_boundary": (
        "the optimizer only proposed two dyadic weight lists; the verifier "
        "reconstructs every basis vector and proves the inequality afresh"
    ),
    "first_shell": (
        "derive the ideal shell on indices 9..16 from the frozen exact "
        "N=8 source, prove unique nearest 2^-9 bins, and impose final balance"
    ),
    "p1": (
        "aggregate eight direct rho coordinates and eight dilates of the "
        "first shell using the frozen 2^-16 weight list"
    ),
    "second_shell": (
        "derive the ideal shell on indices 129..256 from that exact sparse "
        "p1, prove unique nearest 2^-9 bins, and impose final balance"
    ),
    "p2": (
        "reoptimize the original sixteen coordinates and add eight dilates "
        "of the p1-derived second shell using the frozen 2^-16 weight list"
    ),
    "old_energy": (
        "direct canonical Arb Gram evaluation of the exact p1 coefficient "
        "vector; no optimizer objective value is reused"
    ),
    "new_prefix": (
        "P^2*(T+1)+sum(q_M^2/(M(M+1)))-2*P*L_T, including (0,1), "
        "with exact int64 divisor recurrence, certified binary64 reduction, "
        "and Arb Abel compression of L_T"
    ),
    "new_tail": (
        "rho/(T+1)+Q(Q-1)rho/((T+1)(T+2))+"
        "|P|V/(T+1)+P^2/(4(T+1))"
    ),
    "decision": (
        "subtract the certified complete-energy upper bound for p2 from "
        "the direct full-energy enclosure for p1"
    ),
}
REFERENCE = {
    "authors": "H. L. Montgomery and R. C. Vaughan",
    "title": "The large sieve",
    "journal": "Mathematika 20 (1973), 119-134",
    "doi": "10.1112/S0025579300004708",
    "url": "https://doi.org/10.1112/S0025579300004708",
    "used_result": (
        "Theorem 1 in dual form, applied to the Farey frequencies of the "
        "centered periodic residual."
    ),
}
TRUSTED_COMPUTING_BASE = {
    **BASE_TCB,
    "gram": "riemann_lab.nyman canonical Arb natural-dilate Gram builder",
    "formula_scope": (
        "arbitrary-slope absolute-energy prefix and tail helpers are local "
        "to this certificate pending broader reuse"
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
VECTOR_NAMES = {
    "p0_source",
    "shell_1",
    "weights_1",
    "p1",
    "shell_2",
    "weights_2",
    "p2",
    "shell_2_parent_p1_payload_sha256",
}
EXPECTED_CHECKS = {
    "source_and_kernel_commitments_match": True,
    "first_shell_is_exact_source_derived_2^-9_rule": True,
    "first_weights_are_exact_2^-16_list": True,
    "p1_has_55_nonzero_coefficients_through_128": True,
    "second_shell_is_exact_p1_derived_2^-9_rule": True,
    "second_weights_are_exact_2^-16_list": True,
    "p2_has_915_nonzero_coefficients_through_2048": True,
    "optimizer_used_as_evidence": False,
    "complete_gain_strictly_above_claim": True,
    "resolves_rh": False,
}


@dataclass(frozen=True)
class RebasedWitness:
    p0: Mapping[int, Fraction]
    shell1: Mapping[int, Fraction]
    weights1: Mapping[int, Fraction]
    p1: Mapping[int, Fraction]
    shell2: Mapping[int, Fraction]
    weights2: Mapping[int, Fraction]
    p2: Mapping[int, Fraction]
    source: Mapping[str, object]


@dataclass(frozen=True)
class AbsoluteEnergyPrefix:
    cutoff: int
    block_size: int
    block_count: int
    coefficient_denominator: int
    rational_midpoint: Fraction
    rational_error: Fraction
    rational_lower: Fraction
    rational_upper: Fraction
    computed_absolute_term_sum_upper: Fraction
    logarithmic_sum: Any
    logarithmic_cross_term: Any
    slope_square_term: Fraction
    prefix_energy: Any
    last_q_numerator: int
    maximum_q_numerator: int
    maximum_square_numerator: int


@dataclass(frozen=True)
class AbsoluteEnergyTail:
    cutoff: int
    support_limit: int
    farey_spacing_reciprocal: int
    coefficient_sum: Fraction
    harmonic_sum: Fraction
    c0: Fraction
    rho: Fraction
    absolute_residual_bound: Fraction
    mean_term: Fraction
    discrepancy_term: Fraction
    cross_term: Fraction
    slope_square_term: Fraction
    upper_bound: Fraction


def _canonical_body_hash(record: Mapping[str, object], hash_field: str) -> str:
    return content_sha256(
        {key: value for key, value in record.items() if key != hash_field}
    )


def _canonical_int_string(value: object, name: str) -> int:
    if not isinstance(value, str):
        raise ValueError(f"{name} is not a string")
    try:
        result = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} is not an integer") from exc
    if str(result) != value:
        raise ValueError(f"{name} is not canonical")
    return result


def _coefficient_sum(values: Mapping[int, Fraction]) -> Fraction:
    return sum((Fraction(value) for value in values.values()), start=Fraction())


def _harmonic_sum(values: Mapping[int, Fraction]) -> Fraction:
    return sum(
        (
            Fraction(value) / int(index)
            for index, value in values.items()
        ),
        start=Fraction(),
    )


def _clean(values: Mapping[int, Fraction]) -> dict[int, Fraction]:
    result: dict[int, Fraction] = {}
    for raw_index, raw_value in values.items():
        index = int(raw_index)
        value = Fraction(raw_value)
        if index < 1:
            raise ValueError("coefficient indices must be positive")
        if value:
            result[index] = value
    return result


def _add_coefficient(
    target: dict[int, Fraction],
    index: int,
    value: Fraction,
) -> None:
    target[index] = target.get(index, Fraction()) + Fraction(value)
    if not target[index]:
        del target[index]


def _aggregate_basis(
    weights: Sequence[Fraction],
    shells: Sequence[Mapping[int, Fraction]],
) -> dict[int, Fraction]:
    expected = DIRECT_DIMENSION + MULTIPLIER_LIMIT * len(shells)
    if len(weights) != expected:
        raise ValueError("weight count does not match the declared basis")
    result: dict[int, Fraction] = {}
    for index, weight in enumerate(weights[:DIRECT_DIMENSION], start=1):
        _add_coefficient(result, index, Fraction(weight))
    offset = DIRECT_DIMENSION
    for shell in shells:
        for multiplier, weight in enumerate(
            weights[offset : offset + MULTIPLIER_LIMIT],
            start=1,
        ):
            for shell_index, shell_value in shell.items():
                _add_coefficient(
                    result,
                    int(shell_index) * multiplier,
                    Fraction(weight) * Fraction(shell_value),
                )
        offset += MULTIPLIER_LIMIT
    return _clean(result)


def _load_source(root: Path) -> tuple[dict[int, Fraction], dict[str, object]]:
    source_path = root / SOURCE_PATH
    raw_source = source_path.read_bytes()
    if hashlib.sha256(raw_source).hexdigest() != SOURCE_RAW_SHA256:
        raise ValueError("frozen N=8 source raw SHA-256 mismatch")
    source = strict_json_loads(raw_source.decode("utf-8"))
    if not isinstance(source, Mapping) or source.get("schema") != SOURCE_SCHEMA:
        raise ValueError("frozen N=8 source schema changed")
    if source.get("payload_sha256") != SOURCE_PAYLOAD_SHA256:
        raise ValueError("frozen N=8 source payload commitment changed")
    if _canonical_body_hash(source, "payload_sha256") != SOURCE_PAYLOAD_SHA256:
        raise ValueError("frozen N=8 source payload hash mismatch")

    candidate = source.get("candidate")
    if not isinstance(candidate, Mapping) or set(candidate) != {
        "candidate_sha256",
        "coefficients",
        "lower_bound",
        "upper_bound",
    }:
        raise ValueError("frozen N=8 candidate record changed")
    if candidate.get("candidate_sha256") != SOURCE_CANDIDATE_SHA256:
        raise ValueError("frozen N=8 candidate commitment changed")
    if (
        _canonical_body_hash(candidate, "candidate_sha256")
        != SOURCE_CANDIDATE_SHA256
    ):
        raise ValueError("frozen N=8 candidate hash mismatch")
    coefficient_record = candidate.get("coefficients")
    if not isinstance(coefficient_record, Mapping) or set(coefficient_record) != {
        "denominator_exponent",
        "numerators",
    }:
        raise ValueError("frozen N=8 coefficient record changed")
    exponent = _canonical_int_string(
        coefficient_record["denominator_exponent"],
        "source denominator exponent",
    )
    numerators = coefficient_record["numerators"]
    if not isinstance(numerators, list) or len(numerators) != DIRECT_DIMENSION:
        raise ValueError("frozen N=8 coefficient count changed")
    exact_numerators = [
        _canonical_int_string(value, "source coefficient numerator")
        for value in numerators
    ]
    p0 = {
        index: Fraction(numerator, 1 << exponent)
        for index, numerator in enumerate(exact_numerators, start=1)
        if numerator
    }

    contract = source.get("cell_contract")
    if not isinstance(contract, Mapping) or set(contract) != {
        "cell_contract_sha256",
        "cell_id",
        "dilates",
        "dimension",
        "n",
        "schema",
        "shared_kernel_binding",
    }:
        raise ValueError("frozen N=8 cell contract changed")
    if contract.get("cell_contract_sha256") != SOURCE_CELL_CONTRACT_SHA256:
        raise ValueError("frozen N=8 cell-contract commitment changed")
    if (
        _canonical_body_hash(contract, "cell_contract_sha256")
        != SOURCE_CELL_CONTRACT_SHA256
    ):
        raise ValueError("frozen N=8 cell-contract hash mismatch")
    expected_binding = {
        "kernel_contract_sha256": KERNEL_CONTRACT_SHA256,
        "prefix_length": "8",
        "shared_kernel_n": "256",
    }
    if (
        contract.get("n") != "8"
        or contract.get("dimension") != "8"
        or contract.get("dilates")
        != [str(index) for index in range(1, DIRECT_DIMENSION + 1)]
        or contract.get("shared_kernel_binding") != expected_binding
    ):
        raise ValueError("frozen N=8 cell contract semantics changed")

    kernel_path = root / KERNEL_PATH
    raw_kernel = kernel_path.read_bytes()
    if hashlib.sha256(raw_kernel).hexdigest() != KERNEL_RAW_SHA256:
        raise ValueError("frozen shared-kernel raw SHA-256 mismatch")
    kernel = strict_json_loads(raw_kernel.decode("utf-8"))
    if not isinstance(kernel, Mapping) or kernel.get("schema") != KERNEL_SCHEMA:
        raise ValueError("frozen shared-kernel schema changed")
    if kernel.get("payload_sha256") != KERNEL_PAYLOAD_SHA256:
        raise ValueError("frozen shared-kernel payload commitment changed")
    if _canonical_body_hash(kernel, "payload_sha256") != KERNEL_PAYLOAD_SHA256:
        raise ValueError("frozen shared-kernel payload hash mismatch")
    if kernel.get("kernel_contract_sha256") != KERNEL_CONTRACT_SHA256:
        raise ValueError("frozen shared-kernel contract changed")
    prefixes = kernel.get("prefixes")
    if not isinstance(prefixes, list):
        raise ValueError("frozen shared-kernel prefix list changed")
    n8_prefixes = [
        item
        for item in prefixes
        if isinstance(item, Mapping) and item.get("n") == "8"
    ]
    if (
        len(n8_prefixes) != 1
        or n8_prefixes[0].get("prefix_kernel_sha256")
        != KERNEL_N8_PREFIX_SHA256
    ):
        raise ValueError("frozen N=8 kernel prefix commitment changed")

    source_record = {
        "natural_candidate": {
            "path": SOURCE_PATH.as_posix(),
            "raw_sha256": SOURCE_RAW_SHA256,
            "payload_sha256": SOURCE_PAYLOAD_SHA256,
            "candidate_sha256": SOURCE_CANDIDATE_SHA256,
            "cell_contract_sha256": SOURCE_CELL_CONTRACT_SHA256,
        },
        "shared_kernel": {
            "path": KERNEL_PATH.as_posix(),
            "raw_sha256": KERNEL_RAW_SHA256,
            "payload_sha256": KERNEL_PAYLOAD_SHA256,
            "kernel_contract_sha256": KERNEL_CONTRACT_SHA256,
            "n8_prefix_kernel_sha256": KERNEL_N8_PREFIX_SHA256,
        },
    }
    return p0, source_record


def _construct_witness(root: Path) -> RebasedWitness:
    p0, source = _load_source(root)
    shell1 = rounded_ideal_shell(
        p0,
        SHELL_BITS,
        support_cutoff=DIRECT_DIMENSION,
    )
    weights1 = {
        index: value for index, value in enumerate(WEIGHTS_1, start=1) if value
    }
    p1 = _aggregate_basis(WEIGHTS_1, (shell1,))
    shell2 = rounded_ideal_shell(p1, SHELL_BITS, support_cutoff=128)
    weights2 = {
        index: value for index, value in enumerate(WEIGHTS_2, start=1) if value
    }
    p2 = _aggregate_basis(WEIGHTS_2, (shell1, shell2))
    if len(p1) != 55 or max(p1) != 128:
        raise ArithmeticError("rebased p1 support contract changed")
    if len(shell2) != 128 or max(shell2) != 256 or _coefficient_sum(shell2):
        raise ArithmeticError("p1-derived shell contract changed")
    if len(p2) != 915 or max(p2) != 2048:
        raise ArithmeticError("rebased p2 support contract changed")
    return RebasedWitness(
        p0=p0,
        shell1=shell1,
        weights1=weights1,
        p1=p1,
        shell2=shell2,
        weights2=weights2,
        p2=p2,
        source=source,
    )


def _gamma(operation_count: int) -> Fraction:
    if operation_count < 0:
        raise ValueError("operation count must be nonnegative")
    product = operation_count * BINARY64_UNIT_ROUNDOFF
    if product >= 1:
        raise ArithmeticError("binary64 reduction bound exhausted")
    return product / (1 - product)


def _recurrence_operation_bound(
    numerators: Mapping[int, int],
    cutoff: int,
) -> int:
    return sum(
        abs(numerator) * (cutoff // index)
        for index, numerator in numerators.items()
    )


def _compressed_absolute_log_sum(
    coefficients: Mapping[int, Fraction],
    q_t: Fraction,
    cutoff: int,
) -> Any:
    result = _arb_from_fraction(q_t) * arb(cutoff + 1).log()
    for index, coefficient in coefficients.items():
        quotient = cutoff // int(index)
        if not quotient:
            continue
        compressed = (
            quotient * arb(index).log() + arb(quotient + 1).lgamma()
        )
        result -= _arb_from_fraction(Fraction(coefficient)) * compressed
    return result


def absolute_energy_prefix(
    coefficients: Mapping[int, Fraction],
    cutoff: int,
    *,
    block_size: int,
    precision_bits: int,
) -> AbsoluteEnergyPrefix:
    """Enclose the arbitrary-slope absolute energy through interval ``T``.

    On ``[M,M+1)``, put ``P=sum p_n/n`` and
    ``q_M=1+sum p_n floor(M/n)``. Including ``(0,1)``, the prefix is

        (T+1)P^2 + sum q_M^2/[M(M+1)] - 2P L_T,

    with ``L_T`` evaluated by exact Abel compression. This helper is kept
    local until its proof and API receive an independent repository-wide
    audit.
    """

    if cutoff < 1 or block_size < 1:
        raise ValueError("cutoff and block size must be positive")
    if precision_bits < 64:
        raise ValueError("precision must be at least 64 bits")
    require_binary64_tcb()
    if cutoff * (cutoff + 1) >= 1 << 53:
        raise ValueError("cutoff is too large for exact binary64 denominators")
    cleaned = _clean(coefficients)
    denominator, numerators = common_dyadic_numerators(cleaned)
    recurrence_bound = denominator + _recurrence_operation_bound(
        numerators,
        cutoff,
    )
    if recurrence_bound > INT64_MAX:
        raise OverflowError("q divisor recurrence may overflow signed int64")
    if denominator * denominator * cutoff * (cutoff + 1) > 1 << 900:
        raise ValueError("dyadic grid is too fine for the binary64 proof")

    indices = np.asarray(sorted(numerators), dtype=np.int64)
    values = np.asarray(
        [numerators[int(index)] for index in indices],
        dtype=np.int64,
    )
    previous_q = np.int64(denominator)
    block_sums: list[float] = []
    block_absolute_sums: list[Fraction] = []
    block_lengths: list[int] = []
    maximum_q = denominator

    for first in range(1, cutoff + 1, block_size):
        last = min(cutoff, first + block_size - 1)
        length = last - first + 1
        jumps = np.zeros(length, dtype=np.int64)
        for index, numerator in zip(indices, values, strict=True):
            divisor = int(index)
            initial = ((first + divisor - 1) // divisor) * divisor
            if initial <= last:
                jumps[initial - first :: divisor] += numerator
        q_values = previous_q + np.cumsum(jumps, dtype=np.int64)
        maximum_q = max(maximum_q, int(np.max(np.abs(q_values))))

        q_binary64 = q_values.astype(np.float64)
        squares = q_binary64 * q_binary64
        intervals = np.arange(first, last + 1, dtype=np.float64)
        interval_denominators = intervals * (intervals + 1.0)
        terms = (
            squares
            / (denominator * denominator)
            / interval_denominators
        )
        if not np.all(np.isfinite(terms)) or np.any(terms < 0):
            raise ArithmeticError("nonfinite absolute-energy prefix term")
        block_sum = float(np.sum(terms, dtype=np.float64))
        block_absolute = float(np.sum(np.abs(terms), dtype=np.float64))
        if not math.isfinite(block_sum) or not math.isfinite(block_absolute):
            raise ArithmeticError("nonfinite absolute-energy block reduction")
        block_sums.append(block_sum)
        block_absolute_sums.append(Fraction.from_float(block_absolute))
        block_lengths.append(length)
        previous_q = q_values[-1]

    midpoint_float = 0.0
    for block_sum in block_sums:
        midpoint_float += block_sum
    if not math.isfinite(midpoint_float):
        raise ArithmeticError("nonfinite absolute-energy prefix midpoint")
    midpoint = Fraction.from_float(midpoint_float)

    direct_last = denominator + sum(
        numerator * (cutoff // index)
        for index, numerator in numerators.items()
    )
    if int(previous_q) != direct_last:
        raise ArithmeticError("q divisor recurrence endpoint mismatch")

    # q -> binary64, its self-product, and division by M(M+1) contribute
    # at most (1+u)^4-1 relative error. Scaling by the power-of-two
    # coefficient denominator is exact in the guarded normal range.
    term_relative_error = (
        (1 + BINARY64_UNIT_ROUNDOFF) ** 4 - 1
    )
    computed_absolute_upper = Fraction()
    block_reduction_error = Fraction()
    for length, block_absolute in zip(
        block_lengths,
        block_absolute_sums,
        strict=True,
    ):
        reduction_gamma = _gamma(max(0, length - 1))
        exact_computed_absolute = block_absolute / (1 - reduction_gamma)
        computed_absolute_upper += exact_computed_absolute
        block_reduction_error += reduction_gamma * exact_computed_absolute
    term_evaluation_error = (
        term_relative_error
        / (1 - term_relative_error)
        * computed_absolute_upper
    )
    final_reduction_error = _gamma(len(block_sums)) * sum(
        (Fraction.from_float(abs(value)) for value in block_sums),
        start=Fraction(),
    )
    rational_error = (
        term_evaluation_error
        + block_reduction_error
        + final_reduction_error
    )

    p_value = _harmonic_sum(cleaned)
    slope_square = (cutoff + 1) * p_value * p_value
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        logarithmic_sum = _compressed_absolute_log_sum(
            cleaned,
            Fraction(int(previous_q), denominator),
            cutoff,
        )
        logarithmic_cross = (
            -2 * _arb_from_fraction(p_value) * logarithmic_sum
        )
        rational_ball = _arb_from_fraction(midpoint) + _radius_ball(
            rational_error,
            precision_bits,
        )
        prefix_energy = (
            rational_ball
            + logarithmic_cross
            + _arb_from_fraction(slope_square)
        )
    finally:
        ctx.prec = previous_precision
    return AbsoluteEnergyPrefix(
        cutoff=cutoff,
        block_size=block_size,
        block_count=len(block_sums),
        coefficient_denominator=denominator,
        rational_midpoint=midpoint,
        rational_error=rational_error,
        rational_lower=midpoint - rational_error,
        rational_upper=midpoint + rational_error,
        computed_absolute_term_sum_upper=computed_absolute_upper,
        logarithmic_sum=logarithmic_sum,
        logarithmic_cross_term=logarithmic_cross,
        slope_square_term=slope_square,
        prefix_energy=prefix_energy,
        last_q_numerator=int(previous_q),
        maximum_q_numerator=maximum_q,
        maximum_square_numerator=maximum_q * maximum_q,
    )


def _jordan_j2_sieve(limit: int) -> tuple[int, ...]:
    values = [index * index for index in range(limit + 1)]
    if limit >= 1:
        values[1] = 1
    for prime in range(2, limit + 1):
        if values[prime] != prime * prime:
            continue
        square = prime * prime
        for multiple in range(prime, limit + 1, prime):
            values[multiple] = values[multiple] // square * (square - 1)
    return tuple(values)


def _divisor_harmonic_sums(
    coefficients: Mapping[int, Fraction],
    limit: int,
) -> tuple[Fraction, ...]:
    weighted = [Fraction() for _ in range(limit + 1)]
    for index, coefficient in coefficients.items():
        weighted[int(index)] = Fraction(coefficient) / int(index)
    result = [Fraction() for _ in range(limit + 1)]
    for divisor in range(1, limit + 1):
        result[divisor] = sum(
            (
                weighted[multiple]
                for multiple in range(divisor, limit + 1, divisor)
            ),
            start=Fraction(),
        )
    return tuple(result)


def absolute_energy_tail_upper(
    coefficients: Mapping[int, Fraction],
    cutoff: int,
) -> AbsoluteEnergyTail:
    """Return the local arbitrary-slope absolute tail upper bound."""

    if cutoff < 1:
        raise ValueError("tail cutoff must be positive")
    cleaned = _clean(coefficients)
    support_limit = max(cleaned, default=0)
    total = _coefficient_sum(cleaned)
    p_value = _harmonic_sum(cleaned)
    c0 = Fraction(1) - total / 2
    if support_limit:
        divisor_sums = _divisor_harmonic_sums(cleaned, support_limit)
        jordan = _jordan_j2_sieve(support_limit)
        rho = c0 * c0 + sum(
            (
                jordan[divisor]
                * divisor_sums[divisor]
                * divisor_sums[divisor]
                for divisor in range(2, support_limit + 1)
            ),
            start=Fraction(),
        ) / 12
    else:
        rho = c0 * c0
    if rho < 0:
        raise ArithmeticError("periodic residual mean square became negative")
    spacing = (
        support_limit * (support_limit - 1)
        if support_limit >= 2
        else 0
    )
    absolute_residual = abs(c0) + sum(
        (abs(Fraction(value)) for value in cleaned.values()),
        start=Fraction(),
    ) / 2
    mean_term = rho / (cutoff + 1)
    discrepancy_term = (
        spacing * rho / ((cutoff + 1) * (cutoff + 2))
    )
    cross_term = abs(p_value) * absolute_residual / (cutoff + 1)
    slope_square_term = p_value * p_value / (4 * (cutoff + 1))
    upper = mean_term + discrepancy_term + cross_term + slope_square_term
    return AbsoluteEnergyTail(
        cutoff=cutoff,
        support_limit=support_limit,
        farey_spacing_reciprocal=spacing,
        coefficient_sum=total,
        harmonic_sum=p_value,
        c0=c0,
        rho=rho,
        absolute_residual_bound=absolute_residual,
        mean_term=mean_term,
        discrepancy_term=discrepancy_term,
        cross_term=cross_term,
        slope_square_term=slope_square_term,
        upper_bound=upper,
    )


def _direct_full_energy(
    coefficients: Mapping[int, Fraction],
    precision_bits: int,
) -> tuple[Any, str]:
    support_limit = max(coefficients, default=0)
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        system = build_natural_system(tuple(range(1, support_limit + 1)))
        vector = tuple(
            Fraction(coefficients.get(index, Fraction()))
            for index in range(1, support_limit + 1)
        )
        energy = evaluate_natural_distance(system, vector)
        commitment = natural_system_content_sha256(system)
    finally:
        ctx.prec = previous_precision
    return energy, commitment


def _combine_proof(
    old_energy: Any,
    new_prefix: AbsoluteEnergyPrefix,
    new_tail: AbsoluteEnergyTail,
    precision_bits: int,
) -> tuple[Any, Any]:
    previous_precision = ctx.prec
    ctx.prec = precision_bits
    try:
        new_upper = (
            new_prefix.prefix_energy
            + _arb_from_fraction(new_tail.upper_bound)
        )
        gain_lower_computation = old_energy - new_upper
    finally:
        ctx.prec = previous_precision
    return new_upper, gain_lower_computation


def _vector_commitment(values: Mapping[int, Fraction]) -> dict[str, str]:
    record = _dyadic_vector_record(values)
    return {
        "denominator_exponent": str(record["denominator_exponent"]),
        "support_count": str(record["support_count"]),
        "maximum_index": str(record["maximum_index"]),
        "payload_sha256": content_sha256(record),
    }


def _vectors_record(witness: RebasedWitness) -> dict[str, object]:
    p1_commitment = _vector_commitment(witness.p1)
    return {
        "p0_source": _vector_commitment(witness.p0),
        "shell_1": _vector_commitment(witness.shell1),
        "weights_1": _vector_commitment(witness.weights1),
        "p1": p1_commitment,
        "shell_2": _vector_commitment(witness.shell2),
        "weights_2": _vector_commitment(witness.weights2),
        "p2": _vector_commitment(witness.p2),
        "shell_2_parent_p1_payload_sha256": p1_commitment["payload_sha256"],
    }


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
        "weights_1": [_fraction_record(value) for value in WEIGHTS_1],
        "weights_2": [_fraction_record(value) for value in WEIGHTS_2],
    }


def _prefix_record(prefix: AbsoluteEnergyPrefix) -> dict[str, object]:
    return {
        field.name: (
            arb_record(getattr(prefix, field.name))
            if field.name
            in {"logarithmic_sum", "logarithmic_cross_term", "prefix_energy"}
            else (
                str(getattr(prefix, field.name))
                if field.name
                in {
                    "cutoff",
                    "block_size",
                    "block_count",
                    "coefficient_denominator",
                    "last_q_numerator",
                    "maximum_q_numerator",
                    "maximum_square_numerator",
                }
                else _fraction_record(getattr(prefix, field.name))
            )
        )
        for field in fields(AbsoluteEnergyPrefix)
    }


def _tail_record(tail: AbsoluteEnergyTail) -> dict[str, object]:
    return {
        field.name: (
            str(getattr(tail, field.name))
            if field.name
            in {"cutoff", "support_limit", "farey_spacing_reciprocal"}
            else _fraction_record(getattr(tail, field.name))
        )
        for field in fields(AbsoluteEnergyTail)
    }


def _proof_record(
    old_energy: Any,
    old_system_commitment: str,
    new_prefix: AbsoluteEnergyPrefix,
    new_tail: AbsoluteEnergyTail,
    new_upper: Any,
    gain_lower_computation: Any,
) -> dict[str, object]:
    return {
        "old_complete_energy": {
            "support_limit": "128",
            "system_content_sha256": old_system_commitment,
            "enclosure": arb_record(old_energy),
        },
        "new_prefix": _prefix_record(new_prefix),
        "new_tail_upper": _tail_record(new_tail),
        "new_complete_energy_upper_computation": {
            "enclosure": arb_record(new_upper)
        },
        "gain_lower_bound_computation": {
            "enclosure": arb_record(gain_lower_computation)
        },
    }


def _checks(
    witness: RebasedWitness,
    gain_lower_computation: Any,
) -> dict[str, bool]:
    return {
        "source_and_kernel_commitments_match": True,
        "first_shell_is_exact_source_derived_2^-9_rule": (
            witness.shell1
            == rounded_ideal_shell(
                witness.p0,
                SHELL_BITS,
                support_cutoff=DIRECT_DIMENSION,
            )
        ),
        "first_weights_are_exact_2^-16_list": tuple(
            witness.weights1[index] for index in range(1, 17)
        )
        == WEIGHTS_1,
        "p1_has_55_nonzero_coefficients_through_128": (
            len(witness.p1) == 55 and max(witness.p1) == 128
        ),
        "second_shell_is_exact_p1_derived_2^-9_rule": (
            witness.shell2
            == rounded_ideal_shell(
                witness.p1,
                SHELL_BITS,
                support_cutoff=128,
            )
        ),
        "second_weights_are_exact_2^-16_list": tuple(
            witness.weights2[index] for index in range(1, 25)
        )
        == WEIGHTS_2,
        "p2_has_915_nonzero_coefficients_through_2048": (
            len(witness.p2) == 915 and max(witness.p2) == 2048
        ),
        "optimizer_used_as_evidence": False,
        "complete_gain_strictly_above_claim": _strict_lower_proved(
            gain_lower_computation,
            _arb_from_fraction(CLAIMED_GAIN_LOWER_BOUND),
        ),
        "resolves_rh": False,
    }


def _generate_components(
    root: Path,
    *,
    block_size: int,
    precision_bits: int,
) -> tuple[
    RebasedWitness,
    Any,
    str,
    AbsoluteEnergyPrefix,
    AbsoluteEnergyTail,
    Any,
    Any,
]:
    witness = _construct_witness(root)
    old_energy, old_system_commitment = _direct_full_energy(
        witness.p1,
        precision_bits,
    )
    new_prefix = absolute_energy_prefix(
        witness.p2,
        CUTOFF,
        block_size=block_size,
        precision_bits=precision_bits,
    )
    new_tail = absolute_energy_tail_upper(witness.p2, CUTOFF)
    new_upper, gain_lower_computation = _combine_proof(
        old_energy,
        new_prefix,
        new_tail,
        precision_bits,
    )
    return (
        witness,
        old_energy,
        old_system_commitment,
        new_prefix,
        new_tail,
        new_upper,
        gain_lower_computation,
    )


def generate(root: Path) -> dict[str, object]:
    _require_environment()
    (
        witness,
        old_energy,
        old_system_commitment,
        new_prefix,
        new_tail,
        new_upper,
        gain_lower_computation,
    ) = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    checks = _checks(witness, gain_lower_computation)
    if checks != EXPECTED_CHECKS:
        raise ArithmeticError("rebased certificate checks did not close")
    body = {
        "schema": SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "as_of": AS_OF,
        "statement": STATEMENT,
        "configuration": _configuration_record(),
        "claimed_gain_lower_bound": _fraction_record(
            CLAIMED_GAIN_LOWER_BOUND
        ),
        "source": witness.source,
        "vectors": _vectors_record(witness),
        "proof": _proof_record(
            old_energy,
            old_system_commitment,
            new_prefix,
            new_tail,
            new_upper,
            gain_lower_computation,
        ),
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
        raise ValueError("stored generation proof record changed")
    assert isinstance(stored, Mapping)
    old = stored["old_complete_energy"]
    assert isinstance(old, Mapping)
    _ball_from_record(old["enclosure"], "old complete energy")
    new_prefix = stored["new_prefix"]
    assert isinstance(new_prefix, Mapping)
    for name in (
        "logarithmic_sum",
        "logarithmic_cross_term",
        "prefix_energy",
    ):
        _ball_from_record(new_prefix[name], f"new prefix {name}")
    for name in (
        "new_complete_energy_upper_computation",
        "gain_lower_bound_computation",
    ):
        record = stored[name]
        assert isinstance(record, Mapping)
        if set(record) != {"enclosure"}:
            raise ValueError(f"{name} fields changed")
        _ball_from_record(record["enclosure"], name)


def _verify_check_ledger(stored: object, generated: Mapping[str, bool]) -> None:
    """Require a type-canonical JSON boolean ledger before comparing values."""

    if not isinstance(stored, Mapping) or set(stored) != set(generated):
        raise ValueError("rebased check ledger fields changed")
    if any(type(value) is not bool for value in stored.values()):
        raise ValueError("rebased check ledger values are not booleans")
    if stored != generated:
        raise ValueError("rebased check ledger changed")


def verify(root: Path, artifact_path: Path) -> dict[str, object]:
    _require_environment()
    supplied = strict_json_loads(artifact_path.read_text(encoding="utf-8"))
    if not isinstance(supplied, Mapping) or set(supplied) != TOP_LEVEL_FIELDS:
        raise ValueError("rebased certificate top-level fields changed")
    body = {
        key: value for key, value in supplied.items() if key != "payload_sha256"
    }
    if supplied.get("payload_sha256") != content_sha256(body):
        raise ValueError("rebased certificate payload SHA-256 mismatch")
    if supplied.get("payload_sha256") != FROZEN_ARTIFACT_PAYLOAD_SHA256:
        raise ValueError("rebased certificate is not the frozen v1 payload")
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
            raise ValueError(f"rebased certificate {name} changed")
    if supplied.get("configuration") != _configuration_record():
        raise ValueError("rebased certificate configuration changed")
    claim = _parse_fraction(
        supplied.get("claimed_gain_lower_bound"),
        "rebased claimed gain",
    )
    if claim != CLAIMED_GAIN_LOWER_BOUND:
        raise ValueError("rebased claimed gain changed")
    if supplied.get("method") != METHOD or supplied.get("reference") != REFERENCE:
        raise ValueError("rebased method or reference changed")
    if supplied.get("trusted_computing_base") != TRUSTED_COMPUTING_BASE:
        raise ValueError("rebased trusted-computing-base changed")

    generation = _generate_components(
        root,
        block_size=GENERATION_BLOCK_SIZE,
        precision_bits=GENERATION_PRECISION_BITS,
    )
    generation_witness = generation[0]
    if supplied.get("source") != generation_witness.source:
        raise ValueError("rebased source binding changed")
    vectors = supplied.get("vectors")
    if not isinstance(vectors, Mapping) or set(vectors) != VECTOR_NAMES:
        raise ValueError("rebased vector commitment fields changed")
    if vectors != _vectors_record(generation_witness):
        raise ValueError("rebased vector commitments changed")
    generated_proof = _proof_record(*generation[1:])
    _verify_generated_proof(supplied.get("proof"), generated_proof)
    generated_checks = _checks(generation_witness, generation[-1])
    if generated_checks != EXPECTED_CHECKS:
        raise ArithmeticError("generation replay checks did not close")
    _verify_check_ledger(supplied.get("checks"), generated_checks)

    replay = _generate_components(
        root,
        block_size=REPLAY_BLOCK_SIZE,
        precision_bits=REPLAY_PRECISION_BITS,
    )
    replay_witness = replay[0]
    if _vectors_record(replay_witness) != _vectors_record(generation_witness):
        raise ArithmeticError("independent replay reconstructed new vectors")
    replay_checks = _checks(replay_witness, replay[-1])
    if replay_checks != EXPECTED_CHECKS:
        raise ArithmeticError("independent rebased replay did not close")
    return {
        "artifact_payload_sha256": supplied["payload_sha256"],
        "classification": supplied["classification"],
        "hypothesis_status": supplied["hypothesis_status"],
        "replay_old_complete_energy": {
            "lower": replay[1].lower().str(30),
            "upper": replay[1].upper().str(30),
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
