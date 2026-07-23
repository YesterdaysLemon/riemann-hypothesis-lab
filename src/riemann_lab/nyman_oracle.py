"""Independent interval oracles for the Nyman--Beurling normalization.

This module intentionally does not use the finite Vasyunin-sum formula.  It
computes the autocorrelation

    A(lambda) = integral_0^infinity {t} {lambda t} / t^2 dt

directly on a finite interval by partitioning at every exact rational jump of
the two fractional-part functions.  The omitted positive tail is enclosed by
the deliberately generic bound ``[0, 1 / T]``.  Consequently the frozen
``T = 4096`` audit has only about twelve bits of resolving power, regardless
of the working precision used for Arb arithmetic.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Any, Callable, TypeAlias

import flint
from flint import arb, ctx

from .artifacts import content_sha256
from .balls import arb_record


AUTOCORRELATION_SOURCE = "https://arxiv.org/abs/math/0306251"
HARMONIC_BOUND_SOURCE = "https://arxiv.org/abs/math/0306233"

FROZEN_TRUNCATION = 4096
FROZEN_HARMONIC_M = 4096
FROZEN_ORACLE_PRECISION_BITS = 128
FROZEN_RATIO_PAIRS: tuple[tuple[int, int], ...] = (
    (1, 1),
    (1, 2),
    (2, 1),
    (2, 3),
    (3, 2),
    (4, 6),
    (5, 7),
    (7, 5),
    (31, 32),
    (32, 31),
    (127, 128),
    (128, 127),
    (251, 256),
    (256, 251),
)

RatioInput: TypeAlias = Fraction | tuple[int, int]


class NymanOracleError(ValueError):
    """Raised when an oracle input or comparison value is invalid."""


def _canonical_integer(value: Any, name: str, *, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return value


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _ratio_and_metadata(value: RatioInput) -> tuple[Fraction, dict[str, Any]]:
    if isinstance(value, Fraction):
        supplied_numerator = value.numerator
        supplied_denominator = value.denominator
    elif isinstance(value, tuple) and len(value) == 2:
        supplied_numerator = _canonical_integer(value[0], "ratio numerator")
        supplied_denominator = _canonical_integer(
            value[1], "ratio denominator", minimum=1
        )
    else:
        raise TypeError("ratio must be a Fraction or an integer (numerator, denominator) pair")

    if supplied_numerator <= 0:
        raise ValueError("ratio numerator must be positive")
    reduced = Fraction(supplied_numerator, supplied_denominator)
    supplied = {
        "numerator": str(supplied_numerator),
        "denominator": str(supplied_denominator),
    }
    canonical = _fraction_record(reduced)
    return reduced, {
        "supplied": supplied,
        "reduced": canonical,
        "supplied_was_reduced": supplied == canonical,
        "canonicalization": "fractions.Fraction gcd reduction with positive denominator",
    }


def _arb_fraction(value: Fraction) -> Any:
    """Convert an exact rational to an outward-rounded Arb enclosure."""

    return arb(value.numerator) / value.denominator


def _fraction_interval(lower: Fraction, upper: Fraction, bits: int) -> Any:
    """Return a dyadic Arb ball enclosing exact rational endpoints."""

    if lower > upper:
        raise ValueError("interval lower endpoint exceeds upper endpoint")
    scale = 1 << bits
    lower_mantissa = (lower.numerator * scale) // lower.denominator
    upper_mantissa = -((-upper.numerator * scale) // upper.denominator)
    midpoint_mantissa = lower_mantissa + upper_mantissa
    radius_mantissa = upper_mantissa - lower_mantissa
    exponent = -(bits + 1)
    return arb(
        (midpoint_mantissa, exponent),
        (radius_mantissa, exponent),
    )


def _validate_precision(precision_bits: int) -> int:
    precision_bits = _canonical_integer(
        precision_bits, "precision_bits", minimum=64
    )
    return precision_bits


def rational_breakpoints(
    ratio: RatioInput, truncation: int
) -> tuple[Fraction, ...]:
    """Enumerate all exact jumps of ``{t}`` and ``{lambda t}`` on ``[0,T]``."""

    reduced, _ = _ratio_and_metadata(ratio)
    truncation = _canonical_integer(truncation, "truncation", minimum=1)
    endpoint = Fraction(truncation, 1)

    points = {Fraction(0), endpoint}
    points.update(Fraction(index, 1) for index in range(1, truncation + 1))
    scaled_endpoint = reduced * endpoint
    last_scaled_jump = scaled_endpoint.numerator // scaled_endpoint.denominator
    points.update(
        Fraction(index, 1) / reduced
        for index in range(1, last_scaled_jump + 1)
        if Fraction(index, 1) / reduced <= endpoint
    )
    return tuple(sorted(points))


def _integrate_floor_constant_piece(
    ratio: Fraction,
    lower: Fraction,
    upper: Fraction,
    floor_t: int,
    floor_ratio_t: int,
) -> Any:
    """Integrate one exact floor-constant piece using Arb only for logs.

    With ``m = floor(t)`` and ``n = floor(lambda*t)``, the antiderivative
    difference is

      lambda(v-u) - (n + lambda*m) log(v/u)
          + m*n(1/u - 1/v).

    The initial piece begins at zero and has ``m = n = 0``, so only its first
    (nonsingular) term is evaluated.
    """

    if not (Fraction(0) <= lower < upper):
        raise ValueError("piece endpoints must satisfy 0 <= lower < upper")
    if floor_t < 0 or floor_ratio_t < 0:
        raise ValueError("piece floor values must be nonnegative")

    linear = _arb_fraction(ratio * (upper - lower))
    if lower == 0:
        if floor_t != 0 or floor_ratio_t != 0:
            raise ValueError("the piece meeting zero must have both floors zero")
        return linear

    log_coefficient = Fraction(floor_ratio_t, 1) + ratio * floor_t
    log_ratio = _arb_fraction(upper / lower).log()
    reciprocal = Fraction(floor_t * floor_ratio_t, 1) * (
        Fraction(1, 1) / lower - Fraction(1, 1) / upper
    )
    return linear - _arb_fraction(log_coefficient) * log_ratio + _arb_fraction(
        reciprocal
    )


@dataclass(frozen=True)
class TruncatedAutocorrelation:
    """Arb enclosures and exact metadata for one direct autocorrelation run."""

    ratio: Fraction
    ratio_metadata: dict[str, Any]
    truncation: int
    precision_bits: int
    breakpoint_count: int
    truncated: Any
    tail: Any
    total: Any

    @property
    def piece_count(self) -> int:
        return self.breakpoint_count - 1

    def to_record(self) -> dict[str, Any]:
        payload = {
            "schema": "rh-lab/nyman-truncated-autocorrelation-oracle/v1",
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "ratio": self.ratio_metadata,
            "truncation": str(self.truncation),
            "precision_bits": str(self.precision_bits),
            "breakpoint_count": str(self.breakpoint_count),
            "piece_count": str(self.piece_count),
            "integral": {
                "definition": "integral_0^T {t}{lambda*t}/t^2 dt",
                "partition": (
                    "all exact rational jumps of floor(t) and floor(lambda*t)"
                ),
                "piece_formula": (
                    "lambda(v-u)-(n+lambda*m)log(v/u)+m*n(1/u-1/v), "
                    "m=floor(t), n=floor(lambda*t)"
                ),
                "truncated_enclosure": arb_record(self.truncated),
            },
            "tail": {
                "mathematical_bounds": {
                    "lower": _fraction_record(Fraction(0)),
                    "upper": _fraction_record(Fraction(1, self.truncation)),
                },
                "justification": (
                    "0 <= {t}{lambda*t}/t^2 <= 1/t^2 for t >= T"
                ),
                "enclosure": arb_record(self.tail),
            },
            "full_autocorrelation_enclosure": arb_record(self.total),
            "backend": {
                "python_flint": str(flint.__version__),
                "flint": str(flint.__FLINT_VERSION__),
                "arithmetic": (
                    "exact Fraction breakpoints and floor values; Arb-enclosed logs"
                ),
            },
            "source": AUTOCORRELATION_SOURCE,
            "limitation": (
                "The generic tail dominates the resolution; this is a finite "
                "normalization oracle, not evidence resolving RH."
            ),
        }
        return {**payload, "payload_sha256": content_sha256(payload)}


def truncated_inner_product_oracle(
    ratio: RatioInput,
    truncation: int = FROZEN_TRUNCATION,
    precision_bits: int = FROZEN_ORACLE_PRECISION_BITS,
) -> TruncatedAutocorrelation:
    """Rigorously enclose ``A(lambda)`` by direct finite integration plus tail."""

    reduced, ratio_metadata = _ratio_and_metadata(ratio)
    truncation = _canonical_integer(truncation, "truncation", minimum=1)
    precision_bits = _validate_precision(precision_bits)
    breakpoints = rational_breakpoints(reduced, truncation)

    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        truncated = arb(0)
        for lower, upper in zip(breakpoints, breakpoints[1:]):
            midpoint = (lower + upper) / 2
            floor_t = midpoint.numerator // midpoint.denominator
            scaled_midpoint = reduced * midpoint
            floor_ratio_t = (
                scaled_midpoint.numerator // scaled_midpoint.denominator
            )
            truncated += _integrate_floor_constant_piece(
                reduced,
                lower,
                upper,
                floor_t,
                floor_ratio_t,
            )

        tail = _fraction_interval(
            Fraction(0),
            Fraction(1, truncation),
            precision_bits + 16,
        )
        total = truncated + tail
        return TruncatedAutocorrelation(
            ratio=reduced,
            ratio_metadata=ratio_metadata,
            truncation=truncation,
            precision_bits=precision_bits,
            breakpoint_count=len(breakpoints),
            truncated=truncated,
            tail=tail,
            total=total,
        )
    finally:
        ctx.prec = previous_precision


@lru_cache(maxsize=None)
def _exact_harmonic_sum(last_index: int) -> Fraction:
    total = Fraction(0)
    for index in range(1, last_index + 1):
        total += Fraction(1, index)
    return total


@dataclass(frozen=True)
class OneMinusEulerGamma:
    """Independent harmonic-sum enclosure for ``1 - EulerGamma``."""

    last_index: int
    precision_bits: int
    harmonic_sum: Fraction
    enclosure: Any

    def to_record(self) -> dict[str, Any]:
        lower_correction = Fraction(1, 2 * (self.last_index + 1))
        upper_correction = Fraction(1, 2 * self.last_index)
        payload = {
            "schema": "rh-lab/nyman-one-minus-euler-gamma-oracle/v1",
            "classification": "EXPLORATORY",
            "hypothesis_status": "UNRESOLVED",
            "last_harmonic_index": str(self.last_index),
            "precision_bits": str(self.precision_bits),
            "exact_harmonic_sum": _fraction_record(self.harmonic_sum),
            "enclosure": arb_record(self.enclosure),
            "bound": {
                "base": "1-H_M+log(M)",
                "correction_lower": _fraction_record(lower_correction),
                "correction_upper": _fraction_record(upper_correction),
                "inequality": (
                    "1/(2(M+1)) < H_M-log(M)-EulerGamma < 1/(2M)"
                ),
            },
            "backend": {
                "python_flint": str(flint.__version__),
                "flint": str(flint.__FLINT_VERSION__),
                "arithmetic": "exact Fraction harmonic sum; Arb-enclosed log",
            },
            "source": HARMONIC_BOUND_SOURCE,
            "limitation": (
                "This independently checks one normalization constant and has "
                "no implication for the truth of RH."
            ),
        }
        return {**payload, "payload_sha256": content_sha256(payload)}


def one_minus_euler_gamma_oracle(
    last_index: int = FROZEN_HARMONIC_M,
    precision_bits: int = FROZEN_ORACLE_PRECISION_BITS,
) -> OneMinusEulerGamma:
    """Enclose ``1-gamma`` without calling Arb's Euler-constant primitive."""

    last_index = _canonical_integer(last_index, "last_index", minimum=1)
    precision_bits = _validate_precision(precision_bits)
    harmonic_sum = _exact_harmonic_sum(last_index)
    correction_lower = Fraction(1, 2 * (last_index + 1))
    correction_upper = Fraction(1, 2 * last_index)

    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        base = arb(1) - _arb_fraction(harmonic_sum) + arb(last_index).log()
        correction = _fraction_interval(
            correction_lower,
            correction_upper,
            precision_bits + 16,
        )
        enclosure = base + correction
        return OneMinusEulerGamma(
            last_index=last_index,
            precision_bits=precision_bits,
            harmonic_sum=harmonic_sum,
            enclosure=enclosure,
        )
    finally:
        ctx.prec = previous_precision


def generate_frozen_ratio_oracle_table(
    precision_bits: int = FROZEN_ORACLE_PRECISION_BITS,
) -> dict[str, Any]:
    """Generate the deterministic fourteen-ratio, twelve-bit-scope table."""

    precision_bits = _validate_precision(precision_bits)
    entries = [
        truncated_inner_product_oracle(
            supplied,
            truncation=FROZEN_TRUNCATION,
            precision_bits=precision_bits,
        ).to_record()
        for supplied in FROZEN_RATIO_PAIRS
    ]
    payload = {
        "schema": "rh-lab/nyman-frozen-normalization-oracle/v1",
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "precision_bits": str(precision_bits),
        "truncation": str(FROZEN_TRUNCATION),
        "generic_tail_upper_bound": _fraction_record(
            Fraction(1, FROZEN_TRUNCATION)
        ),
        "resolution_scope": (
            "about 12 bits at best because the generic tail is 1/4096; "
            "the 128-bit default controls arithmetic roundoff only"
        ),
        "ordered_supplied_ratio_pairs": [
            {"numerator": str(numerator), "denominator": str(denominator)}
            for numerator, denominator in FROZEN_RATIO_PAIRS
        ],
        "entries": entries,
        "limitation": (
            "This table checks normalization against direct integration. It is "
            "not a high-precision reproduction of a closed formula and does not "
            "resolve RH."
        ),
    }
    return {**payload, "payload_sha256": content_sha256(payload)}


def audit_frozen_ratio_normalization(
    evaluator: Callable[[Fraction], Any],
    *,
    evaluator_label: str,
    precision_bits: int = FROZEN_ORACLE_PRECISION_BITS,
) -> dict[str, Any]:
    """Compare a candidate ``A(lambda)`` evaluator with the frozen oracle.

    The evaluator receives only the reduced ``Fraction``.  The audit retains
    the original supplied pairs separately, so the explicit ``4/6 -> 2/3``
    canonicalization remains visible without treating reduction as a mismatch.
    """

    if not callable(evaluator):
        raise TypeError("evaluator must be callable")
    if not isinstance(evaluator_label, str) or not evaluator_label:
        raise ValueError("evaluator_label must be a nonempty string")
    precision_bits = _validate_precision(precision_bits)

    previous_precision = ctx.prec
    comparisons: list[dict[str, Any]] = []
    try:
        ctx.prec = precision_bits
        for supplied in FROZEN_RATIO_PAIRS:
            reduced, metadata = _ratio_and_metadata(supplied)
            candidate = arb(evaluator(reduced))
            if not candidate.is_finite():
                raise NymanOracleError(
                    f"candidate evaluator returned a nonfinite value for {reduced}"
                )
            oracle = truncated_inner_product_oracle(
                supplied,
                truncation=FROZEN_TRUNCATION,
                precision_bits=precision_bits,
            )
            consistent = bool(candidate.overlaps(oracle.total))
            comparisons.append(
                {
                    "ratio": metadata,
                    "candidate_enclosure": arb_record(candidate),
                    "oracle_full_enclosure": arb_record(oracle.total),
                    "outcome": (
                        "CONSISTENT_WITH_TRUNCATION_ORACLE"
                        if consistent
                        else "DISJOINT_FROM_TRUNCATION_ORACLE"
                    ),
                }
            )
    finally:
        ctx.prec = previous_precision

    all_consistent = all(
        item["outcome"] == "CONSISTENT_WITH_TRUNCATION_ORACLE"
        for item in comparisons
    )
    payload = {
        "schema": "rh-lab/nyman-frozen-normalization-audit/v1",
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": (
            "ALL_14_CONSISTENT_WITH_TRUNCATION_ORACLE"
            if all_consistent
            else "AT_LEAST_ONE_RATIO_DISJOINT_FROM_TRUNCATION_ORACLE"
        ),
        "evaluator_label": evaluator_label,
        "precision_bits": str(precision_bits),
        "truncation": str(FROZEN_TRUNCATION),
        "resolution_scope": (
            "about 12 bits at best; overlap is a normalization check, not an "
            "identity proof"
        ),
        "comparisons": comparisons,
        "limitation": (
            "A passing finite normalization audit has no implication for RH."
        ),
    }
    return {**payload, "payload_sha256": content_sha256(payload)}


__all__ = [
    "AUTOCORRELATION_SOURCE",
    "FROZEN_HARMONIC_M",
    "FROZEN_ORACLE_PRECISION_BITS",
    "FROZEN_RATIO_PAIRS",
    "FROZEN_TRUNCATION",
    "HARMONIC_BOUND_SOURCE",
    "NymanOracleError",
    "OneMinusEulerGamma",
    "TruncatedAutocorrelation",
    "audit_frozen_ratio_normalization",
    "generate_frozen_ratio_oracle_table",
    "one_minus_euler_gamma_oracle",
    "rational_breakpoints",
    "truncated_inner_product_oracle",
]
