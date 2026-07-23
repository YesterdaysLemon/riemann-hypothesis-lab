"""Canonical finite evidence for the Nyman scaled-distance rebound obstruction.

This module binds the frozen natural-distance summary to two kinds of finite
arithmetic:

* exact rational comparisons proving five declines of
  ``k * d_(2**k)^2`` for ``k = 3, ..., 7``; and
* a clean 256-bit Arb comparison proving
  ``d_256^2 * log(256) < 23/500 < 2 + EulerGamma - log(4*pi)``.

The statement that a later scaled rebound must occur additionally uses the
published strong Nyman--Beurling theorem of Baez-Duarte and Burnol's
asymptotic lower bound.  Those analytic theorems are dependencies recorded in
the artifact, not results reproved by this finite audit.
"""

from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import flint
from flint import arb, ctx

from .artifacts import content_sha256
from .balls import arb_record
from .nyman_summary import (
    NymanSummaryError,
    SUMMARY_SCHEMA,
    verify_nyman_summary,
)


REBOUND_SCHEMA = "rh-lab/nyman-forced-rebound-audit/v1"
FROZEN_AUDIT_ID = "nyman-forced-rebound-v1"
FROZEN_PRECISION_BITS = 256
FROZEN_SUMMARY_PAYLOAD_SHA256 = (
    "35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf"
)
FROZEN_INDEX_PAYLOAD_SHA256 = (
    "281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23"
)
FROZEN_PLAN_SHA256 = (
    "27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a"
)
FROZEN_CELL_SIZES = (8, 16, 32, 64, 128, 256)
RATIONAL_SEPARATOR = Fraction(23, 500)

REBOUND_LIMITATION = (
    "This EXPLORATORY artifact certifies finite rational and 256-bit Arb "
    "inequalities derived from the frozen six-cell Nyman batch. The statement "
    "that a later scaled rebound must exist relies on the cited published "
    "Baez-Duarte and Burnol theorems, which this audit does not reprove. It "
    "does not locate a rebound, prove an all-N recurrence, or prove or "
    "disprove the Riemann Hypothesis; the global status remains UNRESOLVED."
)


class NymanReboundError(ValueError):
    """Raised when the frozen rebound audit cannot be generated safely."""


class NymanReboundVerificationError(NymanReboundError):
    """Raised when a supplied rebound artifact does not exactly regenerate."""


def _backend_record() -> dict[str, str]:
    return {
        "python_flint": str(flint.__version__),
        "flint": str(flint.__FLINT_VERSION__),
    }


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _snapshot_json_tree(
    value: Any,
    error_type: type[NymanReboundError],
    *,
    location: str,
) -> Any:
    """Convert one possibly custom Mapping tree to detached JSON built-ins."""

    if value is None or type(value) in {bool, int, str}:
        return value
    if type(value) is float:
        if not math.isfinite(value):
            raise error_type(f"non-finite JSON number at {location}")
        raise error_type(f"floating-point JSON number is noncanonical at {location}")
    if isinstance(value, list):
        return [
            _snapshot_json_tree(
                item,
                error_type,
                location=f"{location}[{index}]",
            )
            for index, item in enumerate(list(value))
        ]
    if isinstance(value, Mapping):
        snapshot: dict[str, Any] = {}
        for key, item in list(value.items()):
            if type(key) is not str:
                raise error_type(f"non-string JSON object key at {location}")
            if key in snapshot:
                raise error_type(f"duplicate JSON object key: {key}")
            snapshot[key] = _snapshot_json_tree(
                item,
                error_type,
                location=f"{location}.{key}",
            )
        return snapshot
    raise error_type(f"non-JSON value at {location}")


def _load_object(
    source: Mapping[str, Any] | Path,
    error_type: type[NymanReboundError],
    *,
    label: str,
) -> dict[str, Any]:
    if isinstance(source, Path):

        def reject_duplicate_keys(
            pairs: list[tuple[str, Any]],
        ) -> dict[str, Any]:
            value: dict[str, Any] = {}
            for key, item in pairs:
                if key in value:
                    raise error_type(f"duplicate JSON object key: {key}")
                value[key] = item
            return value

        def reject_nonstandard_constant(value: str) -> Any:
            raise error_type(f"nonstandard JSON constant: {value}")

        try:
            loaded = json.loads(
                source.read_text(encoding="utf-8"),
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=reject_nonstandard_constant,
            )
        except error_type:
            raise
        except (OSError, json.JSONDecodeError) as exc:
            raise error_type(f"cannot read {label}: {source}") from exc
    elif isinstance(source, Mapping):
        loaded = source
    else:
        raise TypeError(f"{label} must be a mapping or Path")

    snapshot = _snapshot_json_tree(loaded, error_type, location=label)
    if not isinstance(snapshot, dict):
        raise error_type(f"{label} must be a JSON object")
    return snapshot


def _require_payload_hash(
    artifact: Mapping[str, Any],
    error_type: type[NymanReboundError],
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type("rebound audit payload hash mismatch")


def _canonical_integer(
    value: Any,
    name: str,
    *,
    minimum: int | None = None,
) -> int:
    if not isinstance(value, str):
        raise NymanReboundError(f"{name} must be a canonical integer string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise NymanReboundError(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise NymanReboundError(f"{name} is not canonical")
    if minimum is not None and parsed < minimum:
        raise NymanReboundError(f"{name} must be at least {minimum}")
    return parsed


def _source_fraction(value: Any, name: str) -> Fraction:
    if not isinstance(value, Mapping) or set(value) != {"numerator", "denominator"}:
        raise NymanReboundError(f"{name} has noncanonical fraction fields")
    numerator = _canonical_integer(value["numerator"], f"{name} numerator")
    denominator = _canonical_integer(
        value["denominator"],
        f"{name} denominator",
        minimum=1,
    )
    result = Fraction(numerator, denominator)
    if dict(value) != _fraction_record(result):
        raise NymanReboundError(f"{name} is not a reduced canonical fraction")
    return result


@contextmanager
def _clean_precision() -> Iterable[None]:
    """Run the numerical core from a clean pinned FLINT boundary."""

    previous_precision = ctx.prec
    ctx.cleanup()
    ctx.prec = FROZEN_PRECISION_BITS
    try:
        yield
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()


def _validate_frozen_summary(summary: Mapping[str, Any]) -> dict[int, Mapping[str, Any]]:
    if summary.get("schema") != SUMMARY_SCHEMA:
        raise NymanReboundError("source summary schema changed")
    if summary.get("payload_sha256") != FROZEN_SUMMARY_PAYLOAD_SHA256:
        raise NymanReboundError("source summary payload hash is not frozen")

    source = summary.get("source")
    if not isinstance(source, Mapping):
        raise NymanReboundError("source summary provenance is missing")
    plan = source.get("plan")
    index = source.get("index")
    if (
        not isinstance(plan, Mapping)
        or plan.get("plan_sha256") != FROZEN_PLAN_SHA256
    ):
        raise NymanReboundError("source summary plan hash is not frozen")
    if (
        not isinstance(index, Mapping)
        or index.get("payload_sha256") != FROZEN_INDEX_PAYLOAD_SHA256
    ):
        raise NymanReboundError("source summary index hash is not frozen")

    counts = summary.get("counts")
    if not isinstance(counts, Mapping):
        raise NymanReboundError("source summary counts are missing")
    if counts.get("cells") != "6" or counts.get("by_status") != {
        "FINITE_DISTANCE_BRACKET_CERTIFIED": "6",
        "INCONCLUSIVE": "0",
    }:
        raise NymanReboundError("source summary does not contain six certified cells")

    cells = summary.get("cells")
    if not isinstance(cells, list) or len(cells) != len(FROZEN_CELL_SIZES):
        raise NymanReboundError("source summary cell list changed")
    by_n: dict[int, Mapping[str, Any]] = {}
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise NymanReboundError("source summary cell must be an object")
        n = _canonical_integer(cell.get("n"), "source cell n", minimum=1)
        if n in by_n:
            raise NymanReboundError("source summary contains a duplicate cell")
        if cell.get("status") != "FINITE_DISTANCE_BRACKET_CERTIFIED":
            raise NymanReboundError("source summary contains an uncertified cell")
        by_n[n] = cell
    if tuple(by_n) != FROZEN_CELL_SIZES:
        raise NymanReboundError("source summary cells are not the frozen ordered six")
    return by_n


def _cell_bound(
    cell: Mapping[str, Any],
    side: str,
) -> Fraction:
    bounds = cell.get("bounds")
    if not isinstance(bounds, Mapping):
        raise NymanReboundError("source cell bounds are missing")
    bound = bounds.get(side)
    if not isinstance(bound, Mapping) or set(bound) != {"dyadic", "exact_fraction"}:
        raise NymanReboundError(f"source cell {side} bound changed")
    return _source_fraction(
        bound["exact_fraction"],
        f"source cell {cell.get('n')} {side} bound",
    )


def _scaled_declines(
    by_n: Mapping[int, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    declines: list[dict[str, Any]] = []
    for k in range(3, 8):
        n = 1 << k
        next_n = 2 * n
        lower = _cell_bound(by_n[n], "lower")
        next_upper = _cell_bound(by_n[next_n], "upper")
        factor = Fraction(k, k + 1)
        margin = factor * lower - next_upper
        if margin <= 0:
            raise NymanReboundError(
                f"scaled decline exact margin is not positive for k={k}"
            )
        declines.append(
            {
                "k": str(k),
                "n": str(n),
                "next_n": str(next_n),
                "source_lower_bound": _fraction_record(lower),
                "source_next_upper_bound": _fraction_record(next_upper),
                "contraction_factor": _fraction_record(factor),
                "exact_margin_expression": (
                    "k/(k+1)*L_(2^k)-U_(2^(k+1))"
                ),
                "exact_positive_margin": _fraction_record(margin),
                "decision": "SCALED_DECLINE_CERTIFIED",
                "derived_strict_inequality": (
                    f"{k + 1}*d_{next_n}^2 < {k}*d_{n}^2"
                ),
            }
        )
    return declines


def _numerical_gap(
    lower_256: Fraction,
    upper_256: Fraction,
) -> dict[str, Any]:
    with _clean_precision():
        exact_lower = arb(lower_256.numerator) / lower_256.denominator
        exact_upper = arb(upper_256.numerator) / upper_256.denominator
        log_256 = arb(256).log()
        scaled_upper = exact_upper * log_256
        separator = arb(RATIONAL_SEPARATOR.numerator) / RATIONAL_SEPARATOR.denominator
        finite_margin = separator - scaled_upper
        c0 = arb(2) + arb.const_euler() - (4 * arb.pi()).log()
        asymptotic_margin = c0 - separator
        beta_two_anchor_floor = 9 * exact_lower * arb(2).log()
        beta_two_compatibility_margin = beta_two_anchor_floor - c0
        finite_positive = bool(finite_margin > 0)
        asymptotic_positive = bool(asymptotic_margin > 0)
        beta_two_compatible = bool(beta_two_compatibility_margin > 0)
        if not finite_positive:
            raise NymanReboundError(
                "Arb did not certify U_256*log(256) < 23/500"
            )
        if not asymptotic_positive:
            raise NymanReboundError(
                "Arb did not certify 23/500 < 2+EulerGamma-log(4*pi)"
            )
        if not beta_two_compatible:
            raise NymanReboundError(
                "Arb did not certify C0 < 9*L_256*log(2)"
            )
        return {
            "precision_bits": str(FROZEN_PRECISION_BITS),
            "source_lower_bound_d256_squared": _fraction_record(lower_256),
            "source_upper_bound_d256_squared": _fraction_record(upper_256),
            "log_256": arb_record(log_256),
            "scaled_upper_bound": arb_record(scaled_upper),
            "rational_separator": _fraction_record(RATIONAL_SEPARATOR),
            "finite_positive_margin": arb_record(finite_margin),
            "c0_formula": "2+EulerGamma-log(4*pi)",
            "c0": arb_record(c0),
            "asymptotic_positive_margin": arb_record(asymptotic_margin),
            "beta_two_anchor_floor": arb_record(beta_two_anchor_floor),
            "beta_two_compatibility_margin": arb_record(
                beta_two_compatibility_margin
            ),
            "checks": {
                "d256_squared_at_most_exact_source_upper": True,
                "scaled_upper_strictly_below_23_over_500": finite_positive,
                "23_over_500_strictly_below_c0": asymptotic_positive,
                "c0_strictly_below_9_l256_log_2": beta_two_compatible,
            },
            "certified_chain": "d_256^2*log(256) < 23/500 < C0",
        }


def _dependency_manifest() -> dict[str, Any]:
    return {
        "finite_inputs": {
            "summary_schema": SUMMARY_SCHEMA,
            "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
            "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
            "plan_sha256": FROZEN_PLAN_SHA256,
            "required_certified_cells": [str(n) for n in FROZEN_CELL_SIZES],
        },
        "numeric_backend": {
            **_backend_record(),
            "arithmetic": "python-flint Arb real balls",
            "precision_bits": str(FROZEN_PRECISION_BITS),
            "cache_boundary": (
                "clean FLINT cache before and after the complete pinned "
                "256-bit numerical core"
            ),
        },
        "published_analytic_dependencies": [
            {
                "key": "baez-duarte-strong-natural-dilate-criterion",
                "title": (
                    "A strengthening of the Nyman-Beurling criterion for the "
                    "Riemann Hypothesis"
                ),
                "url": "https://arxiv.org/abs/math/0202141",
                "role": (
                    "RH is equivalent to d_N tending to zero for the nested "
                    "natural-dilate spaces; if RH is false their decreasing "
                    "distances have a positive limit."
                ),
                "reproved_by_this_artifact": False,
            },
            {
                "key": "burnol-asymptotic-lower-bound",
                "title": (
                    "A lower bound in an approximation problem involving the "
                    "zeros of the Riemann zeta function"
                ),
                "url": "https://arxiv.org/abs/math/0103058",
                "role": (
                    "For the continuum distance D(lambda), the lower bound "
                    "gives a multiplicity-square zero sum. Inclusion of the "
                    "natural span in the lambda=1/N continuum span transfers "
                    "that lower bound to d_N."
                ),
                "reproved_by_this_artifact": False,
            },
            {
                "key": "bettin-conrey-farmer-natural-distance-normalization",
                "title": (
                    "An optimal choice of Dirichlet polynomials for the "
                    "Nyman--Beurling criterion"
                ),
                "url": "https://arxiv.org/abs/1211.5191",
                "role": (
                    "States the unconditional multiplicity-square lower bound "
                    "directly for the natural Dirichlet-polynomial distance "
                    "and records the C0 zero-sum identity under RH."
                ),
                "reproved_by_this_artifact": False,
            },
        ],
        "analytic_inference": {
            "statement": (
                "The published-theorem dichotomy gives "
                "liminf_N d_N^2*log(N) >= C0; together with the finite chain, "
                "a later dyadic scaled value must exceed the N=256 value."
            ),
            "derivation_steps": [
                (
                    "If RH is false, the strong natural-dilate criterion and "
                    "nestedness give d_N decreasing to a positive limit, so "
                    "d_N^2*log(N) tends to positive infinity."
                ),
                (
                    "If RH is true, natural-space inclusion transfers "
                    "Burnol's bound: liminf d_N^2*log(N) is at least the sum "
                    "over critical-line zeros of m(rho)^2/|rho|^2."
                ),
                (
                    "Under RH, m(rho)^2 is at least m(rho) and the classical "
                    "Hadamard-product zero-sum identity becomes "
                    "sum m(rho)/|rho|^2 = 2+EulerGamma-log(4*pi)=C0."
                ),
                (
                    "A dyadic subsequence has liminf at least the full-sequence "
                    "liminf; because its N=256 value is below C0, it must "
                    "eventually exceed that value and have an adjacent rise."
                ),
            ],
            "machine_checked_by_this_artifact": False,
            "locates_rebound": False,
            "resolves_rh": False,
        },
    }


def generate_nyman_rebound_audit(
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Generate the canonical frozen finite rebound-obstruction audit."""

    source_summary = _load_object(
        summary,
        NymanReboundError,
        label="source Nyman summary",
    )
    try:
        source_verification = verify_nyman_summary(
            source_summary,
            checkpoint_dir,
        )
    except NymanSummaryError as exc:
        raise NymanReboundError("source Nyman summary verification failed") from exc
    by_n = _validate_frozen_summary(source_summary)
    if source_verification != {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
        "plan_sha256": FROZEN_PLAN_SHA256,
        "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
        "verified_cells": "6",
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }:
        raise NymanReboundError("source summary verification result changed")

    lower_128 = _cell_bound(by_n[128], "lower")
    lower_256 = _cell_bound(by_n[256], "lower")
    upper_256 = _cell_bound(by_n[256], "upper")
    declines = _scaled_declines(by_n)
    payload = {
        "schema": REBOUND_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "FINITE_FORCED_REBOUND_PRECONDITIONS_CERTIFIED",
        "source_verification": source_verification,
        "frozen_source_bounds": {
            "n128_lower": _fraction_record(lower_128),
            "n256_lower": _fraction_record(lower_256),
            "n256_upper": _fraction_record(upper_256),
        },
        "five_exact_scaled_declines": {
            "count": str(len(declines)),
            "all_certified": all(
                item["decision"] == "SCALED_DECLINE_CERTIFIED"
                for item in declines
            ),
            "records": declines,
        },
        "n256_scaled_gap": _numerical_gap(lower_256, upper_256),
        "dependency_manifest": _dependency_manifest(),
        "finite_conclusion": (
            "Five exact dyadic scaled declines are certified through N=256, "
            "and the N=256 scaled upper bound lies strictly below both 23/500 "
            "and the published-theorem asymptotic obstruction C0."
        ),
        "limitation": REBOUND_LIMITATION,
    }
    return _with_payload_hash(payload)


def verify_nyman_rebound_audit(
    artifact: Mapping[str, Any] | Path,
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Exactly regenerate the frozen rebound audit and reject every change."""

    supplied = _load_object(
        artifact,
        NymanReboundVerificationError,
        label="Nyman rebound audit",
    )
    _require_payload_hash(supplied, NymanReboundVerificationError)
    expected = generate_nyman_rebound_audit(summary, checkpoint_dir)
    if supplied != expected:
        raise NymanReboundVerificationError(
            "Nyman rebound audit does not canonically regenerate"
        )
    return {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_REBOUND_AUDIT",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": supplied["audit_outcome"],
        "artifact_payload_sha256": supplied["payload_sha256"],
        "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
        "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
        "plan_sha256": FROZEN_PLAN_SHA256,
        "verified_scaled_declines": supplied["five_exact_scaled_declines"]["count"],
        "precision_bits": str(FROZEN_PRECISION_BITS),
    }


__all__ = [
    "FROZEN_CELL_SIZES",
    "FROZEN_AUDIT_ID",
    "FROZEN_INDEX_PAYLOAD_SHA256",
    "FROZEN_PLAN_SHA256",
    "FROZEN_PRECISION_BITS",
    "FROZEN_SUMMARY_PAYLOAD_SHA256",
    "NymanReboundError",
    "NymanReboundVerificationError",
    "RATIONAL_SEPARATOR",
    "REBOUND_LIMITATION",
    "REBOUND_SCHEMA",
    "generate_nyman_rebound_audit",
    "verify_nyman_rebound_audit",
]
