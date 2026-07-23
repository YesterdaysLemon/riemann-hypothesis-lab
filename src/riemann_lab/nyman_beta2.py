"""Finite 256-to-512 Nyman contraction certificate.

The module deliberately separates an untrusted approximate solve from the
sign-bearing audit.  ``propose_nyman_beta2_candidate`` builds the canonical
512-dimensional natural-dilate system and rounds its approximate solution to
one exact dyadic vector.  The audit and verifier use only that frozen vector,
direct interval evaluation, and a replay of the already certified N=256 lower
endpoint.

This is a finite result.  It does not establish the all-scale recurrence that
would imply the Riemann Hypothesis.
"""

from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import flint
from flint import ctx

from . import nyman as core
from .artifacts import content_sha256
from .balls import arb_to_dyadic
from .nyman_summary import NymanSummaryError, SUMMARY_SCHEMA, verify_nyman_summary


CANDIDATE_SCHEMA = "rh-lab/nyman-beta2-n512-candidate/v1"
AUDIT_SCHEMA = "rh-lab/nyman-beta2-n512-audit/v1"
FROZEN_CANDIDATE_ID = "nyman-beta2-n512-candidate-v1"
FROZEN_AUDIT_ID = "nyman-beta2-n512-v1"

FROZEN_OLD_N = 256
FROZEN_NEW_N = 512
FROZEN_GENERATION_BITS = 768
FROZEN_REPLAY_BITS = 1536
FROZEN_COEFFICIENT_GRID_BITS = 256
FROZEN_BOUND_GRID_BITS = 128

FROZEN_SUMMARY_PAYLOAD_SHA256 = (
    "35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf"
)
FROZEN_PLAN_SHA256 = (
    "27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a"
)
FROZEN_INDEX_PAYLOAD_SHA256 = (
    "281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23"
)
FROZEN_N256_CELL_PAYLOAD_SHA256 = (
    "89ea59864d96b41f2527ce537325d62da93b7ef237dfbb576f2b979a4e74aaf9"
)
FROZEN_N256_CELL_CONTRACT_SHA256 = (
    "d9d8fa7702a4023bb0a1e94b35d0a9dbae58742986e10449112dc3fdd5e7f97c"
)
FROZEN_N256_CANDIDATE_SHA256 = (
    "10abf61719970a8e2aa92f5806a9cbd8f2e7f94092ad17d06b1f9fc20f100234"
)
FROZEN_N256_LOWER_CERTIFICATE_SHA256 = (
    "fce2efbaf63ed4ab07eb4268cb413a809c634f02ac5cc61c49bb2476c9a36fcc"
)
FROZEN_N256_UPPER_CERTIFICATE_SHA256 = (
    "0e45ca4a3d1c7bcd4ea07f61eb702b6c549c4ae4d7cf16e045addde2e4309933"
)
FROZEN_N256_PREFIX_KERNEL_SHA256 = (
    "c5b9d836f93b406339ece4c971dde01c45f655d8b58ac2f5ed3a59c7ea285676"
)

FROZEN_L256 = Fraction(
    2801788463381838392697853210438663515,
    1 << FROZEN_BOUND_GRID_BITS,
)
FROZEN_U256 = Fraction(
    2801788463381838392697853210438663771,
    1 << FROZEN_BOUND_GRID_BITS,
)
FROZEN_U512_NUMERATOR = 2513498367794989737282886562436186984
FROZEN_U512 = Fraction(FROZEN_U512_NUMERATOR, 1 << FROZEN_BOUND_GRID_BITS)

STRONG_CONTRACTION = Fraction(449, 500)
BETA_TWO_CONTRACTION = Fraction(9, 10)
STRONG_GAIN = 1 - STRONG_CONTRACTION
BETA_TWO_GAIN = 1 - BETA_TWO_CONTRACTION

FROZEN_CANDIDATE_PAYLOAD_SHA256 = (
    "c6be030ede77f9a6a43d2219620d9cb78d3724f44bf0bb87705e1f0722371343"
)
FROZEN_COEFFICIENT_VECTOR_SHA256 = (
    "b5b2e8c1acb06c9d1424f5cdbb56503f214a05459b7c3f95448fa2549291aaef"
)
FROZEN_GENERATION_KERNEL_SHA256 = (
    "fa1b758f9d0a0ef7457976e11f0598f27650556f39629cd8e8544eae5bd5cef8"
)
FROZEN_GENERATION_PREFIX_512_SHA256 = (
    "d25dec076c37a8d1849521c807ec0253763153c4f66efc333ff24b791826b6a2"
)
FROZEN_AUDIT_PAYLOAD_SHA256 = (
    "2ef885119528faaffb1a67f6d8382de710a8c534d8413b621be80f5f4c4d08b6"
)

CANDIDATE_LIMITATION = (
    "This EXPLORATORY candidate is the exact dyadic output of an approximate "
    "solve and certifies nothing by itself. Only direct interval evaluation "
    "in the separate audit can give it a finite evidentiary role. It gives no "
    "N=512 lower bound, no all-scale recurrence, and no proof or disproof of "
    "the Riemann Hypothesis; the global status remains UNRESOLVED."
)

AUDIT_LIMITATION = (
    "This CERTIFIED_FINITE artifact certifies one finite 256-to-512 "
    "contraction using a stored exact dyadic witness and a fresh replay of "
    "the frozen N=256 lower endpoint. Its verifier additionally performs a "
    "same-backend Arb replay at 1536 bits; that replay is not a clean-room "
    "implementation. The artifact gives no N=512 lower bound, no all-scale "
    "recurrence, and no proof or disproof of the Riemann Hypothesis; the "
    "global status remains UNRESOLVED."
)


class NymanBeta2Error(ValueError):
    """Raised when the finite beta=2 audit cannot be generated safely."""


class NymanBeta2VerificationError(NymanBeta2Error):
    """Raised when a supplied beta=2 artifact fails exact reproduction."""


def _backend_record() -> dict[str, str]:
    return {
        "python_flint": str(flint.__version__),
        "flint": str(flint.__FLINT_VERSION__),
    }


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": content_sha256(body)}


def _require_payload_hash(
    artifact: Mapping[str, Any],
    error_type: type[NymanBeta2Error],
    *,
    label: str,
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type(f"{label} payload hash mismatch")


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _dyadic_record(numerator: int, denominator_exponent: int) -> dict[str, str]:
    return {
        "numerator": str(numerator),
        "denominator_exponent": str(denominator_exponent),
    }


def _canonical_integer(
    value: Any,
    name: str,
    *,
    minimum: int | None = None,
    error_type: type[NymanBeta2Error] = NymanBeta2Error,
) -> int:
    if not isinstance(value, str):
        raise error_type(f"{name} must be a canonical integer string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise error_type(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise error_type(f"{name} is not canonical")
    if minimum is not None and parsed < minimum:
        raise error_type(f"{name} must be at least {minimum}")
    return parsed


def _source_fraction(
    value: Any,
    name: str,
    *,
    error_type: type[NymanBeta2Error] = NymanBeta2Error,
) -> Fraction:
    if not isinstance(value, Mapping) or set(value) != {"numerator", "denominator"}:
        raise error_type(f"{name} has noncanonical fraction fields")
    numerator = _canonical_integer(
        value["numerator"], f"{name} numerator", error_type=error_type
    )
    denominator = _canonical_integer(
        value["denominator"],
        f"{name} denominator",
        minimum=1,
        error_type=error_type,
    )
    result = Fraction(numerator, denominator)
    if dict(value) != _fraction_record(result):
        raise error_type(f"{name} is not a reduced canonical fraction")
    return result


def _snapshot_json_tree(
    value: Any,
    error_type: type[NymanBeta2Error],
    *,
    location: str,
) -> Any:
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
    error_type: type[NymanBeta2Error],
    *,
    label: str,
) -> dict[str, Any]:
    if isinstance(source, Path):

        def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
            result: dict[str, Any] = {}
            for key, item in pairs:
                if key in result:
                    raise error_type(f"duplicate JSON object key: {key}")
                result[key] = item
            return result

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


@contextmanager
def _clean_precision(precision_bits: int) -> Iterable[None]:
    previous_precision = ctx.prec
    ctx.cleanup()
    ctx.prec = precision_bits
    try:
        yield
    finally:
        ctx.prec = previous_precision
        ctx.cleanup()


def _dyadic_fraction(mantissa: Any, exponent: Any) -> Fraction:
    exact_mantissa = int(mantissa)
    exact_exponent = int(exponent)
    if exact_exponent >= 0:
        return Fraction(exact_mantissa << exact_exponent, 1)
    return Fraction(exact_mantissa, 1 << -exact_exponent)


def _prefix_kernel_sha256(system: core.NaturalSystem, n: int) -> str:
    prefix = core.prefix_natural_system(system, n)
    body = {
        "precision_bits": str(prefix.precision_bits),
        "dilates": [str(value) for value in prefix.dilates],
        "target": [arb_to_dyadic(value) for value in prefix.target],
        "gram_upper_triangle": [
            {
                "row": str(row),
                "column": str(column),
                "value": arb_to_dyadic(prefix.gram[row][column]),
            }
            for row in range(n)
            for column in range(row, n)
        ],
    }
    return content_sha256(body)


def _candidate_coefficients(numerators: tuple[int, ...]) -> dict[str, Any]:
    return {
        "encoding": "signed-integers-over-common-power-of-two",
        "numerators": [str(value) for value in numerators],
        "denominator_exponent": str(FROZEN_COEFFICIENT_GRID_BITS),
    }


def _candidate_body(
    numerators: tuple[int, ...],
    system: core.NaturalSystem,
) -> dict[str, Any]:
    coefficients = _candidate_coefficients(numerators)
    vector_sha256 = content_sha256(coefficients)
    prefix_256_sha256 = _prefix_kernel_sha256(system, FROZEN_OLD_N)
    prefix_512_sha256 = _prefix_kernel_sha256(system, FROZEN_NEW_N)
    return {
        "schema": CANDIDATE_SCHEMA,
        "candidate_id": FROZEN_CANDIDATE_ID,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "n": str(FROZEN_NEW_N),
        "policy": {
            "generation_precision_bits": str(FROZEN_GENERATION_BITS),
            "coefficient_grid": f"2^-{FROZEN_COEFFICIENT_GRID_BITS}",
            "rounding": "exact-rational-round-ties-to-even",
            "failed_generation_behavior": "fail closed; no adaptive policy change",
        },
        "backend": _backend_record(),
        "kernel": {
            "core_system_content_sha256": core.natural_system_content_sha256(system),
            "prefix_256_numeric_sha256": prefix_256_sha256,
            "prefix_512_numeric_sha256": prefix_512_sha256,
            "old_prefix_matches_frozen_v1_numeric_kernel": (
                prefix_256_sha256 == FROZEN_N256_PREFIX_KERNEL_SHA256
            ),
        },
        "solver_role": "approximate-untrusted-candidate-generator-only",
        "coefficients": coefficients,
        "coefficient_vector_sha256": vector_sha256,
        "limitation": CANDIDATE_LIMITATION,
    }


def propose_nyman_beta2_candidate() -> dict[str, Any]:
    """Generate the frozen exact candidate; the solve has no certificate role."""

    try:
        with _clean_precision(FROZEN_GENERATION_BITS):
            system = core.build_natural_system(
                tuple(range(1, FROZEN_NEW_N + 1))
            )
            numerators = core.propose_dyadic_coefficients(
                system, FROZEN_COEFFICIENT_GRID_BITS
            )
            proposal = _with_payload_hash(_candidate_body(numerators, system))
    except Exception as exc:  # pragma: no cover - backend-specific failures
        raise NymanBeta2Error("canonical N=512 candidate generation failed") from exc
    _validate_candidate(proposal, NymanBeta2Error, require_frozen=True)
    return proposal


def _parse_numerators(
    candidate: Mapping[str, Any],
    error_type: type[NymanBeta2Error],
) -> tuple[int, ...]:
    coefficients = candidate.get("coefficients")
    if not isinstance(coefficients, Mapping) or set(coefficients) != {
        "encoding",
        "numerators",
        "denominator_exponent",
    }:
        raise error_type("candidate coefficients changed fields")
    if coefficients.get("encoding") != (
        "signed-integers-over-common-power-of-two"
    ):
        raise error_type("candidate coefficient encoding changed")
    exponent = _canonical_integer(
        coefficients.get("denominator_exponent"),
        "candidate coefficient denominator exponent",
        minimum=1,
        error_type=error_type,
    )
    if exponent != FROZEN_COEFFICIENT_GRID_BITS:
        raise error_type("candidate coefficient grid changed")
    raw_numerators = coefficients.get("numerators")
    if not isinstance(raw_numerators, list) or len(raw_numerators) != FROZEN_NEW_N:
        raise error_type("candidate coefficient vector dimension changed")
    return tuple(
        _canonical_integer(
            value,
            f"candidate coefficient numerator {index}",
            error_type=error_type,
        )
        for index, value in enumerate(raw_numerators)
    )


def _validate_candidate(
    candidate: Mapping[str, Any],
    error_type: type[NymanBeta2Error],
    *,
    require_frozen: bool,
) -> tuple[int, ...]:
    _require_payload_hash(candidate, error_type, label="candidate")
    expected_fields = {
        "schema",
        "candidate_id",
        "classification",
        "hypothesis_status",
        "n",
        "policy",
        "backend",
        "kernel",
        "solver_role",
        "coefficients",
        "coefficient_vector_sha256",
        "limitation",
        "payload_sha256",
    }
    if set(candidate) != expected_fields:
        raise error_type("candidate fields changed")
    if candidate.get("schema") != CANDIDATE_SCHEMA:
        raise error_type("candidate schema changed")
    if candidate.get("candidate_id") != FROZEN_CANDIDATE_ID:
        raise error_type("candidate id changed")
    if candidate.get("classification") != "EXPLORATORY":
        raise error_type("candidate classification changed")
    if candidate.get("hypothesis_status") != "UNRESOLVED":
        raise error_type("candidate improperly resolves RH")
    if candidate.get("n") != str(FROZEN_NEW_N):
        raise error_type("candidate dimension changed")
    if candidate.get("policy") != {
        "generation_precision_bits": str(FROZEN_GENERATION_BITS),
        "coefficient_grid": f"2^-{FROZEN_COEFFICIENT_GRID_BITS}",
        "rounding": "exact-rational-round-ties-to-even",
        "failed_generation_behavior": "fail closed; no adaptive policy change",
    }:
        raise error_type("candidate arithmetic policy changed")
    if candidate.get("backend") != _backend_record():
        raise error_type("candidate backend contract changed")
    if candidate.get("solver_role") != (
        "approximate-untrusted-candidate-generator-only"
    ):
        raise error_type("approximate solve was given a certificate role")
    if candidate.get("limitation") != CANDIDATE_LIMITATION:
        raise error_type("candidate limitation changed")

    numerators = _parse_numerators(candidate, error_type)
    vector_sha256 = content_sha256(_candidate_coefficients(numerators))
    if candidate.get("coefficient_vector_sha256") != vector_sha256:
        raise error_type("candidate coefficient vector hash mismatch")
    kernel = candidate.get("kernel")
    if not isinstance(kernel, Mapping) or set(kernel) != {
        "core_system_content_sha256",
        "prefix_256_numeric_sha256",
        "prefix_512_numeric_sha256",
        "old_prefix_matches_frozen_v1_numeric_kernel",
    }:
        raise error_type("candidate kernel binding changed")
    if kernel.get("old_prefix_matches_frozen_v1_numeric_kernel") is not True:
        raise error_type("candidate does not match the frozen N=256 numeric prefix")
    if kernel.get("prefix_256_numeric_sha256") != FROZEN_N256_PREFIX_KERNEL_SHA256:
        raise error_type("candidate N=256 numeric prefix hash changed")

    if require_frozen:
        frozen_values = {
            "candidate payload": FROZEN_CANDIDATE_PAYLOAD_SHA256,
            "coefficient vector": FROZEN_COEFFICIENT_VECTOR_SHA256,
            "generation kernel": FROZEN_GENERATION_KERNEL_SHA256,
            "generation prefix 512": FROZEN_GENERATION_PREFIX_512_SHA256,
        }
        if any(not value for value in frozen_values.values()):
            raise error_type("candidate freeze constants are incomplete")
        expected = {
            "candidate payload": candidate.get("payload_sha256"),
            "coefficient vector": vector_sha256,
            "generation kernel": kernel.get("core_system_content_sha256"),
            "generation prefix 512": kernel.get("prefix_512_numeric_sha256"),
        }
        for label, frozen in frozen_values.items():
            if expected[label] != frozen:
                raise error_type(f"frozen {label} hash changed")
    return numerators


def _validate_source_summary(
    summary: Mapping[str, Any],
    checkpoint_dir: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if summary.get("schema") != SUMMARY_SCHEMA:
        raise NymanBeta2Error("source summary schema changed")
    if summary.get("payload_sha256") != FROZEN_SUMMARY_PAYLOAD_SHA256:
        raise NymanBeta2Error("source summary payload hash changed")
    try:
        verification = verify_nyman_summary(summary, checkpoint_dir)
    except NymanSummaryError as exc:
        raise NymanBeta2Error("source Nyman summary verification failed") from exc
    expected_verification = {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
        "plan_sha256": FROZEN_PLAN_SHA256,
        "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
        "verified_cells": "6",
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }
    if verification != expected_verification:
        raise NymanBeta2Error("source summary verification result changed")

    cells = summary.get("cells")
    if not isinstance(cells, list):
        raise NymanBeta2Error("source summary cells are missing")
    matching = [cell for cell in cells if isinstance(cell, Mapping) and cell.get("n") == "256"]
    if len(matching) != 1:
        raise NymanBeta2Error("source summary must contain exactly one N=256 cell")
    source_cell = dict(matching[0])
    if source_cell.get("cell_id") != "n-0256" or source_cell.get("status") != (
        "FINITE_DISTANCE_BRACKET_CERTIFIED"
    ):
        raise NymanBeta2Error("source N=256 cell is not certified")
    if source_cell.get("decisions") != {
        "lower": "LOWER_BOUND_CERTIFIED",
        "upper": "UPPER_BOUND_CERTIFIED",
    }:
        raise NymanBeta2Error("source N=256 certificate decisions changed")
    hashes = source_cell.get("hashes")
    expected_hashes = {
        "candidate_sha256": FROZEN_N256_CANDIDATE_SHA256,
        "cell_contract_sha256": FROZEN_N256_CELL_CONTRACT_SHA256,
        "cell_payload_sha256": FROZEN_N256_CELL_PAYLOAD_SHA256,
        "lower_certificate_payload_sha256": FROZEN_N256_LOWER_CERTIFICATE_SHA256,
        "prefix_core_system_content_sha256": (
            "1587f9ab711c2b787f0072dcd82476cef3ee1ca0198ca167eccb0c86c95b8726"
        ),
        "prefix_kernel_sha256": FROZEN_N256_PREFIX_KERNEL_SHA256,
        "upper_certificate_payload_sha256": FROZEN_N256_UPPER_CERTIFICATE_SHA256,
    }
    if hashes != expected_hashes:
        raise NymanBeta2Error("source N=256 hash bindings changed")
    bounds = source_cell.get("bounds")
    if not isinstance(bounds, Mapping):
        raise NymanBeta2Error("source N=256 bounds are missing")
    lower = bounds.get("lower")
    upper = bounds.get("upper")
    if not isinstance(lower, Mapping) or not isinstance(upper, Mapping):
        raise NymanBeta2Error("source N=256 endpoint records are missing")
    if _source_fraction(lower.get("exact_fraction"), "source L256") != FROZEN_L256:
        raise NymanBeta2Error("source L256 changed")
    if _source_fraction(upper.get("exact_fraction"), "source U256") != FROZEN_U256:
        raise NymanBeta2Error("source U256 changed")
    if lower.get("dyadic") != _dyadic_record(
        FROZEN_L256.numerator, FROZEN_BOUND_GRID_BITS
    ):
        raise NymanBeta2Error("source L256 dyadic encoding changed")
    if upper.get("dyadic") != _dyadic_record(
        FROZEN_U256.numerator, FROZEN_BOUND_GRID_BITS
    ):
        raise NymanBeta2Error("source U256 dyadic encoding changed")
    return verification, source_cell


def _source_record(
    verification: Mapping[str, Any], source_cell: Mapping[str, Any]
) -> dict[str, Any]:
    return {
        "summary": {
            "schema": SUMMARY_SCHEMA,
            "payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
            "plan_sha256": FROZEN_PLAN_SHA256,
            "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
            "verification": dict(verification),
        },
        "n256": {
            "cell_id": "n-0256",
            "cell_payload_sha256": FROZEN_N256_CELL_PAYLOAD_SHA256,
            "cell_contract_sha256": FROZEN_N256_CELL_CONTRACT_SHA256,
            "candidate_sha256": FROZEN_N256_CANDIDATE_SHA256,
            "lower_certificate_payload_sha256": (
                FROZEN_N256_LOWER_CERTIFICATE_SHA256
            ),
            "upper_certificate_payload_sha256": (
                FROZEN_N256_UPPER_CERTIFICATE_SHA256
            ),
            "numeric_prefix_sha256": FROZEN_N256_PREFIX_KERNEL_SHA256,
            "lower_bound": _fraction_record(FROZEN_L256),
            "upper_bound": _fraction_record(FROZEN_U256),
            "terminal_statement": source_cell.get("terminal_statement"),
        },
        "source_verification_scope": (
            "structural regeneration of frozen v1 plus a fresh numerical "
            "N=256 lower-certificate replay inside the max-512 kernel"
        ),
    }


def _validate_lower_certificate(record: Mapping[str, Any]) -> None:
    if record.get("decision") != "LOWER_BOUND_CERTIFIED":
        raise NymanBeta2Error("fresh N=256 lower certificate is inconclusive")
    checks = record.get("checks")
    if not isinstance(checks, Mapping) or checks.get(
        "fixed_order_interval_ldlt_positive"
    ) is not True:
        raise NymanBeta2Error("fresh N=256 lower certificate check failed")
    if _source_fraction(
        record.get("claimed_lower_bound"), "fresh N=256 lower bound"
    ) != FROZEN_L256:
        raise NymanBeta2Error("fresh N=256 lower certificate changed its endpoint")


def _validate_upper_certificate(record: Mapping[str, Any]) -> None:
    if record.get("decision") != "UPPER_BOUND_CERTIFIED":
        raise NymanBeta2Error("N=512 direct upper certificate is inconclusive")
    checks = record.get("checks")
    if not isinstance(checks, Mapping) or checks.get(
        "energy_at_most_claimed_upper_bound"
    ) is not True or checks.get("energy_strictly_below_claimed_upper_bound") is not True:
        raise NymanBeta2Error("N=512 direct upper certificate is not strict")
    if _source_fraction(
        record.get("claimed_upper_bound"), "N=512 claimed upper bound"
    ) != FROZEN_U512:
        raise NymanBeta2Error("N=512 direct certificate changed its endpoint")


def _comparison_record() -> dict[str, Any]:
    strong_threshold = STRONG_CONTRACTION * FROZEN_L256
    beta_two_threshold = BETA_TWO_CONTRACTION * FROZEN_L256
    schur_safe_threshold = FROZEN_L256 - BETA_TWO_GAIN * FROZEN_U256
    strong_margin = strong_threshold - FROZEN_U512
    beta_two_margin = beta_two_threshold - FROZEN_U512
    schur_safe_margin = schur_safe_threshold - FROZEN_U512
    if min(strong_margin, beta_two_margin, schur_safe_margin) <= 0:
        raise NymanBeta2Error("frozen N=512 endpoint misses a contraction threshold")
    return {
        "stored_upper_bound_u512": {
            "dyadic": _dyadic_record(
                FROZEN_U512_NUMERATOR, FROZEN_BOUND_GRID_BITS
            ),
            "exact_fraction": _fraction_record(FROZEN_U512),
        },
        "source_lower_bound_l256": _fraction_record(FROZEN_L256),
        "source_upper_bound_u256": _fraction_record(FROZEN_U256),
        "strong_contraction_factor": _fraction_record(STRONG_CONTRACTION),
        "strong_gain_fraction": _fraction_record(STRONG_GAIN),
        "strong_threshold_449_over_500_times_l256": _fraction_record(
            strong_threshold
        ),
        "strong_exact_positive_margin": _fraction_record(strong_margin),
        "beta_two_contraction_factor": _fraction_record(BETA_TWO_CONTRACTION),
        "beta_two_gain_fraction": _fraction_record(BETA_TWO_GAIN),
        "beta_two_threshold_9_over_10_times_l256": _fraction_record(
            beta_two_threshold
        ),
        "beta_two_exact_positive_margin": _fraction_record(beta_two_margin),
        "schur_safe_threshold_l256_minus_u256_over_10": _fraction_record(
            schur_safe_threshold
        ),
        "schur_safe_exact_positive_margin": _fraction_record(schur_safe_margin),
        "checks": {
            "u512_strictly_below_449_over_500_l256": True,
            "u512_strictly_below_9_over_10_l256": True,
            "u512_strictly_below_l256_minus_u256_over_10": True,
            "449_over_500_strictly_below_9_over_10": True,
            "all_comparisons_exact_rational": True,
        },
        "certified_finite_chain": (
            "d_512^2 <= E(c_512) < U_512 < (449/500)*L_256 "
            "< (449/500)*d_256^2 < (9/10)*d_256^2"
        ),
        "certified_gain_statement": (
            "d_256^2-d_512^2 > (51/500)*d_256^2 > d_256^2/10"
        ),
    }


def _theorem_bridge() -> dict[str, Any]:
    return {
        "finite_schur_identity": {
            "block_system": "G_512=[[G_256,C],[C^T,D]], b_512=[b_256;e]",
            "definitions": (
                "S=D-C^T*G_256^-1*C; "
                "t=e-C^T*G_256^-1*b_256"
            ),
            "identity": "d_512^2=d_256^2-t^T*S^-1*t",
            "captured_gain": "q=t^T*S^-1*t=d_256^2-d_512^2",
            "machine_gate_role": "analytic identity only; no interval block solve is used",
        },
        "endpoint_implication": [
            "the fresh interval LDL certificate gives L_256 < d_256^2",
            "the stored exact witness and direct evaluation give d_512^2 < U_512",
            "exact Fraction arithmetic gives U_512 < (449/500)*L_256",
            "therefore d_512^2 < (449/500)*d_256^2",
            "therefore q > (51/500)*d_256^2 > d_256^2/10",
        ],
        "all_scale_target_not_proved": (
            "the k=8 step is certified here; still unproved for every k>=9: "
            "d_(2^(k+1))^2 <= ((k+1)/(k+2))*d_(2^k)^2"
        ),
        "conditional_telescoping_consequence": (
            "if the remaining k>=9 target were proved, then for k>=9, "
            "d_(2^k)^2 < 449*d_256^2/(50*(k+1)), hence the distances "
            "would tend to zero"
        ),
        "resolves_rh": False,
    }


def _derive_audit(
    candidate: Mapping[str, Any],
    summary: Mapping[str, Any],
    checkpoint_dir: Path,
) -> dict[str, Any]:
    numerators = _validate_candidate(
        candidate, NymanBeta2Error, require_frozen=True
    )
    verification, source_cell = _validate_source_summary(summary, checkpoint_dir)

    try:
        with _clean_precision(FROZEN_GENERATION_BITS):
            system = core.build_natural_system(
                tuple(range(1, FROZEN_NEW_N + 1))
            )
            kernel_sha256 = core.natural_system_content_sha256(system)
            prefix_256_sha256 = _prefix_kernel_sha256(system, FROZEN_OLD_N)
            prefix_512_sha256 = _prefix_kernel_sha256(system, FROZEN_NEW_N)
            candidate_kernel = candidate["kernel"]
            if kernel_sha256 != candidate_kernel["core_system_content_sha256"]:
                raise NymanBeta2Error(
                    "candidate generation kernel does not reproduce"
                )
            if prefix_256_sha256 != FROZEN_N256_PREFIX_KERNEL_SHA256:
                raise NymanBeta2Error(
                    "max-512 kernel does not reproduce the frozen N=256 prefix"
                )
            if prefix_512_sha256 != candidate_kernel["prefix_512_numeric_sha256"]:
                raise NymanBeta2Error(
                    "candidate N=512 numeric prefix does not reproduce"
                )

            prefix_256 = core.prefix_natural_system(system, FROZEN_OLD_N)
            lower_certificate = core.certify_augmented_lower_bound(
                tuple(range(1, FROZEN_OLD_N + 1)),
                FROZEN_L256,
                precision_bits=FROZEN_GENERATION_BITS,
                system=prefix_256,
            )
            upper_certificate = core.certify_dyadic_upper_bound(
                tuple(range(1, FROZEN_NEW_N + 1)),
                numerators,
                FROZEN_COEFFICIENT_GRID_BITS,
                FROZEN_U512,
                precision_bits=FROZEN_GENERATION_BITS,
                system=system,
            )
    except NymanBeta2Error:
        raise
    except Exception as exc:  # pragma: no cover - backend-specific failures
        raise NymanBeta2Error("finite contraction certification failed") from exc
    _validate_lower_certificate(lower_certificate)
    _validate_upper_certificate(upper_certificate)

    body = {
        "schema": AUDIT_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "source": _source_record(verification, source_cell),
        "candidate": dict(candidate),
        "generation": {
            "precision_bits": str(FROZEN_GENERATION_BITS),
            "backend": _backend_record(),
            "max_512_kernel": {
                "core_system_content_sha256": kernel_sha256,
                "prefix_256_numeric_sha256": prefix_256_sha256,
                "prefix_512_numeric_sha256": prefix_512_sha256,
                "prefix_nesting_check": True,
                "old_prefix_core_hash_equality_required": False,
                "old_prefix_core_hash_equality_note": (
                    "core provenance commits to source dimension; the separate "
                    "numeric-prefix hash is the valid equality check"
                ),
            },
            "source_n256_lower_certificate_replay": lower_certificate,
            "n512_direct_upper_certificate": upper_certificate,
            "exact_comparisons": _comparison_record(),
        },
        "theorem_bridge": _theorem_bridge(),
        "audit_outcome": "N512_STRONG_FINITE_CONTRACTION_CERTIFIED",
        "n512_lower_bound_attempted": False,
        "terminal_statement": (
            "d_512^2 < (449/500)*d_256^2; finite beta=2 step certified"
        ),
        "limitation": AUDIT_LIMITATION,
    }
    return _with_payload_hash(body)


def generate_nyman_beta2_audit(
    candidate: Mapping[str, Any] | Path,
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Generate the canonical finite 256-to-512 contraction artifact."""

    source_candidate = _load_object(
        candidate, NymanBeta2Error, label="Nyman beta=2 candidate"
    )
    source_summary = _load_object(
        summary, NymanBeta2Error, label="source Nyman summary"
    )
    artifact = _derive_audit(source_candidate, source_summary, checkpoint_dir)
    if artifact["payload_sha256"] != FROZEN_AUDIT_PAYLOAD_SHA256:
        raise NymanBeta2Error("frozen audit payload hash changed")
    return artifact


def _dyadic_interval(record: Any, name: str) -> tuple[Fraction, Fraction]:
    if not isinstance(record, Mapping) or set(record) != {
        "mid_mantissa",
        "mid_exponent",
        "radius_mantissa",
        "radius_exponent",
    }:
        raise NymanBeta2VerificationError(f"{name} dyadic enclosure changed")
    midpoint = _dyadic_fraction(record["mid_mantissa"], record["mid_exponent"])
    radius = _dyadic_fraction(record["radius_mantissa"], record["radius_exponent"])
    if radius < 0:
        raise NymanBeta2VerificationError(f"{name} has negative radius")
    return midpoint - radius, midpoint + radius


def _record_interval(record: Any, name: str) -> tuple[Fraction, Fraction]:
    if not isinstance(record, Mapping) or set(record) != {
        "dyadic",
        "display",
        "is_exact",
    }:
        raise NymanBeta2VerificationError(f"{name} Arb record changed")
    return _dyadic_interval(record["dyadic"], name)


def _assert_replay_contained(
    generation_record: Any,
    replay_record: Any,
    name: str,
) -> None:
    generation_lower, generation_upper = _record_interval(
        generation_record, f"generation {name}"
    )
    replay_lower, replay_upper = _record_interval(replay_record, f"replay {name}")
    if not (
        generation_lower <= replay_lower
        and replay_upper <= generation_upper
    ):
        raise NymanBeta2VerificationError(
            f"higher-precision {name} escaped the generation enclosure"
        )


def _preflight_audit(
    supplied: Mapping[str, Any],
    source_candidate: Mapping[str, Any],
) -> None:
    """Reject cheap structural mutations before rebuilding a 512 kernel."""

    _require_payload_hash(
        supplied, NymanBeta2VerificationError, label="Nyman beta=2 audit"
    )
    if not FROZEN_AUDIT_PAYLOAD_SHA256:
        raise NymanBeta2VerificationError("frozen audit payload hash is missing")
    if supplied.get("payload_sha256") != FROZEN_AUDIT_PAYLOAD_SHA256:
        raise NymanBeta2VerificationError("frozen audit payload hash changed")
    expected_fields = {
        "schema",
        "audit_id",
        "classification",
        "hypothesis_status",
        "source",
        "candidate",
        "generation",
        "theorem_bridge",
        "audit_outcome",
        "n512_lower_bound_attempted",
        "terminal_statement",
        "limitation",
        "payload_sha256",
    }
    if set(supplied) != expected_fields:
        raise NymanBeta2VerificationError("audit fields changed")
    expected_preamble = {
        "schema": AUDIT_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "N512_STRONG_FINITE_CONTRACTION_CERTIFIED",
        "n512_lower_bound_attempted": False,
        "terminal_statement": (
            "d_512^2 < (449/500)*d_256^2; finite beta=2 step certified"
        ),
        "limitation": AUDIT_LIMITATION,
    }
    for key, expected in expected_preamble.items():
        if supplied.get(key) != expected:
            raise NymanBeta2VerificationError(f"audit {key} changed")
    _validate_candidate(
        source_candidate,
        NymanBeta2VerificationError,
        require_frozen=True,
    )
    if supplied.get("candidate") != source_candidate:
        raise NymanBeta2VerificationError("embedded candidate changed")
    source = supplied.get("source")
    if not isinstance(source, Mapping):
        raise NymanBeta2VerificationError("audit source binding is missing")
    summary = source.get("summary")
    n256 = source.get("n256")
    if not isinstance(summary, Mapping) or not isinstance(n256, Mapping):
        raise NymanBeta2VerificationError("audit source records changed")
    if (
        summary.get("schema") != SUMMARY_SCHEMA
        or summary.get("payload_sha256") != FROZEN_SUMMARY_PAYLOAD_SHA256
        or summary.get("plan_sha256") != FROZEN_PLAN_SHA256
        or summary.get("index_payload_sha256") != FROZEN_INDEX_PAYLOAD_SHA256
    ):
        raise NymanBeta2VerificationError("audit summary binding changed")
    if (
        n256.get("cell_payload_sha256") != FROZEN_N256_CELL_PAYLOAD_SHA256
        or n256.get("cell_contract_sha256")
        != FROZEN_N256_CELL_CONTRACT_SHA256
        or n256.get("candidate_sha256") != FROZEN_N256_CANDIDATE_SHA256
        or n256.get("lower_certificate_payload_sha256")
        != FROZEN_N256_LOWER_CERTIFICATE_SHA256
        or n256.get("upper_certificate_payload_sha256")
        != FROZEN_N256_UPPER_CERTIFICATE_SHA256
        or n256.get("numeric_prefix_sha256")
        != FROZEN_N256_PREFIX_KERNEL_SHA256
        or n256.get("lower_bound") != _fraction_record(FROZEN_L256)
        or n256.get("upper_bound") != _fraction_record(FROZEN_U256)
    ):
        raise NymanBeta2VerificationError("audit N=256 binding changed")
    generation = supplied.get("generation")
    if not isinstance(generation, Mapping):
        raise NymanBeta2VerificationError("audit generation evidence is missing")
    if generation.get("precision_bits") != str(FROZEN_GENERATION_BITS):
        raise NymanBeta2VerificationError("audit generation precision changed")
    if generation.get("backend") != _backend_record():
        raise NymanBeta2VerificationError("audit backend contract changed")
    kernel = generation.get("max_512_kernel")
    if not isinstance(kernel, Mapping) or (
        kernel.get("core_system_content_sha256")
        != FROZEN_GENERATION_KERNEL_SHA256
        or kernel.get("prefix_256_numeric_sha256")
        != FROZEN_N256_PREFIX_KERNEL_SHA256
        or kernel.get("prefix_512_numeric_sha256")
        != FROZEN_GENERATION_PREFIX_512_SHA256
        or kernel.get("prefix_nesting_check") is not True
    ):
        raise NymanBeta2VerificationError("audit max-512 kernel binding changed")
    if generation.get("exact_comparisons") != _comparison_record():
        raise NymanBeta2VerificationError("audit exact comparisons changed")
    if supplied.get("theorem_bridge") != _theorem_bridge():
        raise NymanBeta2VerificationError("audit theorem bridge changed")


def verify_nyman_beta2_audit(
    artifact: Mapping[str, Any] | Path,
    candidate: Mapping[str, Any] | Path,
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
    *,
    replay_precision_bits: int = FROZEN_REPLAY_BITS,
) -> dict[str, Any]:
    """Regenerate at 768 bits and replay the stored witness at higher precision."""

    if (
        isinstance(replay_precision_bits, bool)
        or not isinstance(replay_precision_bits, int)
        or replay_precision_bits != FROZEN_REPLAY_BITS
    ):
        raise NymanBeta2VerificationError(
            f"replay precision must equal the frozen {FROZEN_REPLAY_BITS} bits"
        )
    supplied = _load_object(
        artifact, NymanBeta2VerificationError, label="Nyman beta=2 audit"
    )
    source_candidate = _load_object(
        candidate, NymanBeta2VerificationError, label="Nyman beta=2 candidate"
    )
    _preflight_audit(supplied, source_candidate)
    source_summary = _load_object(
        summary, NymanBeta2VerificationError, label="source Nyman summary"
    )
    try:
        expected = _derive_audit(source_candidate, source_summary, checkpoint_dir)
    except NymanBeta2Error as exc:
        raise NymanBeta2VerificationError(str(exc)) from exc
    if supplied != expected:
        raise NymanBeta2VerificationError(
            "Nyman beta=2 audit does not canonically regenerate"
        )

    numerators = _validate_candidate(
        source_candidate,
        NymanBeta2VerificationError,
        require_frozen=True,
    )
    try:
        with _clean_precision(replay_precision_bits):
            replay_system = core.build_natural_system(
                tuple(range(1, FROZEN_NEW_N + 1))
            )
            replay_prefix = core.prefix_natural_system(
                replay_system, FROZEN_OLD_N
            )
            replay_lower = core.certify_augmented_lower_bound(
                tuple(range(1, FROZEN_OLD_N + 1)),
                FROZEN_L256,
                precision_bits=replay_precision_bits,
                system=replay_prefix,
            )
            replay_upper = core.certify_dyadic_upper_bound(
                tuple(range(1, FROZEN_NEW_N + 1)),
                numerators,
                FROZEN_COEFFICIENT_GRID_BITS,
                FROZEN_U512,
                precision_bits=replay_precision_bits,
                system=replay_system,
            )
            replay_kernel_sha256 = core.natural_system_content_sha256(
                replay_system
            )
            replay_prefix_256_sha256 = _prefix_kernel_sha256(
                replay_system, FROZEN_OLD_N
            )
            replay_prefix_512_sha256 = _prefix_kernel_sha256(
                replay_system, FROZEN_NEW_N
            )
    except Exception as exc:  # pragma: no cover - backend-specific failures
        raise NymanBeta2VerificationError(
            "higher-precision contraction replay failed"
        ) from exc
    try:
        _validate_lower_certificate(replay_lower)
        _validate_upper_certificate(replay_upper)
    except NymanBeta2Error as exc:
        raise NymanBeta2VerificationError(str(exc)) from exc

    generation_upper = supplied["generation"]["n512_direct_upper_certificate"]
    _assert_replay_contained(
        generation_upper["energy"], replay_upper["energy"], "N=512 energy"
    )
    _assert_replay_contained(
        generation_upper["margin"], replay_upper["margin"], "N=512 upper margin"
    )
    return {
        "classification": "REPRODUCED_CERTIFIED_FINITE_NYMAN_CONTRACTION",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": supplied["audit_outcome"],
        "artifact_payload_sha256": supplied["payload_sha256"],
        "candidate_payload_sha256": source_candidate["payload_sha256"],
        "generation_precision_bits": str(FROZEN_GENERATION_BITS),
        "replay_precision_bits": str(replay_precision_bits),
        "generation_kernel_sha256": supplied["generation"]["max_512_kernel"][
            "core_system_content_sha256"
        ],
        "replay_kernel_sha256": replay_kernel_sha256,
        "replay_prefix_256_numeric_sha256": replay_prefix_256_sha256,
        "replay_prefix_512_numeric_sha256": replay_prefix_512_sha256,
        "stored_candidate_replayed_without_regeneration": True,
        "source_n256_lower_replayed": True,
        "n512_direct_upper_replayed": True,
        "generation_energy_contains_replay": True,
        "generation_margin_contains_replay": True,
        "same_backend_replay_only": True,
    }


__all__ = [
    "AUDIT_SCHEMA",
    "AUDIT_LIMITATION",
    "CANDIDATE_SCHEMA",
    "CANDIDATE_LIMITATION",
    "FROZEN_AUDIT_ID",
    "FROZEN_AUDIT_PAYLOAD_SHA256",
    "FROZEN_BOUND_GRID_BITS",
    "FROZEN_CANDIDATE_ID",
    "FROZEN_COEFFICIENT_GRID_BITS",
    "FROZEN_GENERATION_BITS",
    "FROZEN_L256",
    "FROZEN_NEW_N",
    "FROZEN_OLD_N",
    "FROZEN_REPLAY_BITS",
    "FROZEN_U256",
    "FROZEN_U512",
    "NymanBeta2Error",
    "NymanBeta2VerificationError",
    "generate_nyman_beta2_audit",
    "propose_nyman_beta2_candidate",
    "verify_nyman_beta2_audit",
]
