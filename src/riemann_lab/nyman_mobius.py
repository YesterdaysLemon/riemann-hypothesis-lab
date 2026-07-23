"""Finite core/tail audits for the log-tapered Mobius Nyman family.

For literature coefficients

``a_n = mu(n) * (1 - log(n) / log(N))``

the repository coefficient is ``c_n = -a_n``.  In the reciprocal variable
``t = 1/x`` the residual is

``R_N(t) = 1_[1,infinity)(t) + sum_(n<=N) a_n {t/n}``.

This module certifies finite evaluations of its total Gram energy, its core
integral over ``0 < t < N``, and the complementary tail.  It also checks the
exact divisor identities that turn the core into a weighted prime-counting
error.  Nothing here controls the limit as ``N`` tends to infinity, so every
artifact remains exploratory and RH remains unresolved.
"""

from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping

import flint
from flint import arb, ctx

from .artifacts import content_sha256
from .balls import arb_from_dyadic, arb_record
from .nyman import (
    NaturalSystem,
    build_natural_system,
    natural_system_content_sha256,
    prefix_natural_system,
)
from .nyman_summary import (
    NymanSummaryError,
    SUMMARY_SCHEMA,
    verify_nyman_summary,
)


MOBIUS_SCHEMA = "rh-lab/nyman-mobius-core-tail-audit/v1"
FROZEN_AUDIT_ID = "nyman-mobius-core-tail-v1"
FROZEN_CELL_SIZES = (8, 16, 32, 64, 128, 256)
FROZEN_MAX_N = 256
FROZEN_GENERATION_BITS = 256
FROZEN_REPLAY_BITS = 512
FROZEN_GENERATION_KERNEL_SHA256 = (
    "ad7a07732ba1019584bc0983861a34cca089a54d6046503ae3616c87b44e13c9"
)
FROZEN_SUMMARY_PAYLOAD_SHA256 = (
    "35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf"
)
FROZEN_INDEX_PAYLOAD_SHA256 = (
    "281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23"
)
FROZEN_PLAN_SHA256 = (
    "27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a"
)
FROZEN_INTEGER_RATIO_SEPARATORS = {
    8: 34,
    16: 26,
    32: 21,
    64: 19,
    128: 16,
    256: 15,
}

MOBIUS_LIMITATION = (
    "This EXPLORATORY artifact evaluates one explicit log-tapered Mobius "
    "family at six finite dimensions. It proves exact finite divisor "
    "identities and finite Arb enclosures, not an all-N upper bound. The "
    "prime-counting formula remains valid through floor(t)=N but generally "
    "fails for floor(t)>N; the tail uses truncated divisors and is recovered "
    "from the full Gram energy. "
    "Finite decreases or fitted rates cannot prove or disprove the Riemann "
    "Hypothesis; the global status remains UNRESOLVED."
)


class NymanMobiusError(ValueError):
    """Raised when the frozen Mobius audit cannot be generated safely."""


class NymanMobiusVerificationError(NymanMobiusError):
    """Raised when a supplied Mobius audit does not exactly regenerate."""


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


def _arb_from_fraction(value: Fraction) -> Any:
    return arb(value.numerator) / value.denominator


def _snapshot_json_tree(
    value: Any,
    error_type: type[NymanMobiusError],
    *,
    location: str,
) -> Any:
    """Detach an arbitrary Mapping tree into strict JSON built-ins."""

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
    error_type: type[NymanMobiusError],
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
    error_type: type[NymanMobiusError],
) -> None:
    supplied = artifact.get("payload_sha256")
    body = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied != content_sha256(body):
        raise error_type("Mobius audit payload hash mismatch")


def _canonical_integer(
    value: Any,
    name: str,
    *,
    minimum: int | None = None,
) -> int:
    if not isinstance(value, str):
        raise NymanMobiusError(f"{name} must be a canonical integer string")
    try:
        parsed = int(value)
    except ValueError as exc:
        raise NymanMobiusError(f"{name} is not an integer") from exc
    if value != str(parsed):
        raise NymanMobiusError(f"{name} is not canonical")
    if minimum is not None and parsed < minimum:
        raise NymanMobiusError(f"{name} must be at least {minimum}")
    return parsed


def _source_fraction(value: Any, name: str) -> Fraction:
    if not isinstance(value, Mapping) or set(value) != {"numerator", "denominator"}:
        raise NymanMobiusError(f"{name} has noncanonical fraction fields")
    numerator = _canonical_integer(value["numerator"], f"{name} numerator")
    denominator = _canonical_integer(
        value["denominator"],
        f"{name} denominator",
        minimum=1,
    )
    result = Fraction(numerator, denominator)
    if dict(value) != _fraction_record(result):
        raise NymanMobiusError(f"{name} is not a reduced canonical fraction")
    return result


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


def _validate_frozen_summary(
    summary: Mapping[str, Any],
) -> dict[int, Mapping[str, Any]]:
    if summary.get("schema") != SUMMARY_SCHEMA:
        raise NymanMobiusError("source summary schema changed")
    if summary.get("payload_sha256") != FROZEN_SUMMARY_PAYLOAD_SHA256:
        raise NymanMobiusError("source summary payload hash is not frozen")

    source = summary.get("source")
    if not isinstance(source, Mapping):
        raise NymanMobiusError("source summary provenance is missing")
    plan = source.get("plan")
    index = source.get("index")
    if not isinstance(plan, Mapping) or plan.get("plan_sha256") != FROZEN_PLAN_SHA256:
        raise NymanMobiusError("source summary plan hash is not frozen")
    if (
        not isinstance(index, Mapping)
        or index.get("payload_sha256") != FROZEN_INDEX_PAYLOAD_SHA256
    ):
        raise NymanMobiusError("source summary index hash is not frozen")

    cells = summary.get("cells")
    if not isinstance(cells, list) or len(cells) != len(FROZEN_CELL_SIZES):
        raise NymanMobiusError("source summary cell list changed")
    by_n: dict[int, Mapping[str, Any]] = {}
    for cell in cells:
        if not isinstance(cell, Mapping):
            raise NymanMobiusError("source summary cell must be an object")
        n = _canonical_integer(cell.get("n"), "source cell n", minimum=1)
        if n in by_n:
            raise NymanMobiusError("source summary contains a duplicate cell")
        if cell.get("status") != "FINITE_DISTANCE_BRACKET_CERTIFIED":
            raise NymanMobiusError("source summary contains an uncertified cell")
        by_n[n] = cell
    if tuple(by_n) != FROZEN_CELL_SIZES:
        raise NymanMobiusError("source summary cells are not the frozen ordered six")
    return by_n


def _cell_bound(cell: Mapping[str, Any], side: str) -> Fraction:
    bounds = cell.get("bounds")
    if not isinstance(bounds, Mapping):
        raise NymanMobiusError("source cell bounds are missing")
    bound = bounds.get(side)
    if not isinstance(bound, Mapping) or set(bound) != {"dyadic", "exact_fraction"}:
        raise NymanMobiusError(f"source cell {side} bound changed")
    return _source_fraction(
        bound["exact_fraction"],
        f"source cell {cell.get('n')} {side} bound",
    )


def _mobius_sieve(limit: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes: list[int] = []
    composite = [False] * (limit + 1)
    for value in range(2, limit + 1):
        if not composite[value]:
            primes.append(value)
            mu[value] = -1
        for prime in primes:
            product = value * prime
            if product > limit:
                break
            composite[product] = True
            if value % prime == 0:
                mu[product] = 0
                break
            mu[product] = -mu[value]
    return tuple(mu), tuple(primes)


def _prime_power_bases(limit: int, primes: Iterable[int]) -> tuple[int, ...]:
    bases = [0] * (limit + 1)
    for prime in primes:
        power = prime
        while power <= limit:
            bases[power] = prime
            power *= prime
    return tuple(bases)


def _divisors(value: int) -> Iterable[int]:
    for divisor in range(1, value + 1):
        if value % divisor == 0:
            yield divisor


def _arithmetic_manifest() -> tuple[dict[str, Any], tuple[int, ...], tuple[int, ...]]:
    mu, primes = _mobius_sieve(FROZEN_MAX_N)
    prime_power_base = _prime_power_bases(FROZEN_MAX_N, primes)
    transcript: list[dict[str, Any]] = []
    prime_power_count = 0
    for value in range(1, FROZEN_MAX_N + 1):
        divisor_sum = 0
        log_product = Fraction(1)
        for divisor in _divisors(value):
            coefficient = mu[divisor]
            divisor_sum += coefficient
            if coefficient == 1:
                log_product *= divisor
            elif coefficient == -1:
                log_product /= divisor
        base = prime_power_base[value]
        expected_sum = 1 if value == 1 else 0
        expected_product = Fraction(1, base) if base else Fraction(1)
        if divisor_sum != expected_sum or log_product != expected_product:
            raise NymanMobiusError(
                f"exact Mobius divisor identity failed at k={value}"
            )
        if base:
            prime_power_count += 1
        transcript.append(
            {
                "k": str(value),
                "mobius_divisor_sum": str(divisor_sum),
                "mobius_log_product": _fraction_record(log_product),
                "prime_power_base": str(base) if base else None,
            }
        )

    manifest = {
        "max_n": str(FROZEN_MAX_N),
        "mobius_values_n_1_through_256": [
            str(mu[value]) for value in range(1, FROZEN_MAX_N + 1)
        ],
        "prime_power_base_n_1_through_256": [
            str(prime_power_base[value]) if prime_power_base[value] else None
            for value in range(1, FROZEN_MAX_N + 1)
        ],
        "identity_transcript_sha256": content_sha256(transcript),
        "checks": {
            "checked_k": str(FROZEN_MAX_N),
            "prime_power_count": str(prime_power_count),
            "sum_mu_over_divisors_is_delta_k1": True,
            "product_d_power_mu_d_is_inverse_prime_exactly_at_prime_powers": True,
            "derived_log_identity": (
                "sum_(d|k) mu(d)*log(d) = -Lambda(k)"
            ),
            "derived_bcf_divisor_identity_through_cutoff": (
                "h_N(k)=1_(k=1)+Lambda(k)/log(N) for k<=N"
            ),
        },
    }
    return manifest, mu, prime_power_base


def _bcf_coefficients(mu: tuple[int, ...], n: int) -> tuple[Any, ...]:
    log_n = arb(n).log()
    values = []
    for index in range(1, n + 1):
        if index == n or mu[index] == 0:
            values.append(arb(0))
        else:
            values.append(
                arb(mu[index]) * (1 - arb(index).log() / log_n)
            )
    return tuple(values)


def _candidate_energy(system: NaturalSystem, coefficients: tuple[Any, ...]) -> Any:
    """Evaluate ``1 + 2 a^T b + a^T G a`` with Arb coefficients."""

    size = len(coefficients)
    energy = arb(1)
    for row in range(size):
        value = coefficients[row]
        energy += 2 * system.target[row] * value
        energy += system.gram[row][row] * value * value
        for column in range(row + 1, size):
            energy += (
                2
                * system.gram[row][column]
                * value
                * coefficients[column]
            )
    return energy


def _core_values(
    coefficients: tuple[Any, ...],
    prime_power_base: tuple[int, ...],
    n: int,
) -> tuple[Any, Any, Any, Any, Any]:
    log_n = arb(n).log()
    slope = sum(
        (coefficients[index - 1] / index for index in range(1, n + 1)),
        arb(0),
    )
    kappa = log_n * slope
    general_core = slope * slope
    prime_core = slope * slope
    cumulative_h = arb(0)
    psi = arb(0)

    for interval in range(1, n):
        h_value = sum(
            (
                coefficients[divisor - 1]
                for divisor in _divisors(interval)
            ),
            arb(0),
        )
        cumulative_h += h_value
        general_q = 1 - cumulative_h
        if prime_power_base[interval]:
            psi += arb(prime_power_base[interval]).log()
        prime_q = -psi / log_n
        interval_log = (arb(interval + 1) / interval).log()
        general_core += (
            slope * slope
            + 2 * slope * general_q * interval_log
            + general_q * general_q / (interval * (interval + 1))
        )
        prime_core += (
            slope * slope
            + 2 * slope * prime_q * interval_log
            + prime_q * prime_q / (interval * (interval + 1))
        )

    difference = general_core - prime_core
    if not general_core.overlaps(prime_core) or not difference.contains(0):
        raise NymanMobiusError(
            f"general and prime-counting core paths disagree at N={n}"
        )
    return slope, kappa, general_core, prime_core, general_core.intersection(prime_core)


def _compute_numerical_bundle(precision_bits: int) -> dict[str, Any]:
    arithmetic, mu, prime_power_base = _arithmetic_manifest()
    with _clean_precision(precision_bits):
        full_system = build_natural_system(range(1, FROZEN_MAX_N + 1))
        kernel_hash = natural_system_content_sha256(full_system)
        if (
            precision_bits == FROZEN_GENERATION_BITS
            and kernel_hash != FROZEN_GENERATION_KERNEL_SHA256
        ):
            raise NymanMobiusError("frozen generation kernel commitment changed")

        records: list[dict[str, Any]] = []
        raw: list[dict[str, Any]] = []
        for n in FROZEN_CELL_SIZES:
            system = prefix_natural_system(full_system, n)
            coefficients = _bcf_coefficients(mu, n)
            slope, kappa, general_core, prime_core, core = _core_values(
                coefficients,
                prime_power_base,
                n,
            )
            total = _candidate_energy(system, coefficients)
            tail = total - core
            log_n = arb(n).log()
            scaled_total = total * log_n
            scaled_core = core * log_n
            scaled_tail = tail * log_n
            core_share = core / total
            checks = {
                "candidate_total_strictly_positive": bool(total > 0),
                "core_strictly_positive": bool(core > 0),
                "tail_strictly_positive": bool(tail > 0),
                "core_strictly_below_total": bool(core < total),
                "general_core_overlaps_prime_core": general_core.overlaps(prime_core),
                "general_minus_prime_core_contains_zero": (
                    (general_core - prime_core).contains(0)
                ),
                "endpoint_coefficient_is_exact_zero": coefficients[-1].is_zero(),
            }
            if not all(checks.values()):
                raise NymanMobiusError(f"finite Mobius audit did not close at N={n}")
            records.append(
                {
                    "n": str(n),
                    "log_n": arb_record(log_n),
                    "slope_s_n": arb_record(slope),
                    "kappa_n_equals_log_n_times_s_n": arb_record(kappa),
                    "candidate_total_error_squared": arb_record(total),
                    "general_divisor_core_error_squared": arb_record(general_core),
                    "prime_counting_core_error_squared": arb_record(prime_core),
                    "certified_core_intersection": arb_record(core),
                    "tail_error_squared_by_total_minus_core": arb_record(tail),
                    "scaled_candidate_total_times_log_n": arb_record(scaled_total),
                    "scaled_core_times_log_n": arb_record(scaled_core),
                    "scaled_tail_times_log_n": arb_record(scaled_tail),
                    "core_fraction_of_candidate_total": arb_record(core_share),
                    "checks": checks,
                }
            )
            raw.append(
                {
                    "total": total,
                    "core": core,
                    "tail": tail,
                    "scaled_total": scaled_total,
                    "scaled_core": scaled_core,
                    "scaled_tail": scaled_tail,
                    "core_share": core_share,
                }
            )

        trend_checks = {
            "candidate_total_strictly_decreases_on_grid": all(
                raw[index - 1]["total"] > raw[index]["total"]
                for index in range(1, len(raw))
            ),
            "core_strictly_decreases_on_grid": all(
                raw[index - 1]["core"] > raw[index]["core"]
                for index in range(1, len(raw))
            ),
            "tail_strictly_decreases_on_grid": all(
                raw[index - 1]["tail"] > raw[index]["tail"]
                for index in range(1, len(raw))
            ),
            "scaled_candidate_total_strictly_decreases_on_grid": all(
                raw[index - 1]["scaled_total"] > raw[index]["scaled_total"]
                for index in range(1, len(raw))
            ),
            "scaled_core_strictly_decreases_on_grid": all(
                raw[index - 1]["scaled_core"] > raw[index]["scaled_core"]
                for index in range(1, len(raw))
            ),
            "scaled_tail_strictly_decreases_on_grid": all(
                raw[index - 1]["scaled_tail"] > raw[index]["scaled_tail"]
                for index in range(1, len(raw))
            ),
            "core_fraction_strictly_increases_on_grid": all(
                raw[index - 1]["core_share"] < raw[index]["core_share"]
                for index in range(1, len(raw))
            ),
        }
        if not all(trend_checks.values()):
            raise NymanMobiusError("a frozen finite trend did not certify")

        return {
            "precision_bits": str(precision_bits),
            "backend": _backend_record(),
            "kernel": {
                "construction": (
                    "canonical-rational-vasyunin-with-reciprocity-intersection"
                ),
                "max_n": str(FROZEN_MAX_N),
                "system_content_sha256": kernel_hash,
            },
            "exact_arithmetic": arithmetic,
            "records": records,
            "finite_grid_trends": trend_checks,
        }


@lru_cache(maxsize=2)
def _cached_numerical_json(precision_bits: int) -> str:
    """Cache only immutable JSON text, never mutable artifacts or Arb values."""

    return json.dumps(
        _compute_numerical_bundle(precision_bits),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _numerical_bundle(precision_bits: int) -> dict[str, Any]:
    value = json.loads(_cached_numerical_json(precision_bits))
    if not isinstance(value, dict):
        raise NymanMobiusError("cached numerical bundle is not an object")
    return value


def _source_comparisons(
    records: list[Mapping[str, Any]],
    by_n: Mapping[int, Mapping[str, Any]],
    precision_bits: int,
) -> list[dict[str, Any]]:
    comparisons: list[dict[str, Any]] = []
    with _clean_precision(precision_bits):
        for record in records:
            n = _canonical_integer(record.get("n"), "numerical record n", minimum=1)
            lower = _cell_bound(by_n[n], "lower")
            upper = _cell_bound(by_n[n], "upper")
            if not lower < upper:
                raise NymanMobiusError("source optimized-distance bracket is empty")
            total_record = record.get("candidate_total_error_squared")
            if not isinstance(total_record, Mapping):
                raise NymanMobiusError("candidate total record is missing")
            total = arb_from_dyadic(dict(total_record["dyadic"]))
            separator = FROZEN_INTEGER_RATIO_SEPARATORS[n]
            strict_margin = total - separator * _arb_from_fraction(upper)
            if not strict_margin > 0:
                raise NymanMobiusError(
                    f"candidate/optimum integer separation failed at N={n}"
                )
            comparisons.append(
                {
                    "n": str(n),
                    "source_optimal_lower_bound": _fraction_record(lower),
                    "source_optimal_upper_bound": _fraction_record(upper),
                    "candidate_over_source_upper": arb_record(
                        total / _arb_from_fraction(upper)
                    ),
                    "candidate_over_source_lower": arb_record(
                        total / _arb_from_fraction(lower)
                    ),
                    "certified_strict_integer_factor": str(separator),
                    "strict_factor_margin": arb_record(strict_margin),
                    "derived_inequality": (
                        f"E_BCF({n}) > {separator}*d_{n}^2"
                    ),
                }
            )
    return comparisons


def _dependency_manifest() -> dict[str, Any]:
    return {
        "published_sources": [
            {
                "key": "baez-duarte-strong-natural-dilate-criterion",
                "url": "https://arxiv.org/abs/math/0202141",
                "role": "d_N tending to zero is equivalent to RH",
                "reproved_by_this_artifact": False,
            },
            {
                "key": "bettin-conrey-farmer-log-taper",
                "url": "https://arxiv.org/abs/1211.5191",
                "role": (
                    "source of the log-tapered Mobius family and its "
                    "conditional asymptotic analysis"
                ),
                "reproved_by_this_artifact": False,
            },
            {
                "key": "baez-duarte-balazard-landreau-saias-kernel",
                "url": "https://arxiv.org/abs/math/0306251",
                "role": "rational autocorrelation formula used by the Gram kernel",
                "reproved_by_this_artifact": False,
            },
        ],
        "manual_exact_identity": {
            "general_residual": (
                "R(t)=t*S for 0<t<1 and R(t)=t*S+1-H(floor(t)) "
                "for t>=1"
            ),
            "core_formula": (
                "C_N=S^2+sum_(m=1)^(N-1)[S^2+2*S*q_m*log((m+1)/m)"
                "+q_m^2/(m*(m+1))]"
            ),
            "bcf_core_residual": (
                "R(t)=(kappa_N*t-psi(floor(t)))/log(N) for 1<=t<N"
            ),
            "tail_boundary": (
                "the psi shortcut generally fails for floor(t)>N because "
                "divisors above N are absent"
            ),
            "machine_checked_finitely": True,
        },
        "decisive_all_scale_target": {
            "statement": (
                "There exist K and N0 such that E_BCF(N)<=K/log(N) for every "
                "N>=N0; equivalently log(N)^2*(C_N+T_N)<=K*log(N)."
            ),
            "consequence": (
                "Then d_N^2<=E_BCF(N) tends to zero, so the strong natural-"
                "dilate criterion proves RH. A dyadic all-scale version also "
                "suffices by monotonicity of d_N."
            ),
            "proved_by_this_artifact": False,
        },
    }


def generate_nyman_mobius_audit(
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Generate the canonical six-cell log-tapered Mobius core/tail audit."""

    source_summary = _load_object(
        summary,
        NymanMobiusError,
        label="source Nyman summary",
    )
    try:
        source_verification = verify_nyman_summary(source_summary, checkpoint_dir)
    except NymanSummaryError as exc:
        raise NymanMobiusError("source Nyman summary verification failed") from exc
    expected_source_verification = {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_SUMMARY",
        "hypothesis_status": "UNRESOLVED",
        "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
        "plan_sha256": FROZEN_PLAN_SHA256,
        "index_payload_sha256": FROZEN_INDEX_PAYLOAD_SHA256,
        "verified_cells": "6",
        "numerical_replay_performed": False,
        "generation_kernel_rebuilt": False,
    }
    if source_verification != expected_source_verification:
        raise NymanMobiusError("source Nyman summary verification result changed")
    by_n = _validate_frozen_summary(source_summary)

    numerical = _numerical_bundle(FROZEN_GENERATION_BITS)
    records = numerical.get("records")
    if not isinstance(records, list):
        raise NymanMobiusError("numerical record list is missing")
    payload = {
        "schema": MOBIUS_SCHEMA,
        "audit_id": FROZEN_AUDIT_ID,
        "classification": "EXPLORATORY",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "FINITE_EXPLICIT_FAMILY_ENERGY_CERTIFIED",
        "family": {
            "literature_coefficients": (
                "a_n=mu(n)*(1-log(n)/log(N)), 1<=n<=N"
            ),
            "repository_coefficients": "c_n=-a_n",
            "reciprocal_variable_residual": (
                "R_N(t)=1_[1,infinity)(t)+sum_(n<=N) a_n*{t/n}"
            ),
            "candidate_relation": "d_N^2 <= E_BCF(N)=||R_N||^2",
            "endpoint_coefficient": "a_N=0 exactly",
        },
        "source_verification": source_verification,
        "generation": numerical,
        "candidate_vs_certified_optimum": _source_comparisons(
            records,
            by_n,
            FROZEN_GENERATION_BITS,
        ),
        "dependency_manifest": _dependency_manifest(),
        "finite_conclusion": (
            "At N=8,16,32,64,128,256, two finite core formulas agree, the "
            "Gram total splits into a strictly positive core and tail, and "
            "the total, core, tail, and their log-scaled values all strictly "
            "decrease across this grid. At N=256 the core exceeds 98 percent "
            "of the explicit-family energy and E_BCF(256)>15*d_256^2."
        ),
        "limitation": MOBIUS_LIMITATION,
    }
    return _with_payload_hash(payload)


def _exact_binary_fraction(mantissa: Any, exponent: Any, name: str) -> Fraction:
    coefficient = _canonical_integer(mantissa, f"{name} mantissa")
    power = _canonical_integer(exponent, f"{name} exponent")
    if power >= 0:
        return Fraction(coefficient << power)
    return Fraction(coefficient, 1 << -power)


def _exact_dyadic_endpoints(record: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    if set(record) != {"display", "dyadic", "is_exact"}:
        raise NymanMobiusVerificationError("Arb replay record fields changed")
    encoded = record.get("dyadic")
    if not isinstance(encoded, Mapping) or set(encoded) != {
        "mid_mantissa",
        "mid_exponent",
        "radius_mantissa",
        "radius_exponent",
    }:
        raise NymanMobiusVerificationError("Arb replay dyadic fields changed")
    midpoint = _exact_binary_fraction(
        encoded["mid_mantissa"],
        encoded["mid_exponent"],
        "Arb replay midpoint",
    )
    radius = _exact_binary_fraction(
        encoded["radius_mantissa"],
        encoded["radius_exponent"],
        "Arb replay radius",
    )
    if radius < 0:
        raise NymanMobiusVerificationError("Arb replay radius is negative")
    if record.get("is_exact") is not (radius == 0):
        raise NymanMobiusVerificationError("Arb replay exactness flag changed")
    return midpoint - radius, midpoint + radius


def _assert_replay_contains_generation(
    generation: Mapping[str, Any],
    replay: Mapping[str, Any],
) -> None:
    generation_records = generation.get("records")
    replay_records = replay.get("records")
    if not isinstance(generation_records, list) or not isinstance(replay_records, list):
        raise NymanMobiusVerificationError("replay record list is missing")
    if len(generation_records) != len(replay_records):
        raise NymanMobiusVerificationError("replay record count changed")
    fields = (
        "log_n",
        "slope_s_n",
        "kappa_n_equals_log_n_times_s_n",
        "candidate_total_error_squared",
        "general_divisor_core_error_squared",
        "prime_counting_core_error_squared",
        "certified_core_intersection",
        "tail_error_squared_by_total_minus_core",
        "scaled_candidate_total_times_log_n",
        "scaled_core_times_log_n",
        "scaled_tail_times_log_n",
        "core_fraction_of_candidate_total",
    )
    for low_record, high_record in zip(generation_records, replay_records, strict=True):
        if low_record.get("n") != high_record.get("n"):
            raise NymanMobiusVerificationError("replay N ordering changed")
        for field in fields:
            low_value = low_record.get(field)
            high_value = high_record.get(field)
            if not isinstance(low_value, Mapping) or not isinstance(high_value, Mapping):
                raise NymanMobiusVerificationError(f"replay field {field} is missing")
            low_lower, low_upper = _exact_dyadic_endpoints(low_value)
            high_lower, high_upper = _exact_dyadic_endpoints(high_value)
            if not (low_lower <= high_lower and high_upper <= low_upper):
                raise NymanMobiusVerificationError(
                    f"512-bit replay escaped 256-bit enclosure for N="
                    f"{low_record.get('n')} field {field}"
                )


def _assert_comparison_replay_contains_generation(
    generation: list[Mapping[str, Any]],
    replay: list[Mapping[str, Any]],
) -> None:
    if len(generation) != len(replay):
        raise NymanMobiusVerificationError("comparison replay count changed")
    fields = (
        "candidate_over_source_upper",
        "candidate_over_source_lower",
        "strict_factor_margin",
    )
    for low_record, high_record in zip(generation, replay, strict=True):
        if low_record.get("n") != high_record.get("n"):
            raise NymanMobiusVerificationError("comparison replay N ordering changed")
        for field in fields:
            low_value = low_record.get(field)
            high_value = high_record.get(field)
            if not isinstance(low_value, Mapping) or not isinstance(high_value, Mapping):
                raise NymanMobiusVerificationError(
                    f"comparison replay field {field} is missing"
                )
            low_lower, low_upper = _exact_dyadic_endpoints(low_value)
            high_lower, high_upper = _exact_dyadic_endpoints(high_value)
            if not (low_lower <= high_lower and high_upper <= low_upper):
                raise NymanMobiusVerificationError(
                    f"512-bit comparison replay escaped 256-bit enclosure for "
                    f"N={low_record.get('n')} field {field}"
                )


def verify_nyman_mobius_audit(
    artifact: Mapping[str, Any] | Path,
    summary: Mapping[str, Any] | Path,
    checkpoint_dir: Path,
) -> dict[str, Any]:
    """Exactly regenerate at 256 bits and replay all scalars at 512 bits."""

    supplied = _load_object(
        artifact,
        NymanMobiusVerificationError,
        label="Nyman Mobius audit",
    )
    _require_payload_hash(supplied, NymanMobiusVerificationError)
    source_summary = _load_object(
        summary,
        NymanMobiusVerificationError,
        label="source Nyman summary",
    )
    expected = generate_nyman_mobius_audit(source_summary, checkpoint_dir)
    if supplied != expected:
        raise NymanMobiusVerificationError(
            "Nyman Mobius audit does not canonically regenerate"
        )
    replay = _numerical_bundle(FROZEN_REPLAY_BITS)
    generation = supplied.get("generation")
    if not isinstance(generation, Mapping):
        raise NymanMobiusVerificationError("generation bundle is missing")
    with _clean_precision(FROZEN_REPLAY_BITS):
        _assert_replay_contains_generation(generation, replay)
    by_n = _validate_frozen_summary(source_summary)
    replay_records = replay.get("records")
    generation_comparisons = supplied.get("candidate_vs_certified_optimum")
    if not isinstance(replay_records, list) or not isinstance(
        generation_comparisons,
        list,
    ):
        raise NymanMobiusVerificationError("comparison replay inputs are missing")
    replay_comparisons = _source_comparisons(
        replay_records,
        by_n,
        FROZEN_REPLAY_BITS,
    )
    with _clean_precision(FROZEN_REPLAY_BITS):
        _assert_comparison_replay_contains_generation(
            generation_comparisons,
            replay_comparisons,
        )
    replay_kernel = replay.get("kernel")
    if not isinstance(replay_kernel, Mapping):
        raise NymanMobiusVerificationError("replay kernel record is missing")
    return {
        "classification": "REPRODUCED_EXPLORATORY_NYMAN_MOBIUS_AUDIT",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": supplied["audit_outcome"],
        "artifact_payload_sha256": supplied["payload_sha256"],
        "summary_payload_sha256": FROZEN_SUMMARY_PAYLOAD_SHA256,
        "verified_cells": str(len(FROZEN_CELL_SIZES)),
        "generation_precision_bits": str(FROZEN_GENERATION_BITS),
        "replay_precision_bits": str(FROZEN_REPLAY_BITS),
        "replay_kernel_system_content_sha256": replay_kernel[
            "system_content_sha256"
        ],
        "all_generation_enclosures_contain_replay": True,
        "all_comparison_enclosures_contain_replay": True,
        "same_backend_replay_only": True,
    }


__all__ = [
    "FROZEN_AUDIT_ID",
    "FROZEN_CELL_SIZES",
    "FROZEN_GENERATION_BITS",
    "FROZEN_GENERATION_KERNEL_SHA256",
    "FROZEN_INDEX_PAYLOAD_SHA256",
    "FROZEN_MAX_N",
    "FROZEN_PLAN_SHA256",
    "FROZEN_REPLAY_BITS",
    "FROZEN_SUMMARY_PAYLOAD_SHA256",
    "MOBIUS_LIMITATION",
    "MOBIUS_SCHEMA",
    "NymanMobiusError",
    "NymanMobiusVerificationError",
    "generate_nyman_mobius_audit",
    "verify_nyman_mobius_audit",
]
