"""Exact finite audit of Báez-Duarte's greedy Vasyunin correction.

For positive integers ``n`` and ``m`` define

``h_n(m) = floor(m / n) - 2 floor(m / (2 n))``.

Starting with no coefficients, the greedy correction chooses

``c_n = 1 - sum_(k<n) c_k h_k(n)``.

Báez-Duarte proved the closed form

``c_(2^r m) = 2^max(r-1, 0) mu(m)`` for odd ``m``

and used it to prove divergence of the infinite correction.  This module
checks the recurrence, closed form, and interval interpolation property over
one configurable *finite* prefix using only exact integer and rational
arithmetic.  The theorem in the cited 2005 preprint is recorded as an external
symbolic bridge; it is not claimed as a consequence of a finite Python run.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


VASYUNIN_AUDIT_SCHEMA = "rh-lab/nyman-vasyunin-greedy-prefix-audit/v1"
DEFAULT_AUDIT_LIMIT = 64
MAX_AUDIT_LIMIT = 4096

VASYUNIN_LIMITATION = (
    "This CERTIFIED_FINITE audit checks an exact finite prefix only. The "
    "identity involving log(2) and the divergence of the infinite greedy "
    "correction are recorded from the cited-preprint theorem bridge and are "
    "not proved by finite enumeration. The result does not cover every "
    "possible Vasyunin correction and does not prove or disprove the Riemann "
    "Hypothesis; the global status remains UNRESOLVED."
)

_AUDIT_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "audit_outcome",
    "finite_scope",
    "definitions",
    "coefficient_records",
    "power_of_two_increment_records",
    "checks",
    "external_theorem_bridge",
    "limitation",
    "payload_sha256",
}


class NymanVasyuninError(ValueError):
    """Raised when an exact finite greedy-correction audit cannot close."""


class NymanVasyuninVerificationError(NymanVasyuninError):
    """Raised when a supplied finite audit fails exact regeneration."""


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
    return value


def _audit_limit(value: Any) -> int:
    limit = _positive_integer(value, "limit")
    if limit > MAX_AUDIT_LIMIT:
        raise ValueError(f"limit must be at most {MAX_AUDIT_LIMIT}")
    return limit


def _two_adic_parts(n: int) -> tuple[int, int]:
    exponent = 0
    odd_part = n
    while odd_part % 2 == 0:
        exponent += 1
        odd_part //= 2
    return exponent, odd_part


def mobius(n: int) -> int:
    """Return the exact Möbius value ``mu(n)`` by integer factorization."""

    value = _positive_integer(n, "n")
    remaining = value
    distinct_prime_factors = 0
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime = 3 if prime == 2 else prime + 2
            continue
        remaining //= prime
        distinct_prime_factors += 1
        if remaining % prime == 0:
            return 0
        prime = 3 if prime == 2 else prime + 2
    if remaining > 1:
        distinct_prime_factors += 1
    return -1 if distinct_prime_factors % 2 else 1


def vasyunin_seed(n: int, m: int) -> int:
    """Return ``h_n(m)`` exactly for positive integer arguments."""

    divisor = _positive_integer(n, "n")
    argument = _positive_integer(m, "m")
    return argument // divisor - 2 * (argument // (2 * divisor))


def vasyunin_formula_coefficient(n: int) -> int:
    """Return ``2^max(r-1,0) mu(m)`` for ``n=2^r m`` and odd ``m``."""

    value = _positive_integer(n, "n")
    exponent, odd_part = _two_adic_parts(value)
    return (1 << max(exponent - 1, 0)) * mobius(odd_part)


def _recurrence_coefficient(
    n: int,
    prior_coefficients: Sequence[int],
) -> int:
    if len(prior_coefficients) != n - 1:
        raise NymanVasyuninError(
            "recurrence requires exactly the coefficients before n"
        )
    return 1 - sum(
        coefficient * vasyunin_seed(k, n)
        for k, coefficient in enumerate(prior_coefficients, start=1)
    )


def vasyunin_greedy_coefficients(limit: int) -> tuple[int, ...]:
    """Return ``(c_1, ..., c_limit)`` from the exact greedy recurrence."""

    checked_limit = _audit_limit(limit)
    coefficients: list[int] = []
    for n in range(1, checked_limit + 1):
        coefficients.append(_recurrence_coefficient(n, coefficients))
    return tuple(coefficients)


def _coefficient_snapshot(coefficients: Sequence[int]) -> tuple[int, ...]:
    if isinstance(coefficients, (str, bytes, bytearray)) or not isinstance(
        coefficients,
        Sequence,
    ):
        raise TypeError("coefficients must be a nonempty integer sequence")
    snapshot = tuple(coefficients)
    if not snapshot:
        raise ValueError("coefficients must be nonempty")
    for index, coefficient in enumerate(snapshot, start=1):
        if isinstance(coefficient, bool) or not isinstance(coefficient, int):
            raise TypeError(f"coefficient {index} must be an integer")
    return snapshot


def vasyunin_phi(coefficients: Sequence[int], m: int) -> int:
    """Evaluate ``sum_k c_k h_k(m)`` with exact integer arithmetic."""

    snapshot = _coefficient_snapshot(coefficients)
    argument = _positive_integer(m, "m")
    return sum(
        coefficient * vasyunin_seed(k, argument)
        for k, coefficient in enumerate(snapshot, start=1)
    )


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _content_sha256(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _with_payload_hash(payload: Mapping[str, Any]) -> dict[str, Any]:
    body = dict(payload)
    body.pop("payload_sha256", None)
    return {**body, "payload_sha256": _content_sha256(body)}


def _external_theorem_bridge() -> dict[str, Any]:
    return {
        "citation": {
            "author": "Luis Baez-Duarte",
            "title": "A divergent Vasyunin correction",
            "year": "2005",
            "primary_url": "https://arxiv.org/abs/math/0506318",
        },
        "cited_preprint_statements": {
            "coefficient_closed_form": (
                "for n=2^r*m with m odd, "
                "c_n=2^max(r-1,0)*mu(m)"
            ),
            "greedy_sequence_result": (
                "the infinite greedy correction diverges in "
                "L1((0,infinity),x^-2 dx), hence cannot converge in the "
                "required L2 space"
            ),
        },
        "symbolic_increment_bridge": {
            "finite_difference": "phi_n-phi_(n-1)=c_n*h_n",
            "seed_weighted_integral": (
                "integral_0^infinity h_n(x)*x^-2 dx=log(2)/n"
            ),
            "l1_increment": (
                "norm_1(phi_n-phi_(n-1))=(|c_n|/n)*log(2)"
            ),
            "l2_squared_increment": (
                "norm_2(phi_n-phi_(n-1))^2=(c_n^2/n)*log(2)"
            ),
            "nontrivial_power_of_two_coefficients": (
                "for n=2^r and r>=1, c_n=n/2"
            ),
            "nontrivial_power_of_two_l1_log2_multiplier": {
                "symbolic": "1/2",
                "exact_fraction": {"numerator": "1", "denominator": "2"},
            },
            "nontrivial_power_of_two_l2_squared_log2_multiplier": {
                "symbolic": "n/4",
            },
        },
        "machine_scope": {
            "finite_prefix_recurrence_and_closed_form_checked": True,
            "finite_prefix_interval_interpolation_checked": True,
            "finite_prefix_rational_increment_multipliers_checked": True,
            "log2_integral_identity_reproved_by_finite_enumeration": False,
            "infinite_closed_form_reproved_by_finite_enumeration": False,
            "infinite_l1_divergence_reproved_by_finite_enumeration": False,
            "all_vasyunin_corrections_covered": False,
            "resolves_rh": False,
        },
    }


def audit_vasyunin_greedy(
    limit: int = DEFAULT_AUDIT_LIMIT,
) -> dict[str, Any]:
    """Certify the recurrence and interpolation identities through ``limit``.

    The output contains only canonical JSON values.  ``CERTIFIED_FINITE``
    refers solely to the enumerated prefix; the infinite theorem remains an
    explicitly external dependency from the cited preprint.
    """

    checked_limit = _audit_limit(limit)
    coefficients = vasyunin_greedy_coefficients(checked_limit)
    records: list[dict[str, Any]] = []
    power_records: list[dict[str, Any]] = []

    for n, recurrence_value in enumerate(coefficients, start=1):
        exponent, odd_part = _two_adic_parts(n)
        mu_odd = mobius(odd_part)
        formula_value = vasyunin_formula_coefficient(n)
        if recurrence_value != formula_value:
            raise NymanVasyuninError(
                f"recurrence and closed form disagree at n={n}"
            )

        checked_seed_values = [
            vasyunin_seed(k, n) for k in range(1, n + 1)
        ]
        if any(value not in {0, 1} for value in checked_seed_values):
            raise NymanVasyuninError(f"seed is not binary at argument n={n}")
        if vasyunin_seed(n, n) != 1 or any(
            vasyunin_seed(n, prior_m) != 0
            for prior_m in range(1, n)
        ):
            raise NymanVasyuninError(
                f"new seed does not preserve prior intervals at n={n}"
            )

        phi_on_new_interval = sum(
            coefficients[k - 1] * checked_seed_values[k - 1]
            for k in range(1, n + 1)
        )
        if phi_on_new_interval != 1:
            raise NymanVasyuninError(
                f"greedy interpolation failed on interval n={n}"
            )

        l1_multiplier = Fraction(abs(recurrence_value), n)
        l2_squared_multiplier = Fraction(recurrence_value**2, n)
        records.append(
            {
                "n": str(n),
                "two_adic_exponent_r": str(exponent),
                "odd_part_m": str(odd_part),
                "mobius_of_odd_part": str(mu_odd),
                "recurrence_coefficient": str(recurrence_value),
                "closed_form_coefficient": str(formula_value),
                "phi_n_on_new_integer_interval": str(phi_on_new_interval),
                "l1_increment_log2_multiplier": _fraction_record(
                    l1_multiplier
                ),
                "l2_squared_increment_log2_multiplier": _fraction_record(
                    l2_squared_multiplier
                ),
            }
        )

        if n >= 2 and n & (n - 1) == 0:
            if recurrence_value != n // 2:
                raise NymanVasyuninError(
                    f"power-of-two coefficient specialization failed at n={n}"
                )
            if l1_multiplier != Fraction(1, 2):
                raise NymanVasyuninError(
                    f"power-of-two L1 multiplier failed at n={n}"
                )
            if l2_squared_multiplier != Fraction(n, 4):
                raise NymanVasyuninError(
                    f"power-of-two L2 multiplier failed at n={n}"
                )
            power_records.append(
                {
                    "n": str(n),
                    "two_adic_exponent_r": str(exponent),
                    "coefficient_c_n": str(recurrence_value),
                    "l1_increment": {
                        "expression": "(1/2)*log(2)",
                        "symbolic_rational_multiplier": "1/2",
                        "exact_rational_multiplier": _fraction_record(
                            l1_multiplier
                        ),
                    },
                    "l2_squared_increment": {
                        "expression": f"({n}/4)*log(2)",
                        "symbolic_rational_multiplier": "n/4",
                        "substituted_rational_multiplier": f"{n}/4",
                        "exact_reduced_rational_multiplier": _fraction_record(
                            l2_squared_multiplier
                        ),
                    },
                }
            )

    final_phi_values = [
        vasyunin_phi(coefficients, m)
        for m in range(1, checked_limit + 1)
    ]
    if final_phi_values != [1] * checked_limit:
        raise NymanVasyuninError(
            "final greedy correction is not one on every checked interval"
        )

    payload = {
        "schema": VASYUNIN_AUDIT_SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "GREEDY_VASYUNIN_FINITE_PREFIX_IDENTITIES_CERTIFIED",
        "finite_scope": {
            "limit": str(checked_limit),
            "coefficient_count": str(len(coefficients)),
            "integer_intervals": f"[1,{checked_limit + 1})",
            "interval_convention": (
                "h_n is constant on each [m,m+1), so its value there is "
                "the exact integer seed h_n(m)"
            ),
        },
        "definitions": {
            "seed": "h_n(m)=floor(m/n)-2*floor(m/(2*n))",
            "recurrence": "c_n=1-sum_(k<n)c_k*h_k(n)",
            "closed_form": (
                "c_(2^r*m)=2^max(r-1,0)*mu(m), with m odd"
            ),
            "greedy_correction": "phi_n=sum_(k<=n)c_k*h_k",
        },
        "coefficient_records": records,
        "power_of_two_increment_records": power_records,
        "checks": {
            "recurrence_matches_closed_form_through_limit": True,
            "all_checked_seed_values_are_binary": True,
            "each_new_seed_vanishes_on_all_prior_integer_intervals": True,
            "inductive_phi_n_equals_one_on_intervals_1_through_n": True,
            "final_phi_limit_equals_one_on_every_checked_interval": True,
            "power_of_two_l1_multiplier_is_one_half": True,
            "power_of_two_l2_squared_multiplier_is_n_over_four": True,
            "integer_and_fraction_arithmetic_only": True,
        },
        "external_theorem_bridge": _external_theorem_bridge(),
        "limitation": VASYUNIN_LIMITATION,
    }
    return _with_payload_hash(payload)


def _snapshot_json_tree(value: Any, *, location: str) -> Any:
    if value is None or type(value) in {bool, int, str}:
        return value
    if type(value) is float:
        raise NymanVasyuninVerificationError(
            f"floating-point value is noncanonical at {location}"
        )
    if isinstance(value, list):
        return [
            _snapshot_json_tree(item, location=f"{location}[{index}]")
            for index, item in enumerate(list(value))
        ]
    if isinstance(value, Mapping):
        snapshot: dict[str, Any] = {}
        for key, item in list(value.items()):
            if type(key) is not str:
                raise NymanVasyuninVerificationError(
                    f"non-string key at {location}"
                )
            if key in snapshot:
                raise NymanVasyuninVerificationError(
                    f"duplicate JSON object key: {key}"
                )
            snapshot[key] = _snapshot_json_tree(
                item,
                location=f"{location}.{key}",
            )
        return snapshot
    raise NymanVasyuninVerificationError(
        f"non-JSON value at {location}"
    )


def _load_audit_object(
    source: Mapping[str, Any] | Path,
) -> dict[str, Any]:
    if isinstance(source, Path):

        def reject_duplicate_keys(
            pairs: list[tuple[str, Any]],
        ) -> dict[str, Any]:
            value: dict[str, Any] = {}
            for key, item in pairs:
                if key in value:
                    raise NymanVasyuninVerificationError(
                        f"duplicate JSON object key: {key}"
                    )
                value[key] = item
            return value

        def reject_nonstandard_constant(value: str) -> Any:
            raise NymanVasyuninVerificationError(
                f"nonstandard JSON constant: {value}"
            )

        try:
            loaded = json.loads(
                source.read_text(encoding="utf-8"),
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=reject_nonstandard_constant,
            )
        except NymanVasyuninVerificationError:
            raise
        except (OSError, RecursionError, ValueError) as exc:
            raise NymanVasyuninVerificationError(
                f"cannot read or parse Vasyunin audit: {source}"
            ) from exc
    elif isinstance(source, Mapping):
        loaded = source
    else:
        raise TypeError("artifact must be a mapping or Path")

    try:
        snapshot = _snapshot_json_tree(loaded, location="artifact")
    except NymanVasyuninVerificationError:
        raise
    except (OverflowError, RecursionError, TypeError, ValueError) as exc:
        raise NymanVasyuninVerificationError(
            "artifact is not a finite canonical JSON tree"
        ) from exc
    if not isinstance(snapshot, dict):
        raise NymanVasyuninVerificationError("artifact must be an object")
    return snapshot


def _canonical_limit_string(value: Any) -> int:
    if not isinstance(value, str):
        raise NymanVasyuninVerificationError(
            "finite scope limit must be a canonical integer string"
        )
    try:
        parsed = int(value)
    except ValueError as exc:
        raise NymanVasyuninVerificationError(
            "finite scope limit is not an integer"
        ) from exc
    if value != str(parsed):
        raise NymanVasyuninVerificationError(
            "finite scope limit is not canonical"
        )
    try:
        return _audit_limit(parsed)
    except (TypeError, ValueError) as exc:
        raise NymanVasyuninVerificationError(
            "finite scope limit is outside the audit range"
        ) from exc


def verify_vasyunin_greedy_audit(
    artifact: Mapping[str, Any] | Path,
) -> dict[str, Any]:
    """Strictly load and regenerate a finite audit, rejecting mutations."""

    supplied = _load_audit_object(artifact)
    if set(supplied) != _AUDIT_FIELDS:
        raise NymanVasyuninVerificationError("audit fields changed")
    if supplied.get("schema") != VASYUNIN_AUDIT_SCHEMA:
        raise NymanVasyuninVerificationError("unexpected audit schema")
    if supplied.get("classification") != "CERTIFIED_FINITE":
        raise NymanVasyuninVerificationError(
            "finite audit was improperly reclassified"
        )
    if supplied.get("hypothesis_status") != "UNRESOLVED":
        raise NymanVasyuninVerificationError(
            "finite audit changed the hypothesis status"
        )
    if supplied.get("limitation") != VASYUNIN_LIMITATION:
        raise NymanVasyuninVerificationError("audit limitation changed")

    supplied_hash = supplied.get("payload_sha256")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if supplied_hash != _content_sha256(body):
        raise NymanVasyuninVerificationError("audit payload hash mismatch")

    finite_scope = supplied.get("finite_scope")
    if not isinstance(finite_scope, Mapping):
        raise NymanVasyuninVerificationError("finite scope is missing")
    limit = _canonical_limit_string(finite_scope.get("limit"))
    expected = audit_vasyunin_greedy(limit)
    if supplied != expected:
        raise NymanVasyuninVerificationError(
            "audit does not exactly regenerate"
        )
    return {
        "classification": (
            "REPRODUCED_CERTIFIED_FINITE_VASYUNIN_GREEDY_PREFIX"
        ),
        "hypothesis_status": "UNRESOLVED",
        "verified_limit": str(limit),
        "payload_sha256": supplied["payload_sha256"],
        "infinite_divergence_machine_reproved": False,
    }


__all__ = [
    "DEFAULT_AUDIT_LIMIT",
    "MAX_AUDIT_LIMIT",
    "NymanVasyuninError",
    "NymanVasyuninVerificationError",
    "VASYUNIN_AUDIT_SCHEMA",
    "VASYUNIN_LIMITATION",
    "audit_vasyunin_greedy",
    "mobius",
    "vasyunin_formula_coefficient",
    "vasyunin_greedy_coefficients",
    "vasyunin_phi",
    "vasyunin_seed",
    "verify_vasyunin_greedy_audit",
]
