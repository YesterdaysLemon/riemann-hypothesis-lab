"""Exact finite evidence for the sharp alias--Mobius truncation identity.

For a finitely supported arithmetic function ``y`` put ``a = mu * y`` and

``S_X = sum_(n<=X) a_n / n``.

The associated sharp Nyman section is ``g_X(t) = sum_(n<=X) a_n {t/n}``.
For ``t <= X`` Dirichlet inversion gives the exact identity

``g_X(t) - f_y(t) = t S_X``, where
``f_y(t) = -sum_(j<=floor(t)) y_j``.

This module certifies those *finite* identities for a frozen compact shell.
The analytic theorem that no nonzero finite shell can have
``S_X = o(X**(-1/2))`` is documented as a human-auditable theorem bridge; a
finite Python run is not presented as a proof of that infinite statement.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

from .artifacts import content_sha256


ALIAS_AUDIT_SCHEMA = "rh-lab/nyman-alias-sharp-truncation-audit/v1"
DEFAULT_AUDIT_LIMIT = 4_096
MAX_AUDIT_LIMIT = 4_096
FROZEN_SHELL = ((9, 1), (16, -1))

ALIAS_LIMITATION = (
    "This CERTIFIED_FINITE artifact checks exact arithmetic identities only "
    "through its declared limit for the frozen shell y_9=1, y_16=-1. The "
    "all-shell sharp-truncation obstruction uses the cited critical-line "
    "zero theorem plus the human-auditable Mellin and Abelian argument; it "
    "is not proved by finite enumeration. The obstruction does not cover "
    "scale-dependent shells, non-sharp summability methods, or arbitrary "
    "multiscale constructions, and it does not prove or disprove the "
    "Riemann Hypothesis. The global status remains UNRESOLVED."
)

_AUDIT_FIELDS = {
    "schema",
    "classification",
    "hypothesis_status",
    "audit_outcome",
    "finite_scope",
    "frozen_shell",
    "definitions",
    "vector_bindings",
    "checkpoint_records",
    "checks",
    "external_theorem_bridge",
    "limitation",
    "payload_sha256",
}


class NymanAliasError(ValueError):
    """Raised when an exact alias audit cannot close."""


class NymanAliasVerificationError(NymanAliasError):
    """Raised when a supplied alias audit fails exact regeneration."""


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
    return value


def _audit_limit(value: Any) -> int:
    limit = _positive_integer(value, "limit")
    if limit < max(index for index, _ in FROZEN_SHELL):
        raise ValueError("limit must include the complete frozen shell")
    if limit > MAX_AUDIT_LIMIT:
        raise ValueError(f"limit must be at most {MAX_AUDIT_LIMIT}")
    return limit


def _shell_snapshot(shell: Mapping[int, int]) -> tuple[tuple[int, int], ...]:
    if not isinstance(shell, Mapping):
        raise TypeError("shell must be a mapping")
    snapshot: list[tuple[int, int]] = []
    for index, coefficient in list(shell.items()):
        checked_index = _positive_integer(index, "shell index")
        if isinstance(coefficient, bool) or not isinstance(coefficient, int):
            raise TypeError("shell coefficients must be integers")
        if coefficient:
            snapshot.append((checked_index, coefficient))
    snapshot.sort()
    if not snapshot:
        raise ValueError("shell must be nonzero")
    if len({index for index, _ in snapshot}) != len(snapshot):
        raise ValueError("shell indices must be unique")
    return tuple(snapshot)


def mobius_sieve(limit: int) -> tuple[int, ...]:
    """Return ``mu(0), ..., mu(limit)`` by an exact linear sieve."""

    checked_limit = _positive_integer(limit, "limit")
    mu = [0] * (checked_limit + 1)
    mu[1] = 1
    primes: list[int] = []
    composite = [False] * (checked_limit + 1)
    for value in range(2, checked_limit + 1):
        if not composite[value]:
            primes.append(value)
            mu[value] = -1
        for prime in primes:
            product = value * prime
            if product > checked_limit:
                break
            composite[product] = True
            if value % prime == 0:
                mu[product] = 0
                break
            mu[product] = -mu[value]
    return tuple(mu)


def alias_coefficients(
    shell: Mapping[int, int],
    limit: int,
) -> tuple[int, ...]:
    """Return the exact prefix of ``a = mu * y`` including index zero."""

    checked_limit = _positive_integer(limit, "limit")
    support = _shell_snapshot(shell)
    mu = mobius_sieve(checked_limit)
    coefficients = [0] * (checked_limit + 1)
    for index, shell_coefficient in support:
        if index > checked_limit:
            continue
        for multiplier in range(1, checked_limit // index + 1):
            coefficients[index * multiplier] += (
                shell_coefficient * mu[multiplier]
            )
    return tuple(coefficients)


def divisor_sum_prefix(coefficients: Sequence[int]) -> tuple[int, ...]:
    """Return ``(1*a)(n)`` for a finite coefficient prefix."""

    if isinstance(coefficients, (str, bytes, bytearray)) or not isinstance(
        coefficients, Sequence
    ):
        raise TypeError("coefficients must be an integer sequence")
    snapshot = tuple(coefficients)
    if len(snapshot) < 2:
        raise ValueError("coefficients must include indices zero and one")
    if snapshot[0] != 0:
        raise ValueError("coefficient index zero must be zero")
    for value in snapshot:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError("coefficients must be integers")
    limit = len(snapshot) - 1
    result = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        coefficient = snapshot[divisor]
        if not coefficient:
            continue
        for multiple in range(divisor, limit + 1, divisor):
            result[multiple] += coefficient
    return tuple(result)


def _fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _checkpoints(limit: int) -> tuple[int, ...]:
    values: list[int] = []
    checkpoint = 16
    while checkpoint <= limit:
        values.append(checkpoint)
        checkpoint *= 2
    if not values or values[-1] != limit:
        values.append(limit)
    return tuple(values)


def _external_theorem_bridge() -> dict[str, Any]:
    return {
        "theorem": (
            "For every nonzero finitely supported y, if a=mu*y and "
            "S_X=sum_(n<=X)a_n/n, then S_X is not o(X^(-1/2)). "
            "Consequently the sharp sections g_X=sum_(n<=X)a_n*{t/n} "
            "do not converge to f_y in L2((0,infinity),dt/t^2)."
        ),
        "proof_outline": [
            (
                "For Re(z)>0, D(z)=sum_n a_n*n^(-1-z)="
                "B_y(1+z)/zeta(1+z), where B_y(s)=sum_j y_j*j^(-s)."
            ),
            (
                "If S_X=o(X^(-1/2)), partial summation extends D "
                "holomorphically to Re(z)>-1/2 and a boundary Abelian "
                "lemma forces epsilon*D(-1/2+epsilon+i*gamma) to zero."
            ),
            (
                "A nonzero finite Dirichlet polynomial is an entire "
                "exponential polynomial with only O(T) zeros, counted "
                "with multiplicity, in a bounded vertical strip."
            ),
            (
                "Conrey proved that a positive proportion, hence order "
                "T*log(T), of zeta zeros are simple and on the critical "
                "line. One is not a zero of B_y, giving a boundary pole "
                "that contradicts the Abelian limit."
            ),
            (
                "The exact identity ||g_X-f_y||_2^2 >= X*|S_X|^2 then "
                "excludes sharp-section norm convergence."
            ),
        ],
        "primary_sources": [
            {
                "author": "Luis Baez-Duarte",
                "title": "Moebius-convolutions and the Riemann hypothesis",
                "year": "2005",
                "primary_url": "https://arxiv.org/abs/math/0504402",
                "relevance": (
                    "Mellin-convolution RH criteria and the analogous "
                    "critical x^(-1/2) non-little-o obstruction"
                ),
            },
            {
                "author": "J. B. Conrey",
                "title": (
                    "More than two fifths of the zeros of the Riemann "
                    "zeta function are on the critical line"
                ),
                "year": "1989",
                "primary_url": "https://doi.org/10.1515/crll.1989.399.1",
                "relevance": (
                    "positive-proportion theorem for simple critical-line "
                    "zeros"
                ),
            },
        ],
        "machine_scope": {
            "finite_mobius_convolution_checked": True,
            "finite_divisor_inversion_checked": True,
            "finite_slope_factorization_checked": True,
            "mellin_continuation_machine_proved": False,
            "exponential_polynomial_zero_count_machine_proved": False,
            "conrey_zero_theorem_machine_reproved": False,
            "infinite_non_little_o_conclusion_machine_proved": False,
            "scale_dependent_regularizations_covered": False,
            "resolves_rh": False,
        },
    }


def audit_alias_sharp_truncation(
    limit: int = DEFAULT_AUDIT_LIMIT,
) -> dict[str, Any]:
    """Generate the frozen exact finite alias audit through ``limit``."""

    checked_limit = _audit_limit(limit)
    shell = dict(FROZEN_SHELL)
    mu = mobius_sieve(checked_limit)
    coefficients = alias_coefficients(shell, checked_limit)
    inverted = divisor_sum_prefix(coefficients)

    expected_y = [0] * (checked_limit + 1)
    for index, value in FROZEN_SHELL:
        expected_y[index] = value
    if list(inverted) != expected_y:
        raise NymanAliasError("Dirichlet inversion failed in the finite prefix")

    direct_slopes = [Fraction(0)] * (checked_limit + 1)
    mobius_harmonic = [Fraction(0)] * (checked_limit + 1)
    for index in range(1, checked_limit + 1):
        direct_slopes[index] = direct_slopes[index - 1] + Fraction(
            coefficients[index], index
        )
        mobius_harmonic[index] = mobius_harmonic[index - 1] + Fraction(
            mu[index], index
        )

    for x in range(1, checked_limit + 1):
        factored = sum(
            Fraction(value, index) * mobius_harmonic[x // index]
            for index, value in FROZEN_SHELL
        )
        if direct_slopes[x] != factored:
            raise NymanAliasError(f"slope factorization failed at X={x}")

    cumulative_inverted = 0
    cumulative_y = 0
    for argument in range(1, checked_limit + 1):
        cumulative_inverted += inverted[argument]
        cumulative_y += expected_y[argument]
        if cumulative_inverted != cumulative_y:
            raise NymanAliasError(
                f"integer-interval alias identity failed at m={argument}"
            )

    checkpoint_records: list[dict[str, Any]] = []
    nonzero_count = 0
    next_checkpoint_index = 0
    checkpoints = _checkpoints(checked_limit)
    for index in range(1, checked_limit + 1):
        nonzero_count += int(coefficients[index] != 0)
        if (
            next_checkpoint_index < len(checkpoints)
            and index == checkpoints[next_checkpoint_index]
        ):
            slope = direct_slopes[index]
            checkpoint_records.append(
                {
                    "X": str(index),
                    "nonzero_alias_coefficients_through_X": str(nonzero_count),
                    "S_X": _fraction_record(slope),
                    "X_times_S_X_squared": _fraction_record(
                        index * slope * slope
                    ),
                }
            )
            next_checkpoint_index += 1

    shell_sum = sum(value for _, value in FROZEN_SHELL)
    shell_harmonic_sum = sum(
        Fraction(value, index) for index, value in FROZEN_SHELL
    )
    payload = {
        "schema": ALIAS_AUDIT_SCHEMA,
        "classification": "CERTIFIED_FINITE",
        "hypothesis_status": "UNRESOLVED",
        "audit_outcome": "SHARP_ALIAS_FINITE_IDENTITIES_CERTIFIED",
        "finite_scope": {
            "limit": str(checked_limit),
            "checked_integer_arguments": f"[1,{checked_limit}]",
            "checkpoint_count": str(len(checkpoint_records)),
            "arithmetic": "exact CPython integers and Fraction",
        },
        "frozen_shell": {
            "records": [
                {"index": str(index), "coefficient": str(value)}
                for index, value in FROZEN_SHELL
            ],
            "support_interval": "(8,16]",
            "coefficient_sum": str(shell_sum),
            "harmonic_sum": _fraction_record(shell_harmonic_sum),
            "target_step": "f_y(t)=-1 on [9,16) and 0 elsewhere",
        },
        "definitions": {
            "extension": "a_n=sum_(j|n)y_j*mu(n/j)",
            "slope": "S_X=sum_(n<=X)a_n/n",
            "factored_slope": (
                "S_X=sum_j (y_j/j)*sum_(k<=X/j)mu(k)/k"
            ),
            "target": "f_y(t)=-sum_(j<=floor(t))y_j",
            "sharp_section": "g_X(t)=sum_(n<=X)a_n*{t/n}",
            "finite_error_identity": "g_X(t)-f_y(t)=t*S_X for 0<t<=X",
            "norm_lower_bound": "||g_X-f_y||_2^2>=X*|S_X|^2",
        },
        "vector_bindings": {
            "mobius_prefix_sha256": content_sha256(
                [str(value) for value in mu[1:]]
            ),
            "alias_coefficient_prefix_sha256": content_sha256(
                [str(value) for value in coefficients[1:]]
            ),
            "divisor_inversion_prefix_sha256": content_sha256(
                [str(value) for value in inverted[1:]]
            ),
        },
        "checkpoint_records": checkpoint_records,
        "checks": {
            "mobius_sieve_completed": True,
            "a_equals_mu_convolved_with_y_through_limit": True,
            "one_convolved_with_a_equals_y_through_limit": True,
            "direct_and_factored_S_X_match_for_every_X": True,
            "cumulative_floor_alias_matches_target_on_every_integer_interval": True,
            "frozen_shell_has_zero_coefficient_sum": shell_sum == 0,
            "integer_and_fraction_arithmetic_only": True,
        },
        "external_theorem_bridge": _external_theorem_bridge(),
        "limitation": ALIAS_LIMITATION,
    }
    return {**payload, "payload_sha256": content_sha256(payload)}


def _snapshot_json_tree(value: Any, *, location: str) -> Any:
    if value is None or type(value) in {bool, int, str}:
        return value
    if type(value) is float:
        if not math.isfinite(value):
            raise NymanAliasVerificationError(
                f"non-finite JSON number at {location}"
            )
        raise NymanAliasVerificationError(
            f"floating-point JSON number is noncanonical at {location}"
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
                raise NymanAliasVerificationError(
                    f"non-string JSON object key at {location}"
                )
            if key in snapshot:
                raise NymanAliasVerificationError(
                    f"duplicate JSON object key: {key}"
                )
            snapshot[key] = _snapshot_json_tree(
                item, location=f"{location}.{key}"
            )
        return snapshot
    raise NymanAliasVerificationError(f"non-JSON value at {location}")


def _load_audit(source: Mapping[str, Any] | Path) -> dict[str, Any]:
    if isinstance(source, Path):

        def reject_duplicate_keys(
            pairs: list[tuple[str, Any]],
        ) -> dict[str, Any]:
            value: dict[str, Any] = {}
            for key, item in pairs:
                if key in value:
                    raise NymanAliasVerificationError(
                        f"duplicate JSON object key: {key}"
                    )
                value[key] = item
            return value

        def reject_nonstandard_constant(value: str) -> Any:
            raise NymanAliasVerificationError(
                f"nonstandard JSON constant: {value}"
            )

        try:
            loaded = json.loads(
                source.read_text(encoding="utf-8"),
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=reject_nonstandard_constant,
            )
        except NymanAliasVerificationError:
            raise
        except (OSError, json.JSONDecodeError, RecursionError) as exc:
            raise NymanAliasVerificationError(
                f"cannot read or parse alias audit: {source}"
            ) from exc
    elif isinstance(source, Mapping):
        loaded = source
    else:
        raise TypeError("alias audit must be a mapping or Path")

    if not isinstance(loaded, Mapping):
        raise NymanAliasVerificationError("alias audit must be a JSON object")
    try:
        snapshot = _snapshot_json_tree(loaded, location="alias audit")
    except RecursionError as exc:
        raise NymanAliasVerificationError(
            "alias audit nesting is too deep"
        ) from exc
    if not isinstance(snapshot, dict):
        raise NymanAliasVerificationError("alias audit must be a JSON object")
    return snapshot


def _canonical_positive_integer_text(value: Any, name: str) -> int:
    if not isinstance(value, str) or not value or not value.isascii():
        raise NymanAliasVerificationError(f"{name} is not canonical")
    if not value.isdigit() or (len(value) > 1 and value.startswith("0")):
        raise NymanAliasVerificationError(f"{name} is not canonical")
    parsed = int(value)
    if parsed < 1:
        raise NymanAliasVerificationError(f"{name} must be positive")
    return parsed


def verify_alias_sharp_truncation_audit(
    source: Mapping[str, Any] | Path,
) -> dict[str, Any]:
    """Strictly validate and exactly regenerate a finite alias audit."""

    supplied = _load_audit(source)
    if set(supplied) != _AUDIT_FIELDS:
        raise NymanAliasVerificationError("alias audit schema fields changed")
    if supplied.get("schema") != ALIAS_AUDIT_SCHEMA:
        raise NymanAliasVerificationError("alias audit schema changed")
    if supplied.get("classification") != "CERTIFIED_FINITE":
        raise NymanAliasVerificationError("alias audit classification changed")
    if supplied.get("hypothesis_status") != "UNRESOLVED":
        raise NymanAliasVerificationError("alias audit RH status changed")

    expected_hash = supplied.get("payload_sha256")
    if (
        not isinstance(expected_hash, str)
        or len(expected_hash) != 64
        or any(character not in "0123456789abcdef" for character in expected_hash)
    ):
        raise NymanAliasVerificationError("payload hash is invalid")
    body = {
        key: value
        for key, value in supplied.items()
        if key != "payload_sha256"
    }
    if content_sha256(body) != expected_hash:
        raise NymanAliasVerificationError("payload hash mismatch")

    finite_scope = supplied.get("finite_scope")
    if not isinstance(finite_scope, dict):
        raise NymanAliasVerificationError("finite scope is missing")
    limit = _canonical_positive_integer_text(
        finite_scope.get("limit"), "finite limit"
    )
    try:
        regenerated = audit_alias_sharp_truncation(limit)
    except (NymanAliasError, TypeError, ValueError) as exc:
        raise NymanAliasVerificationError(
            "alias audit cannot be regenerated"
        ) from exc
    if supplied != regenerated:
        raise NymanAliasVerificationError(
            "alias audit does not exactly regenerate"
        )
    return {
        "classification": "REPRODUCED_CERTIFIED_FINITE_ALIAS_IDENTITIES",
        "hypothesis_status": "UNRESOLVED",
        "verified_limit": str(limit),
        "payload_sha256": expected_hash,
        "infinite_sharp_truncation_obstruction_machine_reproved": False,
    }


__all__ = [
    "ALIAS_AUDIT_SCHEMA",
    "ALIAS_LIMITATION",
    "DEFAULT_AUDIT_LIMIT",
    "FROZEN_SHELL",
    "MAX_AUDIT_LIMIT",
    "NymanAliasError",
    "NymanAliasVerificationError",
    "alias_coefficients",
    "audit_alias_sharp_truncation",
    "divisor_sum_prefix",
    "mobius_sieve",
    "verify_alias_sharp_truncation_audit",
]
