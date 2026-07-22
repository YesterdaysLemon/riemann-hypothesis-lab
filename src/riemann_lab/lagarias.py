"""Rigorous finite checks of the Lagarias criterion."""

from __future__ import annotations

from typing import Any

import flint
from flint import arb, ctx

from .arithmetic import divisor_sums_up_to
from .artifacts import content_sha256
from .balls import arb_record


LAGARIAS_SOURCE = "https://arxiv.org/abs/math/0008177"
LAGARIAS_CLAIM = (
    "The Lagarias inequality holds for every integer in the stated finite range, "
    "with equality at n=1 only."
)
LAGARIAS_LIMITATION = (
    "This checks finitely many integers. Lagarias's equivalence has a universal "
    "quantifier, so this result does not prove RH."
)
LAGARIAS_ARITHMETIC = "exact integer sigma; Arb ball transcendental operations"


class LagariasCertificateError(ValueError):
    """Raised when a serialized Lagarias certificate fails replay."""


def certify_lagarias_range(limit: int, precision_bits: int = 192) -> dict[str, Any]:
    """Certify the Lagarias inequality for every integer in ``[1, limit]``.

    Lagarias proved that RH is equivalent to

        sigma(n) <= H_n + exp(H_n) log(H_n)

    for all positive integers, with equality only at one.  This routine uses
    exact integer divisor sums and Arb balls for every transcendental
    operation.  It certifies only the requested finite range.
    """

    if isinstance(limit, bool) or not isinstance(limit, int):
        raise TypeError("limit must be an integer")
    if limit < 1:
        raise ValueError("limit must be at least 1")
    if precision_bits < 64:
        raise ValueError("precision_bits must be at least 64")

    sigma = divisor_sums_up_to(limit)
    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        harmonic = arb(0)
        inconclusive: list[int] = []
        violations: list[int] = []
        equality_at_one = False
        tightest_n: int | None = None
        tightest_margin: Any | None = None

        for n in range(1, limit + 1):
            harmonic += arb(1) / n
            rhs = harmonic + harmonic.exp() * harmonic.log()
            margin = rhs - sigma[n]

            if n == 1:
                equality_at_one = margin.is_zero()
                if not equality_at_one:
                    inconclusive.append(n)
                continue

            if margin > 0:
                if tightest_margin is None or margin.lower() < tightest_margin.lower():
                    tightest_n = n
                    tightest_margin = margin
            elif margin < 0:
                violations.append(n)
            else:
                inconclusive.append(n)

        finite_range_certified = (
            equality_at_one and not inconclusive and not violations
        )
        payload = {
            "schema": "rh-lab/lagarias-certificate/v1",
            "classification": (
                "CERTIFIED_FINITE" if finite_range_certified else "INCONCLUSIVE"
            ),
            "claim": LAGARIAS_CLAIM,
            "scope": {"first_n": "1", "last_n": str(limit)},
            "precision_bits": str(precision_bits),
            "backend": {
                "python_flint": flint.__version__,
                "flint": flint.__FLINT_VERSION__,
                "arithmetic": LAGARIAS_ARITHMETIC,
            },
            "checks": {
                "equality_at_one": equality_at_one,
                "inconclusive_n": [str(n) for n in inconclusive],
                "violating_n": [str(n) for n in violations],
            },
            "smallest_recorded_margin_lower_bound": (
                {
                    "n": str(tightest_n),
                    "sigma_n": str(sigma[tightest_n]),
                    "margin": arb_record(tightest_margin),
                }
                if tightest_n is not None
                else None
            ),
            "source": LAGARIAS_SOURCE,
            "limitation": LAGARIAS_LIMITATION,
        }
        return {**payload, "payload_sha256": content_sha256(payload)}
    finally:
        ctx.prec = previous_precision


def verify_lagarias_certificate(artifact: dict[str, Any]) -> dict[str, Any]:
    """Validate hashes and replay a Lagarias certificate at higher precision."""

    supplied_hash = artifact.get("payload_sha256")
    payload = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied_hash != content_sha256(payload):
        raise LagariasCertificateError("payload hash mismatch")
    if artifact.get("schema") != "rh-lab/lagarias-certificate/v1":
        raise LagariasCertificateError("unexpected schema")
    if artifact.get("classification") != "CERTIFIED_FINITE":
        raise LagariasCertificateError("artifact is not a finite certificate")

    try:
        first_n = int(artifact["scope"]["first_n"])
        last_n = int(artifact["scope"]["last_n"])
        original_bits = int(artifact["precision_bits"])
    except (KeyError, TypeError, ValueError) as exc:
        raise LagariasCertificateError("invalid scope or precision") from exc
    if first_n != 1 or last_n < 1:
        raise LagariasCertificateError("certificate must cover a prefix from one")

    canonical = certify_lagarias_range(last_n, original_bits)
    if artifact != canonical:
        raise LagariasCertificateError(
            "artifact does not match canonical regeneration at its stated precision"
        )

    replay_bits = max(original_bits + 64, original_bits * 2)
    replay = certify_lagarias_range(last_n, replay_bits)
    if replay["classification"] != "CERTIFIED_FINITE":
        raise LagariasCertificateError("higher-precision replay was inconclusive")
    original_minimum = artifact["smallest_recorded_margin_lower_bound"]
    replay_minimum = replay["smallest_recorded_margin_lower_bound"]
    if (original_minimum is None) != (replay_minimum is None):
        raise LagariasCertificateError("margin summary changed on replay")
    if (
        original_minimum is not None
        and replay_minimum["n"] != original_minimum["n"]
    ):
        raise LagariasCertificateError("minimum-lower-bound index changed on replay")
    return {
        "classification": "REPRODUCED",
        "original_precision_bits": str(original_bits),
        "replay_precision_bits": str(replay_bits),
        "last_n": str(last_n),
        "replay_payload_sha256": replay["payload_sha256"],
    }
