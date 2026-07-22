"""Enclosure-preserving JSON encodings for FLINT real balls."""

from __future__ import annotations

from typing import Any


def arb_to_dyadic(ball: Any) -> dict[str, str]:
    """Encode an Arb ball as exact dyadic midpoint and outward radius pairs.

    The represented enclosure is ``mid_mantissa * 2**mid_exponent`` plus or
    minus ``radius_mantissa * 2**radius_exponent``.  Strings avoid JSON number
    truncation in other languages.
    """

    mid_mantissa, mid_exponent = ball.mid().man_exp()
    radius_mantissa, radius_exponent = ball.rad().man_exp()
    return {
        "mid_mantissa": str(mid_mantissa),
        "mid_exponent": str(mid_exponent),
        "radius_mantissa": str(radius_mantissa),
        "radius_exponent": str(radius_exponent),
    }


def arb_record(ball: Any, digits: int = 40) -> dict[str, Any]:
    """Return an exact-dyadic encoding plus a human-readable enclosure."""

    return {
        "dyadic": arb_to_dyadic(ball),
        "display": ball.str(digits),
        "is_exact": ball.is_exact(),
    }


def arb_from_dyadic(encoded: dict[str, str]) -> Any:
    """Reconstruct a safe outward enclosure without decimal or float.

    Arb's public constructor can round the radius outward by one magnitude ulp;
    reconstruction is enclosure-preserving, not bitwise round-trip identity.
    """

    from flint import arb

    required = {
        "mid_mantissa",
        "mid_exponent",
        "radius_mantissa",
        "radius_exponent",
    }
    if set(encoded) != required:
        raise ValueError("dyadic ball fields are invalid")
    midpoint = (int(encoded["mid_mantissa"]), int(encoded["mid_exponent"]))
    radius = (int(encoded["radius_mantissa"]), int(encoded["radius_exponent"]))
    if radius[0] < 0:
        raise ValueError("ball radius cannot be negative")
    return arb(midpoint, radius)
