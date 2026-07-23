"""Small-scale full-Gram cross-check for a balanced-multiplier scout vector.

This deliberately slow helper rounds the supplied floating multiplier and the
derived ideal shell to exact dyadics, forms the complete finite natural-dilate
coefficient vector, and evaluates both old and new energies with the rigorous
Arb autocorrelation kernel.  It is intended for independent small-N checks,
not for the large scout grid and not as a frozen certificate.  The supplied
``c_2`` is reported but replaced by the exactly balanced value implied by the
rounded coefficients ``c_3,...,c_K``.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

from flint import ctx

from riemann_lab.nyman import build_natural_system, evaluate_natural_distance


def _round_dyadic(value: float, exponent: int) -> Fraction:
    return Fraction(round(math.ldexp(value, exponent)), 1 << exponent)


def _candidate(path: Path, n: int) -> tuple[Fraction, ...]:
    artifact: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    record = artifact["candidate"]["coefficients"]
    exponent = int(record["denominator_exponent"])
    values = tuple(
        Fraction(int(value), 1 << exponent) for value in record["numerators"]
    )
    if len(values) != n:
        raise ValueError(f"stored candidate dimension is not N={n}")
    return values


def _ideal_shell(
    coefficients: tuple[Fraction, ...],
    exponent: int,
) -> dict[int, Fraction]:
    n = len(coefficients)
    floats = tuple(float(value) for value in coefficients)
    slope = -sum(value / index for index, value in enumerate(floats, start=1))
    cumulative: dict[int, float] = {}
    for m in range(n + 1, 2 * n):
        intercept = 1.0 + sum(
            value * (m // index)
            for index, value in enumerate(floats, start=1)
        )
        weight = 1.0 / (m * (m + 1.0))
        cumulative[m] = -(
            intercept + slope * math.log1p(1.0 / m) / weight
        )

    shell: dict[int, Fraction] = {
        n + 1: _round_dyadic(cumulative[n + 1], exponent)
    }
    for m in range(n + 2, 2 * n):
        shell[m] = _round_dyadic(
            cumulative[m] - cumulative[m - 1], exponent
        )
    # Enforce sum y=0 exactly after rounding all preceding entries.
    shell[2 * n] = -sum(shell.values(), start=Fraction(0))
    return shell


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--multiplier", type=float, nargs="+", required=True)
    parser.add_argument("--bits", type=int, default=192)
    parser.add_argument("--dyadic-bits", type=int, default=80)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()

    path = (
        args.root
        / "results"
        / "nyman-natural-v1"
        / "cells"
        / f"n-{args.n:04d}.json"
    )
    old = _candidate(path, args.n)
    shell = _ideal_shell(old, args.dyadic_bits)
    if len(args.multiplier) < 2 or args.multiplier[0] != 1.0:
        raise ValueError("multiplier must have c1=1 and K>=2")
    z_values = tuple(
        _round_dyadic(value / index, args.dyadic_bits)
        for index, value in enumerate(args.multiplier[2:], start=3)
    )
    supplied_c2 = _round_dyadic(args.multiplier[1], args.dyadic_bits)
    balanced_c2 = -2 * (1 + sum(z_values, start=Fraction(0)))
    multiplier = (
        Fraction(1),
        balanced_c2,
        *(index * value for index, value in enumerate(z_values, start=3)),
    )

    added: dict[int, Fraction] = {}
    for j, y_value in shell.items():
        for k, c_value in enumerate(multiplier, start=1):
            added[j * k] = added.get(j * k, Fraction(0)) + y_value * c_value

    total = {index: value for index, value in enumerate(old, start=1)}
    for index, value in added.items():
        total[index] = total.get(index, Fraction(0)) + value
    total = {index: value for index, value in total.items() if value}

    ctx.prec = args.bits
    old_system = build_natural_system(tuple(range(1, args.n + 1)))
    old_energy = evaluate_natural_distance(old_system, old)
    dilates = tuple(sorted(total))
    new_system = build_natural_system(dilates)
    new_energy = evaluate_natural_distance(
        new_system,
        tuple(total[index] for index in dilates),
    )
    gain = old_energy - new_energy
    print(
        json.dumps(
            {
                "classification": "EXPLORATORY_FULL_GRAM_CROSS_CHECK",
                "hypothesis_status": "UNRESOLVED",
                "n": args.n,
                "multiplier_limit": len(multiplier),
                "aggregate_dimension": len(dilates),
                "dyadic_bits": args.dyadic_bits,
                "supplied_c2_after_dyadic_rounding": str(supplied_c2),
                "balanced_c2_used": str(balanced_c2),
                "c2_balance_adjustment": str(balanced_c2 - supplied_c2),
                "old_energy": str(old_energy),
                "new_energy": str(new_energy),
                "direct_gain": str(gain),
                "direct_gain_fraction": str(gain / old_energy),
                "harmonic_balance_after_rounding": str(
                    sum(
                        value / index
                        for index, value in enumerate(multiplier, start=1)
                    )
                ),
                "limitation": (
                    "The shell and multiplier are rounded scout values. This "
                    "full finite-Gram evaluation is a cross-check, not a "
                    "frozen claim and not evidence resolving RH."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
