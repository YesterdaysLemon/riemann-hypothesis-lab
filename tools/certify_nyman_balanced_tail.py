"""Exact arithmetic core for the balanced-multiplier tail bound.

This module certifies only the infinite interval tail after a caller-supplied
cutoff ``T``.  It does not evaluate the finite prefix and therefore cannot, by
itself, certify a full direct gain.

For a finitely supported balanced coefficient vector ``a``,

    sum(a_n) = sum(a_n / n) = 0,

put

    phi_n(M) = {M/n} - (n - 1)/(2n),
    D(M) = sum_n a_n phi_n(M).

The old finite Nyman vector is denoted by ``p``.  Exact period means follow
from

    mean(phi_m phi_n) = (gcd(m,n)^2 - 1)/(12mn).

The resulting tail of the direct-gain expression ``2X-E`` is enclosed by the
``TailBound`` returned from :func:`tail_center_radius`.  Every field of that
object is a :class:`fractions.Fraction`.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from flint import arb, ctx


ZERO = Fraction(0)
ONE = Fraction(1)
SCOUT_SCHEMA = "rh-lab/nyman-balanced-multiplier-scout/v1"
CANDIDATE_SCHEMA = "rh-lab/nyman-natural-cell/v1"
FROZEN_SCOUT_PAYLOAD_SHA256 = (
    "f31c330b373b3e8616fb1837193d8fbac6e31ce56dd4986364342dc4369f3144"
)
FROZEN_SCOUT_RAW_SHA256 = (
    "6dbd2db984632cec02722d53354dd3dbf87cbcd20be5198aaec991fb10e8cc74"
)


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"nonstandard JSON constant {value!r} is forbidden")


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r} is forbidden")
        result[key] = value
    return result


def strict_json_loads(
    raw: str,
    *,
    parse_float: type[Decimal] | None = None,
) -> object:
    options: dict[str, object] = {
        "object_pairs_hook": _unique_json_object,
        "parse_constant": _reject_json_constant,
    }
    if parse_float is not None:
        options["parse_float"] = parse_float
    return json.loads(raw, **options)


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _verify_payload_hash(
    record: Mapping[str, object],
    *,
    name: str,
    expected: str | None = None,
) -> str:
    supplied = record.get("payload_sha256")
    if not isinstance(supplied, str) or len(supplied) != 64:
        raise ValueError(f"{name} has no valid payload SHA-256")
    body = {key: value for key, value in record.items() if key != "payload_sha256"}
    computed = _canonical_sha256(body)
    if supplied != computed:
        raise ValueError(f"{name} canonical payload SHA-256 mismatch")
    if expected is not None and supplied != expected:
        raise ValueError(f"{name} is not the frozen payload")
    return supplied


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


def coefficient_sum(values: Mapping[int, Fraction]) -> Fraction:
    return sum((Fraction(value) for value in values.values()), start=ZERO)


def harmonic_sum(values: Mapping[int, Fraction]) -> Fraction:
    return sum(
        (
            Fraction(value) / int(index)
            for index, value in values.items()
        ),
        start=ZERO,
    )


def require_balanced(values: Mapping[int, Fraction], name: str = "a") -> None:
    total = coefficient_sum(values)
    harmonic = harmonic_sum(values)
    if total or harmonic:
        raise ValueError(
            f"{name} is not balanced: sum={total}, harmonic_sum={harmonic}"
        )


def round_fraction_to_dyadic(value: Fraction, bits: int) -> Fraction:
    """Round to the nearest multiple of ``2**-bits``, ties to even."""

    if bits < 0:
        raise ValueError("dyadic bit count must be nonnegative")
    scale = 1 << bits
    numerator = round(Fraction(value) * scale)
    return Fraction(numerator, scale)


def round_decimal_to_dyadic(value: Decimal, bits: int) -> Fraction:
    return round_fraction_to_dyadic(Fraction(value), bits)


def load_exact_candidate(path: Path, n: int | None = None) -> dict[int, Fraction]:
    """Load a stored exact-dyadic natural Nyman candidate."""

    artifact = json.loads(path.read_text(encoding="utf-8"))
    record = artifact["candidate"]["coefficients"]
    exponent = int(record["denominator_exponent"])
    denominator = 1 << exponent
    numerators = tuple(int(value) for value in record["numerators"])
    if n is not None and len(numerators) != n:
        raise ValueError(f"stored candidate dimension is not N={n}")
    return {
        index: Fraction(numerator, denominator)
        for index, numerator in enumerate(numerators, start=1)
    }


def round_coefficients(
    values: Mapping[int, Fraction],
    bits: int | None,
) -> dict[int, Fraction]:
    if bits is None:
        return _clean(values)
    return _clean(
        {
            index: round_fraction_to_dyadic(Fraction(value), bits)
            for index, value in values.items()
        }
    )


def _arb_from_fraction(value: Fraction) -> arb:
    exact = Fraction(value)
    return arb(exact.numerator) / exact.denominator


def _ideal_shell_arb(
    coefficients: Mapping[int, Fraction],
) -> tuple[arb, ...]:
    indices = tuple(sorted(coefficients))
    if indices != tuple(range(1, len(indices) + 1)):
        raise ValueError("the old candidate must have contiguous indices 1..N")
    n = len(indices)
    slope = -harmonic_sum(coefficients)
    slope_ball = _arb_from_fraction(slope)
    cumulative: dict[int, arb] = {}
    for m in range(n + 1, 2 * n):
        intercept = ONE + sum(
            (
                coefficients[index] * (m // index)
                for index in indices
            ),
            start=ZERO,
        )
        log_mass = (arb(m + 1) / m).log()
        averaged_residual = (
            _arb_from_fraction(intercept)
            + slope_ball * log_mass * (m * (m + 1))
        )
        cumulative[m] = -averaged_residual

    values = [cumulative[n + 1]]
    values.extend(
        cumulative[m] - cumulative[m - 1]
        for m in range(n + 2, 2 * n)
    )
    # The last value is imposed exactly after dyadic rounding, so this
    # unrounded diagnostic value is not returned.
    values.append(-cumulative[2 * n - 1])
    return tuple(values)


def _unique_nearest_integer(value: arb) -> int | None:
    """Return the nearest integer only when one open rounding bin contains it."""

    midpoint = value.mid()
    candidate_record = (midpoint + arb(1) / 2).floor().unique_fmpz()
    if candidate_record is None:
        return None
    candidate = int(candidate_record)
    difference = value - candidate
    half = arb(1) / 2
    if difference.lower() > -half and difference.upper() < half:
        return candidate
    return None


def rounded_ideal_shell(
    coefficients: Mapping[int, Fraction],
    bits: int,
) -> dict[int, Fraction]:
    """Construct the ideal first shell and round it to exact dyadics.

    Arb precision is increased until every exact logarithmic expression is
    enclosed strictly inside one nearest-dyadic rounding bin.  The final shell
    coefficient is then replaced by the exact negative sum of the preceding
    coefficients.
    """

    if bits < 0:
        raise ValueError("dyadic bit count must be nonnegative")
    n = len(coefficients)
    if n < 2:
        raise ValueError("the old candidate must have dimension at least two")
    rounded: tuple[Fraction, ...] | None = None
    previous_precision = ctx.prec
    try:
        for precision in (128, 192, 256, 384, 512, 768, 1024):
            ctx.prec = max(precision, bits + 64)
            values = _ideal_shell_arb(coefficients)
            integers = tuple(
                _unique_nearest_integer(value * (1 << bits))
                for value in values[:-1]
            )
            if all(value is not None for value in integers):
                rounded = tuple(
                    Fraction(int(value), 1 << bits)
                    for value in integers
                )
                break
        else:
            raise ArithmeticError(
                "Arb could not prove unique ideal-shell dyadic rounding bins"
            )
    finally:
        ctx.prec = previous_precision
    assert rounded is not None
    result = {
        n + offset: value
        for offset, value in enumerate(rounded, start=1)
        if value
    }
    last = -coefficient_sum(result)
    if last:
        result[2 * n] = last
    elif 2 * n in result:
        del result[2 * n]
    if coefficient_sum(result):
        raise ArithmeticError("failed to enforce the exact shell sum")
    return result


def rounded_balanced_multiplier(
    supplied: Sequence[Decimal | Fraction | float | int],
    bits: int,
) -> dict[int, Fraction]:
    """Round ``z_k=c_k/k`` for ``k>=3`` and impose balance exactly."""

    if len(supplied) < 2:
        raise ValueError("a balanced multiplier needs K >= 2")
    first = Fraction(supplied[0])
    if first != 1:
        raise ValueError("the supplied multiplier must have c_1=1")
    z_values = {
        index: round_fraction_to_dyadic(
            Fraction(value) / index,
            bits,
        )
        for index, value in enumerate(supplied[2:], start=3)
    }
    multiplier: dict[int, Fraction] = {1: ONE}
    multiplier[2] = -2 * (
        ONE + sum(z_values.values(), start=ZERO)
    )
    multiplier.update(
        {
            index: index * value
            for index, value in z_values.items()
            if value
        }
    )
    if harmonic_sum(multiplier):
        raise ArithmeticError("failed to enforce exact multiplier balance")
    return _clean(multiplier)


def convolve_coefficients(
    left: Mapping[int, Fraction],
    right: Mapping[int, Fraction],
) -> dict[int, Fraction]:
    """Return the Dirichlet convolution induced by Nyman dilations."""

    result: dict[int, Fraction] = {}
    for left_index, left_value in left.items():
        for right_index, right_value in right.items():
            index = int(left_index) * int(right_index)
            result[index] = (
                result.get(index, ZERO)
                + Fraction(left_value) * Fraction(right_value)
            )
    return _clean(result)


def _scout_cell(
    artifact: Mapping[str, object],
    n: int,
    multiplier_limit: int,
    section: str = "n256_width_sweep",
) -> Mapping[str, object]:
    if section not in {"n256_width_sweep", "fixed_width_grid"}:
        raise ValueError("unsupported frozen scout section")
    collection = artifact.get(section)
    if not isinstance(collection, Mapping):
        raise ValueError(f"frozen scout has no {section}")
    cells = collection.get("cells")
    if not isinstance(cells, list):
        raise ValueError(f"frozen scout {section} has no cells")
    matches = [
        cell
        for cell in cells
        if isinstance(cell, Mapping)
        and int(cell.get("n", -1)) == n
        and int(cell.get("multiplier_limit", -1)) == multiplier_limit
        and cell.get("objective") == "direct-gain"
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one direct-gain N={n}, K={multiplier_limit} cell"
        )
    return matches[0]


@dataclass(frozen=True)
class ExactTailInputs:
    """Exact vectors derived from one frozen exploratory scout cell."""

    n: int
    multiplier_limit: int
    shell_bits: int
    z_bits: int
    old_bits: int | None
    source_candidate: Mapping[int, Fraction]
    old_coefficients: Mapping[int, Fraction]
    shell: Mapping[int, Fraction]
    multiplier: Mapping[int, Fraction]
    added_coefficients: Mapping[int, Fraction]
    source_candidate_path: str
    source_candidate_raw_sha256: str
    source_candidate_payload_sha256: str
    source_candidate_candidate_sha256: str
    scout_raw_sha256: str
    scout_payload_sha256: str
    scout_cell_payload_sha256: str


def build_exact_vectors(
    root: Path,
    *,
    n: int = 256,
    multiplier_limit: int = 16,
    shell_bits: int = 80,
    z_bits: int | None = None,
    old_bits: int | None = None,
    scout_section: str = "n256_width_sweep",
) -> ExactTailInputs:
    """Load, dyadically round, balance, and convolve the frozen scout vectors."""

    if z_bits is None:
        z_bits = shell_bits
    scout_path = root / "results" / "nyman-balanced-multiplier-scout-v1.json"
    raw_scout_bytes = scout_path.read_bytes()
    scout_raw_sha256 = hashlib.sha256(raw_scout_bytes).hexdigest()
    if scout_raw_sha256 != FROZEN_SCOUT_RAW_SHA256:
        raise ValueError("frozen scout raw SHA-256 mismatch")
    raw_scout = raw_scout_bytes.decode("utf-8")
    hash_artifact = strict_json_loads(raw_scout)
    if not isinstance(hash_artifact, Mapping):
        raise ValueError("frozen scout must be a JSON object")
    if hash_artifact.get("schema") != SCOUT_SCHEMA:
        raise ValueError("unexpected frozen scout schema")
    scout_payload_sha256 = _verify_payload_hash(
        hash_artifact,
        name="frozen scout",
        expected=FROZEN_SCOUT_PAYLOAD_SHA256,
    )
    hash_cell = _scout_cell(
        hash_artifact,
        n,
        multiplier_limit,
        scout_section,
    )
    scout_cell_payload_sha256 = _verify_payload_hash(
        hash_cell,
        name="frozen scout cell",
    )

    # Reparse decimal literals exactly only after authenticating their native
    # JSON representation and the selected cell.
    artifact = strict_json_loads(raw_scout, parse_float=Decimal)
    cell = _scout_cell(artifact, n, multiplier_limit, scout_section)
    source_binding = cell.get("source_binding")
    if not isinstance(source_binding, Mapping):
        raise ValueError("scout cell has no source binding")
    relative_path = str(source_binding["path"])
    expected_relative_path = (
        f"results/nyman-natural-v1/cells/n-{n:04d}.json"
    )
    if relative_path != expected_relative_path:
        raise ValueError("scout candidate path does not match requested N")
    candidate_path = root / Path(relative_path)
    raw_candidate = candidate_path.read_bytes()
    raw_sha256 = hashlib.sha256(raw_candidate).hexdigest()
    if raw_sha256 != source_binding["raw_sha256"]:
        raise ValueError("candidate raw SHA-256 does not match frozen scout")
    candidate_artifact = strict_json_loads(raw_candidate.decode("utf-8"))
    if not isinstance(candidate_artifact, Mapping):
        raise ValueError("bound candidate artifact must be a JSON object")
    if candidate_artifact.get("schema") != CANDIDATE_SCHEMA:
        raise ValueError("unexpected candidate artifact schema")
    if int(candidate_artifact.get("cell_contract", {}).get("n", -1)) != n:
        raise ValueError("candidate artifact dimension does not match request")
    candidate_payload_sha256 = _verify_payload_hash(
        candidate_artifact,
        name="candidate artifact",
    )
    candidate_record = candidate_artifact.get("candidate")
    if not isinstance(candidate_record, Mapping):
        raise ValueError("candidate artifact has no candidate record")
    coefficient_record = candidate_record.get("coefficients")
    if not isinstance(coefficient_record, Mapping):
        raise ValueError("candidate artifact has no coefficient record")
    if str(coefficient_record.get("denominator_exponent")) != str(
        source_binding.get("coefficient_denominator_exponent")
    ):
        raise ValueError("candidate coefficient denominator binding changed")
    supplied_candidate_sha256 = candidate_record.get("candidate_sha256")
    candidate_body = {
        key: value
        for key, value in candidate_record.items()
        if key != "candidate_sha256"
    }
    candidate_sha256 = _canonical_sha256(candidate_body)
    if supplied_candidate_sha256 != candidate_sha256:
        raise ValueError("candidate payload SHA-256 mismatch")
    if candidate_sha256 != source_binding["candidate_sha256"]:
        raise ValueError("candidate payload SHA-256 does not match frozen scout")

    source_candidate = load_exact_candidate(candidate_path, n)
    old_coefficients = round_coefficients(source_candidate, old_bits)
    shell = rounded_ideal_shell(source_candidate, shell_bits)
    supplied = cell.get("coefficients")
    if not isinstance(supplied, list):
        raise ValueError("scout cell has no multiplier coefficients")
    if len(supplied) != multiplier_limit:
        raise ValueError("scout cell multiplier coefficient count mismatch")
    multiplier = rounded_balanced_multiplier(supplied, z_bits)
    added = convolve_coefficients(shell, multiplier)
    require_balanced(added, "convolved added vector")

    return ExactTailInputs(
        n=n,
        multiplier_limit=multiplier_limit,
        shell_bits=shell_bits,
        z_bits=z_bits,
        old_bits=old_bits,
        source_candidate=source_candidate,
        old_coefficients=old_coefficients,
        shell=shell,
        multiplier=multiplier,
        added_coefficients=added,
        source_candidate_path=relative_path,
        source_candidate_raw_sha256=raw_sha256,
        source_candidate_payload_sha256=candidate_payload_sha256,
        source_candidate_candidate_sha256=candidate_sha256,
        scout_raw_sha256=scout_raw_sha256,
        scout_payload_sha256=scout_payload_sha256,
        scout_cell_payload_sha256=scout_cell_payload_sha256,
    )


def load_frozen_n256_k16(
    root: Path,
    *,
    dyadic_bits: int = 80,
    old_bits: int | None = None,
) -> ExactTailInputs:
    """Convenience wrapper for the milestone's frozen ``N=256, K=16`` cell."""

    return build_exact_vectors(
        root,
        n=256,
        multiplier_limit=16,
        shell_bits=dyadic_bits,
        z_bits=dyadic_bits,
        old_bits=old_bits,
    )


def centered_fractional_part(interval: int, index: int) -> Fraction:
    """Return ``phi_index(interval)`` exactly for integer ``interval``."""

    if interval < 0:
        raise ValueError("interval must be nonnegative")
    if index < 1:
        raise ValueError("index must be positive")
    residue = interval % index
    return Fraction(2 * residue - (index - 1), 2 * index)


def phi_covariance(left: int, right: int) -> Fraction:
    if left < 1 or right < 1:
        raise ValueError("indices must be positive")
    gcd = math.gcd(left, right)
    return Fraction(gcd * gcd - 1, 12 * left * right)


def alias_step(
    values: Mapping[int, Fraction],
    interval: int,
) -> Fraction:
    return sum(
        (
            Fraction(value) * centered_fractional_part(interval, index)
            for index, value in values.items()
        ),
        start=ZERO,
    )


def old_periodic_center(
    old: Mapping[int, Fraction],
    interval: int,
) -> Fraction:
    c0 = ONE - coefficient_sum(old) / 2
    return c0 - alias_step(old, interval)


def jordan_j2_sieve(limit: int) -> tuple[int, ...]:
    """Return ``J_2(0),...,J_2(limit)`` in sieve time."""

    if limit < 0:
        raise ValueError("limit must be nonnegative")
    values = [index * index for index in range(limit + 1)]
    if limit >= 1:
        values[1] = 1
    for prime in range(2, limit + 1):
        if values[prime] != prime * prime:
            continue
        square = prime * prime
        factor = square - 1
        for multiple in range(prime, limit + 1, prime):
            values[multiple] = values[multiple] // square * factor
    return tuple(values)


def _coefficient_array(
    values: Mapping[int, Fraction],
    limit: int,
    *,
    harmonic: bool,
) -> list[Fraction]:
    result = [ZERO for _ in range(limit + 1)]
    for index, value in values.items():
        if index > limit:
            raise ValueError("coefficient index exceeds requested limit")
        result[index] = (
            Fraction(value) / index if harmonic else Fraction(value)
        )
    return result


def divisor_harmonic_sums(
    values: Mapping[int, Fraction],
    limit: int | None = None,
) -> tuple[Fraction, ...]:
    """Return ``S_d=sum_{d|n} values[n]/n`` in ``O(X log X)``."""

    cleaned = _clean(values)
    support_limit = max(cleaned, default=0)
    if limit is None:
        limit = support_limit
    if limit < support_limit:
        raise ValueError("limit does not cover the coefficient support")
    weighted = _coefficient_array(cleaned, limit, harmonic=True)
    sums = [ZERO for _ in range(limit + 1)]
    for divisor in range(1, limit + 1):
        total = ZERO
        for multiple in range(divisor, limit + 1, divisor):
            total += weighted[multiple]
        sums[divisor] = total
    return tuple(sums)


def jordan_gcd_means(
    old: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
) -> tuple[Fraction, Fraction]:
    """Return exact ``(mu, nu)`` using Jordan-``J_2`` divisor aggregation.

    ``mu`` is the period mean of ``D(M)^2``.  ``nu`` is the period mean of
    ``D(M) * (C0-sum p_m phi_m(M))``.  The added vector must be balanced.
    """

    cleaned_old = _clean(old)
    cleaned_added = _clean(added)
    require_balanced(cleaned_added, "added vector")
    limit = max(
        max(cleaned_old, default=0),
        max(cleaned_added, default=0),
    )
    if limit == 0:
        return ZERO, ZERO
    old_sums = divisor_harmonic_sums(cleaned_old, limit)
    added_sums = divisor_harmonic_sums(cleaned_added, limit)
    jordan = jordan_j2_sieve(limit)
    mu_numerator = ZERO
    nu_numerator = ZERO
    for divisor in range(1, limit + 1):
        weight = jordan[divisor]
        added_value = added_sums[divisor]
        if added_value:
            mu_numerator += weight * added_value * added_value
            nu_numerator += (
                weight * old_sums[divisor] * added_value
            )
    mu = mu_numerator / 12
    nu = -nu_numerator / 12
    if mu < 0:
        raise ArithmeticError("exact Jordan mean of D^2 became negative")
    return mu, nu


def old_periodic_mean_square(
    old: Mapping[int, Fraction],
) -> Fraction:
    """Return the exact period mean of ``old_periodic_center(old, M)^2``.

    If ``P_d = sum_{d|n} old[n]/n`` and
    ``C0 = 1-sum_n old[n]/2``, the covariance identity gives

        mean(r_bar^2) = C0^2 + (1/12) sum_{d>=2} J_2(d) P_d^2.

    The missing ``d=1`` term is exactly the ``-1`` in
    ``gcd(m,n)^2-1``.
    """

    cleaned_old = _clean(old)
    c0 = ONE - coefficient_sum(cleaned_old) / 2
    limit = max(cleaned_old, default=0)
    if limit < 2:
        return c0 * c0
    old_sums = divisor_harmonic_sums(cleaned_old, limit)
    jordan = jordan_j2_sieve(limit)
    variance = sum(
        (
            jordan[divisor]
            * old_sums[divisor]
            * old_sums[divisor]
            for divisor in range(2, limit + 1)
        ),
        start=ZERO,
    ) / 12
    result = c0 * c0 + variance
    if result < 0:
        raise ArithmeticError("exact old periodic mean square became negative")
    return result


def farey_spacing_reciprocal(maximum_denominator: int) -> int:
    """Return a safe reciprocal spacing for reduced Farey frequencies.

    Distinct reduced fractions whose denominators are at most ``Q`` are
    separated modulo one by at least ``1/(Q*(Q-1))`` for ``Q >= 2``.  The
    same bound remains valid after adjoining frequency zero.
    """

    if maximum_denominator < 0:
        raise ValueError("maximum denominator must be nonnegative")
    if maximum_denominator < 2:
        return 0
    return maximum_denominator * (maximum_denominator - 1)


def absolute_moments(
    values: Mapping[int, Fraction],
) -> tuple[Fraction, Fraction]:
    """Return ``A0=sum|a_n|`` and ``A1=sum n|a_n|``."""

    a0 = sum((abs(Fraction(value)) for value in values.values()), start=ZERO)
    a1 = sum(
        (
            int(index) * abs(Fraction(value))
            for index, value in values.items()
        ),
        start=ZERO,
    )
    return a0, a1


def lcm_mobius_weights(limit: int) -> tuple[int, ...]:
    """Return integer weights ``d*prod_{p|d}(1-p)`` through ``limit``."""

    if limit < 0:
        raise ValueError("limit must be nonnegative")
    radical_transform = [1 for _ in range(limit + 1)]
    for prime in range(2, limit + 1):
        if radical_transform[prime] != 1:
            continue
        factor = 1 - prime
        for multiple in range(prime, limit + 1, prime):
            radical_transform[multiple] *= factor
    result = [0 for _ in range(limit + 1)]
    if limit >= 1:
        result[1] = 1
    for index in range(2, limit + 1):
        result[index] = index * radical_transform[index]
    return tuple(result)


def scaled_divisor_absolute_sums(
    values: Mapping[int, Fraction],
    limit: int | None = None,
) -> tuple[Fraction, ...]:
    """Return ``B_d=sum_{d|n}(n/d)|values[n]|`` in ``O(X log X)``."""

    cleaned = _clean(values)
    support_limit = max(cleaned, default=0)
    if limit is None:
        limit = support_limit
    if limit < support_limit:
        raise ValueError("limit does not cover the coefficient support")
    absolute = _coefficient_array(cleaned, limit, harmonic=False)
    absolute = [abs(value) for value in absolute]
    sums = [ZERO for _ in range(limit + 1)]
    for divisor in range(1, limit + 1):
        total = ZERO
        for quotient, multiple in enumerate(
            range(divisor, limit + 1, divisor),
            start=1,
        ):
            total += quotient * absolute[multiple]
        sums[divisor] = total
    return tuple(sums)


def lcm_absolute_self(values: Mapping[int, Fraction]) -> Fraction:
    """Return ``sum_{m,n}|v_m v_n| lcm(m,n)`` by divisor aggregation."""

    cleaned = _clean(values)
    limit = max(cleaned, default=0)
    if limit == 0:
        return ZERO
    divisor_sums = scaled_divisor_absolute_sums(cleaned, limit)
    weights = lcm_mobius_weights(limit)
    result = sum(
        (
            weights[divisor]
            * divisor_sums[divisor]
            * divisor_sums[divisor]
            for divisor in range(1, limit + 1)
        ),
        start=ZERO,
    )
    if result < 0:
        raise ArithmeticError("LCM absolute self-sum became negative")
    return result


def lcm_absolute_cross(
    left: Mapping[int, Fraction],
    right: Mapping[int, Fraction],
) -> Fraction:
    """Return ``sum_{m,n}|left_m right_n| lcm(m,n)``."""

    cleaned_left = _clean(left)
    cleaned_right = _clean(right)
    limit = max(
        max(cleaned_left, default=0),
        max(cleaned_right, default=0),
    )
    if limit == 0:
        return ZERO
    left_sums = scaled_divisor_absolute_sums(cleaned_left, limit)
    right_sums = scaled_divisor_absolute_sums(cleaned_right, limit)
    weights = lcm_mobius_weights(limit)
    result = sum(
        (
            weights[divisor]
            * left_sums[divisor]
            * right_sums[divisor]
            for divisor in range(1, limit + 1)
        ),
        start=ZERO,
    )
    if result < 0:
        raise ArithmeticError("LCM absolute cross-sum became negative")
    return result


@dataclass(frozen=True)
class TailBound:
    """Exact enclosure for the direct-gain tail after intervals ``1..T``."""

    cutoff: int
    mu: Fraction
    nu: Fraction
    p1: Fraction
    c0: Fraction
    a0: Fraction
    a1: Fraction
    lambda_a: Fraction
    lambda_cross: Fraction
    periodic_partial_sum_bound: Fraction
    center: Fraction
    periodic_radius: Fraction
    delta_radius: Fraction
    radius: Fraction
    lower: Fraction
    upper: Fraction


@dataclass(frozen=True)
class LargeSieveTailBound:
    """Exact tail enclosure using the Fourier/Farey large-sieve bound."""

    cutoff: int
    support_limit: int
    farey_spacing_reciprocal: int
    mu: Fraction
    nu: Fraction
    old_periodic_mean_square: Fraction
    new_periodic_mean_square: Fraction
    p1: Fraction
    c0: Fraction
    a0: Fraction
    old_square_discrepancy_bound: Fraction
    new_square_discrepancy_bound: Fraction
    periodic_partial_sum_bound: Fraction
    center: Fraction
    periodic_radius: Fraction
    delta_radius: Fraction
    radius: Fraction
    lower: Fraction
    upper: Fraction


def tail_center_radius(
    old: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    cutoff: int,
) -> TailBound:
    """Enclose the complete ``2X-E`` tail over integer intervals ``M>T``.

    The center is ``(2*nu-mu)/(T+1)``.  The radius contains both the
    zero-mean periodic fluctuation and the logarithmic interval correction
    ``delta_M`` with ``-1/(6M^3) <= delta_M <= 0``.
    """

    if cutoff < 1:
        raise ValueError("tail cutoff must be positive")
    cleaned_old = _clean(old)
    cleaned_added = _clean(added)
    require_balanced(cleaned_added, "added vector")
    mu, nu = jordan_gcd_means(cleaned_old, cleaned_added)
    p1 = harmonic_sum(cleaned_old)
    c0 = ONE - coefficient_sum(cleaned_old) / 2
    a0, a1 = absolute_moments(cleaned_added)
    lambda_a = lcm_absolute_self(cleaned_added)
    lambda_cross = lcm_absolute_cross(cleaned_old, cleaned_added)
    periodic_bound = (
        lambda_a / 3
        + 2 * lambda_cross / 3
        + abs(c0) * a1
    )
    center = (2 * nu - mu) / (cutoff + 1)
    periodic_radius = periodic_bound / (
        (cutoff + 1) * (cutoff + 2)
    )
    delta_radius = abs(p1) * a0 / (12 * cutoff * cutoff)
    radius = periodic_radius + delta_radius
    return TailBound(
        cutoff=cutoff,
        mu=mu,
        nu=nu,
        p1=p1,
        c0=c0,
        a0=a0,
        a1=a1,
        lambda_a=lambda_a,
        lambda_cross=lambda_cross,
        periodic_partial_sum_bound=periodic_bound,
        center=center,
        periodic_radius=periodic_radius,
        delta_radius=delta_radius,
        radius=radius,
        lower=center - radius,
        upper=center + radius,
    )


def large_sieve_tail_center_radius(
    old: Mapping[int, Fraction],
    added: Mapping[int, Fraction],
    cutoff: int,
) -> LargeSieveTailBound:
    """Enclose the complete tail with a Farey large-sieve discrepancy bound.

    For ``g_M=sum a_n phi_n(M)``, every nonzero Fourier frequency has a
    reduced denominator at most the coefficient support limit ``Q``.  The
    Montgomery--Vaughan large sieve and periodic complementation therefore
    bound every interval square discrepancy by its period mean times
    ``D=Q*(Q-1)``.  If ``v`` is the old periodic residual center, then the
    direct-gain summand has the difference-of-squares identity

        2*g_M*v_M-g_M^2 = v_M^2-(v_M-g_M)^2.

    Its zero-mean interval sums are therefore bounded by
    ``D*(rho+tau)``, where ``rho=mean(v^2)`` and
    ``tau=mean((v-g)^2)=rho+mu-2*nu``.  The implemented bound is entirely
    rational and needs no coefficientwise absolute values.
    """

    if cutoff < 1:
        raise ValueError("tail cutoff must be positive")
    cleaned_old = _clean(old)
    cleaned_added = _clean(added)
    require_balanced(cleaned_added, "added vector")
    support_limit = max(
        max(cleaned_old, default=0),
        max(cleaned_added, default=0),
    )
    spacing_reciprocal = farey_spacing_reciprocal(support_limit)
    mu, nu = jordan_gcd_means(cleaned_old, cleaned_added)
    rho = old_periodic_mean_square(cleaned_old)
    tau = rho + mu - 2 * nu
    if tau < 0:
        raise ArithmeticError("exact new periodic mean square became negative")
    p1 = harmonic_sum(cleaned_old)
    c0 = ONE - coefficient_sum(cleaned_old) / 2
    a0, _ = absolute_moments(cleaned_added)
    old_square_bound = spacing_reciprocal * rho
    new_square_bound = spacing_reciprocal * tau
    periodic_bound = old_square_bound + new_square_bound
    center = (2 * nu - mu) / (cutoff + 1)
    periodic_radius = periodic_bound / (
        (cutoff + 1) * (cutoff + 2)
    )
    delta_radius = abs(p1) * a0 / (12 * cutoff * cutoff)
    radius = periodic_radius + delta_radius
    return LargeSieveTailBound(
        cutoff=cutoff,
        support_limit=support_limit,
        farey_spacing_reciprocal=spacing_reciprocal,
        mu=mu,
        nu=nu,
        old_periodic_mean_square=rho,
        new_periodic_mean_square=tau,
        p1=p1,
        c0=c0,
        a0=a0,
        old_square_discrepancy_bound=old_square_bound,
        new_square_discrepancy_bound=new_square_bound,
        periodic_partial_sum_bound=periodic_bound,
        center=center,
        periodic_radius=periodic_radius,
        delta_radius=delta_radius,
        radius=radius,
        lower=center - radius,
        upper=center + radius,
    )


def pairwise_lcm_absolute(
    left: Mapping[int, Fraction],
    right: Mapping[int, Fraction],
) -> Fraction:
    """Slow reference formula used by small independent tests."""

    return sum(
        (
            abs(Fraction(left_value) * Fraction(right_value))
            * math.lcm(int(left_index), int(right_index))
            for left_index, left_value in left.items()
            for right_index, right_value in right.items()
        ),
        start=ZERO,
    )


def period_mean(values: Iterable[Fraction]) -> Fraction:
    sequence = tuple(Fraction(value) for value in values)
    if not sequence:
        raise ValueError("a period must not be empty")
    return sum(sequence, start=ZERO) / len(sequence)
