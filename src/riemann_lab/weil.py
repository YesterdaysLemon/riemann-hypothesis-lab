"""Rigorous finite Weil quadratic-form matrices and mutation controls.

The normalization follows the multiplicative Fourier basis used in the finite
cutoff construction documented in ``docs/weil-matrix-v1.md``.  The matrix
``A = P - R - S`` represents ``Q`` with the sign convention that RH is
equivalent to nonnegativity of the full infinite-dimensional form.  Positive
finite matrices are therefore finite evidence only.
"""

from __future__ import annotations

from math import gcd, isqrt
from typing import Any, Sequence

import flint
from flint import acb, arb, arb_mat, ctx

from .artifacts import content_sha256
from .balls import arb_from_dyadic, arb_record, arb_to_dyadic


BOMBIERI_SOURCE = (
    "https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0"
)
CCM_SOURCE = "https://arxiv.org/abs/2511.22755v1"
CONNES_CONSANI_CORE_SOURCE = "https://arxiv.org/abs/2106.01715v1"
WEIL_CLAIM = (
    "For c=5/2 and modes -4 through 4, the complete finite Weil matrix "
    "A=P-R-S is positive definite by interval LDL^T; deleting its sole "
    "allowed prime-power contribution produces a certified negative control."
)
WEIL_LIMITATION = (
    "Positive definiteness of one finite compression does not prove RH or "
    "exclude negative directions at higher modes or other cutoffs. The "
    "deliberately corrupted negative control is not a counterexample to RH."
)
CONTROL_WITNESS = (
    -73201,
    -108782,
    -200275,
    -1000000,
    0,
    1000000,
    200275,
    108782,
    73201,
)


class WeilCertificateError(ValueError):
    """Raised when a serialized finite Weil certificate fails replay."""


def _validate_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")


def _reduced_cutoff(numerator: int, denominator: int) -> tuple[int, int]:
    _validate_integer(numerator, "cutoff numerator")
    _validate_integer(denominator, "cutoff denominator")
    if denominator <= 0 or numerator <= denominator:
        raise ValueError("cutoff must be a positive rational strictly above one")
    common = gcd(numerator, denominator)
    return numerator // common, denominator // common


def prime_powers_leq(
    cutoff_numerator: int,
    cutoff_denominator: int = 1,
) -> list[tuple[int, int, int]]:
    """Return ``(prime, exponent, power)`` with exact ``power <= cutoff``.

    Integer cross-multiplication is used for the boundary.  No binary floating
    point value participates in deciding whether a prime power is present.
    """

    numerator, denominator = _reduced_cutoff(
        cutoff_numerator, cutoff_denominator
    )
    limit = numerator // denominator
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for candidate in range(2, isqrt(limit) + 1):
        if sieve[candidate]:
            start = candidate * candidate
            sieve[start : limit + 1 : candidate] = b"\x00" * (
                (limit - start) // candidate + 1
            )

    powers: list[tuple[int, int, int]] = []
    for prime in (value for value in range(2, limit + 1) if sieve[value]):
        exponent = 1
        power = prime
        while power * denominator <= numerator:
            powers.append((prime, exponent, power))
            exponent += 1
            power *= prime
    return sorted(powers, key=lambda item: (item[2], item[0], item[1]))


def weil_kernel(m: int, n: int, y: Any, length: Any) -> Any:
    """Evaluate the nonnegative-half autocorrelation kernel ``q_mn(y)``."""

    _validate_integer(m, "m")
    _validate_integer(n, "n")
    pi = arb.pi()
    if m == n:
        return 2 * (1 - y / length) * (2 * pi * n * y / length).cos()
    return (
        (2 * pi * m * y / length).sin()
        - (2 * pi * n * y / length).sin()
    ) / (pi * (n - m))


def pole_entry(m: int, n: int, length: Any) -> Any:
    """Evaluate the closed-form pole contribution ``P_mn``."""

    pi = arb.pi()
    return (
        32
        * length
        * (length / 4).sinh() ** 2
        * (length**2 - 16 * pi**2 * m * n)
        / (
            (length**2 + 16 * pi**2 * m * m)
            * (length**2 + 16 * pi**2 * n * n)
        )
    )


def _segmented_unit_integral(function: Any) -> Any:
    tolerance = arb((1, -(ctx.prec - 16)))
    cuts = [arb(0), arb(1) / 4, arb(1) / 2, 3 * arb(1) / 4, arb(1)]
    integral = acb(0)
    for left, right in zip(cuts, cuts[1:]):
        integral += acb.integral(
            function,
            left,
            right,
            rel_tol=tolerance,
            abs_tol=tolerance,
            deg_limit=80,
            eval_limit=100000,
            depth_limit=30,
            use_heap=True,
        )
    return integral


def pole_entry_integral(m: int, n: int, length: Any) -> Any:
    """Independently enclose ``P_mn`` from its defining integral."""

    pi = acb.pi()

    def integrand(t: Any, analytic: bool) -> Any:
        del analytic
        if m == n:
            kernel = 2 * (1 - t) * (2 * pi * n * t).cos()
        else:
            kernel = (
                (2 * pi * m * t).sin() - (2 * pi * n * t).sin()
            ) / (pi * (n - m))
        return (
            length
            * kernel
            * ((length * t / 2).exp() + (-length * t / 2).exp())
        )

    integral = _segmented_unit_integral(integrand)
    if not integral.is_finite() or not integral.imag.contains(0):
        raise ArithmeticError("pole integral oracle did not enclose a real value")
    return integral.real


def _spectral_parameter(index: int, length: Any) -> Any:
    return acb(arb(1) / 4, -arb.pi() * index / length)


def archimedean_entry(
    m: int,
    n: int,
    cutoff: Any,
    length: Any,
) -> Any:
    """Evaluate ``R_mn`` by a rigorous digamma-Lerch closed form.

    This is algebraically equal to the desingularized integral in the frozen
    specification.  Integer Fourier frequencies make the finite-endpoint tail
    especially simple; the derivation is recorded in the methodology note.
    """

    pi = arb.pi()
    lerch_base = acb((-2 * length).exp())
    endpoint_factor = (-length / 2).exp()
    z_m = _spectral_parameter(m, length)
    z_n = _spectral_parameter(n, length)

    if m != n:
        infinite_part = (z_n.digamma() - z_m.digamma()).imag / (
            2 * pi * (n - m)
        )
        finite_tail = endpoint_factor * (
            lerch_base.lerch_phi(1, z_m)
            - lerch_base.lerch_phi(1, z_n)
        ).imag / (2 * pi * (n - m))
        return infinite_part - finite_tail

    constant = arb.const_euler() + (
        4 * pi * (cutoff - 1) / (cutoff + 1)
    ).log()
    z = z_n
    infinite_part = (
        acb((arb(1) / 2).digamma())
        - z.digamma()
        - z.polygamma(1) / (2 * length)
    ).real
    finite_tail_correction = (
        endpoint_factor * lerch_base.lerch_phi(2, z).real / (2 * length)
        + 2 * (-length).exp().atanh()
    )
    return constant + infinite_part + finite_tail_correction


def _desingularized_archimedean_integrand(
    m: int,
    n: int,
    t: Any,
    length: Any,
) -> Any:
    """Dimensionless Acb integrand on ``t in [0,1]``."""

    imaginary_unit = acb(0, 1)
    pi = acb.pi()
    denominator = (imaginary_unit * length * t).sinc()
    if m != n:
        numerator = (
            m * (2 * pi * m * t).sinc()
            - n * (2 * pi * n * t).sinc()
        )
        return (
            (length * t / 2).exp()
            * numerator
            / ((n - m) * denominator)
        )

    basis = (1 - t) * (2 * pi * n * t).cos()
    exponential = (
        length
        / 2
        * (length * t / 4).exp()
        * (imaginary_unit * length * t / 4).sinc()
    )
    correction = (
        -2 * pi**2 * n * n * t * (pi * n * t).sinc() ** 2
        - (2 * pi * n * t).cos()
    )
    return (exponential * basis + correction) / denominator


def archimedean_entry_integral(
    m: int,
    n: int,
    cutoff: Any,
    length: Any,
) -> Any:
    """Independently enclose ``R_mn`` by segmented Acb quadrature."""

    integral = _segmented_unit_integral(
        lambda t, analytic: _desingularized_archimedean_integrand(
            m, n, t, length
        )
    )
    if not integral.is_finite() or not integral.imag.contains(0):
        raise ArithmeticError("archimedean integral oracle did not enclose a real value")
    if m == n:
        integral += arb.const_euler() + (
            4 * arb.pi() * (cutoff - 1) / (cutoff + 1)
        ).log()
    return integral.real


def prime_power_entry(
    m: int,
    n: int,
    length: Any,
    prime_powers: Sequence[tuple[int, int, int]],
    cutoff_numerator: int,
    cutoff_denominator: int,
) -> Any:
    """Evaluate the complete prime-power contribution ``S_mn``."""

    numerator, denominator = _reduced_cutoff(
        cutoff_numerator, cutoff_denominator
    )
    total = arb(0)
    for prime, exponent, power in prime_powers:
        if power * denominator > numerator:
            raise ValueError("prime-power transcript exceeds the exact cutoff")
        # q(L)=0 exactly.  Preserve that theorem at an exact prime-power
        # boundary instead of replacing zero by a small trigonometric ball.
        if power * denominator == numerator:
            continue
        log_prime = arb(prime).log()
        location = exponent * log_prime
        total += (
            log_prime
            / arb(power).sqrt()
            * weil_kernel(m, n, location, length)
        )
    return total


def _symmetric_components(
    cutoff_numerator: int,
    cutoff_denominator: int,
    degree: int,
) -> dict[str, Any]:
    numerator, denominator = _reduced_cutoff(
        cutoff_numerator, cutoff_denominator
    )
    _validate_integer(degree, "degree")
    if degree < 0:
        raise ValueError("degree must be nonnegative")

    cutoff = arb(numerator) / denominator
    length = cutoff.log()
    indices = list(range(-degree, degree + 1))
    size = len(indices)
    powers = prime_powers_leq(numerator, denominator)
    endpoint_power = numerator // denominator if numerator % denominator == 0 else None
    oracle_pairs = 0
    pole_oracle_pairs = 0
    oracle_accuracy_bound = arb((1, -(ctx.prec - 32)))
    matrices = {
        name: [[arb(0) for _ in range(size)] for _ in range(size)]
        for name in ("pole", "archimedean", "prime_power", "q")
    }

    for row, m in enumerate(indices):
        for column in range(row, size):
            n = indices[column]
            pole = pole_entry(m, n, length)
            pole_oracle = pole_entry_integral(m, n, length)
            if not pole.overlaps(pole_oracle):
                raise ArithmeticError(
                    f"pole evaluators disagree for matrix entry {(m, n)}"
                )
            if not (pole_oracle.rad() < oracle_accuracy_bound):
                raise ArithmeticError(
                    f"pole integral too wide for matrix entry {(m, n)}"
                )
            pole_oracle_pairs += 1
            archimedean_closed_form = archimedean_entry(
                m, n, cutoff, length
            )
            archimedean_oracle = archimedean_entry_integral(
                m, n, cutoff, length
            )
            if not archimedean_closed_form.overlaps(archimedean_oracle):
                raise ArithmeticError(
                    f"archimedean evaluators disagree for matrix entry {(m, n)}"
                )
            if not (archimedean_oracle.rad() < oracle_accuracy_bound):
                raise ArithmeticError(
                    f"archimedean integral too wide for matrix entry {(m, n)}"
                )
            # The source-integral evaluator supplies the certified matrix entry;
            # the special-function identity is an independent cross-check.
            archimedean = archimedean_oracle
            oracle_pairs += 1
            prime_power = prime_power_entry(
                m,
                n,
                length,
                powers,
                cutoff_numerator=numerator,
                cutoff_denominator=denominator,
            )
            q_value = pole - archimedean - prime_power
            for name, value in (
                ("pole", pole),
                ("archimedean", archimedean),
                ("prime_power", prime_power),
                ("q", q_value),
            ):
                matrices[name][row][column] = value
                matrices[name][column][row] = value

    return {
        "cutoff": cutoff,
        "length": length,
        "indices": indices,
        "prime_powers": powers,
        "endpoint_power": endpoint_power,
        "archimedean_oracle_pairs": oracle_pairs,
        "pole_oracle_pairs": pole_oracle_pairs,
        **matrices,
    }


def certify_positive_ldlt(matrix: Sequence[Sequence[Any]]) -> dict[str, Any]:
    """Prove positive definiteness by interval ``LDL^T`` pivots."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    for row in range(size):
        for column in range(row):
            if arb_to_dyadic(matrix[row][column]) != arb_to_dyadic(
                matrix[column][row]
            ):
                raise ValueError("matrix must use identical symmetric enclosures")

    lower = [[arb(0) for _ in range(size)] for _ in range(size)]
    pivots: list[Any] = []
    for row in range(size):
        lower[row][row] = arb(1)
        for column in range(row):
            value = matrix[row][column]
            for prior in range(column):
                value -= (
                    lower[row][prior]
                    * lower[column][prior]
                    * pivots[prior]
                )
            value /= pivots[column]
            lower[row][column] = value

        pivot = matrix[row][row]
        for prior in range(row):
            pivot -= lower[row][prior] ** 2 * pivots[prior]
        pivots.append(pivot)
        if not (pivot > 0):
            return {
                "classification": "INCONCLUSIVE",
                "failed_pivot_index": str(row),
                "pivots": [
                    {"index": str(index), "value": arb_record(value)}
                    for index, value in enumerate(pivots)
                ],
            }

    return {
        "classification": "POSITIVE_DEFINITE",
        "failed_pivot_index": None,
        "pivots": [
            {"index": str(index), "value": arb_record(value)}
            for index, value in enumerate(pivots)
        ],
    }


def evaluate_integer_rayleigh(
    matrix: Sequence[Sequence[Any]],
    witness: Sequence[int],
) -> dict[str, Any]:
    """Evaluate and classify an exact integer-vector Rayleigh quotient."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    if len(witness) != size:
        raise ValueError("witness length must equal matrix dimension")
    if any(isinstance(value, bool) or not isinstance(value, int) for value in witness):
        raise TypeError("witness entries must be integers")
    for row in range(size):
        for column in range(row):
            if arb_to_dyadic(matrix[row][column]) != arb_to_dyadic(
                matrix[column][row]
            ):
                raise ValueError("matrix must use identical symmetric enclosures")

    numerator = arb(0)
    norm_squared = 0
    for row, coefficient in enumerate(witness):
        norm_squared += coefficient * coefficient
        numerator += coefficient * coefficient * matrix[row][row]
        for column in range(row + 1, size):
            numerator += (
                2
                * coefficient
                * witness[column]
                * matrix[row][column]
            )
    if norm_squared == 0:
        raise ValueError("witness must be nonzero")
    quotient = numerator / norm_squared
    is_negative = quotient < 0
    is_positive = quotient > 0
    return {
        "sign": (
            "NEGATIVE"
            if is_negative
            else "POSITIVE"
            if is_positive
            else "INCONCLUSIVE"
        ),
        "witness": [str(value) for value in witness],
        "norm_squared": str(norm_squared),
        "numerator": arb_record(numerator),
        "rayleigh_quotient": arb_record(quotient),
        "upper_bound_is_negative": is_negative,
        "lower_bound_is_positive": is_positive,
    }


def certify_integer_witness(
    matrix: Sequence[Sequence[Any]],
    witness: Sequence[int],
) -> dict[str, Any]:
    """Wrap a negative Rayleigh result as the frozen mutation control."""

    evaluation = evaluate_integer_rayleigh(matrix, witness)
    return {
        "classification": (
            "CONTROL_NEGATIVE"
            if evaluation["upper_bound_is_negative"]
            else "INCONCLUSIVE"
        ),
        "witness": evaluation["witness"],
        "norm_squared": evaluation["norm_squared"],
        "numerator": evaluation["numerator"],
        "rayleigh_quotient": evaluation["rayleigh_quotient"],
        "upper_bound_is_negative": evaluation["upper_bound_is_negative"],
    }


def certify_rump_spectrum(matrix: Sequence[Sequence[Any]]) -> dict[str, Any]:
    """Run FLINT's Rump eigenvalue enclosure as a secondary oracle."""

    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    unordered_eigenvalues = arb_mat(matrix).eig(algorithm="rump")
    eigenvalues: list[Any] = []
    separated_order = True
    for value in unordered_eigenvalues:
        inserted = False
        for index, existing in enumerate(eigenvalues):
            if value.real < existing.real:
                eigenvalues.insert(index, value)
                inserted = True
                break
            if not (existing.real < value.real):
                separated_order = False
        if not inserted:
            eigenvalues.append(value)
    records: list[dict[str, Any]] = []
    all_positive = True
    for index, value in enumerate(eigenvalues):
        real = value.real
        imaginary = value.imag
        certified_real = imaginary.contains(0)
        positive = certified_real and real > 0
        all_positive = all_positive and positive
        records.append(
            {
                "index": str(index),
                "real": arb_record(real),
                "imaginary": arb_record(imaginary),
                "imaginary_contains_zero": certified_real,
                "real_part_is_positive": positive,
            }
        )

    separated_order = separated_order and all(
        left.real < right.real
        for left, right in zip(eigenvalues, eigenvalues[1:])
    )
    smallest_index = 0
    smallest_is_separated = all(
        index == smallest_index
        or eigenvalues[smallest_index].real < value.real
        for index, value in enumerate(eigenvalues)
    )
    smallest = eigenvalues[smallest_index]
    return {
        "classification": (
            "POSITIVE_SPECTRUM"
            if all_positive and separated_order
            else "INCONCLUSIVE"
        ),
        "algorithm": "flint-arb_mat-eig-rump",
        "ordering": "ascending by pairwise-separated real enclosures",
        "eigenvalues": records,
        "smallest_eigenvalue": {
            "ordered_index": str(smallest_index),
            "real": arb_record(smallest.real),
            "imaginary": arb_record(smallest.imag),
            "separated_from_all_others": smallest_is_separated,
        },
    }


def _serialize_entries(components: dict[str, Any]) -> list[dict[str, Any]]:
    indices = components["indices"]
    entries: list[dict[str, Any]] = []
    for row, m in enumerate(indices):
        for column in range(row, len(indices)):
            entries.append(
                {
                    "m": str(m),
                    "n": str(indices[column]),
                    "pole": arb_record(components["pole"][row][column]),
                    "archimedean": arb_record(
                        components["archimedean"][row][column]
                    ),
                    "prime_power": arb_record(
                        components["prime_power"][row][column]
                    ),
                    "q": arb_record(components["q"][row][column]),
                }
            )
    return entries


def _prime_transcript(components: dict[str, Any]) -> list[dict[str, Any]]:
    transcript: list[dict[str, Any]] = []
    for prime, exponent, power in components["prime_powers"]:
        log_prime = arb(prime).log()
        location = (
            components["length"]
            if power == components["endpoint_power"]
            else exponent * log_prime
        )
        transcript.append(
            {
                "prime": str(prime),
                "exponent": str(exponent),
                "power": str(power),
                "log_location": arb_record(location),
                "weight": arb_record(log_prime / arb(power).sqrt()),
            }
        )
    return transcript


def certify_weil_matrix(
    cutoff_numerator: int = 5,
    cutoff_denominator: int = 2,
    degree: int = 4,
    precision_bits: int = 192,
) -> dict[str, Any]:
    """Create a rigorous finite Weil-matrix certificate and mutation control."""

    numerator, denominator = _reduced_cutoff(
        cutoff_numerator, cutoff_denominator
    )
    _validate_integer(degree, "degree")
    _validate_integer(precision_bits, "precision_bits")
    if degree < 0:
        raise ValueError("degree must be nonnegative")
    if precision_bits < 96:
        raise ValueError("precision_bits must be at least 96")
    if (numerator, denominator, degree) != (5, 2, 4):
        raise ValueError(
            "v1 certificate is frozen to cutoff c=5/2 and degree N=4"
        )

    previous_precision = ctx.prec
    try:
        ctx.prec = precision_bits
        components = _symmetric_components(numerator, denominator, degree)
        positive = certify_positive_ldlt(components["q"])
        spectrum = certify_rump_spectrum(components["q"])

        control: dict[str, Any] | None = None
        if (
            numerator == 5
            and denominator == 2
            and degree == 4
            and components["prime_powers"] == [(2, 1, 2)]
        ):
            mutated = [
                [
                    components["pole"][row][column]
                    - components["archimedean"][row][column]
                    for column in range(len(components["indices"]))
                ]
                for row in range(len(components["indices"]))
            ]
            control = {
                "mutation": {
                    "id": "omit-prime-power-2",
                    "operator": "delete",
                    "target": {"prime": "2", "exponent": "1", "power": "2"},
                    "expected_prime_power_transcript_length": "1",
                },
                **certify_integer_witness(mutated, CONTROL_WITNESS),
            }

        ldlt_certified = positive["classification"] == "POSITIVE_DEFINITE"
        rump_certified = (
            spectrum["classification"] == "POSITIVE_SPECTRUM"
            and spectrum["smallest_eigenvalue"][
                "separated_from_all_others"
            ]
        )
        certified = ldlt_certified and rump_certified
        control_passed = (
            control is not None
            and control["classification"] == "CONTROL_NEGATIVE"
        )
        payload = {
            "schema": "rh-lab/weil-matrix-certificate/v1",
            "classification": (
                "CERTIFIED_FINITE" if certified and control_passed else "INCONCLUSIVE"
            ),
            "claim": WEIL_CLAIM,
            "limitation": WEIL_LIMITATION,
            "scope": {
                "cutoff_c": {
                    "numerator": str(numerator),
                    "denominator": str(denominator),
                },
                "log_cutoff": arb_record(components["length"]),
                "degree": str(degree),
                "indices": [str(value) for value in components["indices"]],
                "dimension": str(len(components["indices"])),
            },
            "precision_bits": str(precision_bits),
            "normalization": {
                "haar_measure": "du/u",
                "basis": "L^(-1/2) exp(2*pi*i*n*x/L) on x in [0,L]",
                "length": "L=log(c)",
                "cutoff_parameter": "c=lambda^2",
                "matrix": "A=P-R-S",
                "rh_sign": "RH iff Q(g)>=0 for every admissible autocorrelation",
                "domain": (
                    "the Fourier basis lies in the Laurent-polynomial core "
                    "C[U,U^-1] of the semilocal Weil form"
                ),
                "prime_power_factor": "(log p)*p^(-r/2), with no extra factor 2",
            },
            "source_equations": {
                "bombieri": BOMBIERI_SOURCE,
                "connes_consani_moscovici": CCM_SOURCE,
                "semilocal_core": CONNES_CONSANI_CORE_SOURCE,
                "core_equation_map": (
                    "Connes-Consani Proposition 2.1, Lemma 2.2, "
                    "Proposition 2.3, and Corollary 2.4"
                ),
                "equation_map": (
                    "CCM equations (2.6), (2.8)-(2.10), (3.5), "
                    "(3.10)-(3.18), and (4.2)-(4.4)"
                ),
                "archimedean_evaluation": (
                    "segmented Acb integral of the desingularized source "
                    "integrand; digamma, trigamma, and Lerch Phi closed form "
                    "used as a cross-check"
                ),
            },
            "backend": {
                "python_flint": flint.__version__,
                "flint": flint.__FLINT_VERSION__,
                "arithmetic": (
                    "exact integer/rational combinatorics; Arb/Acb ball "
                    "special functions and interval LDL^T"
                ),
                "quadrature": (
                    "four equal segments on t in [0,1]; "
                    "rel_tol=abs_tol=2^(-(precision_bits-16)); deg_limit=80; "
                    "eval_limit=100000; depth_limit=30; heap enabled"
                ),
            },
            "prime_power_transcript": _prime_transcript(components),
            "matrix_entries_upper_triangle": _serialize_entries(components),
            "positive_definiteness": positive,
            "secondary_rump_spectrum": spectrum,
            "negative_control": control,
            "checks": {
                "exact_prime_power_cutoff": True,
                "symmetric_enclosures_by_construction": True,
                "archimedean_closed_form_matches_integral_oracle": (
                    components["archimedean_oracle_pairs"]
                    == len(components["indices"])
                    * (len(components["indices"]) + 1)
                    // 2
                ),
                "pole_closed_form_matches_integral_oracle": (
                    components["pole_oracle_pairs"]
                    == len(components["indices"])
                    * (len(components["indices"]) + 1)
                    // 2
                ),
                "all_ldlt_pivots_positive": ldlt_certified,
                "rump_eigenvalue_enclosures_positive": rump_certified,
                "rump_smallest_eigenvalue_separated": spectrum[
                    "smallest_eigenvalue"
                ]["separated_from_all_others"],
                "negative_control_detected": (
                    control is not None
                    and control["classification"] == "CONTROL_NEGATIVE"
                ),
            },
        }
        return {**payload, "payload_sha256": content_sha256(payload)}
    finally:
        ctx.prec = previous_precision


def _check_precision_containment(
    original: dict[str, Any], replay: dict[str, Any]
) -> None:
    original_entries = {
        (entry["m"], entry["n"]): entry
        for entry in original["matrix_entries_upper_triangle"]
    }
    replay_entries = {
        (entry["m"], entry["n"]): entry
        for entry in replay["matrix_entries_upper_triangle"]
    }
    if original_entries.keys() != replay_entries.keys():
        raise WeilCertificateError("matrix index schedule changed on replay")
    for key, original_entry in original_entries.items():
        replay_entry = replay_entries[key]
        for component in ("pole", "archimedean", "prime_power", "q"):
            original_ball = arb_from_dyadic(original_entry[component]["dyadic"])
            replay_ball = arb_from_dyadic(replay_entry[component]["dyadic"])
            if not original_ball.contains(replay_ball):
                raise WeilCertificateError(
                    f"higher-precision {component} entry {key} escaped stored enclosure"
                )


def verify_weil_certificate(
    artifact: dict[str, Any],
    replay_precision_bits: int = 384,
) -> dict[str, Any]:
    """Validate, canonically regenerate, and replay a Weil certificate."""

    if (
        isinstance(replay_precision_bits, bool)
        or not isinstance(replay_precision_bits, int)
        or replay_precision_bits < 96
    ):
        raise WeilCertificateError("replay_precision_bits must be an integer >= 96")

    supplied_hash = artifact.get("payload_sha256")
    payload = {key: value for key, value in artifact.items() if key != "payload_sha256"}
    if supplied_hash != content_sha256(payload):
        raise WeilCertificateError("payload hash mismatch")
    if artifact.get("schema") != "rh-lab/weil-matrix-certificate/v1":
        raise WeilCertificateError("unexpected schema")
    if artifact.get("classification") != "CERTIFIED_FINITE":
        raise WeilCertificateError("artifact is not a finite certificate")

    try:
        cutoff = artifact["scope"]["cutoff_c"]
        numerator = int(cutoff["numerator"])
        denominator = int(cutoff["denominator"])
        degree = int(artifact["scope"]["degree"])
        original_bits = int(artifact["precision_bits"])
    except (KeyError, TypeError, ValueError) as exc:
        raise WeilCertificateError("invalid scope or precision") from exc

    try:
        canonical = certify_weil_matrix(
            numerator, denominator, degree, original_bits
        )
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise WeilCertificateError("artifact parameters cannot be regenerated") from exc
    if artifact != canonical:
        raise WeilCertificateError(
            "artifact does not match canonical regeneration at stated precision"
        )

    replay_bits = max(replay_precision_bits, original_bits + 64)
    try:
        replay = certify_weil_matrix(numerator, denominator, degree, replay_bits)
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise WeilCertificateError("higher-precision replay failed") from exc
    if replay["classification"] != "CERTIFIED_FINITE":
        raise WeilCertificateError("higher-precision replay was inconclusive")
    _check_precision_containment(artifact, replay)
    return {
        "classification": "REPRODUCED",
        "original_precision_bits": str(original_bits),
        "replay_precision_bits": str(replay_bits),
        "dimension": artifact["scope"]["dimension"],
        "prime_power_count": str(len(artifact["prime_power_transcript"])),
        "replay_payload_sha256": replay["payload_sha256"],
    }
