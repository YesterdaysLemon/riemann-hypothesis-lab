from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

import pytest
from flint import arb, ctx
from hypothesis import given, settings, strategies as st

from riemann_lab.balls import arb_to_dyadic
from riemann_lab.nyman import (
    N8_FROZEN_DENOMINATOR_EXPONENT,
    N8_FROZEN_DILATES,
    N8_FROZEN_LOWER_BOUND,
    N8_FROZEN_NUMERATORS,
    N8_FROZEN_UPPER_BOUND,
    NaturalSystem,
    autocorrelation_a,
    autocorrelation_reciprocity_residual,
    build_augmented_lower_matrix,
    build_natural_system,
    certify_augmented_lower_bound,
    certify_dyadic_upper_bound,
    certify_frozen_n8_bracket,
    dyadic_coefficients,
    evaluate_natural_distance,
    fixed_order_interval_ldlt,
    natural_gram_entry,
    natural_system_content_sha256,
    natural_target_entry,
    prefix_natural_system,
    propose_dyadic_coefficients,
    round_arb_to_dyadic_numerator,
    round_arb_vector_to_dyadic,
    vasyunin_sum,
)


@pytest.fixture(autouse=True)
def _high_precision() -> None:
    previous_precision = ctx.prec
    ctx.prec = 192
    try:
        yield
    finally:
        ctx.prec = previous_precision


def test_vasyunin_sum_and_a_have_frozen_signs_and_factors() -> None:
    assert vasyunin_sum(7, 1).is_zero()
    assert vasyunin_sum(1, 2).contains(0)
    assert autocorrelation_a(Fraction(1, 1)).overlaps(
        (2 * arb.pi()).log() - arb.const_euler()
    )
    assert autocorrelation_a(Fraction(1, 2)).overlaps(
        arb("0.77220925599087313986130250668208049449 +/- 1e-38")
    )
    assert autocorrelation_a(Fraction(3, 2)).overlaps(
        arb("1.32331052783820652898200776140810746568 +/- 1e-38")
    )


@given(
    st.integers(min_value=1, max_value=12),
    st.integers(min_value=1, max_value=12),
)
@settings(max_examples=30, deadline=None)
def test_reciprocity_and_gram_symmetry_are_rigorous(a: int, b: int) -> None:
    ratio = Fraction(a, b)
    assert autocorrelation_reciprocity_residual(ratio).contains(0)
    forward = natural_gram_entry(a, b)
    reverse = natural_gram_entry(b, a)
    assert arb_to_dyadic(forward) == arb_to_dyadic(reverse)


def test_gram_and_target_normalization_factors() -> None:
    a_half = autocorrelation_a(Fraction(1, 2))
    assert natural_gram_entry(1, 2).overlaps(a_half)
    assert natural_gram_entry(2, 4).overlaps(a_half / 2)
    assert natural_target_entry(1).overlaps(1 - arb.const_euler())
    assert natural_target_entry(2).overlaps(
        (arb(2).log() + 1 - arb.const_euler()) / 2
    )


def test_build_system_mirrors_identical_enclosures() -> None:
    system = build_natural_system((3, 1, 4, 2))
    assert system.dilates == (3, 1, 4, 2)
    for row in range(4):
        for column in range(4):
            assert arb_to_dyadic(system.gram[row][column]) == arb_to_dyadic(
                system.gram[column][row]
            )


def test_direct_energy_has_minus_two_and_full_quadratic_factor() -> None:
    system = NaturalSystem(
        dilates=(1,),
        gram=((arb(3),),),
        target=(arb(2),),
    )
    # 1 - 2*2*(1/2) + 3*(1/2)^2 = -1/4.
    energy = evaluate_natural_distance(system, (Fraction(1, 2),))
    assert energy.is_exact()
    assert energy == arb(-1) / 4


def test_augmented_matrix_has_the_exact_schur_quadratic_identity() -> None:
    system = NaturalSystem(
        dilates=(1, 2),
        gram=((arb(3), arb(1)), (arb(1), arb(2))),
        target=(arb(2), arb(-1)),
    )
    lower = Fraction(1, 8)
    matrix = build_augmented_lower_matrix(system, lower)
    assert matrix[0][2] == -system.target[0]
    assert matrix[1][2] == -system.target[1]
    coefficients = (Fraction(1, 2), Fraction(-1, 4))
    energy = evaluate_natural_distance(system, coefficients)
    vector = tuple(
        arb(value.numerator) / value.denominator for value in coefficients
    ) + (arb(1),)
    quadratic = arb(0)
    for row in range(3):
        for column in range(3):
            quadratic += vector[row] * matrix[row][column] * vector[column]
    assert quadratic == energy - arb(lower.numerator) / lower.denominator


def test_fixed_order_ldlt_implements_the_schur_lower_gate() -> None:
    # For G=[2], b=[1], min E = 1 - b^2/G = 1/2.
    system = NaturalSystem(
        dilates=(1,),
        gram=((arb(2),),),
        target=(arb(1),),
    )
    below = fixed_order_interval_ldlt(
        build_augmented_lower_matrix(system, Fraction(49, 100))
    )
    above = fixed_order_interval_ldlt(
        build_augmented_lower_matrix(system, Fraction(51, 100))
    )
    assert below["classification"] == "POSITIVE_DEFINITE"
    assert above["classification"] == "INCONCLUSIVE"
    assert above["failed_pivot_index"] == "1"


def test_ldlt_rejects_overlapping_but_nonidentical_symmetry() -> None:
    matrix = (
        (arb(3), arb("1 +/- 0.1")),
        (arb("1.05 +/- 0.1"), arb(3)),
    )
    with pytest.raises(ValueError, match="identical symmetric"):
        fixed_order_interval_ldlt(matrix)


@pytest.mark.parametrize(
    ("value", "exponent", "expected"),
    [
        (arb(5) / 4, 1, 2),
        (arb(7) / 4, 1, 4),
        (-arb(5) / 4, 1, -2),
        (-arb(7) / 4, 1, -4),
    ],
)
def test_dyadic_rounding_is_exact_ties_to_even(
    value: object, exponent: int, expected: int
) -> None:
    assert round_arb_to_dyadic_numerator(value, exponent) == expected


def test_dyadic_rounding_fails_closed_on_an_ambiguous_enclosure() -> None:
    assert round_arb_to_dyadic_numerator(arb("1.24 +/- 0.001"), 1) == 2
    with pytest.raises(ArithmeticError, match="unique dyadic rounding"):
        round_arb_to_dyadic_numerator(arb("1.25 +/- 0.01"), 1)


def test_exact_coefficient_validation_rejects_float_and_nondyadic() -> None:
    system = NaturalSystem(dilates=(1,), gram=((arb(1),),), target=(arb(0),))
    with pytest.raises(TypeError, match="exact Fraction"):
        evaluate_natural_distance(system, (0.5,))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="dyadic"):
        evaluate_natural_distance(system, (Fraction(1, 3),))
    with pytest.raises(TypeError, match="exact Fraction"):
        autocorrelation_a(0.5)  # type: ignore[arg-type]


def test_proposal_is_frozen_before_direct_upper_certification() -> None:
    system = build_natural_system(N8_FROZEN_DILATES)
    proposed = propose_dyadic_coefficients(
        system, N8_FROZEN_DENOMINATOR_EXPONENT
    )
    assert proposed == N8_FROZEN_NUMERATORS
    assert round_arb_vector_to_dyadic(
        [arb(value) for value in (1, -2, 3)], 0
    ) == (1, -2, 3)


def test_n8_calibration_and_frozen_finite_bracket() -> None:
    system = build_natural_system(N8_FROZEN_DILATES)
    coefficients = dyadic_coefficients(
        N8_FROZEN_NUMERATORS, N8_FROZEN_DENOMINATOR_EXPONENT
    )
    energy = evaluate_natural_distance(system, coefficients)
    assert energy.overlaps(
        arb("0.0241614215858966850228089430133681017331 +/- 1e-40")
    )

    upper = certify_dyadic_upper_bound(
        N8_FROZEN_DILATES,
        N8_FROZEN_NUMERATORS,
        N8_FROZEN_DENOMINATOR_EXPONENT,
        N8_FROZEN_UPPER_BOUND,
        precision_bits=192,
    )
    lower = certify_augmented_lower_bound(
        N8_FROZEN_DILATES,
        N8_FROZEN_LOWER_BOUND,
        precision_bits=192,
    )
    bracket = certify_frozen_n8_bracket(precision_bits=192)
    assert upper["decision"] == "UPPER_BOUND_CERTIFIED"
    assert upper["checks"]["approximate_solve_used_as_evidence"] is False
    assert lower["decision"] == "LOWER_BOUND_CERTIFIED"
    assert lower["augmented_ldlt"]["classification"] == "POSITIVE_DEFINITE"
    assert bracket["decision"] == "FINITE_BRACKET_CERTIFIED"
    assert bracket["classification"] == "EXPLORATORY"
    assert bracket["hypothesis_status"] == "UNRESOLVED"


def test_inconclusive_bounds_fail_closed() -> None:
    upper = certify_dyadic_upper_bound(
        (1,),
        (0,),
        0,
        Fraction(1, 2),
        precision_bits=128,
    )
    lower = certify_augmented_lower_bound(
        (1,), Fraction(1, 2), precision_bits=128
    )
    assert upper["decision"] == "INCONCLUSIVE"
    # The one-dilate minimum is approximately 0.858, so L=1/2 is valid.
    assert lower["decision"] == "LOWER_BOUND_CERTIFIED"


@pytest.mark.parametrize(
    "operation",
    [
        lambda: vasyunin_sum(2, 4),
        lambda: vasyunin_sum(True, 1),
        lambda: build_natural_system(()),
        lambda: build_natural_system((1, 1)),
        lambda: dyadic_coefficients((1,), -1),
        lambda: certify_augmented_lower_bound((1,), Fraction(-1, 8)),
        lambda: certify_dyadic_upper_bound((1,), (0,), 0, Fraction(1)),
    ],
)
def test_invalid_exact_inputs_are_rejected(operation: object) -> None:
    with pytest.raises((TypeError, ValueError)):
        operation()  # type: ignore[operator]


def test_certificate_precision_context_is_restored() -> None:
    previous_precision = ctx.prec
    certify_dyadic_upper_bound(
        (1,), (0,), 0, Fraction(3, 4), precision_bits=128
    )
    assert ctx.prec == previous_precision


def test_supplied_certificate_system_is_strictly_bound() -> None:
    system = build_natural_system((1, 2))
    prefix = prefix_natural_system(system, 1)
    assert natural_system_content_sha256(prefix)
    prefix_lower = certify_augmented_lower_bound(
        (1,), Fraction(1, 8), precision_bits=192, system=prefix
    )
    assert prefix_lower["kernel_provenance"]["source_dimension"] == "2"
    assert prefix_lower["kernel_provenance"]["supplied_system_reused"] is True
    upper = certify_dyadic_upper_bound(
        (1, 2),
        (0, 0),
        0,
        Fraction(3, 4),
        precision_bits=192,
        system=system,
    )
    assert upper["decision"] == "INCONCLUSIVE"
    with pytest.raises(ValueError, match="dilates"):
        certify_augmented_lower_bound(
            (1,), Fraction(1, 8), precision_bits=192, system=system
        )
    with pytest.raises(ValueError, match="precision"):
        certify_augmented_lower_bound(
            (1, 2), Fraction(1, 8), precision_bits=256, system=system
        )


@pytest.mark.parametrize("mutation", ["manual", "gram", "target", "dilates"])
def test_supplied_system_mutations_cannot_preserve_provenance(mutation: str) -> None:
    canonical = build_natural_system((1, 2))
    if mutation == "manual":
        forged = NaturalSystem(
            canonical.dilates,
            canonical.gram,
            canonical.target,
            canonical.precision_bits,
        )
    elif mutation == "gram":
        forged = replace(
            canonical,
            gram=(
                (canonical.gram[0][0] + 1, canonical.gram[0][1]),
                (canonical.gram[1][0], canonical.gram[1][1]),
            ),
        )
    elif mutation == "target":
        forged = replace(
            canonical,
            target=(canonical.target[0] + 1, canonical.target[1]),
        )
    else:
        forged = replace(canonical, dilates=(2, 1))

    with pytest.raises(ValueError, match="provenance"):
        certify_augmented_lower_bound(
            forged.dilates,
            Fraction(1, 8),
            precision_bits=192,
            system=forged,
        )


def test_same_precision_artifact_is_stable_after_higher_precision_work() -> None:
    arguments = ((1, 2), (1, 1), 2, Fraction(3, 4))
    first = certify_dyadic_upper_bound(*arguments, precision_bits=128)
    certify_dyadic_upper_bound(*arguments, precision_bits=512)
    second = certify_dyadic_upper_bound(*arguments, precision_bits=128)
    assert first == second
