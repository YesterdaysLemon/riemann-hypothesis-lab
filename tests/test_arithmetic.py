from __future__ import annotations

import pytest
from hypothesis import given, strategies as st

from riemann_lab.arithmetic import divisor_sums_up_to


def _brute_sigma(n: int) -> int:
    return sum(divisor for divisor in range(1, n + 1) if n % divisor == 0)


def test_known_divisor_sums() -> None:
    assert divisor_sums_up_to(12) == [0, 1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28]


@given(st.integers(min_value=1, max_value=500))
def test_linear_sieve_matches_brute_force(n: int) -> None:
    assert divisor_sums_up_to(n)[n] == _brute_sigma(n)


def test_limit_validation() -> None:
    for bad in (0, -1):
        try:
            divisor_sums_up_to(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid limit accepted")


@pytest.mark.slow
def test_release_range_matches_independent_divisor_addition_sieve() -> None:
    limit = 1_000_000
    reference = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            reference[multiple] += divisor
    assert divisor_sums_up_to(limit) == reference
