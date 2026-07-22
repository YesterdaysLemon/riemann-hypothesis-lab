"""Exact integer arithmetic used by the falsifier tracks."""

from __future__ import annotations


def divisor_sums_up_to(limit: int) -> list[int]:
    """Return ``sigma(n)`` for every ``0 <= n <= limit``.

    The implementation is a linear sieve.  Every returned value is a Python
    integer; no floating-point arithmetic participates in the divisor sums.
    Entry zero is a sentinel with value zero.
    """

    if isinstance(limit, bool) or not isinstance(limit, int):
        raise TypeError("limit must be an integer")
    if limit < 1:
        raise ValueError("limit must be at least 1")

    least_prime = [0] * (limit + 1)
    prime_power = [0] * (limit + 1)
    prime_power_sum = [0] * (limit + 1)
    sigma = [0] * (limit + 1)
    primes: list[int] = []
    sigma[1] = 1

    for n in range(2, limit + 1):
        if least_prime[n] == 0:
            least_prime[n] = n
            prime_power[n] = n
            prime_power_sum[n] = n + 1
            sigma[n] = n + 1
            primes.append(n)

        for prime in primes:
            composite = n * prime
            if composite > limit:
                break
            least_prime[composite] = prime
            if prime == least_prime[n]:
                prime_power[composite] = prime_power[n] * prime
                prime_power_sum[composite] = (
                    prime_power_sum[n] + prime_power[composite]
                )
                coprime_part = n // prime_power[n]
                sigma[composite] = sigma[coprime_part] * prime_power_sum[composite]
                break

            prime_power[composite] = prime
            prime_power_sum[composite] = prime + 1
            sigma[composite] = sigma[n] * (prime + 1)

    return sigma
