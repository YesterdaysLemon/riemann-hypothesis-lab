# Adaptive Weil witness search v1

## Classification

This search is `EXPLORATORY`. Its only trusted terminal facts are:

- `FINITE_POSITIVE_CERTIFIED`: interval `LDL^T` closed one scheduled finite
  matrix;
- `NEGATIVE_CANDIDATE_QUARANTINED`: one exact integer vector was strictly
  negative for the complete matrix at both discovery and confirmation
  precision, but still awaits independent audit; or
- `INCONCLUSIVE_MAX_PRECISION`: the explicit precision schedule ended without
  either certificate.

The search never emits `PROVED` or `DISPROVED`. Positive finite cells do not
enter the claim ledger or README headline. A negative candidate halts the run
for independent sign, normalization, source-equation, parity, cutoff, and
backend review.

## Why adaptive precision is necessary

The first sweep remained positive while its smallest reported eigenvalue
collapsed from roughly `9.1e-2` at `c=3/2,N=8` to `6.8e-29` at `c=12,N=12`.
At 128 bits, fixed-order interval elimination became inconclusive for
`c=8,N=8`; the same matrix closed positive at 256 bits. At 384 bits, the
larger cases `c=8,10,12` with `N=16` also closed positive, with minimum Rump
enclosures near `1.33e-27`, `3.66e-31`, and `3.43e-34`, respectively.

This is consistent with the prolate near-null phenomenon described by
[Connes and Consani, Sections 2.2--3](https://arxiv.org/abs/2106.01715v1), not
with a sign change. They report about `2c` small directions and eigenvalues far
below ordinary fixed-precision resolution. The later finite-cutoff analysis of
[Connes, Consani, and Moscovici](https://arxiv.org/abs/2511.22755v1) relates
the conditioning heuristically to an `exp(-4*pi*c)` scale while keeping the
missing rigorous prolate/Weil approximation step explicit. A floating negative
eigenvalue is therefore treated only as a candidate generator.

## Bounded engine-validation grid

Version 1 freezes exact rational cutoffs

```text
3/2, 2, 5/2, 3, 4, 5, 8, 9
```

and degrees

```text
4, 8, 12, 16
```

for 32 ordered cells. Each cell attempts `96`, `192`, then `384` bits. Interval
`LDL^T` is the only positive terminal gate; Rump eigenvalues are diagnostic.
If a candidate integer vector is strictly negative, the exact same vector must
remain negative at the next precision or at the reserved 768-bit confirmation
precision.

Approximate eigenvectors are untrusted hints. Candidate coefficients are
derived deterministically from exact dyadic midpoints, projected into reversal
parity sectors, scaled by powers of two, rounded with ties-to-even, divided by
their gcd, and sign-canonicalized. Only a subsequent full Arb evaluation of
`x^T(P-R-S)x` can establish the sign.

The run is serial because FLINT's global precision context is shared mutable
state. Canonically hashed plans, attempts, terminal cells, and the aggregate
index make interruption and resume auditable. Resume refuses plan, backend, or
hash mismatches.

## Executed v1 result

The complete frozen grid terminated with 32
`FINITE_POSITIVE_CERTIFIED` cells, zero `INCONCLUSIVE_MAX_PRECISION` cells, and
zero `NEGATIVE_CANDIDATE_QUARANTINED` cells. Twenty-three cells closed at 96
bits, six at 192 bits, and three at 384 bits. The 44 total precision attempts
serialized 13,236 upper-triangle entries and independently cross-checked the
same number of pole and archimedean entries against their source-integral
oracles. Candidate extraction triggered 443 exact integer-vector evaluations:
370 were strictly positive and 73 were interval-inconclusive.

The smallest separated Rump enclosure occurred at `c=9,N=16`:

```text
[1.483861410967958711918628391031090118348e-29 +/- 1.47e-69]
```

That spectrum is diagnostic; interval `LDL^T` supplied the terminal
positive-definiteness certificate. The aggregate index is
[`results/weil-search-grid-v1/index.json`](../results/weil-search-grid-v1/index.json),
with payload SHA-256
`f0595f9630b8738d518dfabfeb2d4d9921f4e0fcb0d4d09a9c201fe6f47e7d87`.
A fresh same-backend process replayed all 32 cells. The result remains a bounded
null search and leaves RH `UNRESOLVED`.

## Prime-power transition frontier

The bounded grid validates the engine. The higher-value search uses both sides
of exact prime-power transitions. At `q=p^r`, the new term is present but
exactly zero because `q_mn(L)=0`. The matrix is continuous and its derivative
has a one-sided kink. If `c>q` and `h=log(c/q)`, then

\[
q_{mn}(L-h)=\frac{2h}{L}
+O\!\left((m^2+mn+n^2)\frac{h^3}{L^3}\right),
\]

so the new prime contribution is first-order rank one in the reversal-even
sector and only cubic in the odd sector. This motivates parity-separated
candidate extraction.

For a transition `q`, use the exact log-symmetric rational stencil

\[
c_+(q,j)=q\frac{2^j+1}{2^j},\qquad
c_-(q,j)=q\frac{2^j}{2^j+1},\qquad c_0=q,
\]

with `j` in `4,6,8,10,14,22,38`, excluding points that cross another
prime-power transition. The priority order is:

1. `q=7` as calibration against published figures;
2. `q=8,9,11,13` as the immediate frontier;
3. `q=16,17,19` as the second frontier.

Because the near-null band has dimension about `2c`, degree `N<c` is
structurally under-resolved. The first ladders are `N=16,24,32` at `c=8`,
`N=20,30,40` at `c=10`, and `N=24,36,48` at `c=12`; promising cells extend to
`N=4c`. At fixed cutoff, Rayleigh--Ritz nesting requires the true minimum to be
nonincreasing with `N`. Any certified negative vector must remain negative when
zero-padded into a larger basis.

## Post-quarantine promotion gates

The engine enters `NEGATIVE_CANDIDATE_QUARANTINED` after the first two checks
below. It does not implement checks 3--7. A quarantined finite negative may be
escalated beyond quarantine only after an independent audit establishes all of
the following:

1. a primitive integer vector has an interval upper bound `x^T A x < 0`;
2. the same vector is negative at discovery and confirmation precision;
3. the exact cutoff and independently enumerated complete prime transcript
   agree;
4. direct archimedean quadrature overlaps the special-function evaluator;
5. complex-basis and parity-basis quadratic evaluations agree;
6. the zero-padded vector stays negative at a larger degree; and
7. the fixed vector stays negative on a small exact rational cutoff
   neighborhood, as continuity predicts.

The full mathematical normalization continues to come from
[Bombieri](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0), the semilocal
domain/core bridge from
[Connes--Consani](https://doi.org/10.4171/LEM/1049), and verified numerical
linear algebra from [FLINT](https://flintlib.org/doc/acb_mat.html).
