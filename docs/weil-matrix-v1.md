# Finite Weil matrix certificate v1

## Status and scope

This experiment is `CERTIFIED_FINITE`. It does **not** prove or disprove the
Riemann Hypothesis.

The frozen parameters are

\[
  c=\lambda^2=\frac52,\qquad L=\log c,\qquad -4\le m,n\le4.
\]

On `L^2([0,L])`, with multiplicative Haar measure represented by `dx=du/u`,
the orthonormal basis is

\[
  U_n(x)=L^{-1/2}\exp(2\pi i n x/L).
\]

We use `Q` with the sign convention for which RH corresponds to `Q(g)>=0` on
every admissible autocorrelation. The finite matrix is

\[
  A=P-R-S.
\]

This convention and the component formulas follow Bombieri's discussion of
Weil's quadratic functional and the finite-cutoff normalization of
Connes--Consani--Moscovici (CCM), especially CCM equations (2.6),
(2.8)--(2.10), (3.5), (3.10)--(3.18), and (4.2)--(4.4):

- [Bombieri, *Remarks on Weil's quadratic functional in the theory of prime
  numbers, I*](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0)
- [Connes, Consani, and Moscovici, *Zeta Spectral
  Triples*](https://arxiv.org/abs/2511.22755v1)
- [Connes and Consani, *Spectral Triples and
  Zeta-Cycles*](https://arxiv.org/abs/2106.01715v1), Proposition 2.1, Lemma 2.2,
  Proposition 2.3, and Corollary 2.4

The last source proves that the semilocal form is lower-bounded and lower
semicontinuous, that the Laurent polynomials in
`U(u)=u^(i*pi/log(lambda))` form a core, and that the lower bound is the limit
of the finite-compression minimum eigenvalues. Our basis satisfies
`V_n(u)=(-1)^n L^(-1/2) U(u)^n`, so every vector in the frozen matrix is already
in that core. A genuinely negative vector for the complete `A=P-R-S` matrix
would therefore be admissible without a separate smoothing argument.
Conversely, for one fixed `lambda`, proving every compression `E_N` positive
would imply positivity on the full semilocal form domain by the core and
lower-semicontinuity results. Establishing that for a cofinal unbounded family
of cutoffs would cover every compactly supported Weil test function and hence
reach the global criterion. This certificate proves only the single case
`N=4,c=5/2`.

## Frozen formulas

For `0<=y<=L`, the real symmetric correlation kernel is

\[
q_{mn}(y)=
\begin{cases}
\dfrac{\sin(2\pi m y/L)-\sin(2\pi n y/L)}{\pi(n-m)},&m\ne n,\\[6pt]
2(1-y/L)\cos(2\pi n y/L),&m=n.
\end{cases}
\]

It is extended evenly to `[-L,L]` and is zero outside that interval. In
particular, `q_mn(L)=0` exactly.

The pole term is

\[
P_{mn}=\frac{32L\sinh^2(L/4)(L^2-16\pi^2mn)}
{(L^2+16\pi^2m^2)(L^2+16\pi^2n^2)}.
\]

Writing `q0=2` on the diagonal and `q0=0` off it, the archimedean term is

\[
R_{mn}=\frac{q_0}{2}
\left[\gamma+\log\!\left(\frac{4\pi(c-1)}{c+1}\right)\right]
+\int_0^L\frac{e^{y/2}q_{mn}(y)-q_0}{e^y-e^{-y}}\,dy.
\]

The implementation evaluates this integral directly after analytically
removing the apparent singularity at zero. With `t=y/L`, `sinc(z)=sin(z)/z`,
and `a_n=2*pi*n/L`, the limits at `t=0` are finite. Acb integrates over four
equal subintervals of `[0,1]`, with both tolerances set to
`2^(-(precision_bits-16))`. Every result must have radius below
`2^(-(precision_bits-32))`.

An independently coded closed form using digamma, trigamma, and Lerch Phi is
also evaluated for every entry. All 45 upper-triangle integral enclosures must
overlap the corresponding closed-form enclosures. The pole formula is likewise
checked against its defining integral.

The complete prime-power term is

\[
S_{mn}=\sum_{p^r\le c}(\log p)p^{-r/2}q_{mn}(r\log p).
\]

There is no extra factor of two. The cutoff comparison uses exact integer
cross-multiplication. For `c=5/2`, the entire transcript consists of the single
term `p=2,r=1`. If a prime power lies exactly at `c`, it remains in the
transcript while its kernel value is set to the theorem-level exact value
`q_mn(L)=0`; no rounded logarithm decides the endpoint.

## Certified result

At 192-bit working precision, fixed-order interval `LDL^T` elimination proves
all nine pivots positive. The smallest pivot enclosure is

```text
[0.0004498731584304568131774895632183770911656 +/- 4.11e-44]
```

This is a Sylvester positive-definiteness certificate for the frozen `9x9`
ball matrix. As an algorithmically different same-library regression, the
minimum among FLINT's nine separated Rump eigenvalue enclosures is

```text
[2.306064307831340039096603592485617390355e-5 +/- 4.40e-45]
```

The Rump calculation is secondary; interval `LDL^T` is the certificate.

The artifact is [results/weil-matrix-c5-over-2-n4.json](../results/weil-matrix-c5-over-2-n4.json).
It records exact-dyadic balls for every upper-triangle `P`, `R`, `S`, and `A`
entry, all pivots and eigenvalue enclosures, the prime transcript, backend and
quadrature settings, and a canonical payload hash.

## Negative mutation control

The checker deliberately deletes the sole `p^r=2` contribution and evaluates
the resulting corrupted matrix `P-R` on the exact integer vector

```text
[-73201, -108782, -200275, -1000000, 0,
 1000000, 200275, 108782, 73201]
```

Its exact squared norm is `2114603971100`, and its interval Rayleigh quotient
is

```text
[-0.1908453972823702042194412035495257468903 +/- 1.75e-41]
```

The upper endpoint is strictly negative, so the mutation is classified
`CONTROL_NEGATIVE`. This demonstrates that the pipeline detects a known omitted
prime term. The mutated matrix is not the zeta Weil form and has no implication
for RH.

## Replay and limitations

Reproduce the artifact and its 384-bit consistency replay on PowerShell:

```powershell
.\.venv\Scripts\rh-lab.exe weil --bits 192 `
  --output results\weil-matrix-c5-over-2-n4.json
.\.venv\Scripts\rh-lab.exe verify-weil `
  --artifact results\weil-matrix-c5-over-2-n4.json --bits 384
```

The verifier canonically regenerates every field at 192 bits, recomputes at 384
bits, and requires every higher-precision component to remain inside its stored
enclosure. Both passes use FLINT/Arb and the same repository code, so this is a
same-backend consistency replay, not an independent mathematical reproduction.

Most importantly, positivity on this one finite subspace leaves every higher
mode and every other cutoff untested. The published core theorem supplies the
domain and finite-compression limit bridge; it does not turn any one finite
positive matrix into a universal result. A larger positive matrix would still
be finite evidence. A negative vector in a complete frozen matrix would be a
candidate disproof certificate and would still undergo independent sign,
normalization, source-equation, and backend audits before promotion.
