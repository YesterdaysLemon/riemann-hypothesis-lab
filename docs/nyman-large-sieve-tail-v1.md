# Fourier/Farey large-sieve tail certificate

## Result

For the same explicit exact-dyadic `N=256`, `K=16` balanced Nyman
correction used by the first complete-tail certificate, the complete infinite
weighted-`L^2` direct gain is rigorously greater than

\[
\frac1{5000}=0.0002.
\]

The proof encloses only the first `2^22 = 4,194,304` intervals numerically.
The earlier certificate needed `2^26` intervals and proved the weaker lower
bound `7/50000 = 0.00014`.  The improvement comes from a general theorem:
the periodic tail is controlled by its Fourier energy and Farey spacing,
rather than by summing absolute pairwise LCM bounds.

This is a `CERTIFIED_FINITE` contraction and a reusable finite-support tail
lemma.  It is not an all-scale contraction theorem, and the Riemann
Hypothesis remains `UNRESOLVED`.

The machine-readable certificate is
[nyman-large-sieve-tail-v1.json](../results/nyman-large-sieve-tail-v1.json).

## Numerical headline

| quantity | certified value or enclosure |
|---|---:|
| coefficient support limit `Q` | `8192` |
| reciprocal Farey spacing bound `Q(Q-1)` | `67,100,672` |
| old periodic mean square `rho` | `0.7740050628733907...` |
| added mean square `mu` | `10.584429862109205...` |
| old/added cross mean `nu` | `0.07609141036594362...` |
| new periodic mean square `tau=rho+mu-2nu` | `11.206252104250709...` |
| former pairwise-LCM partial-sum bound | `508,362,537,412.84094...` |
| Fourier/Farey partial-sum bound | `803,883,306.6468434...` |
| rigorous reduction factor | greater than `632` |
| prefix cutoff | `2^22` |
| prefix gain | approximately `0.00025738631` |
| tail center | `-0.000002487240923437...` |
| periodic tail radius | `0.000045695442831005...` |
| logarithmic-correction radius | `8.71e-15` |
| complete gain | `[0.0002092036258016804..., 0.0003005945167313684...]` |
| claimed strict lower bound | `1/5000` |

All proof quantities in the JSON artifact are exact integers, exact rational
numbers, or outward-rounded Arb balls.  The decimals in this note are only
readable summaries.

## 1. Periodic setup

Put

\[
\phi_n(M)=\left\{\frac{M}{n}\right\}-\frac{n-1}{2n}
\]

for integer `M`.  For the finitely supported balanced correction `a`, define

\[
g_M=\sum_n a_n\phi_n(M),
\qquad
\sum_n a_n=\sum_n\frac{a_n}{n}=0.
\]

Let `p` be the old exact Nyman vector and put

\[
C_0=1-\frac12\sum_n p_n,
\qquad
v_M=C_0-\sum_n p_n\phi_n(M).
\]

The periodic part of the direct-gain integrand is

\[
h_M=2g_Mv_M-g_M^2.
\]

The exact interval identity and the small nonperiodic logarithmic correction
are derived in the
[first complete-tail proof](nyman-balanced-full-tail-v1.md).  The only bound
replaced here is the uniform partial-sum bound for the zero-mean part of
`h_M`.

## 2. The Fourier coefficients collapse by denominator

Write `e(x)=exp(2*pi*i*x)`.  The normalized discrete Fourier expansion of one
centered sawtooth is

\[
\phi_n(M)
=\sum_{k=1}^{n-1}
-\frac{e(kM/n)}{n(1-e(-k/n))}.
\]

To verify it, take the DFT over one period and use

\[
\sum_{j=0}^{n-1}j e(-kj/n)
=-\frac{n}{1-e(-k/n)}
\quad (1\le k<n).
\]

If `k/n=r/d` is reduced, then `n=dq` and `k=rq`.  All terms at the same
reduced frequency therefore combine.  With

\[
S_d=\sum_{d\mid n}\frac{a_n}{n},
\]

the coefficient at reduced nonzero frequency `r/d` is exactly

\[
\widehat g(r/d)
=-\frac{S_d}{1-e(-r/d)}.
\]

There is no frequency-zero term because every `phi_n` has mean zero.  More
importantly, every active reduced denominator is at most

\[
Q=\max\bigl(\operatorname{supp}(a)\cup\operatorname{supp}(p)\bigr).
\]

For the old residual, write

\[
P_d=\sum_{d\mid n}\frac{p_n}{n}.
\]

Its zero-frequency coefficient is `C_0`, and its coefficient at `r/d` is

\[
\widehat v(r/d)
=\frac{P_d}{1-e(-r/d)}.
\]

## 3. Parseval recovers the exact Jordan means

The elementary identity

\[
\sum_{\substack{1\le r<d\\(r,d)=1}}
\frac1{|1-e(r/d)|^2}
=\frac{J_2(d)}{12}
\]

follows from the standard finite cosecant-square sum and Möbius inversion.
Parseval then gives

\[
\mu=\operatorname{mean}(g_M^2)
=\frac1{12}\sum_{d\ge2}J_2(d)S_d^2,
\]

\[
\rho=\operatorname{mean}(v_M^2)
=C_0^2+\frac1{12}\sum_{d\ge2}J_2(d)P_d^2,
\]

and

\[
\nu=\operatorname{mean}(g_Mv_M)
=-\frac1{12}\sum_{d\ge2}J_2(d)S_dP_d.
\]

These are the same exact Jordan-`J_2` divisor sums independently implemented
by the certificate code.

## 4. Farey spacing

Take two distinct reduced frequencies with denominators at most `Q`.

- If their denominators `d` and `e` differ, their circular distance is at
  least `1/(de)`, and `de <= Q(Q-1)`.
- If `d=e`, two distinct fractions are at least `1/d` apart.
- Frequency zero is at least `1/d >= 1/Q` from a nonzero reduced fraction of
  denominator `d`.

Thus the union of every frequency needed for `v`, `g`, and `v-g` is separated
modulo one by

\[
\delta\ge\frac1{Q(Q-1)}.
\]

No least common multiple of the coefficient indices is computed or stored.
It exists only as a conceptual common period for the finite rational spectrum.

## 5. Two-sided interval energy from the large sieve

Montgomery and Vaughan's Theorem 1 gives the classical large-sieve
upper bound.  In dual form, if

\[
F(M)=\sum_{\alpha\in\mathcal F}c_\alpha e(\alpha M)
\]

has a `delta`-separated frequency set, then every shifted interval `I` of `H`
consecutive integers satisfies

\[
\sum_{M\in I}|F(M)|^2
\le (H+\delta^{-1})\sum_\alpha|c_\alpha|^2.
\]

The source is H. L. Montgomery and R. C. Vaughan,
[“The large sieve,” Mathematika 20 (1973), 119-134](https://doi.org/10.1112/S0025579300004708).
An openly accessible scan is available from the
[University of Michigan repository](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/152543/mtks0025579300004708.pdf).

For completeness, the one-sided theorem gives exactly the two-sided estimate
needed here.  All frequencies are rational, so choose a common integer period
`P`.  For a residual interval of length at most `P`, apply the upper theorem
once to that interval.  Its cyclic complement is another shifted interval;
applying the same upper theorem to the complement and subtracting from the
exact full-period Parseval identity gives the lower bound.  Longer intervals
are reduced by removing complete periods.  Consequently,

\[
\left|
\sum_{M\in I}|F(M)|^2
-H\operatorname{mean}(|F|^2)
\right|
\le \delta^{-1}\operatorname{mean}(|F|^2).
\]

Using the Farey spacing gives the entirely explicit rational constant

\[
D_Q=Q(Q-1).
\]

## 6. Difference of squares removes the cross-term problem

The key algebraic simplification is

\[
h_M=2g_Mv_M-g_M^2
=v_M^2-(v_M-g_M)^2.
\]

Define

\[
\tau=\operatorname{mean}((v-g)^2)
=\rho+\mu-2\nu\ge0.
\]

Applying the two-sided interval energy estimate separately to `v` and `v-g`
and using the triangle inequality proves, for every integer interval `I`,

\[
\boxed{
\left|
\sum_{M\in I}\bigl(h_M-(2\nu-\mu)\bigr)
\right|
\le Q(Q-1)(\rho+\tau).
}
\]

This is the new periodic partial-sum lemma.  It is unconditional, applies to
every finite balanced correction, includes the constant Fourier mode safely,
and contains no floating-point operation.

For the frozen vector,

\[
C_{\rm LS}=Q(Q-1)(\rho+\tau)
=803{,}883{,}306.6468434\ldots,
\]

while the former absolute pairwise-LCM bound was

\[
C_{\rm LCM}=508{,}362{,}537{,}412.84094\ldots.
\]

The exact ratio is stored in the artifact; `C_LS/C_LCM` is approximately
`0.00158132`.

## 7. Infinite tail

Split

\[
h_M=(2\nu-\mu)+z_M.
\]

The mean term telescopes:

\[
\sum_{M>T}\frac{2\nu-\mu}{M(M+1)}
=\frac{2\nu-\mu}{T+1}.
\]

Every interval partial sum of `z_M` has magnitude at most `C_LS`.  Abel
summation against the decreasing weights `1/(M(M+1))` therefore gives

\[
\left|
\sum_{M>T}\frac{z_M}{M(M+1)}
\right|
\le\frac{C_{\rm LS}}{(T+1)(T+2)}.
\]

The logarithmic interval correction is unchanged from the first certificate.
If

\[
P_1=\sum_n\frac{p_n}{n},
\qquad
A_0=\sum_n|a_n|,
\]

then its omitted tail is at most

\[
\frac{|P_1|A_0}{12T^2}.
\]

Hence the complete direct-gain tail is enclosed by center

\[
\frac{2\nu-\mu}{T+1}
\]

and radius

\[
\boxed{
\frac{Q(Q-1)(\rho+\tau)}{(T+1)(T+2)}
+\frac{|P_1|A_0}{12T^2}.
}
\]

## 8. Independent computation boundary

The finite prefix and analytic tail use different mechanisms.

1. The source candidate, scout cell, shell rounding, multiplier rounding, and
   exact convolution are rebound to frozen hashes.
2. The first `2^22` intervals use exact integer divisor recurrences.
3. NumPy binary64 reductions receive an exact rational forward-error radius.
4. Arb encloses the logarithmic contribution at 256 bits during generation.
5. The replay changes the block partition and raises Arb precision to 384
   bits.
6. Exact Jordan sums recompute `mu`, `nu`, `rho`, and `tau`.
7. Exact rational arithmetic computes the Farey constant and the entire
   omitted-tail radius.
8. The stored generation balls must contain the independent replay balls.

The frozen artifact has canonical payload SHA-256

```text
2ee75c9af3625872c8e7b9a59ea76e353859180d636ec7446cac09db99e25918
```

and raw-file SHA-256

```text
a8fe268d922602d6f66cca128b7b02c275f77d2f7111e0a7153c8add3f517d8f
```

## 9. Reproduction

Generate the artifact:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_large_sieve_tail_certificate.py
```

Replay it with the independent block size and precision:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_large_sieve_tail_certificate.py --verify
```

Run the focused tests:

```powershell
.\.venv\Scripts\python.exe -m pytest `
  tests\test_nyman_balanced_tail.py `
  tests\test_nyman_large_sieve_certificate.py -q
```

## 10. What remains for RH

The lemma eliminates the enormous coefficientwise-LCM loss, but one finite
contraction is not enough.  An iterative proof still needs an explicit family
of corrections with uniform direct gain and controlled periodic energy.

For a sequence with support limit `Q_N`, old periodic energy `rho_N`, new
periodic energy `tau_N`, and desired gain scale `G_N`, the large-sieve tail is
controlled once

\[
\frac{Q_N^2(\rho_N+\tau_N)}{T_N^2}=o(G_N)
\]

together with the corresponding mean and logarithmic-tail conditions.  A
particularly concrete sufficient target for fixed multiplier width is a bound
of the form

\[
\rho_N+\tau_N\le C_KN^2G_N,
\]

paired with a cutoff chosen uniformly on a large enough polynomial scale.
Proving such an all-scale estimate, or finding a different uniform mechanism,
is the next mathematical bottleneck.  Until that is done, the repository's RH
status remains `UNRESOLVED`.
