# Literature map

This is a curated map of primary or official sources, not a claim that every
paper's approach is equally promising.

## Status and canonical formulation

- Clay Mathematics Institute, [Millennium Prize Problems](https://www.claymath.org/millennium-problems/):
  RH remains listed among the unsolved problems as of 2026-07-22.
- Enrico Bombieri, [official problem description](https://www.claymath.org/wp-content/uploads/2022/05/riemann.pdf):
  canonical analytic formulation, explicit formula, and context.
- NIST DLMF, [Riemann zeta zeros](https://dlmf.nist.gov/25.10): definitions,
  symmetry, Hardy Z, and the Riemann-Siegel formula.

## Rigorous frontiers

- Platt and Trudgian,
  [The Riemann hypothesis is true up to `3*10^12`](https://arxiv.org/abs/2004.09765):
  the lowest 12,363,153,437,138 positive-ordinate zeros are on the critical
  line through height 3,000,175,332,800. This finite theorem does not prove RH.
- Pratt, Robles, Zaharescu, and Zeindler,
  [More than five-twelfths of the zeros of zeta are on the critical line](https://doi.org/10.1007/s40687-019-0199-8):
  unconditional lower proportion at least 0.417293962.
- Guth and Maynard,
  [New large value estimates for Dirichlet polynomials](https://arxiv.org/abs/2405.20552):
  a strong zero-density frontier, still compatible with exceptional zeros.
- Mossinghoff, Trudgian, and Yang,
  [Explicit zero-free regions for the Riemann zeta-function](https://arxiv.org/abs/2212.06867):
  explicit exclusions near `Re(s)=1`, not the whole critical strip.
- Bellotti, Trudgian, and Yang,
  [Zero-free regions inspired by work of Heath-Brown](https://arxiv.org/abs/2603.21490):
  a 2026 preprint proving the classical-shaped region
  `sigma >= 1 - 1/(4.896 log t)` for `t >= 3`; still only a neighborhood of
  `Re(s)=1`.

## Equivalent criteria used here

- Bombieri,
  [Remarks on Weil's quadratic functional in the theory of prime numbers, I](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0).
- Connes and Consani,
  [Spectral Triples and Zeta-Cycles](https://arxiv.org/abs/2106.01715): the
  lower-semicontinuous semilocal Weil form, its Laurent-polynomial core, and
  the finite-compression lower-bound limit.
- Lagarias, [An Elementary Problem Equivalent to the Riemann Hypothesis](https://arxiv.org/abs/math/0008177).
- Li, [The positivity of a sequence of numbers and the Riemann hypothesis](https://doi.org/10.1006/jnth.1997.2137).
- Bombieri and Lagarias,
  [Complements to Li's criterion](https://doi.org/10.1006/jnth.1999.2392).
- Baez-Duarte,
  [A strengthening of the Nyman-Beurling criterion](https://arxiv.org/abs/math/0202141).
- Baez-Duarte, Balazard, Landreau, and Saias,
  [Sur l'autocorrelation multiplicative de la fonction "partie fractionnaire"](https://arxiv.org/abs/math/0306251):
  rational autocorrelation, reciprocity, and the Vasyunin-sum formula used by
  the frozen natural-distance experiment.
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058):
  the multiplicity-sensitive continuum lower bound used, through natural-space
  inclusion, in the forced dyadic rebound theorem.
- Burnol,
  [On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte](https://arxiv.org/abs/math/0202166):
  an RH-conditional critical-line estimate for `zeta(s)/zeta(s+A)` and an
  important converse guardrail. Its Theorem 4.2 says that square
  integrability of the associated `f_epsilon` for a sequence
  `epsilon -> 0` implies RH.
- Chen and Qi,
  [The best bounds of harmonic sequence](https://arxiv.org/abs/math/0306233):
  the harmonic-number enclosure used for the independent `1-gamma`
  normalization audit.
- Balazard and de Roton,
  [Sur un critere de Baez-Duarte pour l'hypothese de Riemann](https://arxiv.org/abs/0812.1689):
  an RH-conditional explicit natural-dilate coefficient family and upper rate.
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191):
  a direct statement of the unconditional natural-distance lower bound, plus
  the log-tapered Mobius family and its optimal asymptotic under RH and a
  zero-derivative negative-moment assumption
  `sum_(|Im rho|<=T) 1/|zeta'(rho)|^2 << T^(3/2-delta)` for some
  `delta>0`. This assumption in particular excludes multiple zeros.
- Ehm,
  [On certain Gram matrices and their associated series](https://arxiv.org/abs/2405.06349):
  exact Gram-kernel and quadratic-form decompositions for the `q=1` and `q=2`
  Nyman--Beurling problems, including `S_1(x)=O(x^-2)`. The paper constructs
  modified Levinson--Selberg coefficients, but explicitly sets aside
  estimation of its Mobius-inversion error `E_a^(q)(N)` as a major challenge;
  the kernel decay alone does not prove natural-distance convergence.
- Rodgers and Tao,
  [The de Bruijn-Newman constant is non-negative](https://arxiv.org/abs/1801.05914).
- Griffin, Ono, Rolen, Thorner, Tripp, and Wagner,
  [Jensen Polynomials for the Riemann Xi Function](https://arxiv.org/abs/1910.01227).

Every criterion above retains an infinite or universal quantifier. Finite
positivity, finite approximation, and finite inequality checks do not prove RH.

## Audited non-bridges and claimed proofs

The following papers were checked specifically against the natural-dilate
Schur target used in this repository.  Being posted on arXiv is not itself a
correctness certificate.

- Wong,
  [Inequality and Nyman--Beurling--Baez-Duarte criteria](https://arxiv.org/abs/2310.03972),
  claims a proof of RH, but its central finite inequality step is false.  Its
  system `S(epsilon,n)` asks for `||Aa-c||_infinity <= epsilon`, whereas
  Theorem 2.5 only approximates the projection `AA^+c`.  At `n=3` its own
  matrix is

  \[
  A=\begin{pmatrix}
  1&1\\0&2\\1&0\\0&1\\1&2
  \end{pmatrix}.
  \]

  Rows two and four require both `|2a_3-1|<=epsilon` and
  `|a_3-1|<=epsilon`, which is impossible for `epsilon<1/3`.  The later norm
  claim `||P_n||_infinity <= ||P_n||_2=1` also fails for this exact example:
  one row of `P_3=A(A^T A)^-1 A^T` has absolute row sum `10/7`.
- Carvill,
  [Beurling Nyman Geometry and Gram Matrix Structure](https://arxiv.org/abs/2510.18132),
  studies a smoothed sparse `theta=2^-j 3^-k` ladder, not the residualized natural
  block `(S,t)` needed here.  Gram compressibility alone cannot lower-bound
  `t^T S^-1 t`; even `S=I` can have `t=0`.  The posted proof also uses the
  false uniform separation `|a log 2+b log 3| >= (|a|+|b|)log 2`: the pair
  `(a,b)=(8,-5)` already gives `|log(256/243)|`, far below `13 log 2`.
  Therefore this preprint supplies no valid bridge to the repository's
  recurrence target.
- Alouges, Darses, and Hillion,
  [Polynomial approximations in a generalized Nyman--Beurling criterion](https://arxiv.org/abs/2006.02953),
  obtain an unconditional generalized approximation component, but move the
  remaining difficulty into coefficient and Gram control rather than closing
  it.
- Corvalan,
  [Interpolation and Extrapolation Statements equivalent to the Riemann Hypothesis](https://arxiv.org/abs/2312.00211),
  revised in June 2026, gives equivalences and sufficient conditions, not an
  unconditional natural-distance decay theorem.
- Manzur, Noor, and Quintero,
  [A Hardy space approximation supporting zero-free half-planes for the zeta-function](https://arxiv.org/abs/2606.16097),
  proves a zero-free implication for approximation in shifted Hardy spaces
  and studies the critical range numerically; it does not establish the
  `alpha=1/2` closure required for RH.
- Pyvovarov,
  [A few remarks on the Baez-Duarte Criterion](https://arxiv.org/abs/2607.12084),
  revised on 2026-07-21, reduces an exponentially damped Moebius route to an
  explicit global bilinear cancellation.  Its abstract explicitly leaves
  that boundedness problem unresolved.

No source in this audit supplies the uniform
`Gamma_k >= E_k/(k+2)` inequality or an equivalent unconditional
natural-dilate decay.  Proving such an estimate would itself settle RH, so a
reformulation without the decisive norm or correlation bound is not a proof.

## Exact BCF coefficient audit

For the Bettin--Conrey--Farmer taper, use the repository's sign convention

\[
a_n=\mu(n)\left(1-\frac{\log n}{\log N}\right),\qquad c_n=-a_n,
\]

and set

\[
s_N=\sum_{n\le N}\frac{a_n}{n},\quad
h_N(k)=\sum_{\substack{d\mid k\\d\le N}}a_d,\quad
H_N(m)=\sum_{k\le m}h_N(k),\quad q_N(m)=1-H_N(m).
\]

The residual in the equivalent `t`-space is `ts_N` on `0<t<1` and
`ts_N+q_N(m)` on `m<=t<m+1`. Therefore its squared norm has the exact,
unconditional split

\[
\mathcal E_N^{\rm BCF}=\mathcal C_N+\mathcal T_N,
\]

\[
\mathcal C_N=s_N^2+\sum_{m=1}^{N-1}
\left[s_N^2+2s_Nq_N(m)\log\frac{m+1}{m}
+\frac{q_N(m)^2}{m(m+1)}\right],
\]

\[
\mathcal T_N=\sum_{m=N}^{\infty}
\left[s_N^2+2s_Nq_N(m)\log\frac{m+1}{m}
+\frac{q_N(m)^2}{m(m+1)}\right].
\]

Inside the core, the full divisor identity gives
`h_N(k)=1_(k=1)+Lambda(k)/log N` for `k<=N`, hence
`q_N(m)=-psi(m)/log N` for `m<N`. Beyond `N` the divisor truncation means
that this `psi` formula is no longer available. The three terms in each
bracket represent a single nonnegative interval integral and should not be
estimated independently in a way that loses their cancellation. The full
scaled core formula and normalization checks are recorded in
[the theorem-frontier note](nyman-theorem-frontier.md#exact-bcf-coretail-split).

The next proposed statement is only an `EXPLORATORY`, `UNRESOLVED` target:

\[
\mathcal T_N\le A\mathcal C_N+\frac{B}{\log N}
\]

for fixed constants `A,B` and all sufficiently large `N`. It would compare
the truncated-divisor tail with the prime-error core, but it would not by
itself prove RH. One would also need `mathcal C_N=o(1)`. That core estimate is
zero-sensitive at the scale of `psi(x)-x` and must be treated as RH-strength,
not as a routine unconditional prime-number-theorem bound. No cited source is
claimed to prove either missing estimate. Ehm's inversion error is the closest
audited quadratic-form analogue of the tail obstruction, while Burnol's 2002
theorem explains why seemingly technical `L^2` control can already encode RH.

## Current reproducible directions

- Connes,
  [The Riemann Hypothesis: Past, Present and a Letter Through Time](https://arxiv.org/abs/2602.04022):
  2026 finite approximants using small primes; the unproved convergence bridge
  is explicitly the missing step.
- Connes, Consani, and Moscovici,
  [Zeta zeros and prolate wave operators](https://arxiv.org/abs/2310.18423) and
  [Zeta Spectral Triples](https://arxiv.org/abs/2511.22755): finite-prime
  spectral structures whose missing convergence theorem remains explicit.
- Tao, Trudgian, and Yang,
  [New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach](https://arxiv.org/abs/2501.16779)
  and [ANTEDB source](https://github.com/teorth/expdb): a model for executable
  provenance and machine-checkable bound propagation.
- Groskin,
  [High-Precision Approximation of Riemann Zeros via the Truncated Weil Form](https://arxiv.org/abs/2605.20224)
  and [A finite Guinand-Weil dictionary and archimedean tail order](https://arxiv.org/abs/2607.02828):
  current unreviewed preprints with public artifacts directly relevant to the
  finite-form track. They are exploratory inputs, not validation or an RH claim.
- FLINT,
  [rigorous zeta-zero and Turing APIs](https://flintlib.org/doc/acb_dirichlet.html):
  the trusted backend for baseline finite certificates.
