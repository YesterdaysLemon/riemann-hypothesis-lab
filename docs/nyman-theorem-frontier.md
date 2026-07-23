# Nyman theorem frontier after v1

## Status

This is a theorem-level research roadmap. The finite v1 brackets remain
`EXPLORATORY` and the global hypothesis status remains `UNRESOLVED`. One
unconditional consequence is now proved at the published-theorem,
manual-proof layer: the dyadic scaled distances must eventually rebound above
their certified `N=256` level. The
[full proof and its trust boundaries](nyman-forced-rebound-v1.md) are separate
from the still-unproved targets below.

## Exact normalization map

The v1 distance is

\[
d_N^2=\inf_{c_1,\ldots,c_N\in\mathbb R}
\left\|\chi-\sum_{a=1}^N c_a\rho_a\right\|_2^2,
\qquad
\rho_a(x)=\left\{\frac1{ax}\right\}.
\]

Under the isometry `Uf(x)=f(1/x)` from
`L^2((0,infinity),t^-2 dt)` to `L^2((0,infinity),dx)`, the standard functions
`e_a(t)={t/a}` map exactly to `rho_a`, and `1_[1,infinity)` maps to `chi`.
Mellin--Plancherel therefore identifies the same finite problem with

\[
\frac1{2\pi}\int_{-\infty}^{\infty}
\left|1-\zeta(\tfrac12+it)A_N(\tfrac12+it)\right|^2
\frac{dt}{\tfrac14+t^2},
\qquad
A_N(s)=-\sum_{a\le N}c_a a^{-s}.
\]

Thus a Dirichlet coefficient `a_n` corresponds to primal coefficient
`c_n=-a_n`. This sign and the Hilbert-space measure must be translated before
importing any coefficient family from the literature.

The continuum distance `D(lambda)` in the BBLS/Burnol lower-bound papers is
not equal to this natural distance. At `lambda=1/N`, the natural span is a
subspace of the continuum span, so only the direction

\[
d_N\ge D(1/N)
\]

is available.

## Unconditional asymptotic obstruction

Burnol proves

\[
\liminf_{\lambda\to0}
D(\lambda)^2\log(1/\lambda)
\ge
\sum_{\Re\rho=1/2}\frac{m(\rho)^2}{|\rho|^2}.
\]

Here and below the sum is over distinct nontrivial zeros, including both
positive and negative ordinates, and `m(rho)` is the multiplicity.

The subspace inclusion transfers this lower bound to the natural `d_N`.
Combining it with the strong Nyman--Beurling criterion gives the unconditional
dichotomy:

- if RH is false, `d_N` has a positive limit and `d_N^2 log N` diverges;
- if RH is true, every nontrivial zero is on the line and
  `m(rho)^2 >= m(rho)`, giving

\[
\liminf_{N\to\infty}d_N^2\log N
\ge C_0,
\qquad
C_0=2+\gamma-\log(4\pi)
=0.046191417932242\ldots.
\]

This is asymptotic. It supplies no lower bound for an individual finite `N`,
so a finite scaled value below `C_0` is not contradictory. It does, however,
force every sufficiently late value above any fixed smaller separator.

For the dyadic sequence, define

\[
E_k=d_{2^k}^2,\qquad F_k=k\log(2)E_k.
\]

The exact v1 upper endpoint at `N=256` and an Arb scalar replay prove

\[
F_8<U_{256}\log(256)<\frac{23}{500}<C_0.
\]

The liminf bound therefore forces `F_k > 23/500` for every sufficiently large
`k`. In particular, some adjacent future pair has `F_(k+1)>F_k`, even though
the raw distances `E_k` are nonincreasing. The theorem gives no effective
first rebound index.

One useful design consequence is rigorous: an eventual fixed contraction

\[
d_{2N}^2\le qd_N^2,\qquad q<1,
\]

is impossible, because it would force polynomial decay and hence
`d_N^2 log N -> 0`. A viable dyadic theorem must weaken with scale. For
`N_k=2^k`, the model inequality

\[
E_{k+1}\le\left(1-\frac{\alpha}{k+\beta}\right)E_k
\]

with `0 < alpha/(k+beta) < 1` can be compatible with the obstruction only for
`0 < alpha <= 1`.

The certified `N=256` anchor rules out more at the endpoint. If the recurrence
with `alpha=1` holds at every step `k>=8`, exact telescoping gives

\[
E_k\le E_8\frac{\beta+7}{k+\beta-1},
\qquad
\limsup_{k\to\infty}F_k
\le(\beta+7)\log(2)E_8.
\]

Because the finite certificate is strict, `E_8<U_256`, every such recurrence
with

\[
\beta\le
\beta_*:=\frac{C_0}{U_{256}\log(2)}-7
=1.0935664210964617\ldots
\]

contradicts the universal floor, including at the endpoint. Thus `beta=1` is
ruled out as an all-steps-from-`k=8` law. This does not exclude a `beta=1`
recurrence beginning only after a later index. The convenient next integer
endpoint `beta=2`,

\[
E_{k+1}\le\frac{k+1}{k+2}E_k,
\]

remains compatible because its anchor-based ceiling
is not forced below `C_0`: the stored lower endpoint certifies
`9 log(2) E_8 > C_0`, while
`9 log(2) U_256 = 0.051364595008026...`. Recurrences with `0<alpha<1`
likewise remain compatible with this obstruction.

## Known conditional upper families

- Baez-Duarte proved under RH that a smoothed Mobius family
  `c_{a,n}=-mu(a) exp(-c log(a)/loglog(n))` approaches `chi`, with a stated
  logarithmic rate. Raw unsmoothed Mobius partial sums are not known to work
  in this Hilbert norm.
- Balazard and de Roton improved the RH-conditional upper rate using
  `c_{n,N}=-mu(n)n^-epsilon` with a scale-dependent `epsilon`.
- Bettin, Conrey, and Farmer use the log-tapered family
  `c_{n,N}=-mu(n)(1-log(n)/log(N))`. Under RH plus a zero-derivative negative
  moment estimate
  `sum_(|Im rho|<=T) 1/|zeta'(rho)|^2 << T^(3/2-delta)` for some
  `delta>0`, they prove the expected leading scale `C_0/log N`. The moment
  assumption in particular excludes multiple zeros.
- Burnol's 2002 analytic variant proves a conditional uniform estimate for
  `zeta(s)/zeta(s+A)` on the critical line. It also proves that square
  integrability of the associated functions `f_epsilon` along any sequence
  `epsilon -> 0` already implies RH. This is a useful warning: apparently
  technical `L^2` estimates in this setting can contain the whole hypothesis.
- Ehm's 2024 Gram-kernel decomposition makes a related obstruction explicit.
  Its exact quadratic-form formulas isolate a truncated Mobius-inversion error
  `E_a^(q)(N)`; estimating that error is identified as a major challenge and
  is set aside. The paper supplies useful identities and modified
  Levinson--Selberg coefficients, not an unconditional proof that the natural
  distance tends to zero.

These are conditional analyses, not unconditional certificate families. A
proof that any explicit family has error `epsilon_N -> 0` for all sufficiently
large `N` would itself prove RH.

## Exact BCF core/tail split

The following is an exact arithmetic rewriting of the error of the
Bettin--Conrey--Farmer (BCF) log taper. It is an identity, not an asymptotic
estimate. For `N>1`, put

\[
L=\log N,\qquad
a_n=\mu(n)\left(1-\frac{\log n}{L}\right),\qquad
s_N=\sum_{n\le N}\frac{a_n}{n},\qquad
\kappa_N=Ls_N.
\]

The corresponding primal coefficients are `c_n=-a_n`. In the `t`-space
normalization define the residual

\[
r_N(t)=\mathbf 1_{[1,\infty)}(t)+\sum_{n\le N}a_n\{t/n\}.
\]

For integers `k,m>=1`, define the truncated divisor and cumulative sums

\[
h_N(k)=\sum_{\substack{d\mid k\\d\le N}}a_d,\qquad
H_N(m)=\sum_{k\le m}h_N(k),\qquad
q_N(m)=1-H_N(m).
\]

Then `r_N(t)=ts_N` on `0<t<1`, while on every interval
`m<=t<m+1`,

\[
r_N(t)=ts_N+q_N(m).
\]

Consequently the BCF squared error splits exactly as

\[
\|r_N\|_{L^2(t^{-2}dt)}^2=\mathcal C_N+\mathcal T_N,
\]

where

\[
\mathcal C_N=s_N^2+\sum_{m=1}^{N-1}
\left[
s_N^2+2s_Nq_N(m)\log\frac{m+1}{m}
+\frac{q_N(m)^2}{m(m+1)}
\right]
\]

and

\[
\mathcal T_N=\sum_{m=N}^{\infty}
\left[
s_N^2+2s_Nq_N(m)\log\frac{m+1}{m}
+\frac{q_N(m)^2}{m(m+1)}
\right].
\]

Each bracket is the integral of the nonnegative squared residual over one
unit interval; its three displayed terms must not be bounded separately
without preserving their cancellation.

For `k<=N`, the full divisor sum is available and the elementary Mobius
identities give

\[
h_N(k)=\mathbf 1_{k=1}+\frac{\Lambda(k)}{L}.
\]

Thus `H_N(m)=1+\psi(m)/L` and `q_N(m)=-\psi(m)/L` throughout the core
`m<N`. In particular,

\[
L^2\mathcal C_N
=N\kappa_N^2
-2\kappa_N\sum_{m=1}^{N-1}\psi(m)\log\frac{m+1}{m}
+\sum_{m=1}^{N-1}\frac{\psi(m)^2}{m(m+1)}.
\]

For the tail, put `Q_N(m)=Lq_N(m)`. The exact expression is

\[
L^2\mathcal T_N=\sum_{m=N}^{\infty}
\left[
\kappa_N^2+2\kappa_NQ_N(m)\log\frac{m+1}{m}
+\frac{Q_N(m)^2}{m(m+1)}
\right].
\]

Although `Q_N(N)=-psi(N)`, the identity `Q_N(m)=-psi(m)` does not continue
for `m>N`: the divisor sum is truncated at `d<=N`. That truncation is the
arithmetic tail obstruction.

The recommended next lemma is the explicitly `UNRESOLVED` comparison

\[
\boxed{\quad
\mathcal T_N\le A\mathcal C_N+\frac{B}{\log N}
\quad}
\tag{unproved target}
\]

for absolute constants `A,B` and all sufficiently large `N`. Together with
the separate estimate `mathcal C_N=o(1)`, it would give an explicit family
with total error tending to zero and hence prove RH. The comparison alone
does not prove RH.

The core estimate is not a routine consequence of the prime number theorem.
It is a positive quadratic expression at the scale of fluctuations of
`psi(x)-x`; known estimates strong enough to make it vanish are
zero-sensitive and must be treated as RH-strength. This note does not claim a
published equivalence `mathcal C_N -> 0 iff RH`, but it also does not assume
unconditional core decay. Both that decay and the boxed tail comparison
remain `UNRESOLVED`. Ehm's inversion-error term is the closest audited
quadratic-form analogue of the tail obstruction; no theorem in that paper
closes the boxed estimate.

## Prioritized theorem targets

### 1. Explicit all-N upper family

Start from the exact BCF split above and isolate the truncation obstruction.
The first target is the boxed comparison
`mathcal T_N <= A mathcal C_N+B/log N`, potentially by translating
`mathcal T_N` into Ehm's Mobius-inversion error and then controlling it in
dyadic blocks. This is deliberately narrower than assuming the BCF zero
estimate. A complete proof would still need a valid core estimate strong
enough to imply

\[
E(c^{(N)})\le \frac{K}{(\log N)^\alpha}
\]

for some `alpha>0`. Such a total estimate would prove RH through the exact
Mellin identity above. Neither the tail comparison nor the required core
estimate is currently proved here.

### 2. Harmonic dyadic extension theorem

Use the exact autocorrelation kernel to seek a uniform dyadic distance
inequality, perhaps by an explicit extension of each finite optimizer, with
the still-compatible target

\[
E_{k+1}\le\frac{k+1}{k+2}E_k.
\]

Equivalently, the new block must capture at least `E_k/(k+2)` of the residual
energy. If a larger Gram matrix is split into old and new blocks, block
elimination gives the exact Schur identity

\[
E_{k+1}=E_k-t^\mathsf{T}S^{-1}t,
\]

so the target is `t^T S^-1 t >= E_k/(k+2)`. A recurrence with
`0<alpha<1` is another compatible direction.

The [finite N=512 audit](nyman-beta2-n512-v1.md) now closes the `k=8` step
more strongly than this target:

\[
E_9<\frac{449}{500}E_8<\frac9{10}E_8,
\qquad
E_8-E_9>\frac{51}{500}E_8.
\]

This follows from a direct exact-dyadic `N=512` primal witness and the
certified `N=256` lower endpoint. No `N=512` lower bound or interval Schur
solve is needed. It is one finite step, not evidence that the recurrence
persists.

The first canonical low-dimensional attempt has also been settled, in the
negative. The [eight-column arithmetic trial audit](nyman-trial-subspace-v1.md)
uses log-tapered, raw, and quadratic Moebius columns; squarefree and squareful
mass; and divisibility by two and three. At `k=8`, all 265 augmented interval
`LDL^T` pivots certify

\[
F_V>\frac9{10}U_{256},
\qquad
\Gamma_8(V)<\frac1{10}E_8.
\]

Thus this exact span misses the target even though the unrestricted block
passes it. This does not disfavor every low-dimensional or multiscale rule; it
removes one natural ansatz and prevents fitting its approximate gain from being
mistaken for a theorem.

The remaining proof target can be stated on an explicit trial subspace. Let
`P_k` project onto the old span, let `r_k=(I-P_k)chi`, and put
`w_j=(I-P_k)rho_j` for `2^k<j<=2^(k+1)`. Define

\[
(S_k)_{ij}=\langle w_i,w_j\rangle,
\qquad
(t_k)_i=\langle r_k,w_i\rangle.
\]

For an explicit full-column-rank trial matrix `V_k`, set

\[
A_k=V_k^\mathsf T S_kV_k,
\qquad
u_k=V_k^\mathsf Tt_k.
\]

Then `A_k` is positive definite because `S_k` is positive definite.

The best gain available in that subspace is exactly

\[
\Gamma_k(V_k)=u_k^\mathsf T A_k^{-1}u_k.
\]

It is now sufficient to prove, for one explicit arithmetic rule `V_k` and
every `k>=9`,

\[
\boxed{\Gamma_k(V_k)\ge\frac{E_k}{k+2}.}
\]

The inverse-free sufficient form is to exhibit explicit `z_k` satisfying

\[
2u_k^\mathsf Tz_k-z_k^\mathsf TA_kz_k
\ge\frac{E_k}{k+2}.
\]

Together with the certified base, this would give

\[
E_k<\frac{449}{50(k+1)}E_8\longrightarrow0
\qquad(k\ge9),
\]

and hence RH by the strong natural-dilate criterion. No such uniform trial
rule or lower bound is proved here. The earlier five certified extensions
from `N=8` through `N=256` satisfy the stronger `beta=1` comparison, but the
forced-rebound theorem proves that pattern must eventually fail. Finite
coefficients may suggest an ansatz; fitting them is not evidence for a
uniform inequality.

Two exact coordinate calculations sharpen the remaining target. The dyadic
lift

\[
h_n(t)=2\{t/(2n)\}-\{t/n\}
=\lfloor t/n\rfloor-2\lfloor t/(2n)\rfloor
\]

turns old coefficients into an even new-block direction. At `N=256`, the
exploratory same-backend calculation finds that its ideal infinite Moebius
core is strongly aligned with the residual, but the finite truncation tail
cancels about 67.7 percent of the numerator. Separately, cumulative new-block
coordinates exactly minimize each interval's energy on `N+1<=t<2N`, leaving
`s^2 sum(kappa_m)<s^2/(4N)` across those intervals; divisor aliases beyond
`2N` then create a tail-energy increment with no fixed sign. A viable proof
now needs a uniform alignment or tail-stability bound, or a multiscale
construction that cancels those later aliases; none is presently known here.

### 3. Infinite-tail floor for a disproof

Let `S_N=span(rho_1,...,rho_N)`, `r_N=chi-P_N chi`, and
`\mathscr T_N=closure((I-P_N)span{rho_a:a>N})`. A valid counterexample bridge would
need an explicit `delta>0` with

\[
\|P_{\mathscr T_N}r_N\|^2\le d_N^2-\delta.
\]

Then the limiting squared distance would be at least `delta`. A positive
finite lower bound alone says nothing about this tail and cannot disprove RH.
The exact identity is

\[
d_\infty^2=d_N^2-\|P_{\mathscr T_N}r_N\|^2.
\]

Controlling only a finite sub-block of the tail gives the projection from
below, in the wrong direction. A valid certificate must dominate the entire
infinite tail. An off-critical-line zero would provide an exact annihilating
dual witness and hence a positive all-`N` floor, but producing such a zero
would already disprove RH. The construction and its exact constant are
derived in the [forced-rebound proof note](nyman-forced-rebound-v1.md).

## Role of the v1 data

The six certified energy enclosures give the finite diagnostics

```text
N       8        16       32       64       128      256
d_N^2 log N
        .0502423 .0496128 .0487003 .0473116 .0468636 .0456574
```

They are consistent with a logarithmic scale, but do not establish decay to
zero. The `N=256` scaled value being about 1.2 percent below `C_0` does not
conflict with an asymptotic liminf theorem; combined with that theorem, its
strict certified separator forces a later scaled rebound and at least one
adjacent scaled increase. It still supplies neither an effective rebound
index nor a proof or disproof of RH.

The later `N=512` artifact adds a one-sided primal upper certificate rather
than a two-sided distance bracket. Its role is narrower and theorem-directed:
it proves the `k=8` beta=2 transition and moves the uniform recurrence's first
unproved scale to `k=9`. It does not add a seventh point to the two-sided table
above.

## Primary sources

- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann Hypothesis](https://arxiv.org/abs/math/0202141).
- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann hypothesis, 2](https://arxiv.org/abs/math/0205003).
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058).
- Burnol,
  [On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte](https://arxiv.org/abs/math/0202166).
- Balazard and de Roton,
  [Sur un critere de Baez-Duarte pour l'hypothese de Riemann](https://arxiv.org/abs/0812.1689).
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191).
- Ehm,
  [On certain Gram matrices and their associated series](https://arxiv.org/abs/2405.06349).
