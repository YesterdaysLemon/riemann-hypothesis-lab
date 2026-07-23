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
  moment estimate, they prove the expected leading scale `C_0/log N`.

These are conditional analyses, not unconditional certificate families. A
proof that any explicit family has error `epsilon_N -> 0` for all sufficiently
large `N` would itself prove RH.

## Prioritized theorem targets

### 1. Explicit all-N upper family

Start from the Bettin--Conrey--Farmer log taper and derive a completely
explicit main-term-plus-tail inequality. The immediate objective is not to
assume their zero estimate, but to isolate one exact residue or tail lemma
whose unconditional proof would imply

\[
E(c^{(N)})\le \frac{K}{(\log N)^\alpha}
\]

for some `alpha>0`. This single statement would prove RH through the exact
Mellin identity above.

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
`0<alpha<1` is another compatible direction. The five certified extensions
from `N=8` through `N=256` satisfy the stronger `beta=1` comparison, but the
forced-rebound theorem proves that pattern must eventually fail. Finite v1
coefficients may suggest an ansatz; fitting them is not evidence for a
uniform inequality.

### 3. Infinite-tail floor for a disproof

Let `S_N=span(rho_1,...,rho_N)`, `r_N=chi-P_N chi`, and
`T_N=closure((I-P_N)span{rho_a:a>N})`. A valid counterexample bridge would
need an explicit `delta>0` with

\[
\|P_{T_N}r_N\|^2\le d_N^2-\delta.
\]

Then the limiting squared distance would be at least `delta`. A positive
finite lower bound alone says nothing about this tail and cannot disprove RH.
The exact identity is

\[
d_\infty^2=d_N^2-\|P_{T_N}r_N\|^2.
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

## Primary sources

- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann Hypothesis](https://arxiv.org/abs/math/0202141).
- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann hypothesis, 2](https://arxiv.org/abs/math/0205003).
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058).
- Balazard and de Roton,
  [Sur un critere de Baez-Duarte pour l'hypothese de Riemann](https://arxiv.org/abs/0812.1689).
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191).
