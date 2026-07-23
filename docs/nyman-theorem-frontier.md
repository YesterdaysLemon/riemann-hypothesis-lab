# Nyman theorem frontier after v1

## Status

This is a theorem-level research roadmap, not a result. The finite v1 brackets
remain `EXPLORATORY`, the global hypothesis status remains `UNRESOLVED`, and
every proposed bridge below is explicitly unproved.

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
so a finite scaled value below `C_0` is not contradictory.

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

can be compatible with the obstruction only for `0 < alpha <= 1`.

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

Use the exact autocorrelation kernel to seek an explicit extension
`c^(k) -> c^(k+1)` satisfying

\[
E(c^{(k+1)})
\le
\left(1-\frac{\alpha}{k+\beta}\right)E(c^{(k)}),
\qquad 0<\alpha\le1.
\]

Unlike the impossible fixed-factor contraction, iteration would give a
logarithmic-scale upper family. Finite v1 coefficients may suggest an ansatz,
but fitting those coefficients is not evidence for the uniform inequality.

### 3. Infinite-tail floor for a disproof

Let `S_N=span(rho_1,...,rho_N)`, `r_N=chi-P_N chi`, and
`T_N=closure((I-P_N)span{rho_a:a>N})`. A valid counterexample bridge would
need an explicit `delta>0` with

\[
\|P_{T_N}r_N\|^2\le d_N^2-\delta.
\]

Then the limiting distance would be at least `delta`. A positive finite lower
bound alone says nothing about this tail and cannot disprove RH.

## Role of the v1 data

The six certified midpoint displays give the finite diagnostics

```text
N       8        16       32       64       128      256
d_N^2 log N
        .0502423 .0496128 .0487003 .0473116 .0468636 .0456574
```

They are consistent with a logarithmic scale, but do not establish it. In
particular, the `N=256` scaled value being about 1.2 percent below `C_0` does
not conflict with an asymptotic liminf theorem.

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
