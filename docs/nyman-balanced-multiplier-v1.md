# Scale-dependent balanced multipliers: direct gain and support growth

## Verdict

The fixed-shell sharp cutoff is impossible, but finite balanced multipliers
remain a legitimate scale-dependent construction. This audit makes two
important corrections to the first formulation of that route:

1. minimizing the alias error `||g_c-f_y||` is not the same as maximizing the
   improvement of the old Nyman residual; the missing quantity is a signed
   tail correlation;
2. the coefficient support has cutoff at most `2*N*K`, and its exact top
   index is the product of the shell and multiplier endpoints; polynomial
   effective width therefore grows the approximation dimension too quickly
   for an `O(1/log N)` gain recurrence by itself.

The exact interval identities below are unconditional. The numerical tables
are `EXPLORATORY`: they use the stored exact-dyadic finite candidates, a
binary64 NumPy solve, and finite interval cutoffs. The omitted alias-norm tail
is positive, while the omitted direct-gain tail has no fixed sign. None of the
tables proves a limiting estimate, and the Riemann Hypothesis remains
`UNRESOLVED`.

The finite signal is nevertheless useful. A direct-gain optimized multiplier
of fixed width `K=64` has positive computed gain at all six stored scales. At
`N=256`, the multiplier optimized through four million intervals has prefix
gain decreasing from about `7.44%` of the old candidate's certified
energy-bracket midpoint at 500,000 intervals to about `7.08%` at 4,000,000
intervals. Small instances
agree with an independent full Arb Gram evaluation. This keeps a fixed-width
route alive, but does not supply the missing uniform theorem.

**Certified follow-up.** A separately frozen 9-bit rounding of the
`N=256`, `K=16` direction now has its entire signed tail enclosed. The
[complete-tail certificate](nyman-balanced-full-tail-v1.md) proves direct gain
greater than `7/50000` for that one exact vector. The later
[Fourier/Farey large-sieve certificate](nyman-large-sieve-tail-v1.md) proves
the stronger bound `1/5000` from a prefix sixteen times shorter and replaces
the coefficientwise LCM estimate by a general support-energy lemma. The
[six-scale fixed-width suite](nyman-large-sieve-scaling-v1.md) then certifies
positive complete gain for all six stored `K=64` rows under `T=32768N`.
The later [nested-chain audit](nyman-nested-chain-v1.md) certifies two actual
linked updates, then proves a frozen-prefix obstruction: every strictly
append-only harmonic-balanced continuation from a finite seed has a fixed
positive norm floor. Thus the finite gains are real, but this exact recurrence
architecture cannot prove RH. Overlapping or rebased updates remain open, and
RH remains unresolved.

## Setup

Work in

\[
H=L^2((0,\infty),dt/t^2),
\qquad e_n(t)=\{t/n\},
\qquad \chi(t)=\mathbf 1_{[1,\infty)}(t).
\]

Let a nonzero real shell `y` be supported on `N<j<=2N`, with

\[
\sum_j y_j=0,
\qquad
Y(M)=\sum_{j\le M}y_j,
\qquad
f_y(t)=-Y(\lfloor t\rfloor).
\]

Let `c` be supported on `1<=k<=K` and impose

\[
c_1=1,
\qquad
C(1)=\sum_{k\le K}\frac{c_k}{k}=0.
\]

Put

\[
a=y*c,
\qquad
g_c(t)=\sum_n a_ne_n(t),
\qquad
e_c=g_c-f_y.
\]

The exact dilation identity from the preceding
[regularization audit](nyman-alias-regularization-v1.md) gives

\[
e_c(t)=\sum_jy_jR_c(t/j),
\qquad
R_c(t)=\chi(t)+\sum_{k\le K}c_ke_k(t).
\]

Harmonic balance makes `e_c=0` for `0<t<=2N`. The finite convolution also
shows that the coefficient sequence has

\[
\operatorname{supp}a\subset(N,2NK],
\]

so `g_c` belongs to the span of the natural dilates with indices at most
`2NK`.

The second fact must be carried into every contraction argument.

## Exact integer-interval quadratic

Define

\[
F_k(M)=\sum_jy_j\left\lfloor\frac{M}{jk}\right\rfloor.
\]

If

\[
F(M)=\sum_jy_j\left\lfloor\frac{M}{j}\right\rfloor,
\]

then the useful compression

\[
\boxed{F_k(M)=F(\lfloor M/k\rfloor)}
\]

avoids constructing the much larger convolution Gram matrix. On every unit
interval `[M,M+1)`, harmonic balance removes the linear term and gives the
exact constant value

\[
\boxed{
D_c(M)=Y(M)-\sum_{k\le K}c_kF_k(M).
}
\]

Consequently

\[
\boxed{
\|e_c\|^2
=\sum_{M\ge1}\frac{D_c(M)^2}{M(M+1)}.
}
\]

This identity is independent of the Mellin formula and is the primary
time-domain oracle for the scout.

### Automatic parameterization of both constraints

For `3<=k<=K`, introduce real variables `z_k` and set

\[
c_1=1,
\qquad
c_2=-2\left(1+\sum_{k=3}^Kz_k\right),
\qquad
c_k=kz_k.
\]

Then `C(1)=0` identically. Define

\[
D_0(M)=Y(M)-F_1(M)+2F_2(M),
\]

and

\[
D_k(M)=2F_2(M)-kF_k(M)
\qquad(3\le k\le K).
\]

The interval error becomes

\[
D_c(M)=D_0(M)+\sum_{k=3}^Kz_kD_k(M).
\]

With `w_M=1/(M(M+1))`, put

\[
H_{k\ell}=\sum_Mw_MD_k(M)D_\ell(M),
\qquad
b_k=\sum_Mw_MD_0(M)D_k(M),
\qquad
q=\sum_Mw_MD_0(M)^2.
\]

Thus

\[
\|e_c\|^2=q+2b^Tz+z^THz,
\qquad
z_*=-H^{-1}b.
\]

For every nonzero finite `y` and `K>=3`, the full Gram `H` is positive
definite. Indeed, a null vector defines a balanced perturbation `Delta c`
with `Delta c_1=0`. Its zero-norm error would have critical-line Mellin
transform

\[
-\frac{\zeta(s)}s B_y(s)\Delta C(s),
\qquad
B_y(s)=\sum_jy_jj^{-s},
\qquad
\Delta C(s)=\sum_k\Delta c_kk^{-s},
\]

vanishing almost everywhere. The zeta factor is nonzero almost everywhere on
that line, while `B_y` and `Delta C` are finite entire Dirichlet polynomials.
The identity theorem and the fact that nonzero entire functions have isolated
zeros force `Delta C` to vanish identically, hence `Delta c=0`. The constrained
minimizer is therefore unique. For `K=2` there are no variables and the
constraints force

\[
c=(1,-2).
\]

In that case

\[
R_c(t)=\chi(t)+\{t\}-2\{t/2\}
\]

is exactly the indicator of the union of intervals `[2m,2m+1)` for `m>=1`.

## The decisive quantity is signed direct gain

Let `r_N` be the old residual and let `f=f_(y_N)` be its ideal first-shell
projection. Write

\[
G=\|f\|^2=\langle r_N,f\rangle,
\qquad
e=g_c-f,
\qquad
E=\|e\|^2,
\qquad
X=\langle r_N,e\rangle.
\]

Since balance makes `e` vanish through `2N`, while `f` is supported in the
first shell, `f` and `e` are orthogonal. Direct expansion gives

\[
\boxed{
\|r_N\|^2-\|r_N-g_c\|^2=G+2X-E.
}
\]

The alias-quality objective `E` is only one term. A favorable or unfavorable
signed correlation `X` can change the conclusion.

If the complete direction `g_c` is allowed one optimal real rescaling, then

\[
\langle r_N,g_c\rangle=G+X,
\qquad
\|g_c\|^2=G+E,
\]

and its exact one-direction captured gain is

\[
\boxed{
\frac{(G+X)^2}{G+E}.
}
\]

### Direct interval formula

For the stored old coefficient vector, write on `[M,M+1)`

\[
r_N(t)=s_Nt+q_M.
\]

Its weighted interval mass is

\[
R_N(M)
=s_N\log\frac{M+1}{M}
+\frac{q_M}{M(M+1)}.
\]

Therefore

\[
\boxed{
X=\sum_{M\ge2N}D_c(M)R_N(M).
}
\]

The fixed-addition gain is another concave quadratic in the same variables
`z`, with the same negative Hessian `-H`. It is computationally as cheap to
optimize `G+2X-E` directly as it is to minimize `E`.

## Why a fixed alias-error fraction is not enough

Let

\[
E_N=\|r_N\|^2,
\qquad
\gamma=G/E_N,
\qquad
\delta=E/E_N.
\]

Cauchy--Schwarz alone gives the worst-case fixed-addition bound

\[
\frac{G+2X-E}{E_N}
\ge
\gamma-2\sqrt{(1-\gamma)\delta}-\delta.
\]

To retain a fraction `a` of the local gain using only this norm bound, it is
sufficient to require

\[
\delta\le
\left(\sqrt{1-a\gamma}-\sqrt{1-\gamma}\right)^2.
\]

For small `gamma`, the right side is asymptotic to

\[
\frac{(1-a)^2}{4}\gamma^2.
\]

Thus a norm-only theorem needs `E/G=O(gamma)`, not merely `E/G<theta` for one
fixed `theta<1`. At `N=256`, the stored shell has
`gamma=0.1240179648...`; retaining half of its local gain by this worst-case
argument would require approximately `E/G<0.009`. The observed truncated
ratios are far larger. Direct signed optimization is therefore essential.

## The support-growth constraint

A finite balanced correction constructed from `(N,2N]` and `[1,K]` has
coefficient cutoff

\[
N'\le2NK.
\]

More precisely, let `j_+=max supp(y)` and `k_+=max supp(c)`. The top
convolution coefficient is the unique product of the two top indices, so

\[
N'=j_+k_+,
\qquad
Nk_+<N'\le2Nk_+.
\]

Suppose a future theorem yielded a fractional reduction comparable to
`c/log N` at every step. To force the product of residual factors to zero, the
associated series

\[
\sum_r\frac1{\log N_r}
\]

must diverge.

- If `K` is fixed, then `log N_r` grows linearly in `r`, and the series
  diverges.
- If the effective width is bounded above by a fixed power of `log N`, then
  `log N_r=O(r log r)`, and the series still diverges.
- If the effective width `k_+=N^A` for a fixed `A>0`, then `log N_r` grows
  geometrically, and the series converges.

Therefore a polynomial-effective-width, `O(1/log N)` one-shell recurrence
cannot by itself force the Nyman distance to zero. It would need a stronger
gain sequence whose sum still diverges (a constant fractional gain is one
sufficient example), or it must be embedded in a genuinely multi-shell
argument whose dimension accounting is different.

## Exploratory finite results

The frozen
[exploratory grid](../results/nyman-balanced-multiplier-scout-v1.json) was
generated by [generate_nyman_balanced_scout.py](../tools/generate_nyman_balanced_scout.py).
Its payload SHA-256 is

```text
f31c330b373b3e8616fb1837193d8fbac6e31ce56dd4986364342dc4369f3144
```

The artifact labels this as a retrospective exploratory reporting grid, not a
preregistered experiment. Each cell has its own schema, engine identifier,
raw input binding, finite-value gate, and canonical payload hash. The
large-scale vectors remain binary64 values rather than stored exact dyadics;
only the small full-Gram checks below perform exact dyadic rounding.

The underlying [scout_nyman_balanced.py](../tools/scout_nyman_balanced.py) reads
the stored exact-dyadic candidate at each `N`, reconstructs its ideal shell,
builds `F_k(M)` by divisor-style vectorized accumulation, and solves either
the truncated alias-norm or direct-gain quadratic.

The six shell summaries frozen in the artifact are:

| `N` | `sum abs(y_j)` | `B_y(1)` | `G_loc` |
|---:|---:|---:|---:|
| 8 | 4.73396116127 | -0.00246641627669 | 0.00828945337578 |
| 16 | 8.62933782725 | 0.00186619136716 | 0.00499990666433 |
| 32 | 14.7197123536 | -0.000683334159042 | 0.00301511185048 |
| 64 | 26.2493617871 | 0.0000630473524560 | 0.00267055766521 |
| 128 | 46.7065536761 | -0.000163147062309 | 0.00157128330326 |
| 256 | 84.7746416528 | 0.000212769249838 | 0.00102112873590 |

The source candidates are raw-hash-bound, and the complete table is covered by
the artifact payload hash.

### Fixed width `K=64`

The following direct-gain solve uses `M<=1,000,000`. For this exploratory
table, `gamma` is the ideal local gain divided by the midpoint of the stored
old candidate's certified energy bracket. The next two columns use that same
midpoint and the local gain as denominators.

| `N` | `gamma` | direct gain / energy midpoint | retained local gain | alias norm / `G_loc` |
|---:|---:|---:|---:|---:|
| 8 | 0.3430863 | 0.2824231 | 0.8231838 | 0.4067153 |
| 16 | 0.2794177 | 0.1918003 | 0.6864285 | 0.4380473 |
| 32 | 0.2145690 | 0.1239694 | 0.5777599 | 0.4754955 |
| 64 | 0.2347528 | 0.1343526 | 0.5723150 | 0.4897360 |
| 128 | 0.1626832 | 0.0886732 | 0.5450669 | 0.5394769 |
| 256 | 0.1240180 | 0.0724794 | 0.5844269 | 0.5872513 |

This is a finite signal for fixed width, not a uniform lower bound.

At `N=256`, solving the direct-gain problem through four million intervals and
evaluating that one multiplier on successive prefixes gives the convergence
diagnostic

| last interval | alias norm / `G_loc` | direct gain / energy midpoint | retained local gain |
|---:|---:|---:|---:|
| 500,000 | 0.5707998 | 0.0743611 | 0.5995998 |
| 1,000,000 | 0.5865764 | 0.0724298 | 0.5840268 |
| 2,000,000 | 0.5952251 | 0.0713692 | 0.5754747 |
| 4,000,000 | 0.6002071 | 0.0707575 | 0.5705422 |

The direct gain remains positive and appears to stabilize near seven percent,
but the signed tail prevents interpreting the last row as a certified lower
bound.

### Width sweep at `N=256`

The one-million-interval direct-gain sweep is:

| `K` | coefficient cutoff `2*N*K` | direct gain / energy midpoint |
|---:|---:|---:|
| 2 | 1,024 | -0.2692025 |
| 4 | 2,048 | -0.1062973 |
| 8 | 4,096 | 0.0001371 |
| 16 | 8,192 | 0.0319731 |
| 32 | 16,384 | 0.0580626 |
| 64 | 32,768 | 0.0724794 |
| 128 | 65,536 | 0.0842346 |
| 256 | 131,072 | 0.0941347 |

The gain increases with width in this finite truncation, but the support cost
increases at the same time. The `K=256` value also falls from `0.11055` at
125,000 intervals to `0.09413` at one million intervals, so it is especially
unsafe to extrapolate it as a full-tail result.

## Independent full-Gram cross-checks

The helper
[check_nyman_balanced_gram.py](../tools/check_nyman_balanced_gram.py) rounds a
scout shell and multiplier to exact 80-bit dyadics, enforces both constraints
exactly, constructs the complete finite natural-dilate vector, and directly
evaluates the old and new energies with the 192-bit Arb autocorrelation
kernel. Each multiplier in this table was first optimized by the direct-gain
scout through exactly 2,000,000 intervals.

| `N` | `K` | interval result | full Arb Gram result | aggregate dimension |
|---:|---:|---:|---:|---:|
| 8 | 4 | -0.2126716 | -0.2126861273 | 37 |
| 8 | 8 | 0.0307275 | 0.0307131837 | 60 |
| 16 | 16 | 0.0817755 | 0.0817374701 | 213 |

The interval column divides by the stored certified energy-bracket midpoint;
the Arb column divides by a complete Arb evaluation of the exact old
candidate energy. The interval entries use the original binary64 vectors,
while the Arb entries use their exact 80-bit-dyadic roundings. Their agreement
is an independent consistency diagnostic for the interval signs, dilation
indices, balance parameterization,
and observed tail difference. It does not certify the tail direction for the
unrounded vector or any large-scale table.

## Reproduction

Install the optional scout dependency in addition to the base package:

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[scout]"
```

Run the principal large-scale diagnostic with

```powershell
.\.venv\Scripts\python.exe tools\scout_nyman_balanced.py `
  --n 256 --k 64 --interval-limit 4000000 --objective direct-gain
```

Regenerate the complete frozen grid with

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_balanced_scout.py `
  --output results\nyman-balanced-multiplier-scout-v1-new.json
```

For a small full-Gram cross-check, first obtain a multiplier from the scout,
then pass its coefficient list to

```powershell
.\.venv\Scripts\python.exe tools\scout_nyman_balanced.py `
  --n 8 --k 4 --interval-limit 2000000 --objective direct-gain
.\.venv\Scripts\python.exe tools\check_nyman_balanced_gram.py `
  --n 8 --multiplier 1 <remaining-coefficients> `
  --bits 192 --dyadic-bits 80
```

The scout deliberately emits `EXPLORATORY` and `UNRESOLVED`. Its partial alias
norm is a lower diagnostic for the full norm. Its partial direct gain has a
signed omitted tail and is neither a lower nor an upper bound.

## Surviving theorem target

The most credible one-shell target now has fixed or polylogarithmic width.
For example, it would be enough to prove along a cofinal geometric sequence
that one explicitly constructed balanced direction satisfies

\[
\frac{(G_N+X_{N,K})^2}{G_N+E_{N,K}}
\ge \frac{c}{\log N}\,d_N^2
\]

for one fixed `K` and one constant `c>0`, with all terms evaluated over the
full Hilbert space. Since the next dimension is at most `2KN`, the resulting
reciprocal-log gains would have a divergent sum and force the distances to
zero.

The append-only target formerly posed here is now impossible. Every correction
`a=y*c` is supported strictly above the old core and has
`sum_n a_n/n=0`. It consequently vanishes on the entire initial interval, so
the seed residual there is frozen forever. The
[nested-chain proof](nyman-nested-chain-v1.md) shows that no finite seed can
have zero residual on that whole interval: otherwise Möbius inversion forces
`p_n=-mu(n)`, contradicting an exact Bertrand-prime argument for
`sum_(n<=N)mu(n)/n`.

A repaired, non-append-only construction would instead need both:

1. an overlapping or residualized update that reoptimizes old coefficients
   and removes the frozen-core obstruction;
2. uniform scale bounds and a local-gain estimate of order at least
   `c/log N` for that repaired recurrence.

Neither statement is proved here. The finite-support tail itself remains
controlled exactly by `Q(Q-1)(rho+tau)`, but that tool must now be paired with
a recurrence that changes the earlier residual rather than merely appending
balanced outer support.

## Primary sources

- Luis Baez-Duarte,
  [*A strengthening of the Nyman--Beurling criterion for the Riemann
  hypothesis, 2*](https://arxiv.org/abs/math/0205003).
- Jean-Francois Burnol,
  [*On an analytic estimate in the theory of the Riemann Zeta function and a
  Theorem of Baez-Duarte*](https://arxiv.org/abs/math/0202166).
- Sandro Bettin, J. Brian Conrey, and David W. Farmer,
  [*An optimal choice of Dirichlet polynomials for the Nyman--Beurling
  criterion*](https://arxiv.org/abs/1211.5191).
- Alexandre Pyvovarov,
  [*A few remarks on the Baez-Duarte Criterion*](https://arxiv.org/abs/2607.12084).
- H. L. Montgomery and R. C. Vaughan,
  [*The large sieve*](https://doi.org/10.1112/S0025579300004708),
  *Mathematika* 20 (1973), 119-134.
