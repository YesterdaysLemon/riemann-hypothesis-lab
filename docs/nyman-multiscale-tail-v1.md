# Multiscale Nyman corrections and the divisor-alias barrier

## Verdict

This note audits the most direct continuations of the exact dyadic-lift and
local-cancellation identities in the
[eight-column trial-space rejection](nyman-trial-subspace-v1.md). It reaches
three negative conclusions:

1. the theorem in Báez-Duarte's cited 2005 preprint proves that the canonical
   greedy Vasyunin correction built from exactly the same dyadic seeds diverges
   in weighted `L^1`, hence cannot converge in the `L^2` space required by the
   Nyman criterion;
2. an explicit zero-slope triangular block cancels almost all of its intended
   first-shell error, but its later divisor aliases overwhelm that local gain
   at every tested scale;
3. larger structured `q`-lift spaces pass the desired finite gain threshold at
   several smaller scales and then fail at `N=256`.

These results rule out particular algorithms and finite ansatzes. They do not
rule out every multiscale construction, do not prove an infinite-tail lower
bound, and do not prove or disprove the Riemann Hypothesis. The global status
remains `UNRESOLVED`.

## Exact identification with Vasyunin's first correction

Work in the reciprocal Hilbert space

\[
L^2((0,\infty),dt/t^2),
\qquad e_n(t)=\{t/n\},
\qquad \chi(t)=\mathbf 1_{[1,\infty)}(t).
\]

Every seed and greedy approximant below vanishes on `(0,1)`, while
`[1,infinity)` has weighted measure one. Cauchy--Schwarz therefore makes
weighted `L^2` convergence imply weighted `L^1` convergence in this setting.

The dyadic seed from the preceding note is

\[
h_n(t)=2e_{2n}(t)-e_n(t)
=\left\lfloor\frac tn\right\rfloor
-2\left\lfloor\frac t{2n}\right\rfloor.
\]

It is exactly the first seed in Vasyunin's correction, as analyzed by
Báez-Duarte. Define

\[
\varphi_n(t)=\sum_{k\le n}c_kh_k(t),
\qquad
c_n=1-\sum_{k<n}c_kh_k(n),
\qquad c_1=1.
\]

Because `h_n=0` on `[0,n)` and `h_n=1` on `[n,2n)`, induction gives

\[
\varphi_n(t)=\chi(t)
\qquad(0\le t<n+1).
\]

Thus the greedy construction repairs one new unit interval at every step and
converges pointwise to the target. Pointwise repair is not Hilbert-space
convergence.

## Exact coefficient formula and quantitative divergence

The jump of `h_k` at the integer `m` is

\[
h_k(m)-h_k(m-0)=
\begin{cases}
(-1)^{m/k+1},&k\mid m,\\
0,&k\nmid m.
\end{cases}
\]

Since the jumps of `varphi_n` agree with those of `chi` through `m=n`, the
sequence `c` is the Dirichlet inverse of
`epsilon(m)=(-1)^(m+1)`. For `Re(s)>1`,

\[
\sum_{m\ge1}\frac{\epsilon(m)}{m^s}
=(1-2^{1-s})\zeta(s),
\]

and therefore

\[
\frac1{(1-2^{1-s})\zeta(s)}
=\sum_{a\ge0}\frac{2^a}{(2^a)^s}
 \sum_{m\ge1}\frac{\mu(m)}{m^s}.
\]

Writing `n=2^r m` with `m` odd gives the exact formula

\[
\boxed{c_{2^rm}=2^{\max(r-1,0)}\mu(m).}
\]

In particular, `c_(2^r)=2^(r-1)=n/2` for `r>=1`. The binary seed obeys
`h_n(t)=h_1(t/n)` and `h_n^2=h_n`, so

\[
\int_1^\infty h_n(t)\frac{dt}{t^2}
=\frac1n\sum_{j\ge0}
\left(\frac1{2j+1}-\frac1{2j+2}\right)
=\frac{\log2}{n}.
\]

At powers of two, consecutive greedy approximants therefore satisfy the exact
identities

\[
\|\varphi_n-\varphi_{n-1}\|_1=\frac{\log2}{2},
\qquad
\|\varphi_n-\varphi_{n-1}\|_2^2=\frac{n\log2}{4}.
\]

The increments do not tend to zero in `L^1` and grow in `L^2`. Hence the
sequence is not Cauchy in either norm. This is the quantitative core of the
divergence theorem in Báez-Duarte's cited 2005 preprint.

The theorem is exactly relevant but narrowly scoped. It rejects the
one-coefficient-at-a-time rule above. It does not prove that every Vasyunin
seed, batched block correction, regularized construction, or globally
optimized trial space diverges.

## Exact alias--Möbius extension lemma

There is a clean way to cancel every later divisor alias pointwise, but it
exposes rather than solves the norm-convergence problem. Let `y_j` be supported
on `N<j<=2N`, impose

\[
\sum_{j=N+1}^{2N}y_j=0,
\]

and define the infinite arithmetic extension

\[
a_n=(\mu*y)(n)
=\sum_{\substack{j\mid n\\N<j\le2N}}\mu(n/j)y_j.
\]

Dirichlet convolution gives `1*a=y`. Consequently `a_n=0` for `n<=N`,
`a_n=y_n` on the first block, and for every `n>2N`,

\[
\boxed{
a_n=-\sum_{\substack{d\mid n\\N<d<n}}a_d.
}
\]

Thus each later coefficient cancels the accumulated divisor alias exactly.
No limiting argument is used in this identity.

Define

\[
m(x)=\sum_{1\le k\le x}\frac{\mu(k)}k,
\qquad
S_X=\sum_{n\le X}\frac{a_n}{n}.
\]

For integer `X>=2N`, finite rearrangement gives

\[
S_X=\sum_{j=N+1}^{2N}\frac{y_j}{j}m(X/j).
\]

Let `Y_M=sum_(j<=M)y_j`, put
`f_y(t)=-Y_(floor(t))`, and truncate the natural expansion at `X`:

\[
g_X(t)=\sum_{n\le X}a_n\{t/n\}.
\]

For every `0<t<=X`, another finite rearrangement gives the exact identity

\[
g_X(t)=tS_X-Y_{\lfloor t\rfloor}
=tS_X+f_y(t).
\]

The sum-zero condition makes `f_y` a compact step function supported in
`[N+1,2N)`. The prime number theorem implies `m(x)->0`, hence `S_X->0` and
pointwise convergence `g_X(t)->f_y(t)` for every fixed `t`. But on the entire
growing interval `0<t<=X`,

\[
\boxed{
\|g_X-f_y\|_2^2\ge X|S_X|^2.
}
\]

`L^2` convergence of this natural truncation sequence therefore requires at
least `S_X=o(X^-1/2)`, together with separate control for `t>X`. The prime
number theorem supplies only `S_X=o(1)` in the argument used here. Exact
cancellation of every discrete divisor alias is not enough.

The ideal first-shell projection fits this lemma exactly. If

\[
\alpha_m=sL_m+q_mw_m,
\qquad f_m=\alpha_m/w_m=q_m+sH_m,
\]

choose cumulative values `Y_m=-f_m` on `N<m<2N` and the endpoint
`y_(2N)=-Y_(2N-1)`. Then `f_y=f_m` on the first shell and

\[
\langle r_N,f_y\rangle=\|f_y\|^2=G_{\rm loc}.
\]

This construction enforces `sum y_j=0`. It is distinct from the finite
zero-slope block below, whose endpoint instead enforces
`sum y_j/j=0`. The two coincide only under an additional scalar condition.
The associated Dirichlet-series factorization

\[
\sum_{n\ge1}\frac{a_n}{n^z}
=\frac{\sum_jy_jj^{-z}}{\zeta(z)}
\]

is absolutely justified only for `Re(z)>1`. The right-hand side has a
meromorphic continuation, but using this factorization to deduce critical-line
or Hilbert-tail control requires additional analytic estimates not supplied
here.

## The explicit zero-slope triangular block

Let the old optimal residual be `r_N(t)=s t+q_m` on
`I_m=[m,m+1)`. For a new-block vector `y_(N+1),...,y_(2N)`, impose

\[
\sigma=\sum_{j=N+1}^{2N}\frac{y_j}{j}=0.
\]

With

\[
w_m=\frac1{m(m+1)},\quad
L_m=\log\frac{m+1}{m},\quad
H_m=\frac{L_m}{w_m},
\]

choose cumulative coefficients

\[
Y_m=\sum_{j=N+1}^m y_j=-q_m-sH_m
\qquad(N<m<2N),
\]

and choose `y_(2N)` to enforce `sigma=0`. This is the unique intervalwise
minimizer on the first new shell. Its exact local capture is

\[
G_{\rm loc}=\sum_{m=N+1}^{2N-1}
m(m+1)
\left(sL_m+\frac{q_m}{m(m+1)}\right)^2,
\]

and its local remainder is less than `s^2/(4N)`. Beyond `2N`, however, the
divisor aliases

\[
Z_m=\sum_{j=N+1}^{2N}y_j\lfloor m/j\rfloor
\]

add the exact signed tail increment

\[
\Delta T=\sum_{m\ge2N}
\left[2Z_m(sL_m+q_mw_m)+Z_m^2w_m\right].
\]

The following same-backend calculations use the stored exact dyadic old
coefficients and a 192-bit Arb kernel. Write `rhat_N` for that stored-candidate
residual and `Ehat_N=||rhat_N||^2` for its energy. The displayed normalizations
use `Ehat_N`, not the exact optimum `E_N=d_N^2`; the two are separated by less
than the certified `2^-120` bracket width, so this distinction does not change
any displayed digit or exploratory pass/fail label. These are diagnostics, not
interval-certified signs. The fixed-vector gain keeps the old coefficients
unchanged; the one-dimensional gain optimally rescales this triangular vector
and reoptimizes all old coefficients through the residualized Schur formula.

| `N` | Scout threshold | Local capture / `Ehat_N` | Tail increment / `Ehat_N` | Fixed global gain / `Ehat_N` | Best 1D Schur gain / `Ehat_N` |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.200000 | 0.3430863 | 0.4029839 | -0.0598976 | 0.2288541 |
| 16 | 0.166667 | 0.2794177 | 0.4977072 | -0.2182895 | 0.1320096 |
| 32 | 0.142857 | 0.2145690 | 0.4533195 | -0.2387504 | 0.0794496 |
| 64 | 0.125000 | 0.2347528 | 0.5479740 | -0.3132212 | 0.0795756 |
| 128 | 0.111111 | 0.1626832 | 0.3773947 | -0.2147115 | 0.0761680 |
| 256 | 0.100000 | 0.1240180 | 0.3352419 | -0.2112240 | 0.0345434 |

The intended shell cancellation is extremely accurate: its remaining local
energy ranges from roughly `1.3e-5 Ehat_N` at `N=8` to
`1.5e-7 Ehat_N` at `N=256`.
Nevertheless, the later tail is larger than the local capture at every scale.
At `N=256`, even the optimally scaled residualized direction captures only
about 3.45 percent, far below the required 10 percent.

## Möbius dilation does not repair the tail

For the stored dyadic candidate, let

\[
\widehat p_N=\sum_{n\le N}\widehat c_ne_n,
\qquad \widehat r_N=\chi-\widehat p_N,
\qquad T f(t)=f(t/2),
\]

and put

\[
M_N=\sum_{n\le N}\mu(n)h_n,
\qquad
M_\infty=\chi-2T\chi,
\qquad
\tau_N=M_N-M_\infty.
\]

The explicit `2N`-dimensional candidate

\[
g_{2N}=2T\widehat p_N+M_N
\]

obeys the exact error identity

\[
\boxed{\chi-g_{2N}=2T\widehat r_N-\tau_N.}
\]

The identity was numerically cross-checked to better than `1e-51` in the
192-bit scout. Its computed squared-error ratios are

| `N` | 8 | 16 | 32 | 64 | 128 | 256 |
|---:|---:|---:|---:|---:|---:|---:|
| `||chi-g_(2N)||^2/Ehat_N` | 1.5899 | 4.6257 | 9.4541 | 6.1753 | 6.8395 | 9.9079 |

Thus this natural exact decomposition makes the approximation worse at every
tested scale. At `N=256`, `||tau_N||^2/Ehat_N` is about 13.129; the large
tail is doing the dominant work.

## Structured `q`-lift spaces

For `q=2,...,Q`, define a pure lift on the new block by

\[
v_q(j)=
\begin{cases}
q\mu(j/q),&q\mid j,\\
0,&q\nmid j,
\end{cases}
\]

and a linearly tapered companion by multiplying the nonzero value by

\[
1-\frac{j/q}{\lfloor2N/q\rfloor}.
\]

Exact rational rank pruning removes zero and dependent columns before the
Schur solve. With `Q=16`, the retained dimensions are 21 at `N=32` and 30 at
the three larger scales.

| `N` | Scout threshold | `Q=16` pure+tapered gain / `Ehat_N` | Pass? | Correlation-normalized Gram condition estimate |
|---:|---:|---:|:---:|---:|
| 32 | 0.142857 | 0.1480271166 | yes | `2.84e4` |
| 64 | 0.125000 | 0.1293248977 | yes | `9.02e3` |
| 128 | 0.111111 | 0.1270520566 | yes | `3.45e3` |
| 256 | 0.100000 | 0.0886823998 | no | `2.36e3` |

The near-threshold `N=64` result was repeated at 384 bits. These spaces are
substantially ill-conditioned, remain same-backend exploratory calculations,
and fail at the only scale where the unrestricted block has already been
certified to clear the 10-percent target.

## Frozen finite audit and reproduction

The tracked [finite artifact](../results/nyman-vasyunin-greedy-v1.json) checks
the greedy recurrence, the two-adic Möbius coefficient formula, the exact
unit-interval interpolation invariant, and all reduced-rational increment
multipliers through `n=4096`. It contains 4096 coefficient records and twelve
nontrivial power-of-two records. Its internal payload SHA-256 is

```text
62ad0afd58e4d213efaa49e07a3fa984f7297f2e237fad280144315c3a4ad985
```

and its whole canonical-JSON SHA-256 is

```text
120494625d804daa5ce613cffda9266e26b2ed6b9dd580ed1650b9ba7ab96bc1
```

Generation and verification use exact CPython integers and `Fraction`; no
floating-point or interval backend enters this artifact. The verifier rejects
duplicate keys, nonstandard constants, floats, unknown fields, malformed
canonical values, and semantic mutations, then exactly regenerates the full
declared prefix.

On PowerShell, from an installed checkout:

```powershell
.\.venv\Scripts\rh-lab.exe nyman-vasyunin-greedy-audit `
  --limit 4096 --output results\nyman-vasyunin-greedy-v1-new.json

.\.venv\Scripts\rh-lab.exe verify-nyman-vasyunin-greedy-audit `
  --artifact results\nyman-vasyunin-greedy-v1.json
```

The successful replay reports
`REPRODUCED_CERTIFIED_FINITE_VASYUNIN_GREEDY_PREFIX`. This finite program does
not prove the `log(2)` integral identity, the coefficient formula for every
integer, or divergence of the infinite sequence. Those are the explicit
human-auditable steps above, grounded in the theorem in Báez-Duarte's cited
2005 preprint. The
triangular, dilation, and `q`-lift tables are unfrozen same-backend exploratory
scouts and are not part of the finite artifact or claim ledger.

## What remains

The evidence now excludes three tempting shortcuts:

- pointwise greedy repair does not imply norm convergence;
- exact first-shell minimization does not control its divisor tail;
- passing the gain threshold at several smaller scales does not produce a
  stable arithmetic trial family.

A viable route must control a signed quadratic tail or finite-section inverse
without replacing cancellation by absolute values. The separate and
absolute-value estimates tested here do not close that bound. The most credible
surviving targets are a regularized block operator with a proved uniform global
norm, or a genuinely bilinear estimate for the alias form. Neither is proved
here.

## Primary sources

- Báez-Duarte,
  [A divergent Vasyunin correction](https://arxiv.org/abs/math/0506318).
- Werner Ehm,
  [On certain Gram matrices and their associated series](https://arxiv.org/abs/2405.06349).
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191).
- Pyvovarov,
  [A few remarks on the Baez-Duarte Criterion](https://arxiv.org/abs/2607.12084).
