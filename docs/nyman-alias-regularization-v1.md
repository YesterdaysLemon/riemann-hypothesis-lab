# Fixed-shell alias inversion: a sharp-truncation no-go theorem

## Verdict

The raw alias--Mobius extension from the preceding multiscale note cannot
converge in the required Hilbert norm for any nonzero fixed finite shell.
More precisely, if `y` is nonzero and finitely supported, `a=mu*y`, and

\[
S(X)=\sum_{n\le X}\frac{a_n}{n},
\]

then unconditionally

\[
\boxed{S(X)\ne o(X^{-1/2}).}
\]

The exact real-variable identity

\[
\left\|\sum_{n\le X}a_n\{t/n\}-f_y\right\|_2^2
\ge X|S(X)|^2
\]

therefore proves that the sharp sections cannot converge to the compact
target `f_y`. This closes the fixed-shell natural-cutoff route, rather than
merely showing that a known estimate is too weak.

The conclusion is narrow. It does not reject scale-dependent shells,
balanced finite multipliers, Cesaro or Riesz summation, or a genuinely
multiscale construction. It does not prove or disprove the Riemann
Hypothesis. The global status remains `UNRESOLVED`.

## Exact finite alias identity

Work in

\[
\mathcal H=L^2((0,\infty),dt/t^2),
\qquad e_n(t)=\{t/n\}.
\]

Let `y=(y_j)` be any nonzero finitely supported real sequence and define

\[
a_n=(\mu*y)(n)
=\sum_{j\mid n}y_j\mu(n/j),
\qquad
f_y(t)=-\sum_{j\le\lfloor t\rfloor}y_j.
\]

This step target is compactly supported precisely when `sum_j y_j=0`; the
non-little-o theorem itself does not require that extra condition.

Dirichlet inversion gives `1*a=y`. For the sharp section

\[
g_X(t)=\sum_{n\le X}a_ne_n(t)
\]

and every `0<t<=X`, finite rearrangement gives

\[
\begin{aligned}
g_X(t)
&=t\sum_{n\le X}\frac{a_n}{n}
 -\sum_{n\le t}a_n\left\lfloor\frac tn\right\rfloor\\
&=tS(X)-\sum_{j\le t}y_j
=tS(X)+f_y(t).
\end{aligned}
\]

Consequently the contribution of the growing interval `(0,X]` is exactly

\[
\int_0^X|g_X(t)-f_y(t)|^2\frac{dt}{t^2}
=X|S(X)|^2.
\]

This is an equality on `(0,X]`; the full-norm statement is the lower bound
obtained by discarding `t>X`.

The slope also has the exact factorization

\[
S(X)=\sum_j\frac{y_j}{j}
\sum_{k\le X/j}\frac{\mu(k)}k.
\]

The prime number theorem makes `S(X)->0`, which explains the pointwise
convergence for each fixed `t`. The theorem below shows why that pointwise
limit can never be upgraded to full sharp-section norm convergence.

## The non-little-o theorem

Define the finite Dirichlet polynomial

\[
B_y(s)=\sum_jy_jj^{-s}
\]

and the absolutely convergent Dirichlet series, initially for `Re(z)>0`,

\[
D(z)=\sum_{n\ge1}\frac{a_n}{n^{1+z}}
=\frac{B_y(1+z)}{\zeta(1+z)}.
\]

Assume for contradiction that `S(X)=o(X^(-1/2))`, and set

\[
u(x)=x^{1/2}S(x),
\qquad u(x)\longrightarrow0.
\]

Partial summation gives

\[
D(z)=z\int_1^\infty S(x)x^{-z-1}\,dx.
\]

The assumed little-o bound makes the integral locally uniformly convergent
for `Re(z)>-1/2`, so it supplies a holomorphic continuation throughout that
open half-plane.

### Boundary Abelian lemma

For every fixed real `gamma`,

\[
\varepsilon\int_1^\infty
u(x)x^{-1-\varepsilon-i\gamma}\,dx
\longrightarrow0
\qquad(\varepsilon\downarrow0).
\]

To see this, put `x=e^v` and split the resulting integral at a fixed `V`.
The compact part is `O(epsilon)`. On the tail choose `V` so that
`|u(e^v)|<delta`; absolute integration then bounds the tail by `delta`.
Letting `delta` tend to zero proves the claim. Hence, for
`z_0=-1/2+i*gamma`,

\[
\varepsilon D(z_0+\varepsilon)\longrightarrow0.
\]

### Too many critical-line zeros to cancel

The nonzero function `B_y` is a finite exponential polynomial,

\[
B_y(s)=\sum_jy_je^{-s\log j}.
\]

It has only `O(T)` zeros, counted with multiplicity, in a fixed-width
vertical rectangle of height `T`. This follows directly from Jensen's
formula: translate the center to a point where `B_y` is nonzero and use
`|B_y(z)|<=C exp(C'|z|)` on disks of radius `O(T)`.

Conrey proved that a positive proportion of all zeta zeros are simple and on
the critical line. Their number through height `T` is therefore of order
`T log T`. For arbitrarily large `T`, at least one simple critical-line zero

\[
\rho=\frac12+i\gamma
\]

is not a zero of `B_y`. Near `z_0=rho-1`, the quotient

\[
D(z)=\frac{B_y(1+z)}{\zeta(1+z)}
\]

then has a simple pole. It follows that
`epsilon*D(z_0+epsilon)` tends to a nonzero residue, contradicting the
boundary Abelian lemma.

This proves

\[
\boxed{S(X)\ne o(X^{-1/2})}
\]

for every nonzero finite `y`. Combining it with the exact growing-interval
identity proves that `g_X` cannot converge to `f_y` in `H`. Since the prime
number theorem still gives pointwise convergence to `f_y`, the sharp sections
cannot converge in `H` to any other function either.

Baez-Duarte's Proposition 4.1 in *Moebius-convolutions and the Riemann
hypothesis* proves the analogous critical non-little-o statement for his
Mellin-proper function kernels. A finite atomic shell is not literally one of
those kernels; the argument above is the required discrete extension rather
than an unsupported direct invocation of that proposition.

The shell constraints used previously do not evade the theorem:

\[
\sum_jy_j=0\iff B_y(0)=0,
\qquad
\sum_j\frac{y_j}{j}=0\iff B_y(1)=0.
\]

Neither finite condition cancels the order `T log T` supply of nontrivial
critical-line zeros.

## The exact regularized operator

The failed cutoff suggests replacing `mu` by a finite or summable multiplier
`c=(c_k)`. Assume `c_1=1` and either finite support or the explicit sufficient
condition `sum_k |c_k|/sqrt(k)<infinity`. Put

\[
C(s)=\sum_kc_kk^{-s},
\qquad a^{(c)}=y*c,
\qquad g_c=\sum_na_n^{(c)}e_n,
\]

and define the base residual

\[
R_c(t)=\chi(t)+\sum_kc_ke_k(t),
\qquad \chi=\mathbf1_{[1,\infty)}.
\]

There is an exact dilation factorization:

\[
\boxed{
g_c-f_y=\sum_jy_jR_c(t/j).
}
\]

Indeed, `e_k(t/j)=e_(jk)(t)`, so the double sum has coefficients `y*c`,
while the indicator terms equal `-f_y`.

If `y` is supported on `N<j<=2N`, then for `0<t<=2N`,

\[
\boxed{
g_c(t)-f_y(t)=tB_y(1)C(1).
}
\]

Thus a generic shell preserves its exact first-shell cancellation precisely
when

\[
C(1)=\sum_k\frac{c_k}{k}=0.
\]

This harmonic balance is cheap. Given a finite `d` with `d_1=1`, write
`h=sum d_k/k` and set

\[
c=d-2h\,\delta_2.
\]

Then `c_1=1` and `C(1)=0`. On `(0,1)`, the unbalanced residual is exactly
`R_d(t)=t h`, so `|h|<=||R_d||_2`. Any already-convergent family can
therefore be balanced at vanishing norm cost. The hard part is global
critical-line control, not this scalar constraint.

With `C(1)=0`, put `q=1*c` and

\[
d=y*(\delta_1-q),
\qquad D(M)=\sum_{m\le M}d_m.
\]

The exact real-variable norm is

\[
\boxed{
\|g_c-f_y\|_2^2
=\sum_{M\ge1}\frac{|D(M)|^2}{M(M+1)}.
}
\]

Mellin--Plancherel gives the equivalent spectral target

\[
\boxed{
\|g_c-f_y\|_2^2
=\frac1{2\pi}\int_{-\infty}^{\infty}
|B_y(\tfrac12+i\tau)|^2
\frac{|1-\zeta(\tfrac12+i\tau)C(\tfrac12+i\tau)|^2}
{\tfrac14+\tau^2}\,d\tau.
}
\]

These two boxed identities are the correct targets for finite smooth,
Cesaro, Riesz, Abel, or block multipliers.

## A two-point shell exposes the RH barrier

For distinct positive integers `p,q`, take

\[
y=\delta_p-\delta_q,
\qquad B_y(s)=p^{-s}-q^{-s}.
\]

On the critical line,

\[
|p^{-1/2}-q^{-1/2}|
\le |B_y(1/2+i\tau)|
\le p^{-1/2}+q^{-1/2}.
\]

Therefore

\[
|p^{-1/2}-q^{-1/2}|^2\|R_c\|_2^2
\le\|g_c-f_y\|_2^2
\le(p^{-1/2}+q^{-1/2})^2\|R_c\|_2^2.
\]

For this smallest compact pulse, regularized alias convergence is equivalent
to convergence of the base Nyman residual. Existence of a finite balanced
family with `c_1=1` and `||R_c||_2->0` is equivalent to RH: the forward
direction is the natural Nyman criterion, while under RH Baez-Duarte's
smoothed Mobius approximants converge and the balancing argument above adds
the scalar constraint at vanishing cost.

This rules out a universal regularization shortcut. A viable estimate for the
actual ideal shell must exploit special cancellation in its own `B_y`, not a
uniform bound that would already solve the full criterion.

## Audit of candidate summability methods

### Power damping is RH-strength

Baez-Duarte studies

\[
f_{\varepsilon,n}
=\sum_{k\le n}\mu(k)k^{-\varepsilon}\rho_k.
\]

His Corollary 3.1 states that RH is equivalent to Hilbert convergence as
`n->infinity` for every sufficiently small fixed positive `epsilon`. Under RH,
his Proposition 3.1 chooses `epsilon=c/log log n` and obtains the explicit
global upper bound `O((log log n)^(-1/3))`.

The analytically regularized pointwise value can be written without assuming
RH and tends to the target as `epsilon->0`; what is RH-sensitive is its
realization as a Hilbert limit of the natural coefficient sections. Analytic
continuation alone must not be substituted for that membership statement.

Burnol's Mellin--Hardy estimates and later quantitative work refine this same
route. They do not provide an unconditional critical-line norm bound.

### Fixed prime-factor damping fails even after moment cancellation

Let `0<r<1` and

\[
\mu_r(n)=\mu(n)r^{\Omega(n)},
\qquad a^{(r)}=\mu_r*y.
\]

When `supp(y)` is contained in `(N,2N]`, this normalization preserves the
original first-shell coefficients on that interval. Its exact alias kernel is

\[
q_r=1*\mu_r,
\qquad q_r(n)=(1-r)^{\omega(n)}.
\]

The Dirichlet-series and Euler-product identities, initially for `Re(s)>1`,
are

\[
\sum_n\frac{\mu_r(n)}{n^s}
=\prod_p(1-rp^{-s})
=\zeta(s)^{-r}F_r(s),
\]

\[
\sum_n\frac{q_r(n)}{n^s}
=\zeta(s)^{1-r}F_r(s),
\]

where

\[
F_r(s)=\prod_p(1-rp^{-s})(1-p^{-s})^{-r}
\]

is analytic and nonzero for `Re(s)>1/2`. Indeed,

\[
\log F_r(s)
=\sum_p\sum_{k\ge2}\frac{r-r^k}{k}p^{-ks},
\]

which converges absolutely there.

The Selberg--Delange estimate for `mu_r` itself gives

\[
\sum_{n\le x}\mu_r(n)
=O\!\left(\frac{x}{(\log x)^{1+r}}\right).
\]

Partial summation therefore makes the natural-order series
`sum mu_r(n)/n` convergent. Its value is zero: for `s>1` its Dirichlet
series is `zeta(s)^(-r)F_r(s)`, which tends to zero as `s` decreases to one,
and Abel's theorem identifies that limit with the convergent boundary sum.
Thus the fractional-part series has zero harmonic slope. Its floor part is
finite at each fixed `t`, so the pointwise alias defect is legitimately the
cumulative convolution with `q_r`, rather than only a formal Euler-product
calculation.

Let `ell` be the order of the zero of `B_y` at `s=1`; equivalently, the
moments

\[
B_y^{(h)}(1)=\sum_j\frac{y_j(-\log j)^h}{j}
\]

vanish for `0<=h<ell` and not for `h=ell`. Apply Selberg--Delange first to
the multiplicative function `q_r`, then take the finite shifted combination
defined by `y`. Equivalently, factor its Dirichlet series as

\[
B_y(s)\zeta(s)^{1-r}F_r(s)
=\zeta(s)^{1-r-\ell}
\bigl[B_y(s)\zeta(s)^\ell F_r(s)\bigr]
\]

to obtain

\[
\sum_{n\le x}(y*q_r)(n)
=K_{r,y}\frac{x}{(\log x)^{r+\ell}}
\left(1+O(1/\log x)\right).
\]

Here

\[
\boxed{
K_{r,y}
=\frac{F_r(1)B_y^{(\ell)}(1)}
{\ell!\,\Gamma(1-r-\ell)}\ne0.
}
\]

The constant is nonzero because the first nonzero derivative of `B_y` at one
is nonzero, `F_r(1)>0`, and `Gamma(1-r-ell)` is finite and nonzero.

The natural fractional-part sections converge pointwise to a step function
`g_r`. Its exact error relative to `f_y` has interval values

\[
D_r(M)=\sum_{m\le M}\bigl(y-(y*q_r)\bigr)(m)
\sim-K_{r,y}\frac{M}{(\log M)^{r+\ell}}.
\]

Consequently

\[
\|g_r-f_y\|_2^2
=\sum_{M\ge1}\frac{|D_r(M)|^2}{M(M+1)}
\asymp
\sum_M\frac{1}{(\log M)^{2(r+\ell)}},
\]

which diverges. Every fixed `r` therefore has a pointwise limit outside the
Hilbert space for every fixed nonzero finite `y`, even after any finite number
of logarithmic moment cancellations. If its natural finite sections
converged in `H`, an almost-everywhere convergent subsequence would have to
equal that pointwise limit, which is impossible. The constants may depend on
the fixed pair `(r,y)`; no uniform statement as `r->1` or for a
scale-dependent `y_N` is claimed. A coupled limit `r=r_X->1` with an
additional finite cutoff is not covered and returns to the balanced finite
operator above.

### Exponential damping isolates a live bilinear problem

Pyvovarov's July 2026 preprint studies the absolutely convergent family

\[
\sum_{n\ge1}\mu(n)e^{-nu}\gamma_n(t),
\qquad
\gamma_n(t)=\lfloor t/n\rfloor-\lfloor t\rfloor/n.
\]

It proves convergence of the target inner product and reduces the remaining
norm question to explicit global bilinear cancellation. It does not prove
the needed boundedness and does not resolve RH. This is useful confirmation
that damping the coefficients makes the series well-defined without making
the critical global norm automatic.

## Unfrozen finite regularization scouts

For each stored dyadic candidate residual, the preceding note constructs its
ideal compact first-shell projection `f_(y_N)`. A same-backend 128-bit scout
then applied the plateau--Riesz weight

\[
w_n=
\begin{cases}
1,&n\le2N,\\
(RN-n)/((R-2)N),&2N<n<RN,\\
0,&n\ge RN.
\end{cases}
\]

The table reports the full weighted norm
`||g_w-f_(y_N)||^2/||f_(y_N)||^2`. Its denominator is exactly the stored
candidate's ideal local capture, not the old residual energy or the exact
optimum. Values are floating midpoints, not certified intervals.

| `N` | `R=3` | `R=4` | `R=6` | `R=8` |
|---:|---:|---:|---:|---:|
| 8 | 1.7009 | 1.3225 | 0.8386 | 0.7737 |
| 16 | 1.7530 | 1.4216 | 0.9772 | 0.8838 |
| 32 | 1.7646 | 1.4808 | 1.1391 | 0.9759 |
| 64 | 1.7767 | 1.2975 | 1.0269 | 0.9815 |
| 128 | 2.0711 | 1.5284 | 1.1257 | 1.0414 |

The widest taper is best at every tested scale, but its ratio worsens from
about `0.774` to `1.041`; this is not a convergence signal. Truncated dyadic
shell damping looked somewhat better at selected parameters, but its omitted
infinite tail was material and makes those numbers non-promotable.

The coefficient conditioning also deteriorates: `sum_j |y_j|` grows through

```text
4.73, 8.63, 14.72, 26.25, 46.71, 84.77
```

for `N=8,16,32,64,128,256`. Small weight perturbations therefore act on an
increasingly cancellation-heavy shell. The slopes were small but nonmonotone
while the full errors stayed large, confirming that later alias defects, not
only harmonic slope leakage, dominate these finite candidates. These scouts
are deliberately absent from the claim ledger and finite artifact.

## Frozen finite audit

The tracked
[artifact](../results/nyman-alias-sharp-truncation-v1.json) uses the compact
two-point shell

\[
y_9=1,\qquad y_{16}=-1,
\]

so `f_y=-1` on `[9,16)` and vanishes elsewhere. Through `X=4096`, exact
CPython integer and `Fraction` arithmetic verifies:

1. `a=mu*y`;
2. `1*a=y` at every integer through the limit;
3. the direct and factored formulas for `S_X` at every `X`;
4. the cumulative floor-alias identity on every checked integer interval;
5. exact checkpoint values of `S_X` and `X*S_X^2` at all powers of two from
   16 through 4096.

The artifact's internal payload SHA-256 is

```text
c19d949db86c372916c2cbe583c3e2ff7ea96b3399c0937483eef795e8b56575
```

and its whole canonical-JSON SHA-256 is

```text
5cfd2859838d620c7de5375d1d464e35147f16246e89aec09a4e4bf5da416f94
```

Replay with

```powershell
.\.venv\Scripts\rh-lab.exe verify-nyman-alias-sharp-truncation-audit `
  --artifact results\nyman-alias-sharp-truncation-v1.json
```

The successful replay reports
`REPRODUCED_CERTIFIED_FINITE_ALIAS_IDENTITIES`. The artifact explicitly marks
the Mellin continuation, exponential-polynomial zero count, Conrey theorem,
and infinite non-little-o conclusion as human-auditable external steps rather
than checker-proved consequences of a finite prefix.

## Surviving theorem target

The hard cutoff and fixed prime-factor damping are retired. For an ideal
scale-dependent shell `y_N`, the next finite problem is

\[
\mathcal E_{N,K}
=\inf_{\substack{c_1=1,\ C(1)=0\\
\operatorname{supp}c\le K}}
\frac1{2\pi}\int_{-\infty}^{\infty}
|B_{y_N}(1/2+i\tau)|^2
\frac{|1-\zeta(1/2+i\tau)C(1/2+i\tau)|^2}
{1/4+\tau^2}\,d\tau.
\]

This is a constrained finite Gram problem and can use the repository's
existing exact autocorrelation machinery. The analytic target is a
gain-normalized estimate such as

\[
\mathcal E_{N,N^A}\le\theta\,G_{\rm loc}(N)
\]

uniformly over dyadic `N`, with `theta<1` strong enough to close a contraction
argument. The exact two-point comparison above warns that a bound uniform over
all shells would already be RH-equivalent. Any tractable proof must exploit
the special structure of the actual ideal-shell polynomial `B_(y_N)`.

## Primary sources

- Luis Baez-Duarte,
  [*Moebius-convolutions and the Riemann hypothesis*](https://arxiv.org/abs/math/0504402).
- J. B. Conrey,
  [*More than two fifths of the zeros of the Riemann zeta function are on the critical line*](https://doi.org/10.1515/crll.1989.399.1).
- Luis Baez-Duarte,
  [*A strengthening of the Nyman--Beurling criterion for the Riemann hypothesis, 2*](https://arxiv.org/abs/math/0205003).
- Jean-Francois Burnol,
  [*On an analytic estimate in the theory of the Riemann Zeta function and a Theorem of Baez-Duarte*](https://arxiv.org/abs/math/0202166).
- Regis de la Breteche and Gerald Tenenbaum,
  [*Remarks on the Selberg--Delange method*](https://arxiv.org/abs/2010.12929).
- Alexandre Pyvovarov,
  [*A few remarks on the Baez-Duarte Criterion*](https://arxiv.org/abs/2607.12084).
