# Forced dyadic rebound after Nyman v1

## Classification and scope

This note proves an unconditional consequence of the certified `N=256`
natural-dilate upper bound and published Nyman--Beurling lower-bound theory:
the scaled dyadic distances must eventually rise above their certified
`N=256` level.

The conclusion is deliberately split into distinct trust layers:

| Layer | What is established |
|---|---|
| finite numerical input | The v1 certificate machine-checks a strict upper bound for `d_256^2` |
| scalar separator | Arb ball arithmetic checks `U_256 log(256) < 23/500 < C_0`; combined with the finite gate this gives the displayed `d_256` inequality |
| infinite bridge | Published theorems plus the manual argument below prove the universal asymptotic lower bound; the bridge is human-auditable, not checker-proved |
| global status | RH remains `UNRESOLVED`; the theorem proves neither RH nor its negation |

The [finite rebound audit](../results/nyman-forced-rebound-v1.json) has
canonical payload SHA-256
`4f5ecef5798cd273774aa4e2d05e7ac3d6b31cf241cadfe632742ea879bff825`.
It deterministically regenerates the exact rational and 256-bit Arb
preconditions. Its own dependency manifest marks the analytic inference as
not machine-checked: the infinite bridge has not been formally proved by the
repository's Python checker, and no independent theorem-prover claim is made.
The immutable source evidence is published in the
[Nyman v1 release](https://github.com/YesterdaysLemon/riemann-hypothesis-lab/releases/tag/nyman-natural-v1).

## Definitions

Let

\[
\rho_a(x)=\left\{\frac1{ax}\right\},\qquad
\chi=\mathbf 1_{(0,1]},
\]

and let

\[
d_N^2=\inf_{c_1,\ldots,c_N\in\mathbb R}
\left\|\chi-\sum_{a=1}^N c_a\rho_a\right\|_2^2.
\]

For dyadic indices `k >= 0`, define

\[
E_k=d_{2^k}^2,\qquad
F_k=E_k\log(2^k)=k\log(2)\,E_k.
\]

The spaces are nested, so the raw distances satisfy

\[
E_{k+1}\le E_k.
\]

There is no corresponding monotonicity theorem for `F_k`. The rebound in this
note concerns `F_k`, not the nonincreasing raw distances `E_k`.

## Universal asymptotic floor

Let

\[
S_N=\operatorname{span}(\rho_1,\ldots,\rho_N),\qquad
S=\overline{\bigcup_N S_N}.
\]

Because the `S_N` are nested,

\[
d_N\downarrow d_\infty:=\operatorname{dist}(\chi,S).
\]

Baez-Duarte's strong natural-dilate Nyman--Beurling theorem gives

\[
\mathrm{RH}\quad\Longleftrightarrow\quad d_\infty=0.
\]

We now split on the truth value of RH.

### If RH is false

Then `d_infinity > 0`, and consequently

\[
d_N^2\log N\longrightarrow+\infty.
\]

### If RH is true

Let `D(lambda)` denote the continuum Nyman--Beurling distance in Burnol's
lower-bound theorem. At `lambda=1/N`, every natural generator with `a <= N`
belongs to the continuum family, so its approximation space contains `S_N`
and

\[
d_N\ge D(1/N).
\]

Burnol's lower bound therefore transfers in the required direction to the
natural distances (and is also stated in the Dirichlet-polynomial
normalization by Bettin, Conrey, and Farmer):

\[
\liminf_{N\to\infty}d_N^2\log N
\ge
\sum_{\Re\rho=1/2}\frac{m(\rho)^2}{|\rho|^2}.
\]

The sum is over distinct nontrivial zeros, including both positive and
negative ordinates; `m(rho)` denotes multiplicity.

Under RH every nontrivial zero lies on the critical line. Since
`m(rho)^2 >= m(rho)`, where `m(rho)` is the multiplicity,

\[
\sum_{\Re\rho=1/2}\frac{m(\rho)^2}{|\rho|^2}
\ge
\sum_{\Re\rho=1/2}\frac{m(\rho)}{|\rho|^2}
=2+\gamma-\log(4\pi).
\]

Define

\[
C_0=2+\gamma-\log(4\pi)
=0.0461914179322420676286204958\ldots.
\]

### Unconditional conclusion

Both exhaustive cases imply

\[
\boxed{\displaystyle
\liminf_{N\to\infty}d_N^2\log N\ge C_0.}
\]

This use of the RH/false-RH dichotomy is classical logic, not an assumption
that RH is true. In the false-RH branch the left side is infinite.

## Exact `N=256` separator

The v1 upper certificate stores

\[
U_{256}=
\frac{
2801788463381838392697853210438663771
}{2^{128}}
\]

and its direct-energy gate is strict:

\[
E_8=d_{256}^2<U_{256}.
\]

At 256-bit precision, Arb gives the sign-separated enclosures

\[
\frac{23}{500}-U_{256}\log(256)
=
[0.000342582215087952442497295688064985\ldots
\ \mathord{+/-}\ 2.24\mathord{\cdot}10^{-78}]
>0
\]

and

\[
C_0-\frac{23}{500}
=
[0.000191417932242067628620495812990583\ldots
\ \mathord{+/-}\ 5.35\mathord{\cdot}10^{-76}]
>0.
\]

The scalar comparison can be replayed without regenerating a Nyman kernel:

```powershell
.\.venv\Scripts\python.exe -c "from flint import arb,ctx; ctx.prec=256; u=arb(2801788463381838392697853210438663771)/(arb(2)**128); q=arb(23)/500; c0=arb(2)+arb.const_euler()-(arb(4)*arb.pi()).log(); print(q-u*arb(256).log()); print(c0-q)"
```

The canonical audit also verifies all five exact dyadic scaled declines,
source bindings, and separator signs:

```powershell
.\.venv\Scripts\rh-lab.exe verify-nyman-rebound-audit `
  --artifact results\nyman-forced-rebound-v1.json `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1
```

Therefore

\[
\boxed{\displaystyle
F_8=d_{256}^2\log(256)<\frac{23}{500}<C_0.}
\]

The first strict inequality is a machine-certified finite numerical
statement. The passage from `C_0` to the behavior of every sufficiently large
index is the published-theorem and human-audit layer.

## Forced rebound theorem

Since

\[
\liminf_{k\to\infty}F_k
\ge C_0>\frac{23}{500},
\]

there is an index `K` such that

\[
F_k>\frac{23}{500}
\qquad(k\ge K).
\]

Combining this with `F_8 < 23/500` proves:

\[
\boxed{\displaystyle
\text{Every sufficiently large dyadic scaled distance exceeds }F_8.}
\]

It also proves that an adjacent scaled increase must occur. If
`F_(k+1) <= F_k` held at every step from `k=8` onward, the sequence could
never cross from below `23/500` to above it. Hence

\[
\boxed{\displaystyle
\exists\,j\ge8:\quad F_{j+1}>F_j.}
\]

Equivalently,

\[
(j+1)d_{2^{j+1}}^2>j\,d_{2^j}^2.
\]

No effective value of `K` or `j` follows from the liminf theorem. The next
dyadic value need not rise, and the scaled sequence may fluctuate before its
eventual lower floor takes effect.

## Consequence for harmonic contraction recurrences

First consider a recurrence holding at every dyadic step from the certified
`N=256` anchor:

\[
E_{k+1}\le
\left(1-\frac{\alpha}{k+\beta}\right)E_k,
\qquad k\ge8,
\]

Assume these are genuine positive contractions:
`0 < alpha/(k+beta) < 1` for every `k >= 8`.

### The case `alpha > 1`

Iteration gives `E_k=O(k^-alpha)`, hence

\[
F_k=O(k^{1-\alpha})\longrightarrow0,
\]

contradicting the universal positive floor. Thus `alpha > 1` is impossible.

### The endpoint `alpha = 1`

For `alpha=1`, the recurrence telescopes exactly:

\[
E_k
\le
E_8\prod_{j=8}^{k-1}\frac{j+\beta-1}{j+\beta}
=
E_8\frac{\beta+7}{k+\beta-1}.
\]

It follows that

\[
\limsup_{k\to\infty}F_k
\le
(\beta+7)\log(2)E_8
<
(\beta+7)\log(2)U_{256}.
\]

Define the certified-anchor threshold

\[
\beta_*=
\frac{C_0}{U_{256}\log(2)}-7
=1.09356642109646170749\ldots.
\]

If `beta <= beta_*`, the last strict upper bound is below `C_0`, including at
the endpoint because `E_8 < U_256`. This contradicts the universal liminf.
Therefore every all-steps-from-`k=8` recurrence with

\[
\boxed{\alpha=1,\qquad \beta\le\beta_*}
\]

is ruled out. This anchor-dependent threshold does not exclude the same
formula if it begins only after a later index. By contrast, the
`alpha > 1` argument applies to a recurrence beginning at any fixed eventual
index.

In particular, `beta=1` would say

\[
E_{k+1}\le\frac{k}{k+1}E_k,
\]

which is exactly `F_(k+1) <= F_k`; the forced rebound proves that an
all-steps-from-`k=8` law must fail at some future dyadic block.

The first convenient integer endpoint still compatible with the asymptotic
floor is `beta=2`:

\[
E_{k+1}\le\frac{k+1}{k+2}E_k,
\qquad
E_k-E_{k+1}\ge\frac{E_k}{k+2}.
\]

The stored lower endpoint

\[
L_{256}=
\frac{
2801788463381838392697853210438663515
}{2^{128}}
<E_8
\]

also satisfies `9 log(2) L_256 > C_0`; the positive Arb margin is greater
than `0.00517`. Thus the recurrence's anchor-based ceiling is not forced below
the universal floor. For comparison,

\[
\limsup F_k
\le9\log(2)E_8
<9\log(2)U_{256}
=0.05136459500802605350\ldots>C_0.
\]

That ceiling does not contradict the floor, so `beta=2` remains a viable
proof-side target. Recurrences with `0 < alpha < 1` also remain compatible
with this obstruction; their resulting upper scale is `O(k^(1-alpha))`,
which does not force `F_k` to zero.

The five certified blocks from `N=8` through `N=256` satisfy the stronger
finite comparisons

\[
E_{k+1}<\frac{k}{k+1}E_k,\qquad k=3,4,5,6,7,
\]

as exact rational consequences of their stored dyadic endpoints. The theorem
above proves only that this pattern must eventually break, not where.

## Exact finite block gain

The recurrence target has an exact Schur-complement form. Split a larger Gram
system as

\[
G_M=
\begin{pmatrix}
G_N&C\\
C^\mathsf T&D
\end{pmatrix},
\qquad
b_M=
\begin{pmatrix}
b_N\\ e
\end{pmatrix},
\]

and define

\[
S=D-C^\mathsf T G_N^{-1}C,\qquad
t=e-C^\mathsf T G_N^{-1}b_N.
\]

Positive definiteness of the Gram matrix makes `S` positive definite. Block
elimination gives the exact identity

\[
\boxed{\displaystyle
d_M^2=d_N^2-t^\mathsf T S^{-1}t.}
\]

For `N=2^k` and `M=2^(k+1)`, a harmonic extension theorem is therefore
equivalent to a lower bound on the captured block energy:

\[
t^\mathsf T S^{-1}t
\ge\frac{\alpha}{k+\beta}E_k.
\]

For the compatible `alpha=1, beta=2` target, the required gain is
`E_k/(k+2)`. Finite interval Schur calculations may test an ansatz, but a
proof-side result must establish the inequality uniformly at every future
scale.

## Disproof-side dual witness

The same Hilbert-space framework gives a precise conditional route from an
off-critical-line zero to a positive limiting distance.

Suppose

\[
\zeta(\rho)=0,\qquad
\rho=\sigma+i\tau,\qquad
\frac12<\sigma<1.
\]

The functional-analytic use of Beurling--Nyman approximation to obtain
zero-free regions is developed more generally by Delaunay, Fricain, Mosaki,
and Robert; the exact zeta witness needed here is derived below.

In the complexified Hilbert space define

\[
h_\rho(x)=
\begin{cases}
x^{\bar\rho-1},&0<x\le1,\\[1mm]
\dfrac{x^{-1}}{1-\bar\rho},&x>1.
\end{cases}
\]

For `0 < Re(s) < 1`, the Mellin identity is

\[
\int_0^\infty \rho_a(x)x^{s-1}\,dx
=-\frac{\zeta(s)}{s\,a^s}.
\]

At `s=rho` the right side is zero. On `x>1`,
`rho_a(x)=1/(ax)`, and the tail of the Mellin integral is
`1/(a(1-rho))`. With the convention
`\langle f,g\rangle=\int\bar f g`, the tail chosen in `h_rho` contributes
exactly the same quantity. Hence

\[
\langle h_\rho,\rho_a\rangle=0
\qquad(a\ge1).
\]

The overlap and norm are

\[
\langle h_\rho,\chi\rangle=\frac1\rho,
\qquad
\|h_\rho\|^2=
\frac1{2\sigma-1}+\frac1{|1-\rho|^2}.
\]

Cauchy--Schwarz therefore gives, for every `N` and in the infinite closure,

\[
d_N^2\ge d_\infty^2
\ge
\frac{|\langle h_\rho,\chi\rangle|^2}{\|h_\rho\|^2}
=
\boxed{\displaystyle
\frac{(2\sigma-1)|1-\rho|^2}{|\rho|^4}>0.}
\]

A rigorous argument-principle box lying strictly to the right of the critical
line and containing a zeta zero would make the last constant
interval-computable. Such a box would already disprove RH; the Nyman floor
would be a corollary. At `sigma=1/2`, the first part of `h_rho` is not in
`L^2`, so a critical-line zero produces no positive constant floor. Burnol's
scale-dependent cutoff vectors instead give the logarithmic obstruction used
above.

## Why finite tail blocks cannot disprove RH

Let `P_N` be the orthogonal projector onto `S_N`, set

\[
r_N=\chi-P_N\chi,
\]

and define the residualized infinite tail

\[
T_N=
\overline{(I-P_N)
\operatorname{span}\{\rho_a:a>N\}}.
\]

Then

\[
S=S_N\mathbin{\oplus}T_N
\]

and Pythagoras gives the exact identity

\[
\boxed{\displaystyle
d_\infty^2=d_N^2-\|P_{T_N}r_N\|^2.}
\]

Thus a disproof-side tail certificate needs an upper bound

\[
\|P_{T_N}r_N\|^2\le B_N<L_N,
\]

where `L_N < d_N^2` is a certified finite lower endpoint. It would imply
`d_infinity^2 > L_N-B_N > 0`.

Writing

\[
u_a=(I-P_N)\rho_a,\qquad
q_a=\langle r_N,u_a\rangle,\qquad
K_{ab}=\langle u_a,u_b\rangle,
\]

the required statement is equivalent to the all-tail rank-one domination

\[
\left|\sum_{a>N}c_aq_a\right|^2
\le
B_N\sum_{a,b>N}\bar c_aK_{ab}c_b
\]

for every finitely supported tail vector. Checking only
`N < a <= M` bounds a finite projection from below and cannot control the
omitted infinite tail. An analytic all-index majorant, an exact annihilating
dual vector such as the off-line-zero witness, or equivalent global
information is unavoidable.

## What this theorem does not say

- It does not identify the first rebounding dyadic index.
- It does not say that `d_N` itself rises; `d_N` is nonincreasing.
- It does not prove the scaled sequence converges to `C_0`.
- It does not prove any harmonic contraction recurrence.
- It does not provide an off-line zero or a positive infinite-tail floor.
- It does not prove or disprove RH.

## Primary sources

- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann Hypothesis](https://arxiv.org/abs/math/0202141).
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058).
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191).
- Delaunay, Fricain, Mosaki, and Robert,
  [Zero free regions for Dirichlet series](https://arxiv.org/abs/1101.1199).
