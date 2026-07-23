# Nyman natural-distance certificates v1

## Classification

This experiment is `EXPLORATORY`. It computes rigorous two-sided bounds for
six finite-dimensional approximation distances. Finite decay does not prove
that the limiting distance is zero, and a positive finite lower bound does not
disprove RH. Every plan, cell, audit, and index keeps the global hypothesis
status `UNRESOLVED`.

## Criterion and normalization

Let

\[
H=L^2((0,\infty),dx),\qquad
\rho_a(x)=\left\{\frac{1}{ax}\right\},\qquad
\chi=\mathbf 1_{(0,1]}.
\]

For a positive integer `N`, define

\[
d_N^2=\inf_{c\in\mathbb R^N}
\left\|\chi-\sum_{a=1}^N c_a\rho_a\right\|_2^2.
\]

Báez-Duarte's strong Nyman--Beurling theorem proves

\[
\mathrm{RH}\quad\Longleftrightarrow\quad\lim_{N\to\infty}d_N=0.
\]

This is the natural-dilate formulation in `L^2(0,infinity)`. It must not be
mixed with the classical `L^2(0,1)` coefficient-constrained formulation.
Endpoint conventions for `chi` are immaterial in `L^2`.

Define the autocorrelation

\[
A(\lambda)=\int_0^\infty
\{t\}\{\lambda t\}\frac{dt}{t^2}.
\]

The substitution `t=1/(ax)` and an elementary harmonic-sum calculation give

\[
G_{ab}=\langle\rho_a,\rho_b\rangle
=\frac1a A\!\left(\frac ab\right),\qquad
b_a=\langle\chi,\rho_a\rangle
=\frac{\log a+1-\gamma}{a}.
\]

Therefore, for real coefficients `c`,

\[
E(c)=\left\|\chi-\sum_{a=1}^N c_a\rho_a\right\|_2^2
=1-2b^\mathsf T c+c^\mathsf T Gc.
\]

The factor `1/a` in `G_ab` is essential. Reciprocity
`A(lambda)=lambda*A(1/lambda)` implies exact symmetry
`A(a/b)/a=A(b/a)/b`.

## Rational autocorrelation formula

Every positive rational is reduced before arithmetic and hashing. For
coprime `p,q>0`, set

\[
V(p,q)=\sum_{k=1}^{q-1}
\left\{\frac{kp}{q}\right\}\cot\!\left(\frac{\pi k}{q}\right),
\qquad V(p,1)=0.
\]

The Báez-Duarte--Balazard--Landreau--Saias formula is

\[
A(p/q)=
\frac{1-p/q}{2}\log(p/q)
+\frac{1+p/q}{2}\bigl(\log(2\pi)-\gamma\bigr)
-\frac{\pi}{2q}\bigl(V(p,q)+V(q,p)\bigr).
\]

All logarithms, cotangents, `pi`, and Euler's constant are evaluated with Arb
balls. Exact rational reduction, summation indices, and metadata use integer
arithmetic.

## Frozen finite experiment

V1 contains the six natural prefixes

```text
N = 8, 16, 32, 64, 128, 256.
```

The frozen arithmetic policy is:

| Quantity | Policy |
|---|---:|
| generation precision | 768 bits |
| same-backend replay precision | 1536 bits |
| coefficient grid | `2^-256` |
| bound grid | `2^-128` |
| target lower slack | `2^-120` |

An approximate solve may propose coefficients, but it is never sign-bearing.
The proposed vector is rounded ties-to-even to the exact dyadic coefficient
grid and serialized as integers plus a common power-of-two denominator. Every
certificate uses only that stored exact vector.

All six cells run. If either side of a bracket cannot be certified at the
frozen precision, the cell is `INCONCLUSIVE`; the implementation must not
widen the lower slack or change the grid after seeing a result.

## Primal upper certificate

For the stored exact dyadic vector `c_N`, Arb directly encloses `E(c_N)`. Let
`U_N` be the least `2^-128` grid point above the certified upper endpoint.
The checker requires

\[
E(c_N)\le U_N,
\]

which proves `d_N^2 <= U_N`. The approximate linear solve is not replayed or
trusted by this gate.

## Dual lower certificate

For an exact dyadic target `L`, form the augmented matrix

\[
K_N(L)=
\begin{pmatrix}
G_N&-b_N\\
-b_N^\mathsf T&1-L
\end{pmatrix}.
\]

A fixed-order interval `LDL^T` factorization with every pivot strictly
positive proves `K_N(L)` positive definite. Since

\[
(c,1)^\mathsf T K_N(L)(c,1)=E(c)-L,
\]

this proves `d_N^2>L`. V1 freezes `L_N=U_N-2^-120`. A failed or
zero-containing pivot is `INCONCLUSIVE`, not evidence that `d_N^2<=L_N`.

A terminal `FINITE_DISTANCE_BRACKET_CERTIFIED` cell therefore proves the
unconditional finite statement

\[
L_N<d_N^2\le U_N.
\]

## Independent truncated-integral oracle

The formula audit does not import or call the Vasyunin implementation. For
`lambda=p/q`, it splits `[0,T]` at the exact rational breakpoints from
`Z union (q/p)Z`. On an interval `[u,v]` where
`m=floor(t)` and `n=floor(pt/q)`, it evaluates

\[
\int_u^v\frac{\{t\}\{pt/q\}}{t^2}\,dt
=\frac pq(v-u)-\left(n+\frac pqm\right)\log(v/u)
+mn(1/u-1/v).
\]

Nonnegativity gives the rigorous generic tail bound

\[
0\le\int_T^\infty
\frac{\{t\}\{pt/q\}}{t^2}\,dt\le\frac1T.
\]

V1 freezes `T=4096` and the ordered inputs

```text
(1,1), (1,2), (2,1), (2,3), (3,2), (4,6), (5,7), (7,5),
(31,32), (32,31), (127,128), (128,127), (251,256), (256,251).
```

The supplied `4/6` input must record canonical reduction to `2/3`. It is not a
numerical sign control: the identity `V(dp,dq)=d*V(p,q)` makes the unreduced
outer denominator cancel. Reduction is checked through exact metadata and
cache keys.

The generic tail is only `1/4096=2^-12`. This is a rigorous truncated-integral
normalization oracle, not an independent high-precision reproduction of the
autocorrelation formula. A separate finite harmonic-sum interval checks the
`1-gamma` constant in `b_a`.

## Structural audits

- Compute both orientations of every `G_ab` and require reciprocity overlap
  before storing a shared symmetric enclosure.
- Bind every smaller `(G_N,b_N)` to the exact prefix of the shared `N=256`
  kernel.
- Zero-pad `c_N` into the `2N` space and reproduce the same primal value.
- Check `U_2N<L_N` separately. Passing proves strict finite improvement;
  failure does not invalidate either bracket and has no implication for RH.
- Reject a reversed Vasyunin sign, `p` replacing the outer `q`, either omitted
  Vasyunin term, a reversed logarithmic term, a missing `1/a`, and corrupted
  `b_a`.
- Reject mutations to coefficient, bound, plan, kernel, pivot, dimension,
  dyadic encoding, payload hash, or classification.

## Completed frozen run

The 2026-07-22 run completed all six cells. Every cell has status
`FINITE_DISTANCE_BRACKET_CERTIFIED`, no cell is inconclusive, and all five
separate comparisons `U_2N<L_N` are certified. The evidence contains 504
stored exact dyadic coefficients and 510 positive augmented-`LDL^T` pivots.

- plan payload SHA-256:
  `27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a`;
- 768-bit generation-kernel payload SHA-256:
  `6fb54e7f63208f1bd81fb7ef2e592a8a655458dd7744081748590159b868cc26`;
- index payload SHA-256:
  `281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23`;
- compact summary payload SHA-256:
  `35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf`;
- independently rebuilt 1536-bit replay-kernel payload SHA-256:
  `0b914bd5902f119ace638e74947225bce3bc541f719afee7a313d8390d68dae8`;
- normalization-audit payload SHA-256:
  `021a060fb3eb0121c325a0af13cc39095c0c6f26bc5a4b235269611a66a9442c`.

The tracked [summary](../results/nyman-natural-v1-summary.json) gives the exact
dyadic endpoints and binds every cell, certificate, candidate, and kernel
prefix. The separate verifier rebuilt the full kernel at 1536 bits and
re-certified all six stored candidates; it did not regenerate or trust the
approximate solves.

## Interpretation

- Tight or decreasing finite brackets are unconditional finite mathematics,
  but do not prove RH.
- A positive finite lower bound is expected and does not count against RH.
- Numerical or interval failure is inconclusive.
- A proved family of upper bounds tending to zero along an unbounded
  subsequence would prove RH.
- A uniform positive lower bound for all `N`, or along an unbounded
  subsequence, would disprove RH.
- A hypothetical exact finite identity `d_N=0` would prove RH.

No finite prefix can establish either required limiting statement. The
aggregate remains `EXPLORATORY`, stays out of the global claim ledger, and
leaves RH `UNRESOLVED`.

## Primary sources

- Báez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann Hypothesis](https://arxiv.org/abs/math/0202141).
- Báez-Duarte, Balazard, Landreau, and Saias,
  [Sur l'autocorrélation multiplicative de la fonction "partie fractionnaire"](https://arxiv.org/abs/math/0306251).
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058).
- Chen and Qi,
  [The best bounds of harmonic sequence](https://arxiv.org/abs/math/0306233),
  supplying a stronger published enclosure than the deliberately looser
  `1/(2(M+1)) < H_M-log(M)-gamma < 1/(2M)` oracle bound.

The conjectural asymptotic size of `d_N`, RH-conditional Möbius coefficient
families, and unverified Dirichlet-polynomial normalizations are deliberately
excluded from every v1 acceptance gate.
