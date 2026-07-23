# The exact dyadic Schur-product criterion

## Status

This note records a human-auditable reformulation of the strong
Nyman--Beurling criterion. It is an exact equivalence:

\[
\mathrm{RH}
\quad\Longleftrightarrow\quad
\sum_{k\ge 0}\eta_k=+\infty,
\]

where \(\eta_k\) is the *optimal* relative gain obtained by enlarging the
natural-dilate space from dimension \(2^k\) to dimension \(2^{k+1}\).
The equivalence identifies a precise all-scale target. It does not establish
that the series diverges.

In particular, the repository's explicit rebased vectors \(p_1,\ldots,p_4\)
are not claimed to be finite-dimensional orthogonal projections, and their
certified energy decreases are not the \(\eta_k\) in this note. The global
status remains **RH `UNRESOLVED`**.

## 1. Hilbert-space normalization

Work in

\[
H=L^2\!\left((0,\infty),\frac{dt}{t^2}\right),
\qquad
e_n(t)=\left\{\frac{t}{n}\right\},
\qquad
\chi(t)=\mathbf 1_{[1,\infty)}(t).
\]

This is the reciprocal-coordinate image of the repository's
\(L^2((0,\infty),dx)\) normalization under \(t=1/x\). For \(k\ge0\), put

\[
\mathcal V_k=\operatorname{span}(e_1,\ldots,e_{2^k}),
\qquad
P_k=P_{\mathcal V_k},
\qquad
r_k=(I-P_k)\chi,
\qquad
E_k=\lVert r_k\rVert_H^2.
\]

Every \(\mathcal V_k\) is finite-dimensional and hence closed, and

\[
\mathcal V_0\subset\mathcal V_1\subset\cdots.
\]

Let

\[
\mathcal V_\infty=
\overline{\bigcup_{k\ge0}\mathcal V_k}.
\]

Because the dyadic cutoffs are cofinal in the positive integers,
\(\mathcal V_\infty\) is the closed span of all natural dilates \(e_n\).
For increasing closed subspaces, the distances to the subspaces decrease to
the distance to the closure of their union. Indeed, the lower inequality is
immediate from inclusion; for the reverse inequality, approximate any vector
in \(\mathcal V_\infty\) by one in the union, which belongs to some
\(\mathcal V_k\). Consequently

\[
E_k\downarrow \operatorname{dist}(\chi,\mathcal V_\infty)^2.
\]

Baez-Duarte's strong natural-dilate theorem therefore gives

\[
\boxed{\mathrm{RH}\Longleftrightarrow E_k\longrightarrow0.}
\]

This use of only dyadic dimensions loses nothing: monotonicity then extends a
dyadic limit of zero to the full sequence of natural cutoffs.

## 2. Finite independence and strict noncontainment

Two elementary facts ensure that every finite Schur system below is
nondegenerate.

### Linear independence of the natural dilates

Suppose

\[
\sum_{n=1}^{M}a_ne_n=0
\quad\text{almost everywhere.}
\]

On every open unit interval the left side is affine and continuous, so
almost-everywhere vanishing makes it identically zero there. At the integer
\(m\), the jump of \(e_n\) is \(-1\) when \(n\mid m\) and zero otherwise.
Thus, for \(1\le m\le M\),

\[
\sum_{n\mid m}a_n=0.
\]

At \(m=1\) this gives \(a_1=0\); induction on \(m\) then gives \(a_m=0\)
for every \(m\le M\). Hence every finite family
\((e_1,\ldots,e_M)\) is linearly independent and its Gram matrix is positive
definite.

### The target is in no finite natural span

Fix \(Q\ge1\) and suppose, for contradiction, that

\[
\chi=\sum_{n\le Q}c_ne_n
\quad\text{almost everywhere.}
\]

On \(0<t<1\), all \(e_n(t)=t/n\), so

\[
\sum_{n\le Q}\frac{c_n}{n}=0. \tag{1}
\]

On \(m<t<m+1\), for \(1\le m\le Q\), the slope term vanishes by (1), and
\(\chi=1\) gives

\[
1+\sum_{n\le Q}c_n\left\lfloor\frac{m}{n}\right\rfloor=0.
\]

Taking successive differences yields

\[
\sum_{n\mid1}c_n=-1,
\qquad
\sum_{n\mid m}c_n=0\quad(2\le m\le Q).
\]

Moebius inversion forces \(c_n=-\mu(n)\). Equation (1) would then require

\[
\sum_{n\le Q}\frac{\mu(n)}{n}=0. \tag{2}
\]

For \(Q=1\), (2) is plainly false. For \(Q\ge2\), choose by Bertrand's
postulate a prime \(q\) with \(Q/2<q\le Q\), and let
\(L=\operatorname{lcm}(1,\ldots,Q)\). In the integer

\[
L\sum_{n\le Q}\frac{\mu(n)}n
=\sum_{n\le Q}\mu(n)\frac Ln,
\]

every term except \(n=q\) vanishes modulo \(q\). The remaining term is
\(-L/q\), which is nonzero modulo \(q\): \(2q>Q\) makes \(q\) the only
multiple of \(q\) in the range, while \(q^2>Q\) makes the \(q\)-adic
valuation of \(L\) exactly one. This contradicts (2).

Therefore

\[
\chi\notin\mathcal V_k,\qquad E_k>0
\quad\text{for every finite }k.
\]

The prime-existence input is the classical theorem for which Erdos's 1932
paper is linked in the sources below.

## 3. The exact Schur decrement

Write \(N=2^k\), and introduce the new block

\[
\mathcal U_k=\operatorname{span}(e_{N+1},\ldots,e_{2N}),
\qquad
w_j=(I-P_k)e_j\quad(N<j\le2N).
\]

The residualized block

\[
\mathcal W_k=(I-P_k)\mathcal U_k
\]

is orthogonal to \(\mathcal V_k\), and

\[
\mathcal V_{k+1}=\mathcal V_k\mathbin{\mathop{\oplus}^{\perp}}\mathcal W_k.
\]

The vectors \(w_{N+1},\ldots,w_{2N}\) are linearly independent. Otherwise a
nonzero combination of the new \(e_j\)'s would lie in \(\mathcal V_k\),
contradicting the finite independence just proved.

Define the \(N\times N\) Gram matrix and correlation vector

\[
(S_k)_{ij}=\langle w_i,w_j\rangle_H,
\qquad
(t_k)_j=\langle r_k,w_j\rangle_H
              =\langle r_k,e_j\rangle_H,
\quad N<i,j\le2N.
\]

The second equality uses \(r_k\perp\mathcal V_k\). The matrix \(S_k\) is
positive definite. For a coefficient vector \(a\in\mathbb R^N\),

\[
\left\lVert r_k-\sum_{j=N+1}^{2N}a_jw_j\right\rVert_H^2
=E_k-2a^{\mathsf T}t_k+a^{\mathsf T}S_ka.
\]

The unique minimizer is \(a=S_k^{-1}t_k\). Hence the exact optimal decrement
is

\[
\boxed{
\Delta_k:=E_k-E_{k+1}
=t_k^{\mathsf T}S_k^{-1}t_k.
} \tag{3}
\]

This is the Hilbert-space form of the finite block Schur-complement identity;
it is an analytic identity, not a numerical inversion prescription.

Since \(E_{k+1}>0\), the relative optimal gain

\[
\eta_k=\frac{\Delta_k}{E_k}
\]

satisfies \(0\le\eta_k<1\), and (3) gives

\[
E_{k+1}=E_k(1-\eta_k). \tag{4}
\]

## 4. Product and series equivalences

Iterating (4), for any \(K>k_0\),

\[
E_K
=E_{k_0}\prod_{k=k_0}^{K-1}(1-\eta_k). \tag{5}
\]

The finite noncontainment result guarantees that no factor is zero. Taking
logarithms in (5) gives

\[
E_K\longrightarrow0
\quad\Longleftrightarrow\quad
\sum_{k\ge k_0}-\log(1-\eta_k)=+\infty. \tag{6}
\]

For \(0\le x<1\), \(x\le-\log(1-x)\). Conversely, if infinitely many
\(\eta_k\ge1/2\), then \(\sum\eta_k\) already diverges; otherwise,
eventually \(\eta_k<1/2\), where
\(-\log(1-\eta_k)\le2\eta_k\). Thus

\[
\sum_k-\log(1-\eta_k)=+\infty
\quad\Longleftrightarrow\quad
\sum_k\eta_k=+\infty. \tag{7}
\]

Combining the strong criterion, (5), (6), and (7) proves the promised exact
equivalence:

\[
\boxed{
\begin{aligned}
\mathrm{RH}
&\Longleftrightarrow E_k\to0\\
&\Longleftrightarrow
\prod_{k\ge k_0}(1-\eta_k)=0\\
&\Longleftrightarrow
\sum_{k\ge k_0}-\log(1-\eta_k)=+\infty\\
&\Longleftrightarrow
\sum_{k\ge k_0}\eta_k=+\infty .
\end{aligned}
} \tag{8}
\]

Any finite starting index \(k_0\) gives the same conclusion.

## 5. Explicit trial subspaces and an inverse-free certificate form

Equation (3) uses the entire \(N\)-dimensional new block. A proposed
arithmetic rule can instead be represented by a full-column-rank matrix
\(B_k\in\mathbb R^{N\times r}\). Its columns select an \(r\)-dimensional
coefficient subspace. Put

\[
A_k=B_k^{\mathsf T}S_kB_k,
\qquad
u_k=B_k^{\mathsf T}t_k.
\]

Then \(A_k\) is positive definite, and the best gain available in this trial
subspace is

\[
\Gamma_k(B_k)=u_k^{\mathsf T}A_k^{-1}u_k\le\Delta_k. \tag{9}
\]

The variational identity

\[
\Gamma_k(B_k)
=\max_{z\in\mathbb R^r}
\left(2u_k^{\mathsf T}z-z^{\mathsf T}A_kz\right) \tag{10}
\]

turns every explicit \(z\) into a rigorous lower witness:

\[
\Delta_k\ge\Gamma_k(B_k)
\ge2u_k^{\mathsf T}z-z^{\mathsf T}A_kz. \tag{11}
\]

Thus a proof need not invert \(A_k\). For example, it would suffice to give
one explicit rule \(B_k,z_k\) for every sufficiently large \(k\) such that

\[
2u_k^{\mathsf T}z_k-z_k^{\mathsf T}A_kz_k
\ge\frac{E_k}{k+2}. \tag{12}
\]

Then \(E_{k+1}\le((k+1)/(k+2))E_k\); the factors telescope to zero, and
(8) proves RH. No family satisfying (12) is known here.

## 6. An inverse-free full-block correlation target

There is also a universal spectral upper bound that removes the inverse from
the full block. Set

\[
\kappa=\int_0^\infty\{u\}^2\frac{du}{u^2}
=\log(2\pi)-\gamma.
\]

Scaling gives

\[
\lVert e_j\rVert_H^2=\frac{\kappa}{j}.
\]

Since orthogonal projection is contractive,

\[
\begin{aligned}
\lambda_{\max}(S_k)
&\le\operatorname{tr}(S_k)
=\sum_{j=N+1}^{2N}\lVert w_j\rVert_H^2\\
&\le\kappa\sum_{j=N+1}^{2N}\frac1j
=\kappa(H_{2N}-H_N)
<\kappa\log2.
\end{aligned} \tag{13}
\]

Let \(L_k=\kappa(H_{2N}-H_N)\). Applying the variational form to the full
block with the explicit choice \(a=t_k/L_k\), and using
\(t_k^{\mathsf T}S_kt_k\le L_k\lVert t_k\rVert_2^2\), yields

\[
\Delta_k\ge\frac{\lVert t_k\rVert_2^2}{L_k}
=
\frac{
\displaystyle\sum_{j=N+1}^{2N}
\left|\langle r_k,e_j\rangle_H\right|^2
}{
\kappa(H_{2N}-H_N)
}. \tag{14}
\]

Consequently the inverse-free block-correlation estimate

\[
\boxed{
\sum_{j=N+1}^{2N}
\left|\langle r_k,e_j\rangle_H\right|^2
\ge c\,\frac{E_k}{\log N}
} \tag{15}
\]

for one fixed \(c>0\) and all sufficiently large dyadic \(N=2^k\) would give
\(\eta_k\gg1/k\). Its series would diverge, so (8) would prove RH.
The constant in (15) need not be optimized. The unresolved content is the
arithmetic lower bound on the correlations, not the Schur algebra or the
trace estimate.

## 7. Divisor aliases: why local cancellation is insufficient

The reciprocal coordinate makes the main obstruction to (12) or (15)
explicit. Let a proposed balanced block \(y=(y_j)_{N<j\le2N}\) satisfy

\[
\sum_{j=N+1}^{2N}\frac{y_j}{j}=0.
\]

On an integer interval \(I_m=[m,m+1)\), define the divisor-alias value

\[
Z_m=\sum_{j=N+1}^{2N}y_j\left\lfloor\frac{m}{j}\right\rfloor.
\]

If the old residual is \(r_k(t)=st+q_m\) on \(I_m\), then adding the block to
the approximant changes the residual on the tail to \(r_k+Z_m\). With

\[
w_m=\frac1{m(m+1)},
\qquad
L_m=\log\frac{m+1}{m},
\]

the exact tail-energy increment is

\[
\sum_{m\ge2N}
\left[
2Z_m\left(sL_m+q_mw_m\right)+Z_m^2w_m
\right]. \tag{16}
\]

Moreover,

\[
Z_m-Z_{m-1}
=\sum_{\substack{j\mid m\\N<j\le2N}}y_j. \tag{17}
\]

On the first new shell, the cumulative coordinates are triangular and can be
chosen to cancel the interval residual very accurately. Beyond \(2N\),
several divisors from the block can alias at the same integer. There are no
later coefficients in a one-block construction to undo those collisions.
The cross term in (16) is signed, so a small first-shell remainder does not
imply a positive global gain.

The quantitative burden of discarding that sign is severe. In the exact
local-plus-alias decomposition used by the repository, write the ideal local
correction as \(f\), the alias error as \(e\), and let

- \(G=\lVert f\rVert^2=\langle r_k,f\rangle\) be the ideal local capture,
- \(A=\lVert e\rVert^2\) be the squared norm of the alias error, and
- \(X=\langle r_k-f,e\rangle=\langle r_k,e\rangle\) be its signed
  correlation with the residual left by the local correction.

Here \(f\perp e\). These identities give the fixed-addition and optimally
rescaled formulas below directly.

The fixed-addition gain is \(G+2X-A\), while optimal scalar rescaling captures

\[
\frac{(G+X)^2}{G+A}. \tag{18}
\]

Writing \(\gamma=G/E_k\) and \(\delta=A/E_k\), Cauchy--Schwarz alone gives

\[
\frac{G+2X-A}{E_k}
\ge
\gamma-2\sqrt{(1-\gamma)\delta}-\delta. \tag{19}
\]

To retain a fixed fraction \(a\), with \(0\le a\le1\), of \(G\) from this
norm-only inequality, it is sufficient to require

\[
\delta\le
\left(\sqrt{1-a\gamma}-\sqrt{1-\gamma}\right)^2
\sim\frac{(1-a)^2}{4}\gamma^2. \tag{20}
\]

At the harmonic target \(\gamma\asymp1/k\), (20) asks for

\[
\frac{A}{G}=O(1/k),
\]

not merely a fixed fractional bound \(A/G<\theta<1\). This does not prove
that a norm-only route is impossible; it explains why a genuinely signed
bilinear estimate for (16), or a multiscale construction that cancels its
aliases, is the more credible target.

## 8. Burnol's lower bound fixes the sharp scale

Define the natural-dilate distance

\[
d_N^2:=
\operatorname{dist}\!\left(
\chi,\operatorname{span}(e_1,\ldots,e_N)
\right)^2,
\qquad E_k=d_{2^k}^2.
\]

Burnol states his Nyman--Beurling lower bound for the continuum distance
\(D(\lambda)\). Under the reciprocal isometry, the natural functions
\(e_1,\ldots,e_N\) correspond to parameters \(1,1/2,\ldots,1/N\) in
Burnol's continuum space at \(\lambda=1/N\). The natural space is therefore a
subspace of that continuum space, so
\(d_N\ge D(1/N)\). Burnol's bound consequently implies under RH

\[
\liminf_{N\to\infty}d_N^2\log N
\ge
\sum_{\Re\rho=1/2}\frac{m(\rho)^2}{|\rho|^2}.
\]

Under RH the zero sum is at least

\[
C_0=2+\gamma-\log(4\pi)>0.
\]

If RH is false, the strong criterion gives
\(d_N\downarrow d_\infty>0\), so \(d_N^2\log N\to\infty\). The two cases
together give the unconditional guardrail

\[
\boxed{
\liminf_{N\to\infty}d_N^2\log N\ge C_0.
} \tag{21}
\]

Suppose an all-scale proof tried to establish, for some \(c>1\),

\[
\eta_k\ge\frac{c}{k}
\quad\text{for every sufficiently large }k.
\]

Then (5) and \(\log(1-x)\le-x\) would give
\(E_k=O(k^{-c})\), hence

\[
E_k\log(2^k)=O(k^{1-c})\longrightarrow0,
\]

contradicting (21). Thus an eventual harmonic lower bound with leading
constant strictly larger than one is impossible.

The conditional theorem of Bettin, Conrey, and Farmer -- assuming RH and a
strong negative-moment bound for \(\zeta'(\rho)\) -- supplies an explicit
upper witness with Burnol's constant. Combined with Burnol's lower bound, it
gives the optimal-distance scale \(d_N^2\asymp1/\log N\) with that constant.
Substituting the model \(E_k\sim C/(k\log2)\) formally gives
\(\eta_k\sim1/(k+1)\). This identifies the coefficient-one harmonic behavior
as the sharp expected scale. The published asymptotic is not, by itself, a
proved pointwise asymptotic for every \(\eta_k\), and it supplies no
unconditional lower bound of the form (15).

## 9. What the finite certificates do and do not say

The natural-space \(256\to512\) certificate proves the one-sided bound

\[
\eta_8>\frac{51}{500}.
\]

It does so by combining a rigorous lower endpoint for the optimal
\(d_{256}^2\) with a rigorous upper endpoint supplied by an explicit
512-coordinate witness. It does not compute the exact value of \(\eta_8\).

The later rebased certificates for \(p_1,\ldots,p_4\) prove strict complete
energy decreases for those explicit vectors. They are valuable finite
existence results, but they do not assert

\[
p_j=P_k\chi,\qquad
\mathcal E(p_j)=E_k,\qquad\text{or}\qquad
\frac{\mathcal E(p_j)-\mathcal E(p_{j+1})}{\mathcal E(p_j)}=\eta_k
\]

for any \(k\). Finite witness gains therefore cannot be multiplied as the
optimal factors in (5) without a separate theorem identifying the relevant
projections and covering every later dyadic scale.

The exact open problem isolated by this note is to prove enough all-scale
optimal gain -- for example (12) or (15) -- to make
\(\sum_k\eta_k\) diverge. No such proof, and no positive limiting-distance
counterexample, is contained here. **RH remains `UNRESOLVED`.**

## Primary sources

- Luis Baez-Duarte,
  [*A strengthening of the Nyman--Beurling criterion for the Riemann
  hypothesis, 2*](https://arxiv.org/abs/math/0205003), for the strong
  natural-dilate equivalence.
- Jean-Francois Burnol,
  [*A lower bound in an approximation problem involving the zeros of the
  Riemann zeta function*](https://arxiv.org/abs/math/0103058), for the
  zero-sensitive approximation lower bound.
- Sandro Bettin, J. Brian Conrey, and David W. Farmer,
  [*An optimal choice of Dirichlet polynomials for the Nyman--Beurling
  criterion*](https://arxiv.org/abs/1211.5191), for the conditional optimal
  asymptotic and its hypotheses.
- Luis Baez-Duarte, Michel Balazard, Bernard Landreau, and Eric Saias,
  [*Sur l'autocorrelation multiplicative de la fonction "partie
  fractionnaire"*](https://arxiv.org/abs/math/0306251), for the fractional-part
  autocorrelation underlying the exact Gram entries and
  \(\kappa=\log(2\pi)-\gamma\).
- Werner Ehm,
  [*On certain Gram matrices and their associated
  series*](https://arxiv.org/abs/2405.06349), for exact Nyman Gram-form
  decompositions and the Moebius-inversion error perspective.
- Paul Erdos,
  [*Beweis eines Satzes von
  Tschebyschef*](https://users.renyi.hu/~p_erdos/1932-01.pdf) (1932), for an
  elementary proof of Bertrand's postulate used in the finite
  noncontainment argument.
