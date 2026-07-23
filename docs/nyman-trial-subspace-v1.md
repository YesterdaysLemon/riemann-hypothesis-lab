# Arithmetic trial subspaces at the first dyadic bridge

## Verdict

This milestone separates a successful full-block theorem from a failed
low-dimensional ansatz at the same finite scale.  For the natural-dilate
Nyman problem, the existing certificate proves

\[
d_{512}^2<\frac9{10}d_{256}^2.
\]

The audit in this note considers the eight-dimensional arithmetic trial
space defined below.  Its best possible gain is strictly smaller than
`d_256^2/10`.  Thus this particular trial rule cannot explain the certified
`256 -> 512` contraction and cannot serve as the missing base instance of the
proposed all-scale trial-subspace lemma.

This is a finite rejection of one explicit ansatz.  It does not reject other
trial spaces, does not decide the all-scale recurrence, and does not prove or
disprove the Riemann Hypothesis.  The global status remains `UNRESOLVED`.

## The frozen eight-column rule

Let `J={257,...,512}`, put

\[
\ell_j=\log(512/j),
\]

and let `mu` be the Moebius function.  The columns, in frozen order, are

\[
\begin{aligned}
v_1(j)&=\mu(j)\ell_j,&
v_2(j)&=\mu(j),&
v_3(j)&=\mu(j)\ell_j^2,\\
v_4(j)&=|\mu(j)|\ell_j,&
v_5(j)&=|\mu(j)|,&
v_6(j)&=1-|\mu(j)|,\\
v_7(j)&=\mathbf 1_{2\mid j},&
v_8(j)&=\mathbf 1_{3\mid j}.&&
\end{aligned}
\]

Multiplying any column by a nonzero scalar would leave the trial span
unchanged.  In particular, `mu(j)*ell_j` is span-equivalent to the usual
BCF normalization `mu(j)*(1-log(j)/log(512))`.

Let `P_256` be orthogonal projection onto the old natural-dilate span and
write

\[
w_j=(I-P_{256})\rho_j,
\qquad
r_{256}=(I-P_{256})\chi.
\]

With `V` the `256 x 8` matrix above, define

\[
A=V^{\mathsf T}SV,
\qquad
u=V^{\mathsf T}t,
\qquad
\Gamma(V)=u^{\mathsf T}A^{-1}u,
\]

where `S_ij=<w_i,w_j>` and `t_i=<r_256,w_i>`.  The value `Gamma(V)` is the
largest reduction in squared error available from these eight directions
after the old 256 coefficients are reoptimized.

## Certificate route

The sign-bearing calculation does not certify an approximate Schur solve.
Instead, it forms the aggregate space

\[
\mathcal W=
\operatorname{span}\{\rho_1,\ldots,\rho_{256},
\sum_{j\in J}v_1(j)\rho_j,\ldots,
\sum_{j\in J}v_8(j)\rho_j\}
\]

and certifies a lower bound for its exact optimal squared distance

\[
F_V=\operatorname{dist}(\chi,\mathcal W)^2.
\]

The frozen `N=256` upper endpoint is

\[
U_{256}=
\frac{2801788463381838392697853210438663771}{2^{128}}.
\]

Fixed-order interval `LDL^T` positivity of the aggregate augmented matrix
certifies

\[
F_V>\frac9{10}U_{256}
=
\frac{25216096170436545534280678893947973939}
{3402823669209384634633746074317682114560}.
\]

Because `d_256^2 <= U_256` and analytically
`F_V=d_256^2-Gamma(V)`, it follows that

\[
\boxed{\Gamma(V)<\frac1{10}d_{256}^2.}
\]

The earlier full-block certificate gives a gain greater than
`0.1028950255... d_256^2`.  Consequently this eight-column span misses a
strictly positive part of the gain that the full block achieves.  The proof
above only needs the simpler rational thresholds `1/10` and `9/10`.

## Exploratory scout that selected the audit

A same-backend 384-bit Schur scout, independently repeated at 256 bits, gave
the following cumulative diagnostics.  These decimals select and explain the
certificate target; they carry no sign-bearing role.

| Number of frozen columns | Estimated `Gamma(V)/d_256^2` |
|---:|---:|
| 1 | 0.061612 |
| 2 | 0.075333 |
| 3 | 0.078537 |
| 4 | 0.082668 |
| 5 | 0.086941 |
| 6 | 0.086945 |
| 7 | 0.086956 |
| 8 | 0.087021 |

The raw untapered Moebius column alone gives approximately
`0.07499094784 d_256^2`.  This resolves the shorthand in the preceding
`N=512` note: its `0.074991` scout was the raw `mu(j)` direction, not the
log-tapered direction.

Adding rational tapers, squarefree strata, small-prime strata, and the
independent Ehm-inspired `1/j` corrections raised the best small-space scout
only to about `0.09232 d_256^2`; the largest versions were also severely
ill-conditioned.  None reached `0.1 d_256^2`.  Those larger searches remain
exploratory and are not claims of the frozen audit.

## Exact dyadic lift

For this section, pass to the reciprocal coordinate `t=1/x`. The Hilbert
space is `L^2((0,infinity),dt/t^2)`; write `e_n(t)={t/n}` and
`chi_t(t)=1_[1,infinity)(t)`. Under the reciprocal isometry these represent
`rho_n` and `chi`, and `r_N` below denotes the transported old residual.

The failed basis search nevertheless exposes a useful exact coordinate.  Put

\[
h_n(t)=2\left\{\frac{t}{2n}\right\}-\left\{\frac{t}{n}\right\}.
\]

Direct cancellation of the linear terms gives

\[
h_n(t)=
\left\lfloor\frac tn\right\rfloor
-2\left\lfloor\frac{t}{2n}\right\rfloor\in\{0,1\},
\qquad
\|h_n\|^2=\frac{\log 2}{n}.
\]

For coefficients `a_n`, define

\[
F_{a,N}(t)=\sum_{n\le N}a_nh_n(t),
\qquad
D_{a,N}(x)=\sum_{n\le N}a_n\left\lfloor\frac xn\right\rfloor.
\]

Then

\[
F_{a,N}(t)=D_{a,N}(t)-2D_{a,N}(t/2).
\]

Modulo the old span, this is the even new-block vector
`v_(2n)=2a_n` for `N/2<n<=N`.  Therefore, with `r_N` the old optimal
residual,

\[
\langle r_N,F_{a,N}\rangle=t^{\mathsf T}v,
\qquad
v^{\mathsf T}Sv\le\|F_{a,N}\|^2.
\]

For `a_n=mu(n)`, Moebius inversion gives `D_(mu,N)(x)=1` throughout
`1<=x<N+1`, because then `floor(x)<=N`, hence

\[
F_{\mu,N}(t)=
\begin{cases}
1,&1\le t<2,\\
-1,&2\le t<N+1.
\end{cases}
\]

This yields the sharper exact norm floor `||F_(mu,N)||^2 >= 1-1/(N+1)`.
It does not
yield the lower correlation bound needed for a Schur gain: Cauchy--Schwarz
goes in the opposite direction.  At `N=256`, a direct 384-bit scout for this
lifted direction gave only about `0.03114862 d_256^2`.

The reason is visible by comparing it with the formal full Moebius lift

\[
F_\infty=2\mathbf 1_{[1,2)}-\chi_t.
\]

At `N=256`, projecting this ideal function against the old span gives an
exploratory gain ratio of about `0.53930`, far above the target.  But for
`tau_N=F_(mu,N)-F_\infty`, the finite-cutoff correlation is

\[
E_N=d_N^2,
\qquad
s=-\sum_{n\le N}\frac{c_n}{n},
\qquad q_1=1+c_1,
\]

where `c_n` are the old optimal coefficients, and

\[
\langle r_N,F_{\mu,N}\rangle
=-E_N+2s\log2+q_1+\langle r_N,\tau_N\rangle.
\]

The ideal part contributes approximately `0.01198270`, while the cutoff tail
contributes approximately `-0.00811531`; it cancels about 67.7% of the
promising numerator.  The projected tail norm is about `0.06348`, larger
than both the ideal projected norm (`0.03234`) and the final finite-lift norm
(`0.05832`).  Thus the tail is neither small nor harmlessly aligned.  These
figures are exploratory, but the core/tail identity itself is exact.

Polynomial Moebius tapers also have an exact core reduction.  If

\[
a_n=\mu(n)P(n/N),
\qquad P(z)=\sum_r p_rz^r,
\]

then for `x<=N`,

\[
D_{a,N}(x)=
\sum_{m\le x}\sum_r p_rN^{-r}
\prod_{p\mid m}(1-p^r).
\]

The Dirichlet series for the inner arithmetic function is

\[
\sum_{m\ge1}
\frac{\prod_{p\mid m}(1-p^r)}{m^s}
=\frac{\zeta(s)}{\zeta(s-r)}.
\]

Thus the core is explicit, but useful cancellation for positive `r` again
meets the shifted zeros of zeta rather than bypassing them.

## Exact local cancellation and the tail obstruction

There is a second exact coordinate system for an arbitrary new-block vector
`y_(N+1),...,y_(2N)`.  Write

\[
\sigma=\sum_{j=N+1}^{2N}\frac{y_j}{j},
\qquad
Y_m=\sum_{j=N+1}^{m}y_j.
\]

On `m<=t<m+1` with `N+1<=m<2N`, write the old residual as
`r_N(t)=st+q_m`.  Then

\[
r_N(t)-\sum_{j=N+1}^{2N}y_j\{t/j\}
=(s-\sigma)t+q_m+Y_m.
\]

Put

\[
w_m=\frac1{m(m+1)},\quad
L_m=\log\frac{m+1}{m},\quad
H_m=\frac{L_m}{w_m},\quad
\kappa_m=1-\frac{L_m^2}{w_m}.
\]

Completing the square gives the exact interval identity

\[
\int_m^{m+1}
\frac{[(s-\sigma)t+q_m+Y_m]^2}{t^2}\,dt
=w_m[q_m+Y_m+(s-\sigma)H_m]^2
+(s-\sigma)^2\kappa_m.
\]

Moreover,

\[
0<\kappa_m\le\frac1{(2m+1)^2},
\qquad
\sum_{m=N+1}^{2N-1}\kappa_m<\frac1{4N}.
\]

Imposing `sigma=0` leaves the residual unchanged below `N+1`.  The values

\[
Y_m^*=-q_m-sH_m
\]

can be realized simultaneously on the first block, with `y_(2N)` then chosen
to enforce `sigma=0`.  They capture the exact local amount

\[
G_{\mathrm{loc}}=
\sum_{m=N+1}^{2N-1}m(m+1)
\left(
s\log\frac{m+1}{m}
+\frac{q_m}{m(m+1)}
\right)^2,
\]

leaving less than `s^2/(4N)` on those intervals.

Still assuming `sigma=0`, the obstruction begins at `2N`.  Define

\[
Z_m=\sum_{j=N+1}^{2N}y_j\left\lfloor\frac mj\right\rfloor.
\]

The exact tail-energy increment `Delta T`, meaning the new residual's tail
energy minus the old residual's tail energy, is

\[
\Delta T=
\sum_{m\ge2N}
\left[
2Z_m\left(sL_m+q_mw_m\right)+Z_m^2w_m
\right],
\]

which has no fixed sign.  Its jumps satisfy

\[
Z_m-Z_{m-1}=
\sum_{\substack{j\mid m\\N<j\le2N}}y_j.
\]

Inside the first block the only possible block divisor is `j=m`, which makes
triangular cancellation exact.  Beyond `2N`, several block divisors alias
into the same integer and no coefficients beyond `2N` remain to undo them.
Local cancellation therefore does not imply global gain.

This is the same structural type of omitted-divisor error that Ehm isolates
as a major unresolved estimation problem in his Gram-form decomposition.
The conditional BCF analysis controls related Moebius families only after
assuming strong information about zeta zeros.

## Frozen bindings and reproduction

The tracked [artifact](../results/nyman-trial-subspace-v1.json) has canonical
payload SHA-256
`467d6d819700a87f656e917bcb3e63008d7243fa5f412f75b6eec4ade5a67956`.
Its frozen coefficient-matrix, aggregate-system, and lower-certificate hashes
are, respectively,

```text
ec9ea0b1a3cb7439f2e30c552c38e2e5a4620bf5a5ec1285d68901a2a29a7f99
dd9a571997cc10a32fe332ad2127d58913c021054f058aee90253fdb85422fec
a58186360f701a26fede8872e463e28d57f4e2fd6f06a001137495216b29eea0
```

Generation uses a clean 768-bit kernel. The standalone verifier first
regenerates that frozen artifact, then independently rebuilds the kernel,
basis, aggregate system, and 265-pivot certificate at exactly 1536 bits. It
returned
`REPRODUCED_CERTIFIED_FINITE_NYMAN_TRIAL_SUBSPACE_REJECTION`. Both passes use
the same implementation and FLINT/Arb backend, so this is a high-precision
consistency replay rather than an independent clean-room implementation.

On PowerShell, from an installed checkout:

```powershell
.\.venv\Scripts\rh-lab.exe nyman-trial-subspace-audit `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1 `
  --output results\nyman-trial-subspace-v1-new.json

.\.venv\Scripts\rh-lab.exe verify-nyman-trial-subspace-audit `
  --artifact results\nyman-trial-subspace-v1.json `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1 `
  --bits 1536
```

The full replay is intentionally expensive: it reconstructs two large Arb
systems and repeats fixed-order interval elimination. The artifact is
classified `CERTIFIED_FINITE`, explicitly sets `resolves_rh=false`, and binds
the limitation that no other trial subspace or later scale is covered.

## Remaining theorem target

For `N=2^k`, a successful explicit trial family must still prove

\[
\Gamma_k(V_k)\ge\frac{E_k}{k+2}
\qquad(k\ge9).
\]

The finite rejection here says that the frozen eight-column arithmetic rule
does not even reach the corresponding `k=8` threshold.  The exact identities
above narrow the next credible routes to one of the following:

1. a lower alignment bound between `r_N` and an explicit dyadic lift;
2. a tail-stability estimate that controls divisor aliases after `2N`;
3. a multiscale construction in which later blocks cancel aliases left by
   earlier blocks.

No one of these estimates is proved here.

The [multiscale follow-up](nyman-multiscale-tail-v1.md) tests all three
directions. It identifies the dyadic seed with Báez-Duarte's divergent first
Vasyunin correction, proves an exact alias--Möbius extension identity, and
shows why pointwise cancellation still leaves an uncontrolled weighted-norm
tail. It also records exploratory block scouts; none supplies the uniform
estimate above.

## Sources

- Werner Ehm,
  [On certain Gram matrices and their associated series](https://arxiv.org/abs/2405.06349).
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191).
- Baez-Duarte,
  [A strengthening of the Nyman--Beurling criterion for the Riemann Hypothesis](https://arxiv.org/abs/math/0202141).
- Baez-Duarte,
  [A divergent Vasyunin correction](https://arxiv.org/abs/math/0506318).
