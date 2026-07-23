# Log-tapered Möbius core/tail audit

## Classification and scope

This note documents one explicit Nyman--Beurling candidate family at the six
finite dimensions

```text
N = 8, 16, 32, 64, 128, 256.
```

The canonical
[core/tail artifact](../results/nyman-mobius-core-tail-v1.json) has:

| Field | Value |
|---|---|
| `schema` | `rh-lab/nyman-mobius-core-tail-audit/v1` |
| `audit_id` | `nyman-mobius-core-tail-v1` |
| `classification` | `EXPLORATORY` |
| `hypothesis_status` | `UNRESOLVED` |
| `audit_outcome` | `FINITE_EXPLICIT_FAMILY_ENERGY_CERTIFIED` |
| `payload_sha256` | `b6579bf5774e5413ee5f05b2d3f125d1afcd6f8e519a3942b7352924d01acb8e` |

The result has four deliberately separate trust layers:

| Layer | What is established |
|---|---|
| exact identities | The coefficient sign, divisor decomposition, and core formulas below are exact algebraic identities; the artifact checks their finite arithmetic transcript through `k=256` |
| finite interval certificates | Arb encloses the candidate total, two separately assembled same-backend core formulas, their intersection, the complementary tail, scaled values, and comparisons with the certified optimum at all six `N` |
| published results | The strong natural-dilate criterion, the source of the log taper and its conditional asymptotic, and the rational Gram-kernel formula are cited dependencies, not theorems reproved by the artifact |
| unproved bridge | No bound controls every sufficiently large `N`; the decisive all-scale statement remains a research target |

Nothing in this note proves or disproves the Riemann Hypothesis. Finite
decreases and fitted rates cannot replace the missing all-scale theorem.

## Normalization and the coefficient sign

Let

\[
\rho_n(x)=\left\{\frac1{nx}\right\},\qquad
\chi(x)=\mathbf 1_{(0,1]}(x),
\]

and

\[
d_N^2=\inf_{c_1,\ldots,c_N\in\mathbb R}
\left\|\chi-\sum_{n=1}^N c_n\rho_n\right\|_2^2.
\]

For literature Dirichlet coefficients `a_n`, the repository coefficients are

\[
\boxed{c_n=-a_n.}
\]

Under the isometry obtained from `t=1/x`, the repository residual becomes

\[
R_N(t)
=
\mathbf 1_{[1,\infty)}(t)
+\sum_{n\le N}a_n\left\{\frac tn\right\},
\]

and

\[
E_N(a)
:=
\int_0^\infty |R_N(t)|^2\frac{dt}{t^2}.
\]

The plus sign in this reciprocal-variable residual is essential. For
`0<Re(s)<1`,

\[
\int_0^\infty
\left\{\frac tn\right\}t^{-s-1}\,dt
=-\frac{\zeta(s)}s n^{-s},
\]

so the Mellin transform of `R_N` is

\[
\frac{1-\zeta(s)A_N(s)}s,
\qquad
A_N(s)=\sum_{n\le N}a_n n^{-s}.
\]

Thus the real-variable energy here is the same candidate energy used in the
Dirichlet-polynomial normalization. Reversing `c_n=-a_n` would evaluate the
wrong Mellin residual.

## Exact core/tail decomposition

For arbitrary real coefficients supported on `1,...,N`, define

\[
S_N=\sum_{n\le N}\frac{a_n}{n},
\qquad
h_N(k)=\sum_{\substack{d\mid k\\d\le N}}a_d,
\qquad
H_N(m)=\sum_{k\le m}h_N(k),
\qquad
q_m=1-H_N(m).
\]

The floor identity

\[
\begin{aligned}
\sum_{n\le N}a_n\left\lfloor\frac tn\right\rfloor
&=
\sum_{n\le N}a_n
\sum_{\substack{j\ge1\\nj\le t}}1\\
&=
\sum_{k\le t}
\sum_{\substack{n\mid k\\n\le N}}a_n
=H_N(\lfloor t\rfloor)
\end{aligned}
\]

gives the exact piecewise residual

\[
R_N(t)=tS_N,\qquad 0<t<1,
\]

and

\[
R_N(t)=tS_N+q_m,\qquad m\le t<m+1,\quad m\ge1.
\]

Values at integer endpoints do not affect the integral. Splitting at `t=N`,
with the endpoint assigned to either side, therefore gives

\[
C_N:=\int_0^N |R_N(t)|^2\frac{dt}{t^2}
\]

and

\[
\boxed{
C_N
=
S_N^2+
\sum_{m=1}^{N-1}
\left[
S_N^2
+2S_Nq_m\log\frac{m+1}{m}
+\frac{q_m^2}{m(m+1)}
\right].
}
\]

The exact tail is

\[
T_N:=\int_N^\infty |R_N(t)|^2\frac{dt}{t^2}
\]

or, equivalently,

\[
T_N
=
\sum_{m=N}^{\infty}
\left[
S_N^2
+2S_Nq_m\log\frac{m+1}{m}
+\frac{q_m^2}{m(m+1)}
\right].
\]

Every summand in the last display is the nonnegative integral over one unit
interval, even though its expanded form contains cancellation.

In the repository Gram normalization,

\[
b_n=\langle\chi,\rho_n\rangle
=\frac{\log n+1-\gamma}{n},
\qquad
G_{mn}=\langle\rho_m,\rho_n\rangle.
\]

Substituting `c=-a` into the usual quadratic energy gives

\[
\boxed{
E_N(a)=1+2a^{\mathsf T}b+a^{\mathsf T}Ga.
}
\]

The implemented artifact obtains the tail from the exact identity

\[
\boxed{T_N=E_N(a)-C_N.}
\]

It does not replace the omitted intervals by a prime-counting formula.

## Bettin--Conrey--Farmer log taper

For `N>=2`, put `L=log(N)` and define

\[
a_n
=
\mu(n)\left(1-\frac{\log n}{L}\right)
=
\mu(n)\frac{\log(N/n)}L,
\qquad 1\le n\le N.
\]

The endpoint coefficient satisfies `a_N=0` exactly. Define

\[
\kappa_N=LS_N
=
\sum_{n\le N}\frac{\mu(n)\log(N/n)}n.
\]

The sign-bearing divisor identity is

\[
\boxed{
\sum_{d\mid k}\mu(d)\log d=-\Lambda(k).
}
\]

Together with

\[
\sum_{d\mid k}\mu(d)=\mathbf 1_{k=1},
\]

it gives, for `k<=N`,

\[
\boxed{
h_N(k)=\mathbf 1_{k=1}+\frac{\Lambda(k)}L.
}
\]

Consequently, for `1<=m<=N`,

\[
H_N(m)=1+\frac{\psi(m)}L,
\qquad
q_m=-\frac{\psi(m)}L,
\]

and the core residual is

\[
\boxed{
R_N(t)
=
\frac{\kappa_Nt-\psi(\lfloor t\rfloor)}{\log N},
\qquad 1\le t<N.
}
\]

This produces a second exact core formula:

\[
\boxed{
(\log N)^2 C_N
=
N\kappa_N^2
-2\kappa_N
\sum_{m=1}^{N-1}\psi(m)\log\frac{m+1}{m}
+\sum_{m=1}^{N-1}\frac{\psi(m)^2}{m(m+1)}.
}
\]

The middle sum may also be written

\[
\sum_{m=1}^{N-1}\psi(m)\log\frac{m+1}{m}
=
\sum_{k=1}^{N-1}\Lambda(k)\log\frac Nk.
\]

### The cutoff boundary

The shortcut

\[
h_N(k)=\mathbf 1_{k=1}+\frac{\Lambda(k)}{\log N}
\]

is valid only while `k<=N`, because then every divisor of `k` lies inside the
coefficient support. For `k>N`, the definition remains

\[
h_N(k)=\sum_{\substack{d\mid k\\d\le N}}a_d,
\]

and divisors above `N` are absent. Therefore the formula
`q_m=-psi(m)/log(N)` generally fails for `m>N`. Extending it through the
tail would compute a different function and invalidate the energy split.

## What the finite artifact certifies

### Exact finite arithmetic

The field `generation.exact_arithmetic` stores:

- all Möbius values for `1<=n<=256`;
- the prime base for every prime power through `256`;
- exact checks of the divisor sum and log-product identities for every
  `1<=k<=256`;
- `70` detected prime powers;
- identity transcript SHA-256
  `dc5b9898a27b3c6460c40d6e2f9fddd147c919a1e3003991e8ba62e13730b995`.

This is exact integer and rational-exponent bookkeeping. It is distinct from
the Arb evaluation of logarithms and Gram entries.

### Finite Arb records

Generation uses 256-bit Arb with `python-flint 0.9.0` and FLINT `3.6.0`.
The canonical `N=256` Gram system uses the rational Vasyunin construction
with reciprocal-orientation intersection and has content SHA-256

```text
ad7a07732ba1019584bc0983861a34cca089a54d6046503ae3616c87b44e13c9
```

At each `N`, the artifact uses two separate same-backend assemblies:

1. `general_divisor_core_error_squared` from `h_N`, `H_N`, and `q_m`;
2. `prime_counting_core_error_squared` from `psi(m)`;
3. `certified_core_intersection`, the intersection of those overlapping Arb
   enclosures;
4. `candidate_total_error_squared` from the full Gram quadratic form;
5. `tail_error_squared_by_total_minus_core`;
6. the three values multiplied by `log(N)`;
7. the core fraction of the total and `kappa_N`.

For every row, the artifact certifies that both core constructions overlap,
their difference contains zero, the endpoint coefficient is exactly zero,
and

\[
E_N^{\mathrm{BCF}}>C_N>0,\qquad T_N>0.
\]

The following table shows compact decimal views. The complete artifact stores
the rigorous dyadic midpoint/radius record for every value.

| `N` | Candidate `E_N^BCF` | Core `C_N` | Tail `T_N` | Core/total | `E_N^BCF log(N)` |
|---:|---:|---:|---:|---:|---:|
| 8 | `0.824866774831474...` | `0.727107141794654...` | `0.0977596330368202...` | `0.881484336598728...` | `1.71526223793603...` |
| 16 | `0.475766237934790...` | `0.443633152806281...` | `0.0321330851285096...` | `0.932460350133307...` | `1.31910410572045...` |
| 32 | `0.308865046777846...` | `0.296948576503038...` | `0.0119164702748083...` | `0.961418521133667...` | `1.07044468173790...` |
| 64 | `0.216590742422220...` | `0.211282265134962...` | `0.00530847728725784...` | `0.975490747074916...` | `0.900775574672082...` |
| 128 | `0.159843770358151...` | `0.157252636510787...` | `0.00259113384736396...` | `0.983789585033197...` | `0.775566811276766...` |
| 256 | `0.124287707716388...` | `0.122749844062796...` | `0.00153786365359184...` | `0.987626582854830...` | `0.689197393454984...` |

Across this frozen six-point grid, interval comparisons certify that:

- candidate total, core, and tail strictly decrease;
- candidate total times `log(N)`, core times `log(N)`, and tail times
  `log(N)` strictly decrease;
- the core fraction strictly increases.

These are six-point finite facts only. In particular, the decrease of
`E_N^BCF log(N)` is not an asymptotic estimate.

### Comparison with the optimized distance

The explicit BCF candidate is not the finite optimizer. Since

\[
d_N^2\le E_N^{\mathrm{BCF}},
\]

the artifact uses the previously certified strict brackets for `d_N^2` to
prove:

| `N` | Certified finite separation |
|---:|---|
| 8 | `E_BCF(8) > 34*d_8^2` |
| 16 | `E_BCF(16) > 26*d_16^2` |
| 32 | `E_BCF(32) > 21*d_32^2` |
| 64 | `E_BCF(64) > 19*d_64^2` |
| 128 | `E_BCF(128) > 16*d_128^2` |
| 256 | `E_BCF(256) > 15*d_256^2` |

At `N=256`, more than `98%` of the explicit-family energy is in the core,
but the total remains more than fifteen times the optimized squared distance.
Neither observation has an all-scale consequence.

## Artifact contract and replay boundary

The major JSON fields have the following roles:

| Field path | Meaning |
|---|---|
| `family` | frozen coefficient formula, repository sign, reciprocal residual, endpoint zero, and the inequality `d_N^2 <= E_BCF(N)` |
| `source_verification` | binding to the six-cell Nyman summary, plan, and index; it explicitly records that this source-summary check did not rebuild or numerically replay the old optimizer run |
| `generation.precision_bits` | frozen 256-bit generation precision |
| `generation.backend` | exact python-flint and FLINT versions |
| `generation.kernel` | construction, dimension, and content commitment for the newly built BCF Gram system |
| `generation.exact_arithmetic` | Möbius values, prime-power bases, exact identity checks, and transcript commitment |
| `generation.records[*]` | all Arb total/core/tail, scaled, share, slope, kappa, and per-cell checks |
| `generation.finite_grid_trends` | the seven strict comparisons on the frozen grid |
| `candidate_vs_certified_optimum` | source brackets, ratio enclosures, strict integer factors, and positive margins |
| `dependency_manifest` | published dependencies, exact manual identity layer, and the explicitly unproved all-scale target |
| `finite_conclusion` | the permitted finite summary |
| `limitation` | the non-promotion and cutoff warning |

The verifier first requires exact canonical regeneration of the artifact at
256 bits. It then rebuilds all six numerical records at 512 bits and requires
every higher-precision enclosure to lie inside its corresponding generation
enclosure. The replay kernel content SHA-256 is

```text
fdb64334dbb6b234d92b38f2b2956d19bf1afe1077c37562171be9aebf57aae5
```

This is a separate higher-precision replay using the same implementation,
python-flint, FLINT, formulas, and kernel construction. It is not an
independent backend or a clean-room derivation.

The parser rejects duplicate object keys, nonstandard JSON constants,
floating-point JSON values in canonical inputs, payload-hash mutations,
source-provenance changes, and semantic changes that fail exact
regeneration.

## Published conditional result

Bettin, Conrey, and Farmer introduced the log-tapered polynomial

\[
V_N(s)=
\sum_{n\le N}
\left(1-\frac{\log n}{\log N}\right)\frac{\mu(n)}{n^s}.
\]

Their asymptotic analysis is conditional: it assumes RH and a negative
moment estimate for zeta derivatives at the nontrivial zeros,

\[
\sum_{|\Im\rho|\le T}\frac1{|\zeta'(\rho)|^2}
\ll T^{3/2-\delta}
\]

for some `delta>0`. Under those hypotheses, they obtain the expected
`1/log(N)` leading scale. That theorem is useful motivation and a source
check for the coefficient family; it cannot be used as an unconditional
upper bound in an RH proof.

The artifact also cites, without reproving:

- [Baez-Duarte's strong natural-dilate criterion](https://arxiv.org/abs/math/0202141),
  which gives `d_N -> 0` if and only if RH;
- [Bettin--Conrey--Farmer](https://arxiv.org/abs/1211.5191), for this
  coefficient family and its conditional analysis;
- [Baez-Duarte--Balazard--Landreau--Saias](https://arxiv.org/abs/math/0306251),
  for the rational autocorrelation formula used by the Gram kernel.

## The unproved all-scale target

The decisive explicit-family lemma would be:

\[
\boxed{
\text{There exist }K,N_0\text{ such that }
E_N^{\mathrm{BCF}}\le\frac K{\log N}
\quad\text{for every }N\ge N_0.
}
\]

Equivalently,

\[
(\log N)^2(C_N+T_N)\le K\log N.
\]

If this were proved, then

\[
0\le d_N^2\le E_N^{\mathrm{BCF}}\longrightarrow0,
\]

and the published strong natural-dilate criterion would imply RH. A theorem
on all sufficiently large dyadic `N` would also suffice because `d_N` is
nonincreasing.

The artifact sets
`dependency_manifest.decisive_all_scale_target.proved_by_this_artifact` to
`false`. It supplies no uniform bound for the weighted prime-counting core,
no valid `psi` replacement in the tail, and no argument extending six finite
values to all scales. The Riemann Hypothesis therefore remains
`UNRESOLVED`.

## Reproduction

Generate the frozen artifact:

```powershell
.\.venv\Scripts\rh-lab.exe nyman-mobius-core-tail-audit `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1 `
  --output results\nyman-mobius-core-tail-v1.json
```

Regenerate at 256 bits and run the 512-bit same-backend containment replay:

```powershell
.\.venv\Scripts\rh-lab.exe verify-nyman-mobius-core-tail-audit `
  --artifact results\nyman-mobius-core-tail-v1.json `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1
```
