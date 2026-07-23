# Certified finite N=256 to N=512 Nyman contraction

## Verdict

The frozen audit certifies the strict finite inequality

\[
\boxed{d_{512}^2<\frac{449}{500}d_{256}^2}
\]

for the natural-dilate Nyman distances. Equivalently, the true optimal gain
in the `256 -> 512` block satisfies

\[
\frac{d_{256}^2-d_{512}^2}{d_{256}^2}>
\frac{51}{500}=0.102.
\]

This is stronger than the single finite `beta=2` target
`d_512^2 < (9/10)d_256^2`. It is not an all-scale recurrence and does not
prove or disprove the Riemann Hypothesis. The global status remains
`UNRESOLVED`.

The sign-bearing artifact is
[nyman-beta2-n512-v1.json](../results/nyman-beta2-n512-v1.json). Its payload
SHA-256 is
`2ef885119528faaffb1a67f6d8382de710a8c534d8413b621be80f5f4c4d08b6`.
The separate [candidate](../results/nyman-beta2-n512-candidate-v1.json) has
payload SHA-256
`c6be030ede77f9a6a43d2219620d9cb78d3724f44bf0bb87705e1f0722371343`
and certifies nothing by itself.

## Exact certificate chain

The frozen v1 `N=256` augmented-matrix certificate gives

\[
L_{256}<d_{256}^2\le U_{256},
\]

where

\[
L_{256}=
\frac{2801788463381838392697853210438663515}{2^{128}},
\qquad
U_{256}=
\frac{2801788463381838392697853210438663771}{2^{128}}.
\]

At 768 bits, an approximate solve proposed 512 coefficients. Each coefficient
was rounded by exact rational ties-to-even rounding onto the common grid
`2^-256`. That solve has no certificate role. Direct Arb evaluation of the
stored exact vector `c_512` certifies

\[
d_{512}^2\le E(c_{512})<U_{512},
\]

with

\[
U_{512}=
\frac{2513498367794989737282886562436186984}{2^{128}}
=0.00738650783036012675905411602281768653229\ldots.
\]

Exact `Fraction` arithmetic then gives

\[
\frac{449}{500}L_{256}-U_{512}
=
\frac{250767232190113935978562053773285247}
{34028236692093846346337460743176821145600}
>0.
\]

Therefore

\[
d_{512}^2\le E(c_{512})<U_{512}
<\frac{449}{500}L_{256}
<\frac{449}{500}d_{256}^2
<\frac9{10}d_{256}^2.
\]

The endpoint ratio supplies the slightly sharper decimal diagnostic

\[
\frac{d_{512}^2}{d_{256}^2}
<\frac{U_{512}}{L_{256}}
=0.8971049744280571747668119742\ldots,
\]

so the normalized gain is greater than
`0.1028950255719428252331880257...`. This decimal is derived from exact
endpoints; the theorem statement remains the simpler rational `51/500`
lower bound.

As a second exact route, the artifact also checks

\[
U_{512}<L_{256}-\frac{U_{256}}{10}.
\]

It follows directly that

\[
d_{256}^2-d_{512}^2
>L_{256}-U_{512}
>\frac{U_{256}}{10}
\ge\frac{d_{256}^2}{10}.
\]

## Trust separation

The computation has three deliberately different roles.

1. **Proposal.** An approximate 768-bit solve of `G_512 c=b_512` produces a
   possible vector. Exact dyadic rounding freezes it. Neither the solve nor
   its apparent energy decides a sign.
2. **Generation certificate.** A clean 768-bit kernel rebuild directly
   evaluates the exact frozen vector. A separate 257-pivot interval
   `LDL^T` replay re-establishes `L_256<d_256^2` inside the leading prefix of
   the max-512 kernel.
3. **Verification replay.** The verifier rebuilds the generation evidence and
   then repeats both decisive certificates at exactly 1536 bits without
   regenerating the coefficients. The narrower replay energy and margin must
   lie inside the generation enclosures.

No lower bound for `d_512^2` is attempted. A 513-pivot augmented lower
certificate would add substantial runtime and conditioning risk while adding
no logical strength to this contraction result.

Generation uses `python-flint==0.9.0` with FLINT 3.6.0. The verifier is a
same-backend replay, not a clean-room independent implementation. The trusted
code includes the repository's autocorrelation kernel, canonical-system
provenance, interval evaluation, fixed-order `LDL^T`, exact JSON parsing, and
hash checks.

## Frozen bindings

| Object | SHA-256 commitment |
|---|---|
| v1 summary payload | `35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf` |
| v1 index payload | `281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23` |
| v1 plan | `27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a` |
| v1 N=256 cell payload | `89ea59864d96b41f2527ce537325d62da93b7ef237dfbb576f2b979a4e74aaf9` |
| candidate coefficient vector | `b5b2e8c1acb06c9d1424f5cdbb56503f214a05459b7c3f95448fa2549291aaef` |
| 768-bit max-512 kernel | `fa1b758f9d0a0ef7457976e11f0598f27650556f39629cd8e8544eae5bd5cef8` |
| 768-bit N=256 numeric prefix | `c5b9d836f93b406339ece4c971dde01c45f655d8b58ac2f5ed3a59c7ea285676` |
| 768-bit N=512 numeric prefix | `d25dec076c37a8d1849521c807ec0253763153c4f66efc333ff24b791826b6a2` |
| fresh N=256 lower certificate | `9cf5409b271d1b7da0df7a00006ac9304455b46883102c4cfc3c51427a12a9a2` |
| direct N=512 upper certificate | `b75268e9fe3d82f3ef603a163cb7b787f2b3e695dc3ff76a3ac29da5bb06bab4` |
| 1536-bit replay max-512 kernel | `fc9b180e58b1e5fc90ba8790967450038a67bb48b6854928eeeca3ade769ba96` |
| 1536-bit replay N=256 numeric prefix | `d955b5780a5eaebf367573ebd9602fdd2d85abafb9e3bb5b3c8f6dd9e41edf41` |
| 1536-bit replay N=512 numeric prefix | `a0b2736227703e47f61f1ee33faaf0f0240aa60e7b72b6817218bc61ff300b06` |

The old max-256 core-system hash is not required to equal the leading-prefix
core-system hash copied from a max-512 system: core provenance intentionally
commits to the source dimension. The separate numeric-prefix commitment omits
that provenance distinction and does reproduce exactly.

The final verifier returned
`REPRODUCED_CERTIFIED_FINITE_NYMAN_CONTRACTION`, with both decisive
certificates replayed, both generation enclosures containing their replay
counterparts, and `same_backend_replay_only=true`.

## Schur interpretation

Split the exact finite system into old and new blocks:

\[
G_{512}=\begin{pmatrix}G&C\\C^\mathsf T&D\end{pmatrix},
\qquad
b_{512}=\binom b e.
\]

Set

\[
S=D-C^\mathsf T G^{-1}C,
\qquad
t=e-C^\mathsf T G^{-1}b.
\]

Positive definiteness gives the analytic identity

\[
d_{512}^2=d_{256}^2-t^\mathsf T S^{-1}t.
\]

Thus the certified normalized gain lower bound can also be read as

\[
t^\mathsf T S^{-1}t>\frac{51}{500}d_{256}^2.
\]

The artifact does not form an interval Schur complement or use a block solve
as evidence. The endpoint proof above is both cheaper and simpler.

## Remaining all-scale lemma

Write `E_k=d_(2^k)^2`. This artifact settles the `k=8` base more strongly
than the harmonic target:

\[
E_9<\frac{449}{500}E_8<\frac9{10}E_8.
\]

For `N=2^k`, let `P_k` project onto the old span, put
`r_k=(I-P_k)chi`, and residualize the new block by
`w_j=(I-P_k)rho_j` for `N<j<=2N`. Then

\[
(S_k)_{ij}=\langle w_i,w_j\rangle,
\qquad
(t_k)_i=\langle r_k,w_i\rangle.
\]

For an explicit full-column-rank `N x r` trial matrix `V_k`, define

\[
A_k=V_k^\mathsf T S_kV_k,
\qquad
u_k=V_k^\mathsf Tt_k.
\]

Because `S_k` is positive definite, full column rank makes `A_k` positive
definite as well.

The best gain in that trial subspace is exactly

\[
\Gamma_k(V_k)=u_k^\mathsf T A_k^{-1}u_k.
\]

It would now suffice to prove, for one explicit rule `V_k` and every
`k>=9`,

\[
\boxed{\Gamma_k(V_k)\ge\frac{E_k}{k+2}.}
\]

Then `E_(k+1)<=((k+1)/(k+2))E_k`, and the certified base gives

\[
E_k<\frac{449}{50(k+1)}E_8\longrightarrow0.
\]

Monotonicity extends the dyadic limit to all `N`, and Baez-Duarte's strong
natural-dilate criterion would imply RH. This uniform lemma is not proved.

An inverse-free sufficient form is to exhibit explicit vectors `z_k` with

\[
2u_k^\mathsf Tz_k-z_k^\mathsf TA_kz_k
\ge\frac{E_k}{k+2}
\]

for every `k>=9`.

The one-direction Mobius scout at `256 -> 512` was insufficient: its estimated
gain was about `0.074991 E_8`. That one-direction number is an exploratory
diagnostic, not a stored interval certificate. In contrast, the full optimal
gain lower bound `>0.1028950255 E_8` follows rigorously from the endpoint chain
above. A proof-oriented next basis should test the Mobius taper together with
untapered/quadratic Mobius terms, squarefree and squareful mass, and small-prime
divisibility strata. Success at this one block would still be finite evidence;
the required object is a uniform arithmetic bound.

## Reproduction

On PowerShell, from an installed checkout:

```powershell
.\.venv\Scripts\rh-lab.exe propose-nyman-beta2-n512 `
  --output results\nyman-beta2-n512-candidate-new.json

.\.venv\Scripts\rh-lab.exe nyman-beta2-n512-audit `
  --candidate results\nyman-beta2-n512-candidate-new.json `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1 `
  --output results\nyman-beta2-n512-new.json

.\.venv\Scripts\rh-lab.exe verify-nyman-beta2-n512-audit `
  --artifact results\nyman-beta2-n512-new.json `
  --candidate results\nyman-beta2-n512-candidate-new.json `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1 `
  --bits 1536
```

The proposal and generation passes each rebuild the 512-dimensional kernel.
The verifier rebuilds it at both 768 and 1536 bits. Runtime is therefore
substantial and interruption-safe only at command boundaries.
