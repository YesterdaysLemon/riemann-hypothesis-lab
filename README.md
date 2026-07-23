# Riemann Hypothesis Lab

[![certificate gates](https://github.com/YesterdaysLemon/riemann-hypothesis-lab/actions/workflows/ci.yml/badge.svg)](https://github.com/YesterdaysLemon/riemann-hypothesis-lab/actions/workflows/ci.yml)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A public, certificate-first research program aimed at proving or disproving the
Riemann Hypothesis (RH). Search is allowed to be speculative; conclusions are
not.

## Headline results

**Global verdict as of 2026-07-22: `UNRESOLVED`. This repository has not proved
or disproved RH.** Clay Mathematics Institute continues to list RH as an
[unsolved Millennium Prize Problem](https://www.claymath.org/millennium/riemann-hypothesis/).

Current finite results:

| Result | Status | Higher-precision consistency replay (same backend) | Meaning |
|---|---:|---:|---|
| Full `9x9` Weil matrix `A=P-R-S` for `c=5/2`, modes `-4..4`, certified positive definite by nine interval `LDL^T` pivots; the minimum among nine separated Rump eigenvalue enclosures is `[2.30606430783134e-5 +/- 4.40e-45]` | `CERTIFIED_FINITE` | All 45 upper-triangle component enclosures replayed at 384 bits; direct archimedean integrals cross-check a separate special-function formula | One positive finite compression cannot prove RH |
| 81 complete Weil matrices around the prime-power transitions `q=7,8,9`, degrees `12..24`; every cell certified positive definite by interval `LDL^T` | `EXPLORATORY` | All 81 cells, 108 attempts, and 905 stored candidate evaluations replayed; 9 parity and 6 degree-nesting audits also passed | A bounded null search: no negative witness was found, and positive finite compressions cannot prove RH |
| Six natural-dilate Nyman distances for `N=8,16,32,64,128,256`; each has an exact bracket `L_N < d_N^2 <= U_N` of width `2^-120`, and all five exact comparisons `U_2N < L_N` pass | `EXPLORATORY` | A separate same-backend replay rebuilt the full kernel at 1536 bits and re-certified all six stored dyadic candidates and all 510 interval `LDL^T` pivots | Finite decay cannot establish the required limit `d_N -> 0` |
| 10,000 ordered critical-line Hardy-Z roots isolated at 192-bit working precision; Turing counts account for every nontrivial zeta zero below the exact final separator `T ~= 9878.2187131956` | `CERTIFIED_FINITE` | All 10,000 stored enclosures and all ten total-count separators replayed at 384 bits | Calibration only; finite checks cannot prove RH |
| Lagarias inequality certified for every `1 <= n <= 1,000,000`, equality only at `n=1`, using exact integer divisor sums and Arb balls | `CERTIFIED_FINITE` | Full range replayed at 384 bits | A counterexample would disprove RH; a positive finite prefix does not prove it |
| Proof or counterexample | `UNRESOLVED` | — | The actual objective remains open |

The Weil pipeline's deliberate corruption test deletes the only allowed prime
term, `p^r=2`. It then certifies the exact integer witness's Rayleigh quotient
as `[-0.190845397282370... +/- 1.75e-41]`. This is a successful
`CONTROL_NEGATIVE`, not evidence against RH: the corrupted matrix is not the
zeta Weil form. See the [full normalization, certificate, and limitations](docs/weil-matrix-v1.md).

Artifacts: [finite Weil certificate](results/weil-matrix-c5-over-2-n4.json),
[transition-v2 summary](results/weil-transition-q7-q9-v2-summary.json),
[transition-v2 release manifest](results/weil-transition-q7-q9-v2.release.json),
[Nyman v1 summary](results/nyman-natural-v1-summary.json),
[Nyman normalization audit](results/nyman-natural-v1-normalization.json),
[zero certificate](results/zeros-1-10000.json),
[Lagarias certificate](results/lagarias-1-1000000.json), and the
[machine-checked claim ledger](claims/registry.json).

### Adaptive finite Weil search (`EXPLORATORY`)

The frozen [v1 search plan](plans/weil-grid-v1.json) evaluated all 32 pairs of
exact cutoffs `c in {3/2, 2, 5/2, 3, 4, 5, 8, 9}` and degrees
`N in {4, 8, 12, 16}`. Every complete matrix closed
`FINITE_POSITIVE_CERTIFIED` by interval `LDL^T`: 23 cells at 96 bits, six at
192 bits, and the three hardest cells at 384 bits. There were no inconclusive
terminal cells and no interval-negative witnesses. At lower precisions where
`LDL^T` had not yet closed, the engine performed 443 primitive integer-vector
candidate evaluations: 370 were strictly positive and 73 were
interval-inconclusive at the attempted precision; none was negative.

The most ill-conditioned terminal cell was `c=9, N=16` (dimension 33). Its
secondary Rump diagnostic enclosed the smallest eigenvalue as
`[1.483861410967958711918628391031090118348e-29 +/- 1.47e-69]`; the actual
terminal certificate is the complete list of 33 positive interval `LDL^T`
pivots. A separate process replayed all 32 cell artifacts successfully. The
[search index](results/weil-search-grid-v1/index.json) has payload SHA-256
`f0595f9630b8738d518dfabfeb2d4d9921f4e0fcb0d4d09a9c201fe6f47e7d87` and
binds plan SHA-256
`68afe80851df02f1db4eb0b8119818bf45e9afef2d1ade7de4f73200f32827fc`.

This is a bounded null search, not evidence sufficient for RH: it says only
that no negative direction was found in these 32 finite Fourier
compressions. The aggregate remains `EXPLORATORY`, does not enter the global
claim ledger, and leaves the headline verdict `UNRESOLVED`. The adaptive
engine, quarantine gates, and prime-power transition frontier are documented
in the [v1 search methodology](docs/weil-search-plan-v1.md).

### Prime-power transition search v2 (`COMPLETED`, `EXPLORATORY`)

The frozen experiment evaluated all 81 explicitly hashed matrices around the
prime-power cutoffs `q=7,8,9`. For each transition it used the exact point
`c=q` and four log-symmetric rational pairs
`c_-=q*2^j/(2^j+1)`, `c_+=q*(2^j+1)/2^j`, with
`j in {4,6,8,10}`. Each cutoff was evaluated at three scheduled degrees; the
largest matrix has degree 24 and dimension 49. The canonical plan hash is
`33185881cc619fda99b7835f5e71422b0c5cf745cb5411162638a8e08dd05336`.

Every cell closed `FINITE_POSITIVE_CERTIFIED`: 9 at 192 bits and 72 at
384 bits. There were no terminally inconclusive cells, quarantined negative
candidates, or conflicting negative observations. The run used 108 ordered
attempts and stored 905 candidate-evaluation records. A separate same-backend
process replayed all 81 cells, all 108 attempts, and all 905 evaluations.
All 9 reversal-parity audits and all 6 adjacent degree-nesting audits also
reproduced strictly. The aggregate conclusion is
`NO_CERTIFIED_NEGATIVE_WITNESS_FOUND_IN_FINITE_TRANSITION_BATCH`.

The secondary Rump diagnostic covered 107 of 108 search attempts and 3,963
stored eigenvalue enclosures. One failed spectrum construction at
`q-9-j-4-above-n-16`, attempt 0 at 192 bits, is retained explicitly; the
cell subsequently closed at higher precision. Rump spectra are diagnostic
only. The interval `LDL^T` pivots, not the eigensolver, certify each positive
finite matrix.

V2 adds exact transition-contract validation, per-cell precision policies,
dedicated 768- or 1536-bit candidate replay, reversal-parity audits, and
degree-nesting/zero-padding audits. It completes the full bounded plan even if
a finite negative candidate appears. Approximate eigenvectors and Rump
spectra remain diagnostic only; interval `LDL^T` and direct evaluation of an
exact primitive integer vector are the sign-bearing computations. The frozen
[plan](plans/weil-transition-q7-q9-v2.json),
[methodology](docs/weil-transition-v2.md), and
[compact summary](results/weil-transition-q7-q9-v2-summary.json) are tracked.
The summary has payload SHA-256
`32a08241f112fddfe496bd10934bf03760b0bde939b7f91938fd08475f701467`;
the index it binds has payload SHA-256
`b1d00edb4dd3115992e1e377c7cf81475824d16d175396a40fd2cadb7c03d3ff`.

The complete 99-file evidence set is published in the
[v2 evidence release](https://github.com/YesterdaysLemon/riemann-hypothesis-lab/releases/tag/weil-transition-q7-q9-v2)
as a deterministic 200,980,480-byte USTAR archive. Its SHA-256 is
`9e6a76c5328ce2e6763da642bf16ad1982152ba21780ca48500c3447b7c961dc`;
the tracked [release manifest](results/weil-transition-q7-q9-v2.release.json)
binds the plan, index, summary, every cell, and all 15 audits.

This result cannot turn 81 positive matrices into a proof. A negative finite
direction would remain quarantined until every promotion gate and an
independent implementation reproduce the bridge back to the full Weil
criterion. No such direction was found, and the global verdict remains
`UNRESOLVED`.

The zero result is intentionally nowhere near a record. Platt and Trudgian
already proved that the lowest **12,363,153,437,138** positive-ordinate zeros
lie on the critical line through height **3,000,175,332,800**
([paper](https://arxiv.org/abs/2004.09765)). Our 10,000-zero run exists to test
the certificate plumbing before new mathematics is trusted.

### Natural-dilate Nyman distance v1 (`COMPLETED`, `EXPLORATORY`)

For `rho_a(x)={1/(a*x)}` and `chi=1_(0,1]`, define

```text
d_N^2 = inf_c ||chi - sum_(a=1)^N c_a rho_a||_2^2.
```

Baez-Duarte's strong Nyman--Beurling theorem says RH is equivalent to
`d_N -> 0`. The frozen [v1 plan](plans/nyman-natural-v1.json) uses one
provenance-bound `N=256` Gram kernel at 768 bits, exact `2^-256` dyadic
coefficients, direct interval evaluation of the primal error, and positive
fixed-order interval `LDL^T` for the augmented lower-bound matrix. Approximate
linear solves propose coefficients but never decide a sign.

All six cells closed with the following shared decimal prefixes. Every exact
bracket has width `2^-120`; the linked summary stores the full dyadic endpoints.

| `N` | Certified `d_N^2` bracket's shared decimal prefix | Terminal status |
|---:|---:|---|
| 8 | `0.02416142158589668502280894301085133...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |
| 16 | `0.01789402347696943509905506314795876...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |
| 32 | `0.01405194369952985983642504072212329...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |
| 64 | `0.01137604029967365814705897002932676...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |
| 128 | `0.00965854927811190991003539257145257...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |
| 256 | `0.00823371627726102144159803810263543...` | `FINITE_DISTANCE_BRACKET_CERTIFIED` |

The run stores 504 exact coefficient numerators and 510 positive lower-bound
pivots. All five separate diagnostics `U_2N < L_N` pass. The canonical plan
payload SHA-256 is
`27ebd77ca7bfae0ff898bd33fcd739c1c069a79b23b9bff4f03f84ac91f1540a`;
the generation-kernel payload SHA-256 is
`6fb54e7f63208f1bd81fb7ef2e592a8a655458dd7744081748590159b868cc26`;
the index payload SHA-256 is
`281f12c122897407d169e9830871ebd4764f86f03347e28dbc0cfa5b39876d23`;
and the compact summary payload SHA-256 is
`35cf625bd2ff70def7c440065aae20e375691c483666b2411098a4136ec399cf`.

A separate same-backend replay rebuilt the complete kernel at 1536 bits and
re-certified all six stored candidates without regenerating the approximate
solves. Its replay-kernel payload SHA-256 is
`0b914bd5902f119ace638e74947225bce3bc541f719afee7a313d8390d68dae8`.
The independent piecewise-integral/harmonic normalization bundle passed all 14
frozen rational checks and has payload SHA-256
`021a060fb3eb0121c325a0af13cc39095c0c6f26bc5a4b235269611a66a9442c`.
That oracle shares Arb and its generic tail resolves only about 12 bits; it is
not a clean-room backend.

These finite values do **not** prove RH. Burnol's unconditional asymptotic
lower obstruction also rules out an eventual fixed-factor contraction under
`N -> 2N`; the next proof-side target must be an explicit all-`N` upper family
with a weakening, logarithmic-scale contraction. The exact normalization,
known conditional families, and three valid infinite bridges are recorded in
the [post-v1 theorem frontier](docs/nyman-theorem-frontier.md). Full formulas,
certificate gates, oracle scope, and limitations are in the
[v1 methodology](docs/nyman-natural-v1.md).

## Research strategy

The main bet is **Weil explicit-formula positivity**. We define `Q` as the
negative of the explicit-formula side in Bombieri's Clay convention, so RH is
equivalent to `Q(g) >= 0` for every admissible autocorrelation. We will search
nested finite test-function spaces for a rational, interval-certified `Q(g)<0`
witness. Such a witness is decisive only when the full archimedean term and
every prime-power contribution allowed by the support are included. A negative
truncation is not a disproof. Positive finite matrices remain finite evidence;
the proof-side target is a theorem controlling every Fourier mode and every
cutoff. The published semilocal core theorem identifies the finite-compression
limit, but no single computed matrix controls it.

The first frozen matrix now implements that program at `c=5/2`, degree four.
It uses exact prime-power cutoff arithmetic, segmented Acb integration of the
desingularized archimedean term, a separate digamma/Lerch closed-form oracle,
interval `LDL^T`, and a known-negative omission control. Its positive result is
strictly finite and leaves all higher modes and other cutoffs untouched.

Multiple tracks are intended to reduce the risk that one hidden assumption
contaminates every result:

1. Weil quadratic forms and rational negative witnesses;
2. Nyman-Beurling-Baez-Duarte primal/dual approximation certificates;
3. Li-coefficient signs computed by two formulas with rigorous tail bounds;
4. Robin/Lagarias structural counterexample search;
5. de Bruijn-Newman heat-flow reproduction and Connes's 2026 finite spectral
   approximants, with their missing bridges stated rather than hand-waved; and
6. modest direct zero certification as a regression oracle.

The full portfolio, failure shields, and source map live in
[the methodology](docs/methodology.md) and
[the literature map](docs/literature-map.md). The exact shared dependencies are
listed in the [trusted computing base](docs/trusted-computing-base.md).

## Confirmation standard

Under this repository's internal standard, a global `PROVED` or `DISPROVED`
headline is blocked unless the result has:

- a complete theorem or off-line-zero certificate;
- an explicit assumption and axiom manifest;
- a pinned machine checker with no placeholders or unproved custom axioms;
- at least two independent clean-room reproductions; and
- adversarial review of the exact bridge back to RH.

The schema/hash gate enforces file binding, checker manifests, cross-linked
resolution claims, and distinct reproduction backends. It cannot establish
mathematical independence or correctness by itself. See the
[research charter](docs/research-charter.md) for the complete internal standard.
Official Clay recognition is separate and governed by the
[CMI rules](https://www.claymath.org/millennium-problems/rules/), including
publication, a waiting period, and general acceptance.

## Reproduce the baseline

Python 3.11+ and `python-flint==0.9.0` are required. On PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\rh-lab.exe verify-claims
.\.venv\Scripts\rh-lab.exe verify-zeros `
  --artifact results\zeros-1-10000.json --bits 384
.\.venv\Scripts\rh-lab.exe verify-lagarias `
  --artifact results\lagarias-1-1000000.json
.\.venv\Scripts\rh-lab.exe verify-weil `
  --artifact results\weil-matrix-c5-over-2-n4.json --bits 384
.\.venv\Scripts\rh-lab.exe verify-weil-search `
  --index results\weil-search-grid-v1\index.json `
  --checkpoint-dir results\weil-search-grid-v1
.\.venv\Scripts\rh-lab.exe verify-nyman-search `
  --index results\nyman-natural-v1\index.json `
  --checkpoint-dir results\nyman-natural-v1 --bits 1536
.\.venv\Scripts\rh-lab.exe verify-nyman-summary `
  --summary results\nyman-natural-v1-summary.json `
  --checkpoint-dir results\nyman-natural-v1
.\.venv\Scripts\rh-lab.exe verify-nyman-normalization-audit `
  --artifact results\nyman-natural-v1-normalization.json
```

Generate fresh artifacts:

```powershell
.\.venv\Scripts\rh-lab.exe zeros --first 1 --count 10000 --bits 192 `
  --block-size 1000 --output results\zeros-new.json
.\.venv\Scripts\rh-lab.exe lagarias --limit 1000000 --bits 192 `
  --output results\lagarias-new.json
.\.venv\Scripts\rh-lab.exe weil --bits 192 `
  --output results\weil-matrix-new.json
.\.venv\Scripts\rh-lab.exe weil-search `
  --plan plans\weil-grid-v1.json `
  --checkpoint-dir results\weil-search-grid-v1-new
.\.venv\Scripts\rh-lab.exe weil-transition-search `
  --plan plans\weil-transition-q7-q9-v2.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2
.\.venv\Scripts\rh-lab.exe nyman-search `
  --plan plans\nyman-natural-v1.json `
  --checkpoint-dir results\nyman-natural-v1-new
.\.venv\Scripts\rh-lab.exe nyman-normalization-audit `
  --output results\nyman-natural-v1-normalization-new.json
```

Both search engines write atomic per-attempt checkpoints. Resume the v2 batch
without changing its frozen plan:

```powershell
.\.venv\Scripts\rh-lab.exe weil-transition-search `
  --plan plans\weil-transition-q7-q9-v2.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2 --resume
.\.venv\Scripts\rh-lab.exe verify-weil-transition-search `
  --index results\weil-transition-q7-q9-v2\index.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2
.\.venv\Scripts\rh-lab.exe verify-weil-transition-summary `
  --summary results\weil-transition-q7-q9-v2-summary.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2
```

Resume and summarize the frozen Nyman batch without changing its plan:

```powershell
.\.venv\Scripts\rh-lab.exe nyman-search `
  --plan plans\nyman-natural-v1.json `
  --checkpoint-dir results\nyman-natural-v1 --resume
.\.venv\Scripts\rh-lab.exe verify-nyman-search `
  --index results\nyman-natural-v1\index.json `
  --checkpoint-dir results\nyman-natural-v1 --bits 1536
.\.venv\Scripts\rh-lab.exe summarize-nyman `
  --checkpoint-dir results\nyman-natural-v1 `
  --output results\nyman-natural-v1-summary-new.json
```

Weil terminal cell files embed their ordered attempts; the redundant
`attempts/` working directories are not published as primary evidence. The
v2 Weil and v1 Nyman batches are substantial serial computations because
FLINT's precision context is shared.
Parity and nesting audits have their own `weil-parity-audit`,
`weil-nesting-audit`, and corresponding `verify-...` commands; run
`rh-lab <command> --help` for exact witness and precision options.

`python-flint` wraps FLINT/Arb ball arithmetic. FLINT's critical-line routine
isolates Hardy-Z roots and embeds them at real part `1/2`; the decisive finite
evidence is that Turing's total count matches those roots through the separator.
That routine is not an off-line-zero search and may fail to terminate if its
critical-line isolation cannot close. See the
[FLINT documentation](https://flintlib.org/doc/acb_dirichlet.html). Higher-
precision replay uses the same backend, so it is a consistency reproduction,
not an independent implementation.

## Repository map

- `src/riemann_lab/`: certificate generators, exact-dyadic ball encodings,
  same-backend replay verifiers, and schema/hash claim gates;
- `plans/`: canonically hashed exploratory search schedules;
- `results/`: hash-pinned generated evidence;
- `claims/registry.json`: the machine-checked claim ledger;
- `docs/`: research charter, methodology, and primary-source map;
- `tests/`: property tests, higher-precision replay, and adversarial mutations;
- `.github/workflows/ci.yml`: Linux and Windows certificate gates.

Contributions are welcome, especially negative results and adversarial tests.
Read [CONTRIBUTING.md](CONTRIBUTING.md) before proposing a headline claim.
