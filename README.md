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

The initial release establishes a reproducible baseline for future experiments:

| Result | Status | Higher-precision consistency replay (same backend) | Meaning |
|---|---:|---:|---|
| Full `9x9` Weil matrix `A=P-R-S` for `c=5/2`, modes `-4..4`, certified positive definite by nine interval `LDL^T` pivots; the minimum among nine separated Rump eigenvalue enclosures is `[2.30606430783134e-5 +/- 4.40e-45]` | `CERTIFIED_FINITE` | All 45 upper-triangle component enclosures replayed at 384 bits; direct archimedean integrals cross-check a separate special-function formula | One positive finite compression cannot prove RH |
| 10,000 ordered critical-line Hardy-Z roots isolated at 192-bit working precision; Turing counts account for every nontrivial zeta zero below the exact final separator `T ~= 9878.2187131956` | `CERTIFIED_FINITE` | All 10,000 stored enclosures and all ten total-count separators replayed at 384 bits | Calibration only; finite checks cannot prove RH |
| Lagarias inequality certified for every `1 <= n <= 1,000,000`, equality only at `n=1`, using exact integer divisor sums and Arb balls | `CERTIFIED_FINITE` | Full range replayed at 384 bits | A counterexample would disprove RH; a positive finite prefix does not prove it |
| Proof or counterexample | `UNRESOLVED` | — | The actual objective remains open |

The Weil pipeline's deliberate corruption test deletes the only allowed prime
term, `p^r=2`. It then certifies the exact integer witness's Rayleigh quotient
as `[-0.190845397282370... +/- 1.75e-41]`. This is a successful
`CONTROL_NEGATIVE`, not evidence against RH: the corrupted matrix is not the
zeta Weil form. See the [full normalization, certificate, and limitations](docs/weil-matrix-v1.md).

Artifacts: [finite Weil certificate](results/weil-matrix-c5-over-2-n4.json),
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

The zero result is intentionally nowhere near a record. Platt and Trudgian
already proved that the lowest **12,363,153,437,138** positive-ordinate zeros
lie on the critical line through height **3,000,175,332,800**
([paper](https://arxiv.org/abs/2004.09765)). Our 10,000-zero run exists to test
the certificate plumbing before new mathematics is trusted.

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
```

The search writes atomic per-attempt checkpoints so an interrupted cell can be
resumed with `--resume`. Terminal cell files embed those ordered attempts;
therefore the redundant `attempts/` working directory is not committed.

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
