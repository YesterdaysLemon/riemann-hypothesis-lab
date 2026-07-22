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
| 10,000 ordered critical-line Hardy-Z roots isolated at 192-bit working precision; Turing counts account for every nontrivial zeta zero below the exact final separator `T ~= 9878.2187131956` | `CERTIFIED_FINITE` | All 10,000 stored enclosures and all ten total-count separators replayed at 384 bits | Calibration only; finite checks cannot prove RH |
| Lagarias inequality certified for every `1 <= n <= 1,000,000`, equality only at `n=1`, using exact integer divisor sums and Arb balls | `CERTIFIED_FINITE` | Full range replayed at 384 bits | A counterexample would disprove RH; a positive finite prefix does not prove it |
| Proof or counterexample | `UNRESOLVED` | — | The actual objective remains open |

Artifacts: [zero certificate](results/zeros-1-10000.json),
[Lagarias certificate](results/lagarias-1-1000000.json), and the
[machine-checked claim ledger](claims/registry.json).

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
the proof-side target is an explicit density and continuity bound that closes
the infinite-dimensional limit.

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
```

Generate fresh artifacts:

```powershell
.\.venv\Scripts\rh-lab.exe zeros --first 1 --count 10000 --bits 192 `
  --block-size 1000 --output results\zeros-new.json
.\.venv\Scripts\rh-lab.exe lagarias --limit 1000000 --bits 192 `
  --output results\lagarias-new.json
```

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
- `results/`: hash-pinned generated evidence;
- `claims/registry.json`: the machine-checked claim ledger;
- `docs/`: research charter, methodology, and primary-source map;
- `tests/`: property tests, higher-precision replay, and adversarial mutations;
- `.github/workflows/ci.yml`: Linux and Windows certificate gates.

Contributions are welcome, especially negative results and adversarial tests.
Read [CONTRIBUTING.md](CONTRIBUTING.md) before proposing a headline claim.
