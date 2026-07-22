# Methodology

## Research portfolio

The project uses several routes intended to reduce the chance that a
normalization error or hidden assumption contaminates every result. Shared
backends and normalizations are recorded rather than treated as independent.

### 1. Weil explicit-formula positivity (primary proof/disproof track)

Construct finite-dimensional spaces of admissible compactly supported test
functions, evaluate the Weil quadratic form with exact and ball arithmetic,
and search for rational negative witnesses. We normalize `Q` as the negative of
the explicit-formula side in Bombieri's Clay convention, so RH is equivalent to
`Q(g) >= 0` for every admissible autocorrelation. A rigorously negative witness
would disprove RH only if its full archimedean term and every prime-power term
allowed by the support are included; a negative truncated sum is not decisive.
Positive finite matrices are only `CERTIFIED_FINITE`.

The proof target is not a larger matrix. It is a quantitative density and
continuity theorem that upgrades positivity on nested finite spaces to the
entire Weil test-function space. The Davenport-Heilbronn function will be used
as a negative control once this track is implemented. It is not the Riemann
zeta function and its known off-critical zeros do not disprove RH; it only tests
whether the implementation detects failure of the analogous hypothesis.

### 2. Nyman-Beurling-Baez-Duarte approximation

Compute certified primal and dual bounds for finite-dimensional distances,
using exact piecewise integration or interval-enclosed logarithms. Avoid normal
equations in the ill-conditioned Gram systems. The theorem target is an
explicit coefficient family with a proved error bound tending to zero; a
finite decay plot proves nothing.

### 3. Independent falsifiers

- Li coefficients: a single interval-certified negative coefficient disproves
  RH. Each sign must be computed by two formulas, with omitted-zero tails
  bounded without assuming RH.
- Robin/Lagarias inequalities: exact divisor sums plus Arb enclosures make a
  compact same-backend consistency verifier. Search structurally over exponent vectors; every
  pruning rule needs a theorem and a coverage proof.
- Jensen polynomials: non-hyperbolicity is a candidate falsifier only after the
  exact equivalence hypotheses and coefficient enclosures are checked.

### 4. Finite zero certification (calibration, not the main bet)

FLINT/Arb isolates Hardy-Z roots on the critical line and uses Turing's method
to prove total zeta-zero counts. FLINT embeds those roots at real part `1/2`, so
that coordinate is not an independent two-dimensional root enclosure. The
pipeline serializes ordinate enclosures as exact dyadic midpoint/radius pairs,
places separators between blocks, and asks `zeta_nzeros` for the total count
below each one. Equality between isolated line roots and the total count is the
finite-RH evidence. The critical-line API is not an off-line-zero search.

This validates the certificate machinery. It is intentionally far below the
published record and will not be marketed as mathematical progress.

### 5. Current spectral and heat-flow work

Reproduce Alain Connes's 2026 finite-prime spectral approximants, with the
missing convergence theorem displayed as the central gap. Reproduce the
de Bruijn-Newman interval before attempting any improvement. Numerical spectral
agreement or a smaller positive upper bound is not RH.

## Separation of duties

1. **Design:** freeze theorem statements, normalizations, data schemas, and
   promotion gates.
2. **Search:** generate candidates quickly; exploratory arithmetic is allowed
   but cannot emit a headline claim.
3. **Certification:** re-evaluate a candidate with exact or enclosing
   arithmetic and a pinned consistency checker whose trusted computing base is
   documented.
4. **Independent audit:** attack signs, branches, endpoint cases, completeness,
   serialization, and hidden RH assumptions.
5. **Publication:** update the README and claim registry only after the previous
   gate passes.

The baseline's shared dependencies and independence limits are enumerated in
the [trusted computing base](trusted-computing-base.md).

## Failure shields

- Never use the Euler product or Dirichlet series outside its proved
  convergence region.
- The functional equation gives symmetry, not confinement to the critical
  line.
- Never exchange a limit, sum, integral, or derivative without a stated
  domination or uniform-convergence result.
- Never use a zero table or tail estimate that assumes RH to certify RH.
- Arbitrary precision is not interval arithmetic; a displayed decimal is not
  a certificate.
- An asymptotic density of 100% can still leave infinitely many exceptions.
- A self-adjoint finite matrix matching several zero ordinates is not a
  Hilbert-Polya operator.
- Formal verification checks the supplied assumptions, so every headline
  theorem needs an axiom and assumption audit.
