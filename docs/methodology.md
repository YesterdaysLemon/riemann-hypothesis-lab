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

The first frozen implementation uses `c=5/2` and modes `-4..4`. It evaluates
the desingularized archimedean source integral with Acb, cross-checks every
entry against a separately derived digamma/Lerch formula, uses an exact
prime-power transcript, and certifies positive definiteness by interval
`LDL^T`. Deleting the sole `p=2` term produces a rigorously negative mutation
control. See [the v1 specification and result](weil-matrix-v1.md).

The follow-on search engine applies that complete matrix construction to a
canonically hashed grid of exact rational cutoffs and degrees. It checkpoints
each precision attempt atomically, treats approximate eigenvectors only as
candidate generators, and rounds their exact dyadic midpoints to primitive
integer vectors. A negative is retained only when direct Arb evaluation of the
same integer vector has a strictly negative upper bound at two precisions; even
then its status is `NEGATIVE_CANDIDATE_QUARANTINED`, never `DISPROVED`. Interval
`LDL^T` is the sole positive terminal gate. Failed elimination, overlapping
eigenvalue enclosures, and exhausted precision all remain inconclusive. See
[the adaptive search plan](weil-search-plan-v1.md).

The proof target is not one larger matrix. Connes and Consani prove that the
Laurent-polynomial Fourier spaces form a core and that finite-compression
minimum eigenvalues converge to the semilocal lower bound. Our basis vectors
are in that core, so a complete negative finite vector needs no separate
smoothing bridge. The unresolved positive-side target is an analytic theorem
controlling every mode and every cutoff, not a long finite prefix. The
Davenport-Heilbronn function remains a future negative control for the
analogous hypothesis; it is not zeta and cannot disprove RH.

### 2. Nyman-Beurling-Baez-Duarte approximation

Compute certified primal and dual bounds for finite-dimensional distances,
using exact piecewise integration or interval-enclosed logarithms. Avoid normal
equations in the ill-conditioned Gram systems. The theorem target is an
explicit coefficient family with a proved error bound tending to zero; a
finite decay plot proves nothing. The frozen natural-dilate normalization,
two-sided augmented-matrix certificates, independent truncated-integral
oracle, and finite-result semantics are specified in
[Nyman natural-distance certificates v1](nyman-natural-v1.md). The post-v1
[theorem frontier](nyman-theorem-frontier.md) records the exact
Mellin/Dirichlet normalization, Burnol's asymptotic obstruction, and the
all-`N` upper-family or infinite-tail statements that would actually bridge
to a proof or disproof. The separate
[forced-rebound proof](nyman-forced-rebound-v1.md) combines the certified
`N=256` upper endpoint with the published asymptotic lower bound. Its finite
rational and Arb inequalities are reproducible, while its infinite bridge is
human-auditable published-theorem mathematics rather than a Python-checked
proof.

The [log-tapered Mobius core/tail audit](nyman-mobius-core-tail-v1.md)
evaluates the Bettin--Conrey--Farmer candidate on the same six finite
dimensions. It checks the exact divisor identities through `256`, assembles
the core by separate same-backend truncated-divisor and `psi` paths, and obtains
the tail from the complete Gram energy. The `psi` simplification is explicitly
forbidden beyond the cutoff. A 512-bit same-backend replay must lie inside all
stored 256-bit enclosures. The isolated all-scale target
`E_BCF(N)<=K/log(N)` remains unproved; finite decreases are not promoted into
an asymptotic claim.

The [finite `N=512` contraction audit](nyman-beta2-n512-v1.md) separates its
untrusted approximate solve from every sign-bearing step. The solve proposes
one exact `2^-256` dyadic vector. A clean 768-bit rebuild directly evaluates
that vector, freshly re-certifies the frozen `N=256` lower endpoint by
interval `LDL^T`, and checks all comparison thresholds with exact rational
arithmetic. Verification repeats both certificates at exactly 1536 bits
without regenerating coefficients and requires the replay energy and margin
to remain inside the generation enclosures. This certifies only
`d_512^2 < (449/500)d_256^2`. It settles the `k=8` beta=2 step but leaves the
uniform `k>=9` trial-subspace inequality wholly open.

The follow-up [arithmetic trial-subspace audit](nyman-trial-subspace-v1.md)
tests one frozen eight-column rule on the same new block. Rather than trusting
an ill-conditioned approximate Schur solve, it builds the 264-dimensional
aggregate space and certifies the strict restricted-distance lower bound
`F_V>(9/10)U_256` using all 265 positive pivots of an augmented fixed-order
interval `LDL^T`. Combined with the frozen `d_256^2<=U_256` endpoint, this
proves that the best gain in that exact span is less than `d_256^2/10`.
Generation at 768 bits and complete replay at exactly 1536 bits share the
natural-dilate/Arb backend. The result rejects only this finite ansatz; it does
not weaken the full-block contraction or establish anything at later scales.
The accompanying exact dyadic-lift and triangular-cancellation calculations
identify the remaining analytic obstruction for those explicit routes as
control of the truncated divisor tail and its aliases beyond `2N`.

The [multiscale and divisor-alias audit](nyman-multiscale-tail-v1.md) closes one
especially tempting continuation. The dyadic seed is exactly the first
Vasyunin correction studied by Báez-Duarte. The theorem in his cited 2005
preprint proves that the sequential greedy approximants diverge in weighted
`L^1`, despite
pointwise interpolation on every fixed compact interval. A repository audit
checks the exact recurrence, two-adic Möbius coefficient formula, interpolation
invariant, and rational increment multipliers through `4096`; it deliberately
does not recast finite enumeration as a proof of the infinite theorem.

The same follow-up proves an exact alias--Möbius extension lemma: an infinite
arithmetic extension can cancel all divisor aliases pointwise, but its finite
truncations differ from the intended compact step by `t S_X` on `0<t<=X`.
Consequently Hilbert convergence requires at least `S_X=o(X^-1/2)` plus
separate control beyond `X`; the prime-number-theorem limit used there gives
only `S_X=o(1)`. Structured block scouts remain exploratory. The surviving
proof target is a uniform global tail or bilinear operator estimate, not a
pointwise repair rule.

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
