# Trusted computing base

The baseline artifacts are rigorous relative to the following software and
mathematical contracts. They are not clean-room independent proofs.

## Zero certificate

- CPython and this repository's parsing, hashing, schedule, and enclosure code;
- `python-flint==0.9.0` and its bundled FLINT 3.6.0;
- FLINT's Hardy-Z isolation and Turing zero-count implementation; and
- the published theorems behind those FLINT routines.

Generation and replay both use FLINT. Replay at higher precision catches
corrupted transcripts, missing count blocks, inconsistent separators, and many
precision failures, but it cannot expose a bug shared by the same backend. The
critical-line routine constructs Hardy-Z roots and embeds the real coordinate
at `1/2`; it is not a search for off-line roots. A disproof search needs a
separate two-dimensional argument-principle implementation.

## Lagarias certificate

- CPython arbitrary-precision integers and the linear divisor-sum sieve;
- `python-flint==0.9.0` / FLINT 3.6.0 Arb operations for harmonic, exponential,
  logarithmic, and comparison enclosures; and
- Lagarias's published equivalence theorem.

The verifier regenerates the entire artifact at its stated precision and then
reruns the whole range at doubled precision. It therefore rejects contradictory
human-readable fields as well as numerical changes, but both passes share the
same code. The test suite also compares all one million release values against a
separately implemented divisor-addition sieve.

## Finite Weil matrix certificate

- CPython exact integers and rational cutoff comparisons;
- this repository's kernel, pole, archimedean, prime-power, interval `LDL^T`,
  mutation-control, serialization, and replay code;
- `python-flint==0.9.0` / FLINT 3.6.0 Arb and Acb arithmetic, special functions,
  certified quadrature, and the secondary Rump eigenvalue enclosure;
- Bombieri's Weil-functional formulation and the frozen finite-cutoff formulas
  of Connes--Consani--Moscovici.

The 45 archimedean entries are supplied by segmented Acb integration of a
desingularized source integrand and cross-checked against a separately coded
digamma/trigamma/Lerch expression. Pole entries are cross-checked against their
defining integrals. Interval `LDL^T` is the positive-definiteness certificate;
Rump eigenvalues are a regression oracle. Both computations, both precisions,
and both entry formulas still share FLINT and repository code, so they are not
clean-room independent. Connes and Consani's published core theorem supplies
the domain and finite-compression limit bridge, but this one matrix leaves all
higher modes and other cutoffs untested.

## Exploratory finite Weil search

The adaptive grid inherits the entire finite Weil trusted computing base. Its
additional trusted code covers plan canonicalization, exact ties-to-even
rationalization of approximate eigenvectors, primitive integer-vector
evaluation, precision scheduling, atomic checkpoints, and artifact replay.
Approximate eigensolvers and Rump spectra are expressly untrusted hints:
neither can establish a terminal sign. Positive cells close only through
interval `LDL^T`; negative candidates require the same exact integer witness
to have a strictly negative interval upper bound at two precisions.

Arb quadrature is permitted to return slightly different valid enclosures on
repeated calls. Before hashing, the engine therefore outward-rounds every ball
onto a guarded deterministic dyadic grid. Replay requires the same stabilized
evidence, exact input structure and prime-power transcript, and higher-
precision containment where a candidate is confirmed. Generation and
verification still share this repository, CPython, python-flint, FLINT, and
Arb, so a quarantined negative would require a separate implementation and
mathematical audit before it could support any global claim. Every aggregate
search artifact remains `EXPLORATORY` and leaves the hypothesis status
`UNRESOLVED`.

## Exploratory Nyman natural-distance certificates

The frozen natural-dilate experiment trusts CPython exact integers and
`Fraction` arithmetic for rational reduction, breakpoint construction, and
dyadic encodings; this repository's autocorrelation, shared-prefix,
direct-energy, augmented-matrix, interval `LDL^T`, checkpoint, hashing, and
replay code; and `python-flint==0.9.0` / FLINT 3.6.0 Arb arithmetic for
logarithms, cotangents, constants, linear solves, and interval operations.
Baez-Duarte's strong Nyman--Beurling theorem supplies the bridge from the
infinite natural-dilate closure problem to RH.

Approximate linear solves only propose coefficient vectors. The upper gate
directly evaluates `E(c)` on the stored exact dyadic coefficients, and the
lower gate requires every fixed-order interval `LDL^T` pivot of the augmented
matrix to be strictly positive. Canonical shared systems and their prefixes
carry provenance and content hashes, so ordinary public construction,
copying, or serialized-artifact mutation cannot substitute an arbitrary
symmetric matrix into either certificate path. Hostile in-process Python that
calls private constructors or monkeypatches trusted code remains inside the
declared repository-code trust boundary. Generation and replay use clean FLINT
cache boundaries to make same-precision evidence deterministic.

The separate normalization oracle uses exact rational breakpoints and a
piecewise antiderivative rather than the Vasyunin-sum implementation. Its
generic positive tail bound is only `1/T`, so the frozen `T=4096` audit checks
roughly twelve bits of normalization, not the full precision of the closed
formula. It still shares CPython, Arb logarithms, serialization code, and this
repository, and therefore is not a clean-room backend. Higher-precision replay
also shares the same implementation and FLINT. Six certified finite brackets,
even if strictly decreasing, cannot establish the limiting statement required
to prove or disprove RH.

The tracked compact summary is a reporting and integrity layer only. It
rechecks canonical schemas, exact dyadics, semantic cross-links, and all stored
payload commitments, but deliberately does not rebuild a kernel or rerun an
interval certificate. The `verify-nyman-search` command is the numerical
reproduction gate.

### Forced-rebound audit

The rebound artifact first regenerates the frozen Nyman summary from its
checkpoint tree. It then trusts CPython `Fraction` arithmetic for the five
dyadic scaled-decline comparisons and `python-flint==0.9.0` / FLINT 3.6.0 Arb
at a clean 256-bit boundary for three scalar signs. Strict JSON
loading rejects duplicate keys, nonstandard constants, floating-point values,
unknown fields after regeneration, and payload mutations.

That checker certifies only the finite chain
`d_256^2 log(256) < 23/500 < 2+gamma-log(4*pi)` and the five preceding exact
scaled declines, plus the lower-endpoint comparison showing that the
`beta=2` anchor ceiling is not forced below the asymptotic floor. The
conclusion that a later dyadic scaled rebound exists
uses the published strong natural-dilate criterion and asymptotic lower-bound
theory. Those analytic dependencies are recorded in the artifact but are not
formalized or reproved by the Python verifier. The artifact therefore remains
`EXPLORATORY`, and RH remains `UNRESOLVED`.

### Log-tapered Mobius core/tail audit

The Mobius artifact first regenerates the frozen six-cell Nyman summary and
binds its exact optimized-distance brackets. That summary operation is an
integrity check, not a numerical replay of the old optimizer certificates.
The new computation independently rebuilds the canonical natural-dilate Gram
kernel at 256 bits, evaluates the Bettin--Conrey--Farmer coefficient formula,
and recomputes every stored scalar at 512 bits. Generation and replay share
CPython, this repository, python-flint, FLINT, and the rational Vasyunin kernel,
so the replay is not a clean-room backend.

CPython integer and `Fraction` arithmetic checks the Mobius divisor sums and
the rational products encoding
`sum_(d|k) mu(d)log(d)=-Lambda(k)` for every `k<=256`. Arb evaluates the
formula-defined logarithmic coefficients, both finite core expressions, the
full Gram energy, their comparisons, and the tail obtained by subtraction.
The mathematical identity `T_N=E_N-C_N` is part of the manual proof layer;
the checker does not directly sum the infinite tail. Requiring the two core
enclosures to overlap guards the divisor/prime-counting translation, but both
paths still share code and Arb logarithms.

The shortcut `q_m=-psi(m)/log(N)` is used only in the core, where `m<N`.
Beyond the cutoff the divisor sum is truncated and that identity generally
fails. The artifact certifies six finite energies and strict comparisons only.
It neither proves the required all-scale upper estimate nor imports the
conditional Bettin--Conrey--Farmer asymptotic as an unconditional fact; its
classification is `EXPLORATORY` and RH remains `UNRESOLVED`.

### Finite N=512 beta=2 contraction audit

The candidate generator trusts the same canonical natural-dilate kernel and
FLINT/Arb backend as v1. Its 768-bit approximate solve is explicitly untrusted:
it only proposes 512 coefficients, which are rounded by exact rational
ties-to-even arithmetic onto a common `2^-256` grid. The standalone candidate
is `EXPLORATORY` and certifies nothing.

The `CERTIFIED_FINITE` audit rebuilds one max-512 kernel at 768 bits. It
directly evaluates the stored exact vector against a frozen `2^-128` upper
endpoint and separately replays the v1 `N=256` lower endpoint using all 257
positive pivots of the augmented fixed-order interval `LDL^T`. CPython
`Fraction` arithmetic checks the strict comparisons
`U_512 < (449/500)L_256`, `U_512 < (9/10)L_256`, and
`U_512 < L_256-U_256/10`. No approximate solve, Schur-complement solve, or
displayed decimal decides a sign.

Verification first rejects any artifact whose frozen payload, source hashes,
candidate, policy, or theorem bridge has changed. It then regenerates the
768-bit audit and repeats both decisive certificates at exactly 1536 bits
without regenerating coefficients. The 1536-bit energy and upper-margin
enclosures must lie inside their 768-bit counterparts. Both precisions share
CPython, repository code, python-flint, FLINT, and the same analytic kernel;
this is not a clean-room independent reproduction.

The artifact supplies no `N=512` lower bound and proves only one finite
transition. The Schur identity translates it into a normalized gain greater
than `51/500`, but the uniform gain inequality for every later dyadic scale is
unproved. Nothing in this audit resolves RH.

### Arithmetic trial-subspace rejection audit

The trial audit trusts the same CPython, repository code,
`python-flint==0.9.0`, FLINT 3.6.0, Arb logarithms, and canonical
natural-dilate kernel as the preceding Nyman audits. CPython exact integer and
`Fraction` arithmetic constructs the Moebius and divisibility columns, exact
rational thresholds, and the theorem bridge. Arb encloses the logarithmic
columns, aggregate Gram entries, and every interval-elimination operation.

The exploratory Schur solve is not part of the artifact and decides no sign.
The sign-bearing route aggregates the old 256 natural dilates and eight frozen
new-block functions, then proves positivity of the 265-by-265 matrix
`[[G_V,-b_V],[-b_V^T,1-(9/10)U_256]]` by 265 positive fixed-order interval
`LDL^T` pivots. This directly certifies `F_V>(9/10)U_256`; exact rational
arithmetic and the source endpoint `d_256^2<=U_256` then give
`Gamma(V)<d_256^2/10`.

The artifact freezes the ordered coefficient matrix, aggregate system, lower
certificate, source references, and complete payload. Strict loading rejects
duplicate keys, nonstandard constants, floats, cycles, values that cannot be
canonically serialized, and unknown structure. Verification regenerates the
768-bit artifact and then rebuilds the kernel, basis, aggregate system, and all
265 pivots at exactly 1536 bits. Both passes use the same backend and
implementation, so they are consistency checks rather than independent
clean-room reproductions. The certificate rejects one eight-dimensional
finite space only; it does not bound other trial spaces or resolve RH.

### First Vasyunin-correction finite audit

The finite prefix artifact trusts CPython arbitrary-precision integers and
`Fraction`; this repository's greedy recurrence, integer floor seed, Möbius
factorization, two-adic decomposition, canonical JSON, hashing, and exact
regeneration code; and the induction that the binary seed is constant on each
unit interval. It does not use FLINT, Arb, floating-point arithmetic, or an
approximate solve.

The verifier strictly rejects duplicate JSON keys, nonstandard constants,
floating-point values, unknown fields, malformed canonical integers, payload
mutations, and any record that does not exactly regenerate at the declared
limit. Generation and verification nevertheless use the same Python
implementation, so replay is a consistency check rather than a clean-room
formal proof.

Báez-Duarte's coefficient formula and divergence theorem are external
mathematical dependencies from the cited 2005 preprint. The checker confirms
their recurrence and rational increment factors only through `n=4096`; it does
not prove the weighted
integral `||h_n||_1=log(2)/n`, extend the coefficient formula to all integers,
or derive infinite `L^1` divergence from a finite prefix. That cited theorem
rejects the first sequential greedy correction only. It does not reject every
Vasyunin seed or multiscale Nyman construction and does not resolve RH.

### Sharp alias-truncation finite audit

The sharp alias artifact trusts CPython arbitrary-precision integers and
`Fraction`; this repository's linear Möbius sieve, finite Dirichlet
convolution and divisor-sum implementations, canonical JSON, hashing, and
exact regeneration code. It does not use FLINT, Arb, floating point, or an
approximate solve. The frozen shell is `y_9=1, y_16=-1`, and every arithmetic
identity is checked only through `X=4096`.

The verifier rejects duplicate keys, nonstandard constants, floats, unknown
fields, malformed canonical limits, payload mutations, and any tree that does
not exactly regenerate. Generation and verification share the same Python
implementation. The three vector hashes bind the Möbius, alias-coefficient,
and divisor-inversion prefixes, but this remains a consistency replay rather
than a clean-room implementation.

The theorem that `sum_(n<=X)(mu*y)(n)/n` is not `o(X^-1/2)` for every nonzero
finite `y` is a separate human-auditable analytic bridge. The finite checker
does not prove Mellin continuation, Jensen's exponential-polynomial zero
count, Conrey's positive-proportion theorem for simple critical-line zeros,
or the boundary Abelian contradiction. The theorem rejects fixed-shell sharp
truncation only; it does not cover scale-dependent regularization or resolve
RH.

### Balanced-multiplier exploratory scout

The balanced-multiplier scout is not a certificate path. It trusts the stored
exact-dyadic Nyman candidate files, CPython binary64 conversion, NumPy 2.3.5,
the platform BLAS/LAPACK implementation used for the finite normal solve, and
vectorized floating-point divisor accumulation. The reconstructed shell,
quadratic coefficients, solution, interval sums, and displayed convergence
tables are approximate. Reported energy fractions use the midpoint of each
stored candidate's certified energy bracket, not a fresh exact energy
evaluation.

The time-domain identity for the interval error and the formulas for signed
direct gain are human-auditable exact mathematics. Their large-scale numerical
evaluation is truncated. The missing alias-norm tail is positive, but the
missing old-residual cross term has no fixed sign; a truncated direct gain is
therefore neither a certified lower nor a certified upper bound.

The separate small-scale full-Gram helper rounds the shell and multiplier to
exact dyadics, reconstructs harmonic balance exactly, and evaluates a complete
finite coefficient vector with the existing Arb autocorrelation kernel. It is
an independent consistency check for selected rounded vectors, not a verifier
for the large-scale binary64 grid. None of those exploratory grid values is
entered in the claim ledger, and no uniform contraction or implication for RH
is claimed.

### Balanced-multiplier complete-tail certificate

The `CERTIFIED_FINITE` complete-tail path is separate from the exploratory
grid. The grid selects one `N=256`, `K=16` direction, but every sign-bearing
coefficient is then made explicit and exact: the old vector uses a `2^-16`
grid, the shell and multiplier use `2^-9`, and their exact convolution uses
`2^-18`. Python-FLINT 0.9.0 / Arb proves the first 255 shell rounding bins;
the last shell coefficient is imposed by exact zero-sum balance. Arb also
encloses the compressed logarithm and log-gamma evaluations. CPython integers
and `Fraction` perform the balance checks, divisor aggregation, Jordan-`J_2`
means, ordered-pair LCM bounds, tail center, and tail radius exactly.

The first `2^26` interval numerators are computed with signed int64 after
global and blockwise overflow proofs. NumPy 2.3.5 performs only the final
binary64 term divisions and reductions under the required IEEE-754
round-to-nearest-ties-to-even contract. These operations are not assumed
exact: the checker derives an exact rational error radius from `u=2^-53`, the
absolute computed term mass, and the within-block and final-block reduction
bounds. A conservative dyadic-lattice gate rules out subnormal terms or
partial sums. The verifier changes the block size from `2^20` to `2^19` and
raises Arb precision from 256 to 384 bits.

Generation and replay still share the same repository formulas, CPython,
NumPy, Python-FLINT, and FLINT, so this is not a clean-room implementation.
The test suite separately brute-forces small exact covariance, LCM, prefix,
logarithmic-compression, and tail cases and rejects mutations of the source
hash chain. The artifact proves only that one explicit correction has complete
direct gain greater than `7/50000`; it supplies no all-scale contraction and
does not resolve RH.

### Fourier/Farey large-sieve tail certificate

The large-sieve certificate reuses the exact source-vector construction,
integer prefix recurrence, binary64 error analysis, and Arb enclosure boundary
above, but its prefix ends at `2^22`. Generation uses `2^18`-sized blocks and
256-bit Arb; replay changes to `2^17`-sized blocks and 384-bit Arb. The stored
generation balls must contain the replay balls.

The omitted periodic tail is recomputed with exact `Fraction` arithmetic.
Jordan-`J_2` divisor sums give the old, added, cross, and new periodic mean
squares. The coefficient support gives the exact Farey reciprocal-spacing
bound `Q(Q-1)`, and the direct-gain difference-of-squares identity gives the
partial-sum constant `Q(Q-1)(rho+tau)`. No numerical Fourier transform, common
period, square root, or LCM pair sum is trusted by this path. Small tests
enumerate every cyclic interval for independent exact examples.

The analytic bridge is the dual consecutive-interval form of Montgomery and
Vaughan's large-sieve theorem, plus the human-auditable cyclic-complement
argument documented in `docs/nyman-large-sieve-tail-v1.md`. The Python checker
does not formally prove that published theorem. Generation and replay still
share repository code and the same CPython, NumPy, Python-FLINT, and FLINT
families, so this remains reproducible certification rather than independent
formal verification. It proves complete gain greater than `1/5000` for one
explicit correction and a general finite-support tail lemma; it does not prove
the uniform all-scale energy/gain estimate needed for RH.

### Six-scale fixed-width large-sieve suite

The scaling suite applies the same exact-vector, prefix, and Fourier/Farey
tail path to the six frozen `fixed_width_grid` cells with `K=64` and
`N=8,16,32,64,128,256`. Every cutoff is fixed in advance by `T=32768N`.
Generation uses `2^18` blocks and 256-bit Arb; replay changes to `2^17`
blocks and 384-bit Arb. Each source candidate and scout cell is hash-bound,
and each complete sparse dyadic vector is committed by a canonical payload
hash before its balance, prefix, tail, and strict rational gain gate replay.

The artifact stores exact proof records for all six signs but labels its trend
summary `EXPLORATORY`. The rows are independently optimized old candidates,
not states of a nested recurrence, and the checker proves nothing for an
undeclared `N`. The compact vector commitments make the artifact dependent on
the frozen source files and deterministic reconstruction code; it is not a
self-contained formal proof object. The same shared-backend and published-
theorem limitations as the single-cell large-sieve certificate apply.

### Genuine nested-chain certificate and route obstruction

The nested-chain artifact reconstructs the frozen `N=8,K=64` exact update,
forms `p1=p0+a1`, and builds the second ideal shell from that sparse updated
vector with explicit support cutoff 1,024. It then applies the prescribed
exact `K=8` multiplier and outer scale `1/4`, certifies the prefix through
`2^25`, and applies the same Fourier/Farey tail bound. Generation uses
`2^18` blocks and 256-bit Arb; replay changes to the incommensurate block size
196,613 and 448-bit Arb. Exact vector commitments force step 1's updated vector
to be step 2's old vector.

The artifact proves two finite complete gains and their telescoped sum. It
does not prove an indefinitely iterable recurrence. The companion
frozen-prefix theorem is a human-auditable analytic argument: append-only
harmonic-balanced corrections vanish on the initial interval, Möbius inversion
identifies the only hypothetical zero seed prefix, and Bertrand's postulate
rules that seed out. The checker verifies the finite balance, support, update,
and gain identities; it does not formalize Möbius inversion or Bertrand's
postulate. The theorem rejects this recurrence architecture, not RH.

### Rebased Schur certificate with arbitrary harmonic slope

The rebased certificate binds the frozen exact `N=8` source cell and shared
natural-dilate kernel before constructing anything. The first ideal shell is
derived from the source's full `2^-256` coefficient vector, rounded by the
declared unique-nearest `2^-9` rule, and combined with a frozen sixteen-entry
`2^-16` weight list to reconstruct `p1`. A second ideal shell is derived from
that exact sparse `p1`; a new twenty-four-entry weight list then rebuilds
`p2` from the eight direct coordinates and eight dilates of each shell. Exact
vector commitments bind the source, both shells, both weight lists, `p1`, and
`p2`, including the second shell's `p1` parent.

The two weight lists came from an exploratory Schur-complement optimizer. The
optimizer, its objective value, and its numerical conditioning diagnostics are
not trusted evidence. Once the lists are frozen, the verifier reconstructs
every coefficient and decides the sign independently.

`E(p1)` is enclosed by the canonical Arb natural-dilate Gram builder. Because
`p2` has a nonzero harmonic sum, its complete energy uses a local
arbitrary-slope proof path. The prefix includes `(0,1)` and integer intervals
`M=1,...,T`, with

```text
(T+1)P^2 + sum q_M^2/[M(M+1)] - 2P L_T.
```

The `q_M` recurrence is exact signed int64 after a global overflow proof.
NumPy 2.3.5 performs only guarded binary64 conversions, products, divisions,
and fixed reductions under IEEE-754 round-to-nearest-ties-to-even. Their
rounding contribution is covered by an exact rational error radius, including
block and final-reduction gamma bounds. Arb evaluates the Abel-compressed
logarithmic term at 256 bits.

The omitted tail is bounded by exact rational arithmetic:

```text
rho/(T+1) + Q(Q-1)rho/[(T+1)(T+2)]
  + |P|V/(T+1) + P^2/[4(T+1)].
```

Exact Jordan-`J_2` divisor sums compute `rho`, including its constant
coefficient `c0^2`. The `Q(Q-1)` discrepancy term uses the dual
consecutive-interval form of Montgomery and Vaughan's large-sieve theorem,
Farey spacing, a cyclic-complement argument, and Abel summation. The checker
does not formally prove that published theorem.

Generation uses cutoff `T=2^17`, block size `2^16`, and 256-bit Arb. Replay
reconstructs the vectors and complete proof with the incommensurate block size
32,749 and 448-bit Arb. Strict JSON loading, exact schema equality, source and
kernel pins, canonical payload verification, support checks, and the strict
rational threshold reject relevant mutations.

Generation and replay nevertheless share the repository formulas, CPython,
NumPy, Python-FLINT, FLINT, and the same Gram implementation. This is
reproducible finite certification, not clean-room or formal verification. It
proves only `E(p1)-E(p2)>1/5000` for the two explicit vectors. It proves no
third step, all-scale recurrence, convergence theorem, or implication for RH;
the global status remains `UNRESOLVED`.

## Claim ledger

The ledger gate checks schemas, declared SHA-256 binding modes, a pinned checker
manifest, resolution-claim linkage, and distinct declared backend families.
JSON artifacts use a canonical parsed representation so whitespace and checkout
line endings cannot change their identity; arbitrary proof artifacts use raw
bytes. The gate does not
execute arbitrary proof assistants, establish author independence, peer review
the mathematics, or confer Clay recognition. Those are separate promotion
requirements.
