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

## Claim ledger

The ledger gate checks schemas, declared SHA-256 binding modes, a pinned checker
manifest, resolution-claim linkage, and distinct declared backend families.
JSON artifacts use a canonical parsed representation so whitespace and checkout
line endings cannot change their identity; arbitrary proof artifacts use raw
bytes. The gate does not
execute arbitrary proof assistants, establish author independence, peer review
the mathematics, or confer Clay recognition. Those are separate promotion
requirements.
