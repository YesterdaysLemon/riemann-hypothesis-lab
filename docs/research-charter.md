# Research charter and confirmation standard

## Objective

Prove or disprove the Riemann Hypothesis (RH) without weakening its quantifiers,
changing the zeta function, or treating numerical evidence as a theorem. The
project remains `UNRESOLVED` until one of the confirmation gates below closes.

RH says that every nontrivial zero of the meromorphically continued Riemann zeta
function has real part `1/2`. The word *every* is the central difficulty.

## What counts as a confirmed proof

All of the following are required:

1. A complete global argument whose convergence regions, boundary cases,
   multiplicities, branches, limit exchanges, and analytic continuations are
   explicit.
2. A dependency and assumption manifest showing that no lemma equivalent to RH
   has been assumed in disguise.
3. A machine-checked formalization in a pinned proof assistant, with no
   placeholders or new unproved axioms, or an equally inspectable certificate
   checker for every nontrivial analytic estimate.
4. Two independent clean-room reproductions. At least one must not reuse the
   same numerical backend for its trusted core.
5. Adversarial expert review of the exact theorem statement and the bridge from
   any equivalent criterion back to RH.

## What counts as a confirmed disproof

Either of these paths is sufficient after two independent reproductions:

- an interval-certified nontrivial zero in a rectangle whose real projection
  excludes `1/2`, with a rigorous argument-principle or Turing-style count and
  proof that the rectangle contains the asserted number of zeros; or
- a certified violation of a proved criterion genuinely equivalent to RH,
  together with a checked instantiation of every hypothesis of that
  equivalence.

A tiny value of zeta, a floating-point root, a missing sign change, or a zero
count deficit is an anomaly to investigate, not a disproof.

## Classification ladder

`EXPLORATORY` -> `CERTIFIED_FINITE` or `CONDITIONAL` -> independent audit ->
`PROVED`/`DISPROVED`.

`FAILED` is a permanent, useful state. Failed searches and invalidated lemmas
stay in the record so later work cannot silently rediscover them.

The claim-ledger schema/hash gate requires a pinned checker manifest, hash-bound
proof and assumption artifacts, a cross-linked resolution claim, and at least
two artifacts from distinct backend families. This mechanical gate is necessary
but not sufficient: it cannot decide whether the implementations are genuinely
independent or the mathematics is valid.

These are this repository's internal promotion rules. Official Clay recognition
has separate [CMI rules](https://www.claymath.org/millennium-problems/rules/),
including qualifying publication, a waiting period, and general acceptance in
the mathematics community.
