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

## Claim ledger

The ledger gate checks schemas, SHA-256 file binding, a pinned checker manifest,
resolution-claim linkage, and distinct declared backend families. It does not
execute arbitrary proof assistants, establish author independence, peer review
the mathematics, or confer Clay recognition. Those are separate promotion
requirements.
