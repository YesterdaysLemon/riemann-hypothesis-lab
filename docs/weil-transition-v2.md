# Prime-power transition search v2

## Classification

This experiment is `EXPLORATORY`. It evaluates 81 explicitly listed finite
Weil matrices and cannot prove RH. A negative integer-vector evaluation is
only a quarantined candidate until its full matrix, cutoff, prime-power
transcript, normalization, parity identity, degree transport, and independent
backend reproduction have all been checked. Every v2 artifact keeps the
global hypothesis status `UNRESOLVED`.

## Why examine transitions

The finite Weil matrix uses the complete prime-power transcript satisfying
the exact inequality `p^r <= c`. At an exact transition `c=q=p^r`, the new
endpoint kernel is identically zero. Immediately above the transition the new
term participates, so the matrix is continuous but has new one-sided local
structure. This makes close rational neighbors of `q` a controlled place to
look for finite negative directions without ever truncating an allowed prime
term.

For `j in {4,6,8,10}`, v2 uses the exact log-symmetric stencil

\[
c_-(q,j)=q\frac{2^j}{2^j+1},\qquad
c_0(q)=q,\qquad
c_+(q,j)=q\frac{2^j+1}{2^j}.
\]

Thus `c_-(q,j)c_+(q,j)=q^2` exactly. Every rational is reduced before
hashing. Plan validation independently checks `q=p^r`, primality of `p`, the
stencil formula and side label, the absence of an intervening prime-power
transition, unique matrix inputs, four matched below/above pairs plus one
shared exact point for every `(q,N)`, and the canonical plan hash.

## Frozen 81-cell design

There are nine cutoffs and three Fourier degrees at each of the three
transitions:

| Transition | Prime-power identity | Degrees | Attempt bits | Dedicated candidate replay |
|---|---|---|---|---|
| `q=7` | `7^1` | `12` | `192, 384` | `768` |
| `q=7` | `7^1` | `16` | `192, 384, 768` | `1536` |
| `q=7` | `7^1` | `20` | `384, 768` | `1536` |
| `q=8` | `2^3` | `16` | `192, 384, 768` | `1536` |
| `q=8` | `2^3` | `20, 24` | `384, 768` | `1536` |
| `q=9` | `3^2` | `16` | `192, 384, 768` | `1536` |
| `q=9` | `3^2` | `20, 24` | `384, 768` | `1536` |

The explicit cutoffs are:

```text
q=7:  112/17, 448/65, 1792/257, 7168/1025, 7,
      7175/1024, 1799/256, 455/64, 119/16
q=8:  128/17, 512/65, 2048/257, 8192/1025, 8,
      1025/128, 257/32, 65/8, 17/2
q=9:  144/17, 576/65, 2304/257, 9216/1025, 9,
      9225/1024, 2313/256, 585/64, 153/16
```

The canonical execution order is global exact `Fraction(c)` order, then
degree, numerator, and denominator. It is intentionally not grouped by
transition: the outer `q=8` and `q=9` stencil intervals overlap. The engine is
serial because FLINT's precision context is shared mutable state. Atomic
attempt and cell files permit deterministic interruption and resume.

The v1 numerical state machine remains responsible for constructing the
complete `A=P-R-S` matrix, exact source-oracle cross-checks, interval
`LDL^T`, diagnostic Rump enclosures, deterministic integer-candidate
extraction, and direct interval Rayleigh evaluation. V2 binds each invocation
to a separately hashed transition contract and does not stop the remaining
bounded cells when a candidate is found. A v1 next-attempt confirmation is not
enough for v2: the same primitive integer vector must also remain strictly
negative at the cell's dedicated replay precision.

## Structural audits

Reversal centrosymmetry gives an exact orthonormal decomposition into

```text
even: e_0; (e_-k + e_k)/sqrt(2), k=1,...,N
odd:       (e_-k - e_k)/sqrt(2), k=1,...,N.
```

The parity audit regenerates the full matrix, checks all component-wise
centrosymmetry overlaps, checks the orthonormal transform, encloses the
even/odd cross block at zero, and compares every requested integer witness in
the parity and original bases. Approximate block eigenvectors may suggest an
integer vector; only direct interval evaluation in the original basis decides
its sign.

At a fixed rational cutoff, the degree-`N` matrix is the central principal
submatrix of every degree-`M` matrix with `M>N`. The nesting audit regenerates
`P`, `R`, `S`, and `A` at two precisions, checks every shared entry, and
zero-pads any quarantined lower-degree witness. A genuinely negative witness
must remain negative after exact zero-padding. Diagnostic Rump spectra and a
Rayleigh--Ritz consistency check can expose implementation problems, but they
never create a positive or negative classification.

Both audit verifiers reconstruct all certificate-bearing matrix and witness
data. They deliberately do not trust or require byte-for-byte reproduction of
an approximate eigenvector as a proof step.

## Run and replay

From an installed development environment, start the frozen batch with:

```powershell
.\.venv\Scripts\rh-lab.exe weil-transition-search `
  --plan plans\weil-transition-q7-q9-v2.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2
```

An interrupted invocation resumes only when the stored plan, backend, engine,
contracts, and existing checkpoints still agree:

```powershell
.\.venv\Scripts\rh-lab.exe weil-transition-search `
  --plan plans\weil-transition-q7-q9-v2.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2 --resume
.\.venv\Scripts\rh-lab.exe verify-weil-transition-search `
  --index results\weil-transition-q7-q9-v2\index.json `
  --checkpoint-dir results\weil-transition-q7-q9-v2
```

`--max-cells K` bounds only the number of newly completed cells in one
invocation; it does not change the plan or final index. Replay regenerates all
certificate-bearing numerical evidence with the same pinned FLINT backend.
That is a consistency reproduction, not an independent implementation.

## Local expansion is motivation, not a gate

For `c>q`, set `h=log(c/q)>0` and
`epsilon=h/log(c)`. For a sufficiently small displacement and fixed Fourier
indices, the new `q=p^r` endpoint kernel has the local expansion

\[
q_{mn}=2\epsilon-
\frac{4\pi^2}{3}(m^2+mn+n^2)\epsilon^3+O(\epsilon^5),
\]

while the odd block begins at cubic order,

\[
O_{k\ell}=-\frac{8\pi^2}{3}k\ell\epsilon^3+O(\epsilon^5).
\]

This explains the parity emphasis. The expansion is not uniform in growing
degree and is especially unsafe as a numerical approximation for the coarser
stencils. It is never used to prune cells, construct a certificate, or decide
a sign.

## Interpretation

The possible terminal numerical cell states are finite positive, terminally
inconclusive, or a quarantined negative candidate. None resolves RH:

- a positive `LDL^T` certificate controls only one finite compression;
- an inconclusive interval computation says nothing about the true sign; and
- even a finite negative matrix direction needs the full mathematical and
  independent-computation promotion audit before it can bear on RH.

The normalization follows
[Bombieri's official problem description](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0).
The finite-compression/core perspective and prolate conditioning motivation
follow [Connes--Consani](https://doi.org/10.4171/LEM/1049). The broader source
map is maintained in [literature-map.md](literature-map.md).
