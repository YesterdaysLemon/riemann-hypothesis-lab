# Genuine nested Nyman chain and frozen-prefix obstruction v1

Date: 2026-07-23

Status:

- two linked finite contractions: `CERTIFIED_FINITE`;
- strictly append-only harmonic-balanced shell recurrence: `FAILED` as an
  all-scale route;
- Riemann Hypothesis: `UNRESOLVED`.

## 1. Result

The machine certificate
[nyman-nested-chain-v1.json](../results/nyman-nested-chain-v1.json)
constructs three exact coefficient vectors

\[
p_1=p_0+a_1,
\qquad
p_2=p_1+a_2.
\]

Unlike the preceding six-scale table, the second old vector is not loaded from
an independently optimized candidate. The verifier first reconstructs
`p1=p0+a1`, builds the next ideal shell from that exact sparse `p1`, and only
then constructs `a2`.

| update | old support | correction | new support | prefix cutoff | complete-gain claim |
|---|---:|---:|---:|---:|---:|
| `p0 -> p1` | `8` (8 nonzero) | `K=64`, 302 nonzero | `1,024` (310 nonzero) | `2^18` | `> 1/150` |
| `p1 -> p2` | `1,024` (310 nonzero) | `K=8`, 6,293 nonzero | `16,384` (6,603 nonzero) | `2^25` | `> 1/9000` |

The second multiplier is the exact vector

\[
\left(
1,-\frac{317}{256},-\frac{591}{512},\frac5{128},
-\frac{255}{512},\frac{141}{128},-\frac{133}{256},-\frac18
\right),
\]

followed by the outer scale `alpha=1/4`. Its harmonic sum is zero. The
second shell is rounded on the `2^-9` grid, is supported on
`1025,...,2048`, and has exact coefficient sum zero.

At 448-bit replay precision with a block partition different from generation,
the enclosures are

```text
step 1: 0.00677209912099658... <= gain <= 0.00683712021428834...
step 2: 0.000112115654830510... <= gain <= 0.000117329252291817...
total : 0.00688421477572761... <= gain <= 0.00695444946667964...
```

Thus the genuine two-step telescoped gain is strictly greater than
`61/9000`. This remains a finite statement.

## 2. Why the apparent recurrence still cannot reach RH

Work in

\[
H=L^2((0,\infty),dt/t^2),
\qquad e_n(t)=\{t/n\},
\qquad \chi(t)=\mathbf 1_{[1,\infty)}(t).
\]

For a finite vector `p` define

\[
r_p(t)=\chi(t)-\sum_n p_ne_n(t).
\]

The following obstruction is independent of the finite large-sieve tail
estimate.

### Frozen-prefix theorem

Let `Q>=1`, let `p0` be any finite vector supported through `Q`, and suppose
every later correction `a_j` is finite, has minimum support at least `Q+1`,
and obeys

\[
\sum_n\frac{a_{j,n}}n=0.
\]

Then every updated residual equals the seed residual on `0<t<Q+1` and there
is a fixed positive constant

\[
B(p_0)=\int_0^{Q+1}|r_{p_0}(t)|^2\frac{dt}{t^2}>0
\]

such that every residual energy is at least `B(p0)`. Consequently no such
append-only harmonic-balanced recurrence can drive the Nyman distance to
zero.

This covers every fixed-width shell-multiplier chain considered here: a shell
supported on `(Q_j,2Q_j]` convolved with a finite multiplier is supported
strictly above `Q_j`, and harmonic balance of the multiplier gives harmonic
balance of the correction.

### Proof that the prefix is frozen

If `0<t<Q+1` and `n>Q`, then `0<t/n<1`, so

\[
\sum_n a_{j,n}\{t/n\}
=t\sum_n\frac{a_{j,n}}n=0.
\]

Every correction therefore vanishes on the initial interval, and all later
residuals equal `r_p0` there. Endpoints are irrelevant to the `L^2` norm.

### Proof that no finite seed has a zero prefix

Assume for contradiction that `r_p0=0` almost everywhere on `(0,Q+1)`.
Put

\[
P_1=\sum_{n\le Q}\frac{p_{0,n}}n.
\]

On `(0,1)`, the residual is `-t P1`; hence `P1=0`. On each open interval
`(m,m+1)`, `1<=m<=Q`, vanishing now gives

\[
1+\sum_{n\le Q}p_{0,n}\left\lfloor\frac mn\right\rfloor=0.
\]

Successive differences imply

\[
\sum_{n\mid1}p_{0,n}=-1,
\qquad
\sum_{n\mid m}p_{0,n}=0\quad(2\le m\le Q).
\]

Mobius inversion forces \(p_{0,n}=-\mu(n)\) for every `n<=Q`. The already
deduced identity `P1=0` would therefore require

\[
\sum_{n\le Q}\frac{\mu(n)}n=0.
\]

For `Q=1` this is plainly false. For `Q>=2`, Bertrand's postulate supplies a
prime `q` with `Q/2<q<=Q`. Let `L=lcm(1,...,Q)`. Modulo `q`, every term in

\[
L\sum_{n\le Q}\frac{\mu(n)}n
=\sum_{n\le Q}\mu(n)\frac Ln
\]

vanishes except `n=q`: no other `n<=Q` is divisible by `q`. The remaining
term is `-L/q`, which is nonzero modulo `q` because `q^2>Q`. This is a
contradiction, so `B(p0)>0` for every finite seed.

## 3. Exact floor for the certified chain

For this certificate's rounded `N=8` seed,

\[
P_1=\frac{241057}{27525120}.
\]

Every stored correction has harmonic sum zero, so this value is invariant.
Already on `(0,1)`,

\[
\|r_{p_j}\|_H^2
\ge P_1^2
=\frac{58108477249}{757632231014400}
\approx 7.6697472560\times10^{-5}.
\]

The two certified gains are compatible with this floor: they prove genuine
finite contraction, but further contractions in the same architecture can
never make the energy tend to zero.

## 4. What survives

The obstruction rejects append-only harmonic-balanced updates, not the
Nyman--Beurling criterion and not RH. A viable recurrence must alter the
earlier residual. Concrete options are:

1. reoptimize old coefficients when each new shell is introduced;
2. use overlapping multiscale corrections rather than support-disjoint
   append-only ones;
3. allow a controlled harmonic-moment change and generalize the tail proof;
4. study a residualized Schur-complement update in which new directions are
   projected through the current old span.

For a repaired recurrence with energies `E_j`, relative gains
`eta_j=(E_j-E_(j+1))/E_j`, and support growth `Q_(j+1)<=16Q_j`, a sufficient
target remains

\[
\eta_j\ge \frac{c}{\log Q_j}
\]

for some fixed `c>0`: the sum of these relative gains diverges. The present
append-only construction cannot be used to establish that target because of
the frozen-prefix theorem.

## 5. Certificate and trust boundary

Generation uses 256-bit Arb and blocks of `2^18`. Replay uses 448-bit Arb and
the deliberately incommensurate block size `196613`. Each step independently
recomputes its integer prefix and Fourier/Farey large-sieve tail. The verifier
also reconstructs all exact vectors and requires the commitment of step 1's
updated vector to equal the commitment of step 2's old vector.

The proof of the frozen-prefix theorem is human-auditable mathematics. The
JSON verifier checks the exact finite chain and its balance/update identities;
it is not a formal proof assistant for Bertrand's postulate or Mobius
inversion.

The prime-existence input is Bertrand's postulate; an elementary primary
source is Paul Erdos,
[*Beweis eines Satzes von Tschebyschef*](https://users.renyi.hu/~p_erdos/1932-01.pdf)
(1932).

Hashes:

```text
artifact payload SHA-256: 1171b416a8aa949ce380d5687cf8865d5550891ac24ab0d781cfa8261627bcad
raw-file SHA-256:         171325a504be89fd35d380531f9ea3a20a713672beb2fc975a4f5ba1e35855d2
canonical artifact SHA:  c21082e8382f3128fac262c43d81d722a70ae3770be097fe4fafa709b6ee51e6
```

Replay:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_nested_chain_certificate.py --verify
```

The Riemann Hypothesis remains unresolved.
