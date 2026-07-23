# Rebased Schur finite contraction with arbitrary harmonic slope v1

Date: 2026-07-23

Status:

- one explicit rebased finite contraction: `CERTIFIED_FINITE`;
- optimizer output: proposal only, not evidence;
- all-scale recurrence or convergence theorem: not proved;
- Riemann Hypothesis: `UNRESOLVED`.

## 1. Exact public claim

For a finite coefficient vector \(p=(p_n)\), put

\[
\rho_n(t)=\{t/n\},\qquad
\chi(t)=\mathbf 1_{[1,\infty)}(t),
\]

and define its complete Nyman energy by

\[
E(p)=\int_0^\infty
\left|\chi(t)-\sum_n p_n\rho_n(t)\right|^2\frac{dt}{t^2}.
\]

The machine certificate
[nyman-rebased-schur-v1.json](../results/nyman-rebased-schur-v1.json)
constructs two explicit exact-dyadic vectors \(p_1,p_2\) and proves exactly

\[
\boxed{E(p_1)-E(p_2)>\frac1{5000}.}
\]

That is the only public inequality promoted by this certificate. It is one
finite statement about two declared vectors. It neither proves nor disproves
the Riemann Hypothesis.

## 2. How the two vectors are fixed

The construction starts from the frozen exact \(N=8\) natural-dilate source.
Its coefficients retain their full \(2^{-256}\) encoding when the first ideal
shell is derived; they are not first rounded to the later weight grid.

1. Derive the ideal shell on indices \(9,\ldots,16\) from that exact source,
   round it to the \(2^{-9}\) grid with unique-nearest-bin checks, and impose
   the final exact balance condition.
2. Build \(p_1\) from eight direct coefficients and eight dilates of that
   shell, using the frozen list of sixteen \(2^{-16}\)-grid weights.
3. Derive a new ideal shell on indices \(129,\ldots,256\) from the exact
   sparse vector \(p_1\), again with the certified \(2^{-9}\) rounding rule.
4. Rebuild \(p_2\) from eight direct coefficients, eight dilates of the first
   shell, and eight dilates of the second shell, using a new frozen list of
   twenty-four \(2^{-16}\)-grid weights.

Thus the second step is rebased. It reoptimizes the old coordinates instead
of appending a support-disjoint harmonic-balanced correction to \(p_1\).
Consequently the frozen-prefix obstruction for append-only balanced chains
does not apply to this one step.

The vector commitments are:

| object | nonzero entries | largest index | denominator grid | payload SHA-256 |
|---|---:|---:|---:|---|
| source \(p_0\) | 8 | 8 | \(2^{-256}\) | `4200517a85a84817b64be599aa5680abb0a7e45ad35b6ada53692410b56ee2d6` |
| first shell | 7 | 16 | \(2^{-9}\) | `376513fa3ce5dfe30d6b1b738ea0ddf1163e643da2600426bc5c7ca68683eeb8` |
| first weight list | 16 | 16 | \(2^{-16}\) | `c93e6ad2f5d90ee3e19dc86249aa2195d335b763c70c02e5bbfa0c97a74c21f1` |
| \(p_1\) | 55 | 128 | \(2^{-25}\) | `4d1f3bb4359542d5007f560b4d57b40426d858b5831377b0fd4f6d8757ae377b` |
| second shell | 128 | 256 | \(2^{-9}\) | `100d0408424f33cf659c64cdd60895a28dc9c82aedefc40e703d2358a0cadbb6` |
| second weight list | 24 | 24 | \(2^{-16}\) | `65aa8df808e491f02db2e9c20573c430ff3b32200089f849eab6c7a5ebd6110d` |
| \(p_2\) | 915 | 2,048 | \(2^{-25}\) | `b392a0c8d9ea2bf09b5fdf98575906bb8c17b75044c4a01586bdc12fa299520a` |

The second-shell parent commitment is exactly the displayed \(p_1\)
commitment. The verifier reconstructs every vector rather than trusting these
hashes as substitutes for the construction.

## 3. Exact prefix identity for arbitrary harmonic slope

This certificate cannot use the earlier balanced-tail specialization because
the new vector may have nonzero harmonic slope. Let \(p\) be any finite real
coefficient vector supported through \(Q\), and write

\[
P=\sum_n\frac{p_n}{n},\qquad
r_p(t)=\chi(t)-\sum_n p_n\{t/n\}.
\]

For an integer \(M\ge1\), define

\[
q_M=1+\sum_n p_n\left\lfloor\frac Mn\right\rfloor.
\]

Away from interval endpoints, \(t=M+u\) with \(0<u<1\) gives

\[
r_p(t)=q_M-Pt.
\]

On \(0<t<1\), instead, \(r_p(t)=-Pt\). Therefore the energy through the
integer interval \(M=T\), meaning all of \(0<t<T+1\), is exactly

\[
\boxed{
E_{\le T}(p)
=(T+1)P^2
+\sum_{M=1}^{T}\frac{q_M^2}{M(M+1)}
-2P L_T,
}
\]

where

\[
L_T=\sum_{M=1}^{T}q_M\log\frac{M+1}{M}.
\]

Indeed, direct integration on \((M,M+1)\) gives

\[
\int_M^{M+1}\frac{(q_M-Pt)^2}{t^2}\,dt
=\frac{q_M^2}{M(M+1)}
-2Pq_M\log\frac{M+1}{M}+P^2,
\]

and the omitted interval \((0,1)\) contributes the additional \(P^2\).

The logarithmic sum does not require \(T\) separate transcendental
evaluations. Put \(k_n=\lfloor T/n\rfloor\). Abel summation and

\[
q_M-q_{M-1}=\sum_{n\mid M}p_n
\]

give the exact compression

\[
\boxed{
L_T=q_T\log(T+1)
-\sum_n p_n\left(k_n\log n+\log\Gamma(k_n+1)\right).
}
\]

To see this, telescope the logarithms first:

\[
L_T=q_T\log(T+1)
-\sum_{M=1}^{T}(q_M-q_{M-1})\log M,
\]

where the \(M=1\) term vanishes. Reordering the divisor sum turns the last
term into

\[
\sum_n p_n\sum_{j\le T/n}\log(nj)
=\sum_n p_n\left(k_n\log n+\log(k_n!)\right).
\]

The implementation computes the rational square sum by an exact signed-int64
divisor recurrence followed by binary64 divisions and reductions. Binary64
is not treated as exact: the artifact includes a rational error radius for
conversion, multiplication, division, block reduction, and final reduction.
Arb encloses only the compressed logarithmic expression and its combination
with the rational interval.

## 4. Absolute tail theorem for arbitrary harmonic slope

The prefix identity leaves the intervals \(M=T+1,T+2,\ldots\). Center each
finite sawtooth by putting

\[
\phi_n(M)=\left\{\frac Mn\right\}-\frac{n-1}{2n},
\qquad
c_0=1-\frac12\sum_n p_n,
\]

and

\[
v_M=c_0-\sum_n p_n\phi_n(M).
\]

For \(t=M+u\), \(0<u<1\),

\[
\boxed{r_p(M+u)=v_M-P(u-\tfrac12).}
\]

The sequence \(v_M\) is periodic and its exact mean square is

\[
\boxed{
\rho=c_0^2+\frac1{12}\sum_{d=2}^{Q}
J_2(d)\left(\sum_{\substack{n\le Q\\d\mid n}}\frac{p_n}{n}\right)^2,
}
\]

where \(J_2\) is the second Jordan totient. For completeness, over a common
period,

\[
\operatorname{mean}\bigl(\phi_m\phi_n\bigr)
=\frac{\gcd(m,n)^2-1}{12mn}.
\]

Expanding \(v_M^2\) and using
\(\sum_{d\mid g}J_2(d)=g^2\) yields the displayed divisor-sum formula. In
particular, the constant term \(c_0^2\) is part of \(\rho\); it must not be
dropped.

All nonconstant Fourier frequencies of \(v_M\) are reduced fractions whose
denominators are at most \(Q\). Distinct such frequencies have circular
spacing at least

\[
\delta\ge\frac1{Q(Q-1)}
\]

when \(Q\ge2\). The dual consecutive-interval form of Montgomery and
Vaughan's large sieve, together with a cyclic-complement argument over a
common period, therefore gives every integer interval \(I\) the two-sided
bound

\[
\left|
\sum_{M\in I}(v_M^2-\rho)
\right|
\le Q(Q-1)\rho.
\]

For \(Q=1\), \(v_M\) is constant and the same statement holds with the
right-hand side zero. Abel summation against
\(w_M=1/[M(M+1)]\) now proves

\[
\sum_{M>T}\frac{v_M^2}{M(M+1)}
\le
\frac{\rho}{T+1}
+\frac{Q(Q-1)\rho}{(T+1)(T+2)}.
\]

Finally,

\[
V=|c_0|+\frac12\sum_n|p_n|
\]

satisfies \(|v_M|\le V\). Since \(2|u-\tfrac12|\le1\) and
\((u-\tfrac12)^2\le1/4\), summing the absolute cross term and the slope-square
term gives

\[
\sum_{M>T}
\left|\int_0^1
\frac{-2Pv_M(u-\tfrac12)}{(M+u)^2}\,du\right|
\le\frac{|P|V}{T+1},
\]

\[
\sum_{M>T}\int_0^1
\frac{P^2(u-\tfrac12)^2}{(M+u)^2}\,du
\le\frac{P^2}{4(T+1)}.
\]

Combining the three pieces proves the absolute complete-tail bound used by
the certificate:

\[
\boxed{
E_{>T}(p)
\le
\frac{\rho}{T+1}
+\frac{Q(Q-1)\rho}{(T+1)(T+2)}
+\frac{|P|V}{T+1}
+\frac{P^2}{4(T+1)}.
}
\]

Here \(E_{>T}\) starts at \(M=T+1\). This indexing matches the prefix, which
already includes \((0,1)\) and every interval \(M=1,\ldots,T\).

The large-sieve input is H. L. Montgomery and R. C. Vaughan,
[*The large sieve*](https://doi.org/10.1112/S0025579300004708),
Mathematika 20 (1973), 119-134, Theorem 1 in dual form. The repository
verifier checks the finite identities and exact constants; it is not a formal
proof assistant for the published theorem.

## 5. Certificate instantiation

For \(p_2\), the artifact uses

\[
T=131072=2^{17},\qquad Q=2048,\qquad Q(Q-1)=4192256.
\]

Its exact tail ingredients include

\[
c_0=\frac{223}{8192},\qquad
\rho\approx1.36681669635944,\qquad
V=\frac{2234890213}{33554432}.
\]

The four rational upper-bound contributions are:

| tail contribution | value |
|---|---:|
| periodic mean | \(1.04279042698301\times10^{-5}\) |
| large-sieve discrepancy | \(3.33524911444076\times10^{-4}\) |
| harmonic-slope cross term | \(2.95015615847671\times10^{-6}\) |
| slope square | \(6.42881093523137\times10^{-11}\) |
| total tail upper bound | \(0.000346903036160492\) |

The generation proof then records

```text
E(p1) direct Arb enclosure:
    [0.01679752084590578205329038607663539753563 +/- 4.71e-42]

p2 prefix enclosure:
    [0.0161329708 +/- 6.94e-11]

certified p2 complete-energy upper computation:
    [0.0164798739 +/- 5.92e-11]

certified gain-lower computation:
    [0.0003176470 +/- 5.96e-11]
```

The lower endpoint of the last enclosure is strictly greater than
\(1/5000=0.0002\). The old energy is evaluated directly from the canonical
Arb natural-dilate Gram matrix; no optimizer objective value is reused.

Generation uses two binary64 blocks of size \(2^{16}\) and 256-bit Arb.
Verification reconstructs the vectors and repeats the proof with the
incommensurate block size 32,749 and 448-bit Arb.

## 6. Evidence boundary and limitations

The two weight lists were selected by an exploratory, correlation-normalized
Schur-complement search. They are frozen proposals. The certificate's sign
decision comes only after exact reconstruction, direct Gram evaluation of
\(p_1\), and an independently replayed prefix-plus-tail upper bound for
\(p_2\).

The verifier binds the frozen source cell and shared kernel, checks the shell
rounding rules, reconstructs all seven committed vectors, checks the exact
support contracts, regenerates the stored proof, and repeats the strict gain
test at the replay settings. Strict loading and the payload commitment reject
schema, number-encoding, and content changes.

Generation and replay still share the repository implementation, CPython,
NumPy, Python-FLINT, FLINT, and the same mathematical formulas. This is
reproducible finite certification, not independent formal verification.

Most importantly:

- the optimizer is not evidence;
- no third step is certified here;
- no arbitrary-scale continuation is proved;
- no uniform contraction, support-growth, or convergence theorem is proved;
- the finite gain cannot establish the Nyman limit;
- RH remains `UNRESOLVED`.

## 7. Reproduction and commitments

Replay the frozen artifact:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_rebased_schur_certificate.py `
  --verify --artifact results\nyman-rebased-schur-v1.json
```

Relevant files:

- [certificate generator and verifier](../tools/generate_nyman_rebased_schur_certificate.py);
- [focused tests](../tests/test_nyman_rebased_schur_certificate.py);
- [frozen artifact](../results/nyman-rebased-schur-v1.json).

Hashes:

```text
artifact payload SHA-256:  a6e1119fb457e18b521e70203483101a77be7509d52d4f41889951399ce95edd
raw-file SHA-256:          a2e68b81dc37a252b35070f73579cea73705a063060c4dab517d42e2738ef0aa
canonical artifact SHA:   f26192ff4f2646cd74dc3f2e55582580618641722e2576ad0edd4a634891671c
```

The Riemann Hypothesis remains unresolved.
