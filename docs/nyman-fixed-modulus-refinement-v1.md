# Fixed-modulus refinement of the third rebased Nyman step v1

Date: 2026-07-23

Status:

- strengthened finite \(p_2\)-to-\(p_3\) inequality: `CERTIFIED_FINITE`;
- inherited third weight list: optimizer proposal only, not evidence;
- fourth or arbitrary-scale conclusion from this certificate: none;
- Riemann Hypothesis: `UNRESOLVED`.

## 1. Exact public claim

For a finite coefficient vector \(p=(p_n)\), define

\[
E(p)=\int_0^\infty
\left|\mathbf 1_{[1,\infty)}(t)-\sum_n p_n\{t/n\}\right|^2
\frac{dt}{t^2}.
\]

The frozen
[fixed-modulus certificate](../results/nyman-fixed-modulus-refinement-v1.json)
reconstructs the exact \(p_2,p_3\) already published by the
[third-step certificate](../results/nyman-rebased-third-step-v1.json) and
proves the stronger finite inequality

\[
\boxed{E(p_2)-E(p_3)>\frac1{17500}.}
\]

No coefficients were changed. The improvement comes only from a sharper
upper bound on the omitted tail of the same \(p_3\) prefix. This is a finite
statement about two explicit vectors, not an all-scale lemma and not a proof
or disproof of RH.

## 2. Fixed-modulus local spacing

For \(d\ge2\), put

\[
P_d=\sum_{d\mid n}\frac{p_n}{n},
\qquad
\mathcal D=\{d\ge2:P_d\ne0\},
\qquad
Q_A=\max\mathcal D.
\]

The original third-step proof used, for every active denominator \(d\), the
uniform reciprocal-gap charge \(dQ_A\). Fix a positive integer \(L\), and
instead define

\[
m_d=\gcd(d,L),\qquad
R_m=\max_{e\in\mathcal D}\frac{e}{\gcd(m,e)}.
\]

If \(a/d\) and \(b/e\) are distinct reduced frequencies, their circular
distance is at least

\[
\frac{\gcd(d,e)}{de}.
\]

Because \(m_d\mid d\), one has
\(\gcd(d,e)\ge\gcd(m_d,e)\). The reciprocal local gap at denominator \(d\)
is therefore at most \(dR_{m_d}\). Applying Montgomery and Vaughan's
weighted periodic cosecant theorem gives the exact discrepancy constant

\[
\boxed{
C_L=\frac32\left[
Q_Ac_0^2+
\frac1{12}\sum_{d\in\mathcal D}
dR_{m_d}J_2(d)P_d^2
\right].
}
\]

For \(L=1\), every \(m_d=1\) and \(R_1=Q_A\), so this reduces exactly to the
active-\(Q\) constant \((3/2)Q_A\sigma\) used previously. The refinement is
therefore a strict generalization, not a different tail model.

The cited analytic input remains H. L. Montgomery and R. C. Vaughan,
[*Hilbert's inequality*](https://doi.org/10.1112/jlms/s2-8.1.73),
Journal of the London Mathematical Society (2) 8 (1974), 73-82, Theorem 1
in its weighted periodic cosecant form.

## 3. Frozen \(L=840\) computation

For the published \(p_3\), the verifier finds 17,157 active denominators
through \(Q_A=32768\). Choosing \(L=840\) partitions them into 32 active
\(\gcd(d,840)\) classes. All divisor sums, gcd classes, Jordan-\(J_2\)
weights, and rational constants are recomputed exactly.

| quantity | exact-record decimal, approximately |
|---|---:|
| active-\(Q\) spacing constant | \(2146613456.2376066310\) |
| fixed-modulus spacing constant \(C_{840}\) | \(1816711317.6485586848\) |
| active-\(Q\)/fixed-modulus ratio | \(1.18159304419\) |
| periodic mean term at \(T=2^{23}\) | \(7.4824993850\times10^{-7}\) |
| fixed-modulus discrepancy term | \(2.5817011362\times10^{-5}\) |
| harmonic-slope cross term | \(3.6426171847\times10^{-7}\) |
| slope-square term | \(9.9437940684\times10^{-13}\) |
| complete \(p_3\) tail upper bound | \(2.6929524013\times10^{-5}\) |

The new spacing constant is about 15.37% smaller than the active-\(Q\)
constant. The comparison is between two valid worst-case constants for this
one vector; it is not an empirical convergence rate.

At the unchanged cutoff \(T=8388608\), nonnegativity of the omitted
\(p_2\) integral gives the old complete-energy lower bound, while the
\(p_3\) prefix plus the refined tail gives the new complete-energy upper
bound. The generation enclosure has gain lower endpoint

\[
0.00005730562650289445\ldots
>
0.00005714285714285714\ldots
=\frac1{17500}.
\]

Generation uses block size 65,536 and 320-bit Arb. Verification changes to
the incommensurate block size 250,003 and 512-bit Arb and requires the
generation and replay enclosures to overlap.

## 4. Evidence and limitations

Before reconstruction, the verifier pins the parent third-step artifact by
raw LF, canonical, payload, schema, semantic status, and exact \(p_3\)
commitments. It then:

- reconstructs the inherited source lineage and the exact \(p_2,p_3\);
- recomputes the active denominator set and all 32 modulus classes;
- relies on focused exact tests showing that \(L=1\) recovers the prior
  active-\(Q\) formula;
- regenerates the stored 320-bit proof and replays it at 512 bits; and
- applies the strict rational \(1/17500\) gate.

The optimizer that originally proposed the third weight list is inherited
unchanged and is not trusted evidence. Generation and replay share repository
formulas, CPython, NumPy, Python-FLINT, FLINT, and the prefix implementation.
This is reproducible finite certification, not a clean-room implementation or
formal proof of Montgomery and Vaughan's theorem.

This certificate supplies no new vector and no all-scale estimate. A separate
certificate now proves a finite fourth step, but neither result establishes a
fifth step, an arbitrary-scale recurrence, convergence to zero, or RH.

## 5. Reproduction and commitments

Replay the frozen artifact:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_fixed_modulus_refinement_certificate.py `
  --verify --artifact results\nyman-fixed-modulus-refinement-v1.json
```

Relevant files:

- [fixed-modulus tail helper](../tools/certify_nyman_fixed_modulus_tail.py);
- [certificate generator and verifier](../tools/generate_nyman_fixed_modulus_refinement_certificate.py);
- [focused certificate tests](../tests/test_nyman_fixed_modulus_refinement_certificate.py);
- [focused theorem tests](../tests/test_nyman_fixed_modulus_tail.py);
- [frozen artifact](../results/nyman-fixed-modulus-refinement-v1.json);
- [parent third-step proof note](nyman-rebased-third-step-v1.md).

Hashes:

```text
artifact payload SHA-256:  de69623dee4550796471850d62bf78ec3162fdec34620526f21d0ac578f35f92
raw LF file SHA-256:       45571b59b62c613fe5a107db76a1c80a732d223254065d0d23ab08230894db2b
canonical artifact SHA:   5666bfbd314da21f0d286924f2a4b93c3ff97fc025d269dc769f396e627efe54
```

The Riemann Hypothesis remains `UNRESOLVED`.
