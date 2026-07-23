# Certified fourth rebased Nyman step v1

Date: 2026-07-23

Status:

- one explicit fourth rebased finite contraction: `CERTIFIED_FINITE`;
- fourth weight list: optimizer proposal only, not evidence;
- fifth step, all-scale recurrence, or convergence theorem: not proved;
- Riemann Hypothesis: `UNRESOLVED`.

## 1. Exact public claim

For a finite coefficient vector \(p=(p_n)\), define

\[
E(p)=\int_0^\infty
\left|\mathbf 1_{[1,\infty)}(t)-\sum_n p_n\{t/n\}\right|^2
\frac{dt}{t^2}.
\]

The frozen
[fourth-step certificate](../results/nyman-rebased-fourth-step-v1.json)
reconstructs the exact published \(p_3\), derives a fourth balanced shell,
and constructs an exact-dyadic \(p_4\) satisfying

\[
\boxed{E(p_3)-E(p_4)>\frac1{100000}.}
\]

This is one finite inequality. The weight proposal supplies construction data
only; exact reconstruction, interval arithmetic, and the complete-tail bound
decide the sign.

Combining this inequality with the separately certified fixed-modulus
refinement \(E(p_2)-E(p_3)>1/17500\) gives the exact finite corollary

\[
\boxed{
E(p_2)-E(p_4)>
\frac1{17500}+\frac1{100000}
=\frac{47}{700000}.
}
\]

That corollary is only algebraic composition of two finite, hash-bound
certificates. It is not an all-scale claim and does not imply that the process
continues or that \(E(p_k)\) tends to zero.

## 2. Exact construction and lineage

The verifier first pins the third-step parent artifact:

| parent commitment | SHA-256 |
|---|---|
| raw LF bytes | `e5bd17bf94f2ce95611ad09bf75acaecc156db003ce4481b4a6e41c3504340bd` |
| canonical parsed JSON | `811ebaf03089288c99d8155171a07e4b5e5731d53678265354ed2738157fa1ca` |
| parent payload | `bbafc2c92acb954c9c9e550f781612965bc5c76fe59ebd4a1adc5b596c3bee02` |
| exact parent \(p_3\) | `4797b563591fd11321c2ff2ef7adbe6cbe6ff35f7ccf7d639eb780ef50ca9c97` |

It also requires the parent's inherited natural-source and kernel lineage to
match a fresh exact reconstruction.

From exact \(p_3\), the verifier derives the ideal shell on indices
\(32769,\ldots,65536\), rounds it by the declared unique-nearest
\(2^{-9}\) rule, and imposes exact zero-sum balance. The fourth proposal has
33 exact \(2^{-16}\)-grid weights:

- eight direct-coordinate weights;
- eight dilates of each of shells 1, 2, and 3; and
- one undilated fourth-shell weight.

Retaining only the first multiplier of shell 4 is a declared finite basis
choice. It avoids treating a much larger full-dilate candidate as certified
evidence. The resulting update overlaps the old support and changes 13,537
old \(p_3\) coefficients, so it is outside the append-only frozen-prefix
no-go hypotheses.

Exact commitments:

| object | nonzero entries | largest index | denominator grid | payload SHA-256 |
|---|---:|---:|---:|---|
| \(p_3\) | 14,660 | 32,768 | \(2^{-25}\) | `4797b563591fd11321c2ff2ef7adbe6cbe6ff35f7ccf7d639eb780ef50ca9c97` |
| fourth shell | 32,685 | 65,536 | \(2^{-9}\) | `05d1aa0b6a1f7d2f8b8e303f6a699725f632b4211462331dee6ab3c846b88fad` |
| fourth weights | 33 | 33 | \(2^{-16}\) | `eca0b2afabc8ef71888292d21f0d704728957bd740cf532445bebbaf8f7c7436` |
| \(p_4\) | 47,345 | 65,536 | \(2^{-25}\) | `bcc33e26af682975e7c711b74cfd954f246f0fb8756e3d93de6e02e065f801f1` |

The common-denominator numerators of the fourth weight list are:

```text
denominator: 65536
numerators:
-62704, 60496, 60567, 7832, 49590, -32697, 40718, 3700,
52259, -30480, -31835, 931, -17330, 19824, -15002, 847,
41140, -28229, -33848, -3095, -14364, 30892, -10955, 5776,
35212, -28199, -28837, 949, -11048, 28333, -8579, -710, 2390
```

## 3. Complete prefix and tail proof

The certificate uses the arbitrary-harmonic-slope prefix identity and the
denominator-weighted local-spacing tail theorem documented in the
[third-step proof](nyman-rebased-third-step-v1.md). For \(p_4\), with support
limit \(Q=65536\), the exact tail records are approximately:

| quantity | value |
|---|---:|
| harmonic slope \(P\) | \(0.00576792384723181028\) |
| periodic mean square \(\rho\) | \(17.7497922118135126\) |
| denominator-weighted \(\sigma\) | \(222994.538114203232\) |
| local constant \((3/2)Q\sigma\) | \(21921255074.7786345\) |
| global/local constant ratio | \(3.47761173117602\) |
| periodic mean term | \(2.6449251096\times10^{-7}\) |
| local discrepancy term | \(4.8674962049\times10^{-6}\) |
| harmonic-slope cross term | \(9.9345681240\times10^{-8}\) |
| slope-square term | \(1.2393647809\times10^{-13}\) |
| complete tail upper bound | \(5.2313345211\times10^{-6}\) |

At \(T=2^{26}=67108864\), the proof compares bounds in the required
directions:

\[
E(p_3)\ge E_{\le T}(p_3),
\qquad
E(p_4)\le E_{\le T}(p_4)+E_{>T}^{\rm upper}(p_4).
\]

The generation proof has exact old-prefix lower endpoint

\[
0.0160596180957771374198\ldots
\]

and new complete-energy upper endpoint

\[
0.0160490741178841575220\ldots.
\]

Their difference is greater than

\[
0.00001054397789253718\ldots>\frac1{100000}.
\]

Generation uses block size 1,048,576 and 512-bit Arb. Verification rebuilds
the complete construction and proof with the incommensurate block size
250,003 and 640-bit Arb, and it requires the generation and replay
enclosures to overlap.

## 4. Evidence and limitations

The verifier:

- strictly loads the artifact and checks its frozen payload;
- pins the parent in raw, canonical, payload, schema, source-lineage, and
  exact-\(p_3\) forms;
- reconstructs all four shells, all four weight lists, and \(p_4\) exactly;
- checks unique shell rounding, exact balance, dyadic grids, support, overlap,
  and vector commitments;
- regenerates the stored proof at 512 bits and replays at 640 bits; and
- applies the strict rational \(1/100000\) threshold.

The fourth weights were selected by an exploratory optimizer. The optimizer's
objective and diagnostics are not evidence. Generation and replay share the
repository formulas, CPython, NumPy, Python-FLINT, FLINT, and prefix
implementation. This is reproducible finite certification, not a clean-room
implementation or a formal proof of the Montgomery--Vaughan theorem.

The reduced fourth-shell basis is a successful finite engineering choice, not
a theorem that the same choice works again. The certificate proves no fifth
step, arbitrary-scale continuation, uniform contraction, or convergence
theorem. RH remains unresolved.

## 5. Reproduction and commitments

Replay the frozen artifact:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_rebased_fourth_step_certificate.py `
  --verify --artifact results\nyman-rebased-fourth-step-v1.json
```

Relevant files:

- [certificate generator and verifier](../tools/generate_nyman_rebased_fourth_step_certificate.py);
- [focused tests](../tests/test_nyman_rebased_fourth_step_certificate.py);
- [frozen artifact](../results/nyman-rebased-fourth-step-v1.json);
- [fixed-modulus refinement](nyman-fixed-modulus-refinement-v1.md);
- [parent third-step proof note](nyman-rebased-third-step-v1.md).

Hashes:

```text
artifact payload SHA-256:  a807f63928869c6451fa397edcc294ff6b3553b69a54099e8053722bd2b594b0
raw LF file SHA-256:       16bd8e0f0fea73599a8a2fcd0cf45bed60163b5fab0c52314c43cefbb1beb560
canonical artifact SHA:   48b142aa8e21bcc7ac0ec980765f30ca13da831fad54c99a6d9e655db9258f3d
```

The Riemann Hypothesis remains `UNRESOLVED`.
