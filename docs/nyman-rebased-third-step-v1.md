# Certified third rebased Nyman step v1

Date: 2026-07-23

Status:

- one explicit third rebased finite contraction: `CERTIFIED_FINITE`;
- third weight list: optimizer proposal only, not evidence;
- this artifact contains no fourth step; a separate artifact now certifies one;
- all-scale recurrence or convergence theorem: not proved;
- Riemann Hypothesis: `UNRESOLVED`.

## 1. Exact public claim

For a finite coefficient vector \(p=(p_n)\), put

\[
\rho_n(t)=\{t/n\},\qquad
\chi(t)=\mathbf 1_{[1,\infty)}(t),
\]

and define the complete Nyman energy

\[
E(p)=\int_0^\infty
\left|\chi(t)-\sum_n p_n\rho_n(t)\right|^2\frac{dt}{t^2}.
\]

The frozen machine certificate
[nyman-rebased-third-step-v1.json](../results/nyman-rebased-third-step-v1.json)
constructs exact-dyadic \(p_2,p_3\) and proves

\[
\boxed{E(p_2)-E(p_3)>\frac1{20000}.}
\]

This is the only new public inequality. It is a finite statement about two
explicit vectors. It supplies neither a proof nor a disproof of RH.

## 2. Exact lineage and construction

Before reconstructing a coefficient, the third-step verifier pins the
published parent artifact
[nyman-rebased-schur-v1.json](../results/nyman-rebased-schur-v1.json)
four ways:

| parent commitment | SHA-256 |
|---|---|
| raw LF bytes | `a2e68b81dc37a252b35070f73579cea73705a063060c4dab517d42e2738ef0aa` |
| canonical parsed JSON | `f26192ff4f2646cd74dc3f2e55582580618641722e2576ad0edd4a634891671c` |
| parent payload | `a6e1119fb457e18b521e70203483101a77be7509d52d4f41889951399ce95edd` |
| exact parent \(p_2\) | `b392a0c8d9ea2bf09b5fdf98575906bb8c17b75044c4a01586bdc12fa299520a` |

It also requires the inherited natural-source and shared-kernel records to
match a fresh exact reconstruction. The source is therefore a declared
lineage, not an unverified copy of a coefficient list.

From exact \(p_2\), the verifier:

1. derives the ideal shell on indices \(2049,\ldots,4096\);
2. proves that every declared value is in the unique nearest \(2^{-9}\) bin
   and imposes the final exact zero-sum balance;
3. freezes one thirty-two-entry \(2^{-16}\)-grid weight proposal; and
4. rebuilds \(p_3\) from eight direct coordinates and eight dilates of each
   of the first, second, and third shells.

Thus every old coordinate is reoptimized. The difference \(p_3-p_2\) is not
an append-only, support-disjoint, harmonic-balanced correction. The
frozen-prefix no-go theorem for that narrower architecture does not apply.

The exact vector commitments are:

| object | nonzero entries | largest index | denominator grid | payload SHA-256 |
|---|---:|---:|---:|---|
| source \(p_0\) | 8 | 8 | \(2^{-256}\) | `4200517a85a84817b64be599aa5680abb0a7e45ad35b6ada53692410b56ee2d6` |
| first shell | 7 | 16 | \(2^{-9}\) | `376513fa3ce5dfe30d6b1b738ea0ddf1163e643da2600426bc5c7ca68683eeb8` |
| first weights | 16 | 16 | \(2^{-16}\) | `c93e6ad2f5d90ee3e19dc86249aa2195d335b763c70c02e5bbfa0c97a74c21f1` |
| \(p_1\) | 55 | 128 | \(2^{-25}\) | `4d1f3bb4359542d5007f560b4d57b40426d858b5831377b0fd4f6d8757ae377b` |
| second shell | 128 | 256 | \(2^{-9}\) | `100d0408424f33cf659c64cdd60895a28dc9c82aedefc40e703d2358a0cadbb6` |
| second weights | 24 | 24 | \(2^{-16}\) | `65aa8df808e491f02db2e9c20573c430ff3b32200089f849eab6c7a5ebd6110d` |
| \(p_2\) | 915 | 2,048 | \(2^{-25}\) | `b392a0c8d9ea2bf09b5fdf98575906bb8c17b75044c4a01586bdc12fa299520a` |
| third shell | 2,048 | 4,096 | \(2^{-9}\) | `290a8f9dcb6d9d5731a24d97b2e6848ca4a0778bc16b9890e9505ae83f205335` |
| third weights | 32 | 32 | \(2^{-16}\) | `8df8d30ce08214f66aa746204dbb1dadaae95871878d03058f5d632f7bce418f` |
| \(p_3\) | 14,660 | 32,768 | \(2^{-25}\) | `4797b563591fd11321c2ff2ef7adbe6cbe6ff35f7ccf7d639eb780ef50ca9c97` |

The common-denominator numerators of the third weight list are:

```text
denominator: 65536
numerators:
-62698, 60485, 60554, 7849, 49569, -32652, 40703, 3693,
52082, -30277, -31529, 928, -17135, 19550, -14768, 831,
39574, -26774, -31574, -2862, -13563, 28601, -9984, 4575,
23635, -16402, -14389, 949, -4871, 11964, -3181, 19
```

These integers are construction data, not evidence from the optimizer. The
verifier's exact reconstruction and subsequent sign proof are the evidence.

## 3. Arbitrary-slope prefix

For any finite \(p\) supported through \(Q\), define

\[
P=\sum_n\frac{p_n}{n},\qquad
q_M=1+\sum_n p_n\left\lfloor\frac Mn\right\rfloor.
\]

On \(t=M+u\), \(M\ge1\) and \(0<u<1\), the residual is

\[
\chi(t)-\sum_n p_n\{t/n\}=q_M-Pt.
\]

Including the interval \((0,1)\), the exact prefix through the integer
interval \(M=T\) is

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

Writing \(k_n=\lfloor T/n\rfloor\), Abel summation compresses the logarithmic
sum to

\[
\boxed{
L_T=q_T\log(T+1)
-\sum_n p_n\left(k_n\log n+\log\Gamma(k_n+1)\right).
}
\]

The integer recurrence for \(q_M\) is exact. Guarded binary64 operations
evaluate the rational square sum, with conversion, product, division, block
reduction, and final reduction all covered by an exact rational error radius.
Arb encloses the compressed logarithmic expression.

## 4. Denominator-weighted interval discrepancy theorem

The earlier global Farey-spacing bound costs \(Q(Q-1)\rho\). At
\(Q=32768\), that constant is too expensive for the desired cutoff. This
certificate instead retains the denominator of each Fourier frequency.

Center the periodic residual:

\[
\phi_n(M)=\left\{\frac Mn\right\}-\frac{n-1}{2n},\qquad
c_0=1-\frac12\sum_n p_n,
\]

\[
v_M=c_0-\sum_n p_n\phi_n(M).
\]

Then

\[
\chi(M+u)-\sum_n p_n\{(M+u)/n\}
=v_M-P(u-\tfrac12).
\]

For \(d\ge2\), put

\[
P_d=\sum_{\substack{n\le Q\\d\mid n}}\frac{p_n}{n}.
\]

The exact periodic mean square and its denominator-weighted companion are

\[
\rho=c_0^2+\frac1{12}\sum_{d=2}^{Q}J_2(d)P_d^2,
\]

\[
\boxed{
\sigma=c_0^2+\frac1{12}\sum_{d=2}^{Q}d\,J_2(d)P_d^2.
}
\]

Here \(J_2\) is the second Jordan totient. In the Fourier expansion of
\(v_M\), the total squared coefficient mass at exact reduced denominator
\(d\) is \(J_2(d)P_d^2/12\); the zero-frequency mass is \(c_0^2\).

Let a nonzero frequency be \(\alpha=a/d\) in lowest terms. Against any other
reduced frequency \(b/e\) with \(e\le Q\),

\[
\|\alpha-b/e\|\ge\frac1{de}\ge\frac1{dQ}.
\]

For the zero frequency the corresponding gap is at least \(1/Q\).
Montgomery and Vaughan's weighted periodic cosecant inequality therefore
charges the Fourier mass at denominator \(d\) by at most \(dQ\), instead of
charging every frequency by the worst global gap.

More explicitly, for distinct frequencies \(x_r\) and any
\(0<\delta_r\le\min_{s\ne r}\|x_r-x_s\|\), the cited weighted form is

\[
\left|
\sum_{r\ne s}u_r\overline{u_s}\csc\pi(x_r-x_s)
\right|
\le\frac32\sum_r\delta_r^{-1}|u_r|^2.
\]

For clarity, write \(e(x)=\exp(2\pi i x)\). If
\(I=\{A,\ldots,A+H-1\}\), its off-diagonal Fourier kernel is

\[
\sum_{M\in I}e(M(\alpha-\beta))
=e((A+(H-1)/2)(\alpha-\beta))
\frac{\sin(\pi H(\alpha-\beta))}
     {\sin(\pi(\alpha-\beta))}.
\]

Writing the numerator sine as a difference of two phase factors splits this
into two cosecant bilinear forms. The factor \(1/2\) from that difference
averages the two applications of the \(3/2\) weighted theorem; it does not
double the constant. Consequently every consecutive integer interval obeys

\[
\boxed{
\left|\sum_{M\in I}(v_M^2-\rho)\right|
\le C_{\mathrm{loc}},\qquad
C_{\mathrm{loc}}=\frac32Q\sigma.
}
\]

This is the only new analytic bridge used by the third-step certificate.
It is a finite-support theorem, not an all-scale recurrence lemma.

The cited input is H. L. Montgomery and R. C. Vaughan,
[*Hilbert's inequality*](https://doi.org/10.1112/jlms/s2-8.1.73),
Journal of the London Mathematical Society (2) 8 (1974), 73-82, Theorem 1
in its weighted periodic cosecant form. A
[public author-hosted scan](https://personal.science.psu.edu/rcv4/personal/Publications/s2-8-1-73.pdf)
is also available.

## 5. Complete tail

Apply Abel summation to \(v_M^2-\rho\) with
\(w_M=1/[M(M+1)]\). The interval discrepancy theorem gives

\[
\sum_{M>T}\frac{v_M^2}{M(M+1)}
\le
\frac{\rho}{T+1}
+\frac{C_{\mathrm{loc}}}{(T+1)(T+2)}.
\]

Also set

\[
V=|c_0|+\frac12\sum_n|p_n|,
\]

so \(|v_M|\le V\). Bounding the cross and slope-square pieces in
\(v_M-P(u-\tfrac12)\) yields the complete absolute tail bound

\[
\boxed{
E_{>T}(p)
\le
\frac{\rho}{T+1}
+\frac{C_{\mathrm{loc}}}{(T+1)(T+2)}
+\frac{|P|V}{T+1}
+\frac{P^2}{4(T+1)}.
}
\]

For the frozen \(p_3\), \(T=2^{23}\) and \(Q=32768\). The exact rational
records have the following readable decimal values:

| exact tail quantity | value, approximately |
|---|---:|
| \(P\) | \(0.00577631718021394433\) |
| \(c_0\) | \(3569/131072\) |
| \(\rho\) | \(6.27677616837708865\) |
| \(\sigma\) | \(43672.9625699382860\) |
| \(C_{\rm loc}\) | \(2146613456.23760663\) |
| old global \(Q(Q-1)\rho\) | \(6739431414.47146091\) |
| global/local constant ratio | \(3.13956450561142833\) |
| periodic mean term | \(7.48249938503164\times10^{-7}\) |
| local discrepancy term | \(3.05052010467794\times10^{-5}\) |
| harmonic-slope cross term | \(3.64261718470008\times10^{-7}\) |
| slope-square term | \(9.94379406837140\times10^{-13}\) |
| complete tail upper bound | \(3.16177136981319\times10^{-5}\) |

The factor \(3.13956\ldots\) compares two valid worst-case constants for this
one vector. It is not a measured convergence rate.

## 6. Old lower bound, new upper bound

The proof deliberately compares bounds with opposite directions:

\[
E(p_2)\ge E_{\le T}(p_2)
\]

because the omitted weighted square integral is nonnegative, while

\[
E(p_3)\le E_{\le T}(p_3)+E_{>T}^{\rm upper}(p_3).
\]

Thus it neither subtracts two upper bounds nor imports an optimizer objective.
It also avoids requiring a dense direct-Gram evaluation of \(p_2\).

Generation uses cutoff \(T=8388608\), block size 65,536, and 320-bit Arb.
The stored proof reports:

```text
p2 prefix:
    [0.01614340 +/- 6.50e-9]

exact stored p2 prefix lower endpoint:
    [0.01614340236829551328279208632920926902443 +/- 1.71e-42]

p3 prefix:
    [0.01605917 +/- 6.87e-9]

p3 complete-energy upper computation:
    [0.01609078 +/- 4.94e-9]

gain-lower computation:
    [5.262e-5 +/- 5.65e-9]
```

The exact dyadic enclosures imply the conservative outward-readable bounds

\[
E(p_3)<0.016090785,
\]

\[
E(p_2)-E(p_3)>0.0000526174
>0.00005=\frac1{20000}.
\]

The strict margin above the claimed threshold is greater than
\(2.6174\times10^{-6}\).

Verification rebuilds the parent and all three shells, then repeats the
complete proof with the incommensurate block size 250,003 and 512-bit Arb.

## 7. Evidence and trust boundary

The third weight list came from an exploratory numerical optimizer. Its
objective value and conditioning diagnostics are excluded from the public
proof. Once the integers are frozen, the verifier:

- strictly loads JSON and requires exact field sets;
- pins the parent raw, canonical, payload, source-lineage, and \(p_2\) hashes;
- reconstructs every shell, weight list, and rebased vector exactly;
- checks unique shell rounding, exact balance, grid, support, and vector hashes;
- regenerates the stored proof at 320 bits;
- replays with a different block partition at 512 bits; and
- tests the exact strict rational threshold.

The verifier additionally requires the generation and replay enclosures for
the old prefix, new prefix, new complete upper bound, and gain computation to
overlap. A higher-precision sign alone cannot silently replace a disagreement
with the stored proof.

Focused tests exercise the exact vector commitments, source lineage, tail
constants, and strict inequality. Separate mutation tests reject changed
semantic fields, third weights, parent binding, boolean checks, duplicate
keys, and payload mismatches.

The analytic denominator-weighted discrepancy theorem remains a
human-auditable published-theorem bridge. The tests check exact identities,
the declared constant, and exhaustive interval sums for a small exact period;
they are not a formal proof assistant for Montgomery and Vaughan's theorem.

Generation and replay share repository formulas, CPython, NumPy,
Python-FLINT, FLINT, and the prefix implementation. This is reproducible
finite certification, not a clean-room reproduction or formal verification.

## 8. Exploratory trend and engineering outlook

The following observations are explicitly **exploratory**. They are not part
of the certified inequality and must not be extrapolated.

| object | support limit | nonzero count | coefficient \(L^1\), approx. |
|---|---:|---:|---:|
| \(p_1\) | 128 | 55 | 15.0439067 |
| \(p_2\) | 2,048 | 915 | 133.155394 |
| \(p_3\) | 32,768 | 14,660 | 1057.93770 |

The rounded shell \(L^1\) values likewise grow from \(4.73046875\), to
\(55.52734375\), to \(820.80078125\). Between \(p_2\) and \(p_3\), the
periodic mean square \(\rho\) grows from about \(1.36681670\) to
\(6.27677617\), while the harmonic slope \(P\) changes only from about
\(0.00580566\) to \(0.00577632\). The denominator-weighted theorem improves
the discrepancy constant enough for this finite proof, but it does not
control those adverse support and coefficient-growth trends.

If the same full basis shape were attempted once more, its bookkeeping alone
would permit support as high as \(8\cdot65536=524288\). This certificate
contains no fourth shell, weight list, \(p_4\), or fourth contraction. The
separate [fourth-step certificate](nyman-rebased-fourth-step-v1.md) later
certified a reduced basis retaining only the first fourth-shell multiplier,
with support through 65,536. The larger number remains a full-basis
support-budget warning, not a prediction.

The current binary64 prefix path forms \(M(M+1)\) under the exactness guard
\(T(T+1)<2^{53}\). A possible future implementation could divide a scaled
square sequentially by \(M\) and \(M+1\), then enlarge the proved rounding
factor by one operation. That backend change is not implemented or audited
here and cannot be used to support the present claim.

There is no all-scale lemma, no rate estimate, and no justified inference
from the four now-published finite vectors to \(E(p_k)\to0\).

## 9. Reproduction and commitments

Replay the frozen artifact:

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_rebased_third_step_certificate.py `
  --verify --artifact results\nyman-rebased-third-step-v1.json
```

Relevant files:

- [certificate generator and verifier](../tools/generate_nyman_rebased_third_step_certificate.py);
- [focused tests](../tests/test_nyman_rebased_third_step_certificate.py);
- [frozen artifact](../results/nyman-rebased-third-step-v1.json);
- [fixed-modulus refinement](nyman-fixed-modulus-refinement-v1.md);
- [separate fourth-step proof](nyman-rebased-fourth-step-v1.md);
- [parent proof note](nyman-rebased-schur-v1.md).

Hashes:

```text
artifact payload SHA-256:  bbafc2c92acb954c9c9e550f781612965bc5c76fe59ebd4a1adc5b596c3bee02
raw LF file SHA-256:       e5bd17bf94f2ce95611ad09bf75acaecc156db003ce4481b4a6e41c3504340bd
canonical artifact SHA:   811ebaf03089288c99d8155171a07e4b5e5731d53678265354ed2738157fa1ca
```

The Riemann Hypothesis remains `UNRESOLVED`.
