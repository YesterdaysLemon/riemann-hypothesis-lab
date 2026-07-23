# Complete-tail certificate for one balanced Nyman correction

## Verdict

One explicit exact-dyadic balanced correction lowers the complete weighted
`L^2` energy of one explicit `N=256` Nyman vector by more than

\[
\boxed{\frac{7}{50000}=0.00014.}
\]

The old coefficients are rounded to the grid `2^-16`. The correction comes
from an `N=256`, `K=16` shell/multiplier construction, has 3,151 nonzero
coefficients on the grid `2^-18`, and is supported through index 8,192. The
certificate includes every interval from 1 to infinity: a finite prefix
through `T=2^26` is enclosed by exact integer recurrences, an explicit
binary64 forward-error bound, and Arb; the omitted infinite tail is enclosed
analytically by exact periodic means and Abel summation.

The resulting stored enclosure is

\[
0.0001420236276423717967\ldots
\;<\;\text{direct gain}\;<\;
0.0003677819654619422937\ldots.
\]

This is `CERTIFIED_FINITE`. It proves one strict finite-dimensional
contraction, not a uniform contraction theorem. It does not prove or disprove
the Riemann Hypothesis, whose repository status remains `UNRESOLVED`.

## Exact vectors and the claim

Use the natural Nyman functions

\[
e_n(t)=\{t/n\},\qquad
\chi(t)=\mathbf 1_{[1,\infty)}(t),
\]

in `L^2((0,infinity),dt/t^2)`. For the old exact vector `p`, let

\[
r(t)=\chi(t)-\sum_m p_m e_m(t).
\]

The added coefficient vector `a` satisfies the two exact balances

\[
\sum_n a_n=0,
\qquad
\sum_n\frac{a_n}{n}=0.
\]

Put `g=sum a_n e_n`. The new residual is `r-g`, so its complete direct gain
is

\[
\Delta=\lVert r\rVert^2-\lVert r-g\rVert^2
=2\langle r,g\rangle-\lVert g\rVert^2.
\]

The machine-checked claim is exactly

\[
\boxed{\Delta>7/50000.}
\]

Equivalently, the explicit coefficient vector `p+a`, supported through
index 8,192, has squared residual norm more than `7/50000` below that of the
explicit 256-term vector `p`.

The exploratory scout selected this direction but is not part of the proof.
The certificate rederives and stores every exact coefficient. Arb proves that
each of the first 255 ideal-shell values lies strictly inside its selected
9-bit rounding bin; the final shell coefficient and the multiplier coefficient
`c_2` are then imposed by exact balance. The truth of the contraction depends
only on the stored dyadic vectors, not on the heuristic quality of their
selection.

## Exact interval identity

On `[M,M+1)`, harmonic balance makes `g` constant. Define

\[
g_M=-\sum_n a_n\left\lfloor\frac{M}{n}\right\rfloor,
\qquad
P_1=\sum_m\frac{p_m}{m},
\qquad
q_M=1+\sum_m p_m\left\lfloor\frac{M}{m}\right\rfloor.
\]

Then `r(t)=q_M-P_1t`, and direct integration gives

\[
\boxed{
\Delta=\sum_{M\ge1}
\left[
\frac{2g_Mq_M-g_M^2}{M(M+1)}
-2P_1g_M\log\left(1+\frac1M\right)
\right].
}
\]

This is the identity evaluated by the prefix and tail certificates.

## Rigorous prefix through `2^26`

Write the correction over denominator `A=2^18` and the old vector over
denominator `P=2^16`:

\[
a_n=A_n/A,\qquad p_m=P_m/P.
\]

The exact integer recurrences are

\[
G_M=-\sum_n A_n\left\lfloor\frac{M}{n}\right\rfloor,
\qquad
Q_M=P+\sum_mP_m\left\lfloor\frac{M}{m}\right\rfloor.
\]

Since `A/P=4`, the rational interval numerator is

\[
H_M=8G_MQ_M-G_M^2,
\]

and the rational prefix equals

\[
\sum_{M\le T}\frac{H_M}{2^{36}M(M+1)}.
\]

Every divisor jump, cumulative sum, product, and numerator is first evaluated
exactly in signed 64-bit integers with preflight and runtime overflow checks.
At `T=2^26`, the observed maxima are

| Exact integer | Maximum absolute value |
|---|---:|
| `G_M` | 4,568,415 |
| `Q_M` | 12,322,065,530 |
| `H_M` | 421,449,858,668,025,872 |

The integer numerators are converted to pinned IEEE-754 binary64 and reduced
in 64 blocks of size `2^20`. With unit roundoff `u=2^-53`, term evaluation is
bounded by

\[
\tau=(1+u)^2-1,
\]

and a reduction of `k+1` values is bounded by

\[
\gamma_k=\frac{ku}{1-ku}.
\]

The implementation recovers an upper bound for the exact absolute mass of
the computed terms and adds the term, within-block, and final-block errors as
exact rational numbers. For this prefix it obtains

| Quantity | Certified value |
|---|---:|
| rational midpoint | `-0.00011437378224159554` |
| rational forward-error radius | `1.5038252322311281e-11` |
| absolute-term-sum upper bound | `0.12917747154779183` |

The logarithmic sum is not accumulated 67 million times. If
`q_n=floor(T/n)`, Abel summation gives the exact compression

\[
\sum_{M\le T}g_M\log\frac{M+1}{M}
=g_T\log(T+1)
+\sum_na_n\left[q_n\log n+\log\Gamma(q_n+1)\right].
\]

Arb evaluates the resulting 3,152 logarithm/log-gamma terms with outward
rounding. Its contribution to the gain is

\[
0.0003694320313862137410841428507090\ldots.
\]

Combining the rational and logarithmic parts encloses the prefix near
`0.0002550582`, with radius below `6.42e-11`. The verifier repeats this at
384 bits with blocks of size `2^19`, rather than reusing the generation
partition.

## Exact complete-tail enclosure

For an integer interval define the centered periodic sawtooth

\[
\phi_n(M)=\{M/n\}-\frac{n-1}{2n}.
\]

Both balances imply

\[
g_M=\sum_na_n\phi_n(M).
\]

Let

\[
C_0=1-\frac12\sum_mp_m,
\qquad
\bar r_M=C_0-\sum_mp_m\phi_m(M).
\]

The elementary complete-period covariance is

\[
\operatorname{mean}(\phi_m\phi_n)
=\frac{\gcd(m,n)^2-1}{12mn}.
\]

Consequently, with `J_2` denoting the second Jordan totient,

\[
\mu=\operatorname{mean}(g_M^2)
=\frac1{12}\sum_dJ_2(d)
\left(\sum_{d\mid n}\frac{a_n}{n}\right)^2,
\]

and

\[
\nu=\operatorname{mean}(g_M\bar r_M)
=-\frac1{12}\sum_dJ_2(d)
\left(\sum_{d\mid m}\frac{p_m}{m}\right)
\left(\sum_{d\mid n}\frac{a_n}{n}\right).
\]

These are exact rational `O(X log X)` divisor sums; no period of length equal
to a huge least common multiple is enumerated.

The old intercept decomposes as

\[
q_M=\bar r_M+P_1(M+1/2).
\]

Thus the exact interval gain is a periodic main term plus

\[
-2P_1g_M\delta_M,
\qquad
\delta_M=\log(1+1/M)-\frac1{M+1}-\frac1{2M(M+1)}.
\]

The sign and bound used by the checker are

\[
-\frac1{6M^3}\le\delta_M\le0.
\]

For completeness, set `x=1/M` and

\[
d(x)=\log(1+x)-\frac{x(2+x)}{2(1+x)}.
\]

Then `d(0)=0` and

\[
d'(x)=-\frac{x^2}{2(1+x)^2}\le0,
\qquad
\left(d(x)+\frac{x^3}{6}\right)'
=\frac{x^2}{2}\left(1-\frac1{(1+x)^2}\right)\ge0,
\]

which proves both endpoints of the inequality without a numerical
approximation.

For ordered-pair sums define

\[
\Lambda_a=\sum_{m,n}|a_ma_n|\operatorname{lcm}(m,n),
\]

\[
\Lambda_\times=\sum_{m,n}|p_ma_n|\operatorname{lcm}(m,n),
\qquad
A_0=\sum_n|a_n|,
\qquad
A_1=\sum_n n|a_n|.
\]

After subtracting its mean, `phi_m phi_n` has absolute value at most `1/3`.
Its period divides `lcm(m,n)`, so complete periods cancel and every remaining
cyclic partial sum is bounded by `lcm(m,n)/3`. Similarly, centered `phi_n`
has absolute value below `1/2`, period `n`, and partial-sum bound `n/2`.
Therefore the zero-mean periodic part of `2g_M \bar r_M-g_M^2` has arbitrary
partial sums bounded by

\[
C_{\rm per}=\frac{\Lambda_a}{3}
+\frac{2\Lambda_\times}{3}
+|C_0|A_1.
\]

Abel summation of `1/(M(M+1))` and the delta bound now give the full omitted
tail after `T`:

\[
\boxed{
\Delta_{>T}=\frac{2\nu-\mu}{T+1}\;\mathbin{\pm}\;
\left[
\frac{C_{\rm per}}{(T+1)(T+2)}
+\frac{|P_1|A_0}{12T^2}
\right].
}
\]

All quantities in this formula are stored as reduced exact fractions. For the
frozen vector they are approximately

| Quantity | Value |
|---|---:|
| `mu` | `10.584429862109205` |
| `nu` | `0.076091410365943621` |
| `P_1` | `0.0028016899417424974` |
| `C_0` | `-0.0084228515625` |
| `A_0` | `656.08344268798828` |
| `A_1` | `1612448.3183746338` |
| `Lambda_a` | `1510367962502.7656` |
| `Lambda_cross` | `7359804495.7593803` |
| `C_per` | `508362537412.84094` |
| tail center | `-1.5545259246117957e-7` |
| tail radius | `0.00011287915373246328` |

Adding this complete tail to the prefix leaves a lower endpoint above
`0.00014202362`, which clears the simpler exact public claim `7/50000`.

## Artifact, replay, and trust boundary

The canonical artifact is
[nyman-balanced-full-tail-v1.json](../results/nyman-balanced-full-tail-v1.json).
It stores the exact old, shell, multiplier, and convolved correction vectors;
the complete exact tail data; the prefix error budget; outward Arb
enclosures; source hashes; and the unresolved-RH limitation.

Generate it with

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_balanced_tail_certificate.py
```

Replay it independently at the higher precision and different block
partition with

```powershell
.\.venv\Scripts\python.exe tools\generate_nyman_balanced_tail_certificate.py --verify
```

The replay validates the canonical artifact hash, the scout and candidate
hash chain, every exact vector and convolution, both balances, every exact
tail scalar, the overlap of independently reduced prefixes, and the strict
lower bound. The numerical trust base is Python-FLINT 0.9.0 for outward Arb
arithmetic and NumPy 2.3.5 with IEEE-754 round-to-nearest binary64. The
binary64 component has its own exact rational forward-error radius rather
than being trusted as exact arithmetic.

## What remains open

The result certifies one scale only. A proof of RH along this route still
needs an explicit construction at an unbounded sequence of scales and a
uniform lower bound whose accumulated contractions force the Nyman distance
to zero. In particular, a fixed- or polylogarithmic-width recurrence with
gain comparable to `c/log N` remains a viable theorem target. Nothing in this
finite certificate establishes that estimate.
