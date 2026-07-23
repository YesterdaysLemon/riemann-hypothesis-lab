# Six-scale fixed-width large-sieve certificate

## Result

The frozen direct-gain scout contains one `K=64` multiplier at each

```text
N = 8, 16, 32, 64, 128, 256.
```

After exact 9-bit shell and multiplier rounding, exact 16-bit old-vector
rounding, exact balance, and complete infinite-tail certification, every one
of the six corrections has strictly positive direct gain.  All six use the
single cutoff rule

\[
T_N=32768N.
\]

| `N` | `T_N` | support limit `Q` | complete-gain enclosure | certified lower bound |
|---:|---:|---:|---:|---:|
| 8 | 262,144 | 1,024 | `[0.00677209911997..., 0.00683712021532...]` | `1/150` |
| 16 | 524,288 | 2,048 | `[0.00337058812839..., 0.00346095932980...]` | `1/300` |
| 32 | 1,048,576 | 4,032 | `[0.00166751969021..., 0.00180026855100...]` | `1/625` |
| 64 | 2,097,152 | 8,192 | `[0.00141513662112..., 0.00162096679226...]` | `1/750` |
| 128 | 4,194,304 | 16,384 | `[0.00065161169910..., 0.00102882678035...]` | `1/1600` |
| 256 | 8,388,608 | 32,768 | `[0.00023169645419..., 0.00092036206889...]` | `1/5000` |

In particular, all six satisfy the clean finite statement

\[
G_N>\frac1{20N}.
\]

The machine-readable evidence is
[nyman-large-sieve-scaling-v1.json](../results/nyman-large-sieve-scaling-v1.json).
The six signs are `CERTIFIED_FINITE`; interpretation of their trend is
`EXPLORATORY`.  The Riemann Hypothesis remains `UNRESOLVED`.

## Why this test matters

The first complete-tail certificate proved one contraction at `N=256`.  The
[Fourier/Farey theorem](nyman-large-sieve-tail-v1.md) then replaced the
pairwise-LCM tail constant by

\[
C_{\rm LS}=Q(Q-1)(\rho+\tau),
\]

where `rho` and `tau` are the exact old and new periodic residual mean
squares.  That formula is cheap enough to test one fixed-width construction
at every stored scale rather than at a single endpoint.

The suite asks a deliberately uniform finite question:

1. select the already frozen direct-gain `K=64` cell for each `N`;
2. round and balance it by exactly the same rule at every scale;
3. take `T_N=32768N`, without tuning a cutoff per result;
4. enclose the prefix and every omitted interval;
5. require a strict rational gain threshold for every row.

All six rows close.

## Exact construction

For each stored old candidate `p^(N)`:

- the old coefficients are rounded to the `2^-16` grid;
- Arb proves the first `N-1` ideal-shell rounding bins on the `2^-9` grid;
- the final shell coefficient is imposed by exact zero-sum balance;
- `c_k/k` is rounded on the `2^-9` grid for `k>=3`;
- `c_1=1`, and `c_2` is imposed by exact harmonic balance;
- the correction is the exact Dirichlet convolution of shell and multiplier,
  on a grid no finer than `2^-18`;
- both `sum a_n=0` and `sum a_n/n=0` are checked exactly.

The source candidate, scout artifact, selected cell, and all reconstructed
vector commitments are hash-bound.  The compact suite artifact stores a
canonical hash of each full sparse dyadic vector record rather than repeating
more than twenty thousand coefficients.

## Complete-tail calculation

For each row, the prefix evaluator uses exact integer divisor recurrences.
Only the final divisions and reductions use binary64, with the same exact
rational forward-error analysis as the single-cell certificate.  Arb encloses
the logarithmic contribution.

The omitted periodic tail uses the theorem

\[
\left|
\sum_{M\in I}\bigl(h_M-\operatorname{mean}(h)\bigr)
\right|
\le Q(Q-1)(\rho+\tau)
\]

for every integer interval `I`.  Abel summation gives periodic radius

\[
\frac{Q(Q-1)(\rho+\tau)}{(T_N+1)(T_N+2)},
\]

and the small logarithmic remainder receives its separate exact radius.  The
reference theorem is Montgomery and Vaughan,
[“The large sieve,” Mathematika 20 (1973), 119-134](https://doi.org/10.1112/S0025579300004708).

Generation uses blocks of `2^18` and 256-bit Arb.  Replay changes to blocks of
`2^17` and 384-bit Arb, reconstructs every vector, and requires the stored
generation balls to contain the replay balls.

## What the finite scaling data says

The cutoff/support ratio is exactly `256` in five rows and approximately
`260.06` at `N=32`, where the exact top coefficient vanishes and lowers `Q`
from 4096 to 4032.  Thus the computational horizon is essentially a fixed
multiple of support across the whole suite.

The difficult energy, however, grows:

| `N` | `rho+tau` | periodic tail radius | radius / prefix gain, approx. |
|---:|---:|---:|---:|
| 8 | `2.132718...` | `3.2511e-5` | `0.0048` |
| 16 | `2.962747...` | `4.5186e-5` | `0.013` |
| 32 | `4.490229...` | `6.6374e-5` | `0.038` |
| 64 | `6.745476...` | `1.0292e-4` | `0.068` |
| 128 | `12.361347...` | `1.8861e-4` | `0.224` |
| 256 | `22.566891...` | `3.4433e-4` | `0.595` |

The finite evidence therefore has two sides:

- the Fourier/Farey bound is strong enough to turn every stored fixed-width
  signal into a complete proof of positive gain;
- at a fixed `T/Q`, tail headroom visibly shrinks, so these six points do not
  support a constant-headroom extrapolation.

## Why this is not a recurrence

The six rows are cross-sectional and independently optimized.  They are not
successive states of one construction.

A correction built at old dimension `N` may be supported as far as

\[
Q\le2NK=128N.
\]

After applying it, a genuine next recurrence step would have to use that new
residual and its enlarged support.  The stored table instead moves from `N`
to `2N` and reloads a separately optimized old candidate.  Consequently, the
suite proves none of the following:

- that a next multiplier exists for the residual produced by the prior row;
- that the finite lower bounds persist for arbitrary `N`;
- that `rho+tau` has a sufficiently small asymptotic growth rate;
- that the product of repeated contractions forces a Nyman distance to zero.

The clean finite inequality `G_N>1/(20N)` is true only for the six declared
values.  Extending its quantifier would be an unsupported extrapolation.

## Reproduction

Generate the suite:

```powershell
.\.venv\Scripts\python.exe `
  tools\generate_nyman_large_sieve_scaling_certificate.py
```

Replay every cell with the independent settings:

```powershell
.\.venv\Scripts\python.exe `
  tools\generate_nyman_large_sieve_scaling_certificate.py --verify
```

Run focused tests:

```powershell
.\.venv\Scripts\python.exe -m pytest `
  tests\test_nyman_balanced_tail.py `
  tests\test_nyman_large_sieve_certificate.py `
  tests\test_nyman_large_sieve_scaling_certificate.py -q
```

The frozen artifact hashes are:

```text
payload SHA-256: e974817e59e33efd84b18aa6e81d0592e96c7bb5fd5584dfb83f7763de8c0320
raw-file SHA-256: 02886cec84010235b3a6e564b7aba8d42e57bbcf2f4e28783dc91894d01f9cce
whole canonical artifact SHA-256: b8a22d0d419b282b70c705c6df0edad77dcc3a69f91189ae11c2035567d3fcdb
```

## Next theorem target

The suite narrows the live problem to a nested construction.  One needs an
explicit rule that consumes the residual produced at one step and returns the
next balanced correction, together with uniform bounds comparing:

\[
\text{complete gain},\qquad Q^2(\rho+\tau),\qquad
\text{and the current residual norm}.
\]

Without that nested rule and its uniform estimates, the six contractions are
valuable finite evidence and nothing more.  They do not prove or disprove RH.
