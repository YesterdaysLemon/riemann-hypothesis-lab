# Literature map

This is a curated map of primary or official sources, not a claim that every
paper's approach is equally promising.

## Status and canonical formulation

- Clay Mathematics Institute, [Millennium Prize Problems](https://www.claymath.org/millennium-problems/):
  RH remains listed among the unsolved problems as of 2026-07-22.
- Enrico Bombieri, [official problem description](https://www.claymath.org/wp-content/uploads/2022/05/riemann.pdf):
  canonical analytic formulation, explicit formula, and context.
- NIST DLMF, [Riemann zeta zeros](https://dlmf.nist.gov/25.10): definitions,
  symmetry, Hardy Z, and the Riemann-Siegel formula.

## Rigorous frontiers

- Platt and Trudgian,
  [The Riemann hypothesis is true up to `3*10^12`](https://arxiv.org/abs/2004.09765):
  the lowest 12,363,153,437,138 positive-ordinate zeros are on the critical
  line through height 3,000,175,332,800. This finite theorem does not prove RH.
- Pratt, Robles, Zaharescu, and Zeindler,
  [More than five-twelfths of the zeros of zeta are on the critical line](https://doi.org/10.1007/s40687-019-0199-8):
  unconditional lower proportion at least 0.417293962.
- Guth and Maynard,
  [New large value estimates for Dirichlet polynomials](https://arxiv.org/abs/2405.20552):
  a strong zero-density frontier, still compatible with exceptional zeros.
- Mossinghoff, Trudgian, and Yang,
  [Explicit zero-free regions for the Riemann zeta-function](https://arxiv.org/abs/2212.06867):
  explicit exclusions near `Re(s)=1`, not the whole critical strip.
- Bellotti, Trudgian, and Yang,
  [Zero-free regions inspired by work of Heath-Brown](https://arxiv.org/abs/2603.21490):
  a 2026 preprint proving the classical-shaped region
  `sigma >= 1 - 1/(4.896 log t)` for `t >= 3`; still only a neighborhood of
  `Re(s)=1`.

## Equivalent criteria used here

- Bombieri,
  [Remarks on Weil's quadratic functional in the theory of prime numbers, I](https://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0).
- Connes and Consani,
  [Spectral Triples and Zeta-Cycles](https://arxiv.org/abs/2106.01715): the
  lower-semicontinuous semilocal Weil form, its Laurent-polynomial core, and
  the finite-compression lower-bound limit.
- Lagarias, [An Elementary Problem Equivalent to the Riemann Hypothesis](https://arxiv.org/abs/math/0008177).
- Li, [The positivity of a sequence of numbers and the Riemann hypothesis](https://doi.org/10.1006/jnth.1997.2137).
- Bombieri and Lagarias,
  [Complements to Li's criterion](https://doi.org/10.1006/jnth.1999.2392).
- Baez-Duarte,
  [A strengthening of the Nyman-Beurling criterion](https://arxiv.org/abs/math/0202141).
- Baez-Duarte, Balazard, Landreau, and Saias,
  [Sur l'autocorrelation multiplicative de la fonction "partie fractionnaire"](https://arxiv.org/abs/math/0306251):
  rational autocorrelation, reciprocity, and the Vasyunin-sum formula used by
  the frozen natural-distance experiment.
- Burnol,
  [A lower bound in an approximation problem involving the zeros of the Riemann zeta function](https://arxiv.org/abs/math/0103058).
- Chen and Qi,
  [The best bounds of harmonic sequence](https://arxiv.org/abs/math/0306233):
  the harmonic-number enclosure used for the independent `1-gamma`
  normalization audit.
- Balazard and de Roton,
  [Sur un critere de Baez-Duarte pour l'hypothese de Riemann](https://arxiv.org/abs/0812.1689):
  an RH-conditional explicit natural-dilate coefficient family and upper rate.
- Bettin, Conrey, and Farmer,
  [An optimal choice of Dirichlet polynomials for the Nyman--Beurling criterion](https://arxiv.org/abs/1211.5191):
  the log-tapered Mobius family and its optimal asymptotic under RH plus a
  zero-derivative moment assumption.
- Rodgers and Tao,
  [The de Bruijn-Newman constant is non-negative](https://arxiv.org/abs/1801.05914).
- Griffin, Ono, Rolen, Thorner, Tripp, and Wagner,
  [Jensen Polynomials for the Riemann Xi Function](https://arxiv.org/abs/1910.01227).

Every criterion above retains an infinite or universal quantifier. Finite
positivity, finite approximation, and finite inequality checks do not prove RH.

## Current reproducible directions

- Connes,
  [The Riemann Hypothesis: Past, Present and a Letter Through Time](https://arxiv.org/abs/2602.04022):
  2026 finite approximants using small primes; the unproved convergence bridge
  is explicitly the missing step.
- Connes, Consani, and Moscovici,
  [Zeta zeros and prolate wave operators](https://arxiv.org/abs/2310.18423) and
  [Zeta Spectral Triples](https://arxiv.org/abs/2511.22755): finite-prime
  spectral structures whose missing convergence theorem remains explicit.
- Tao, Trudgian, and Yang,
  [New exponent pairs, zero density estimates, and zero additive energy estimates: a systematic approach](https://arxiv.org/abs/2501.16779)
  and [ANTEDB source](https://github.com/teorth/expdb): a model for executable
  provenance and machine-checkable bound propagation.
- Groskin,
  [High-Precision Approximation of Riemann Zeros via the Truncated Weil Form](https://arxiv.org/abs/2605.20224)
  and [A finite Guinand-Weil dictionary and archimedean tail order](https://arxiv.org/abs/2607.02828):
  current unreviewed preprints with public artifacts directly relevant to the
  finite-form track. They are exploratory inputs, not validation or an RH claim.
- FLINT,
  [rigorous zeta-zero and Turing APIs](https://flintlib.org/doc/acb_dirichlet.html):
  the trusted backend for baseline finite certificates.
