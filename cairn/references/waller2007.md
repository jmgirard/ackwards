# waller2007

**Full citation.** Waller, N. (2007). A general method for computing
hierarchical component structures by Goldberg's Bass-Ackwards method.
*Journal of Research in Personality, 41*(4), 745–752.
https://doi.org/10.1016/j.jrp.2006.08.005

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `351a916`) from `cairn/references/sources/waller2007.pdf` (local only;
gitignored). Pagination: journal pages (745–752).
Extraction: unverified — first pass, values not re-read against the source — observed 2026-07-19.

## Why this is a primary source

The mathematical foundation of ackwards' **algebraic edge route** (CLAUDE.md
Invariant 1). Waller proves that the end result of a Bass-Ackwards analysis —
the correlations between component scores at different levels — can be
computed **from the correlation matrix alone**, without ever materializing
scores. Consequences he draws that we inherit:

- works on any dataset for which only `R` is available (raw data missing) —
  exactly why `ackwards()` accepts a correlation matrix as input;
- faster (matters for Monte Carlo);
- for the *factor* (not component) model, avoids the bias of correlating
  **estimated** factor scores.

## The algebra (orthogonal components case)

With `Z` the standardized unrotated PC scores and `T_i` the orthonormal
varimax transformation for the i-level solution (`V_i = Z_i T_i`):

```
Cor(V_i, V_j) = T_i' S T_j        (his Eq. 14)
```

where `S = [I | 0]` is a selection matrix. For the FUPC vs. level-2 case the
correlations are just the first row of `T_2` (Eq. 9–10). §3 extends to oblique
rotations via `T_i^{-1} S T_j'^{-1}` (out of scope for us — DESIGN §2).

**Correspondence to `compute_edges()`.** Our `W'RW / sqrt(diag(...))` route is
the same identity in weight-matrix form: for standardized PCA scores,
`W_i = Q Λ^{-1/2} T_i` (his Eqs. 6–7), so `W_i' R W_j` collapses to
`T_i' S T_j`. We use the `W'RW` form because it also covers EFA/ESEM scoring
weights (tenBerge, Bartlett, …) where the transformation-matrix shortcut does
not apply, and we standardize by the real score SDs because non-PCA scores are
not unit-variance (Invariant 1).

## The §4 caveat — why the scores cross-check exists

Waller notes that if **estimated** factor scores (regression, Bartlett, …) are
substituted into the score-correlation procedure, the results will **not**
agree with the algebraic values (factor-score indeterminacy; Guttman 1955,
McDonald & Mulaik 1979); the algebra is the *correct* between-level
correlation and estimated-score correlations only approximate it. This is the
design rationale behind keeping both routes and testing their agreement for
linear engines (Invariant 2) — the two agree in our implementation because
`compute_edges()` correlates the *same* linear composites the scores route
materializes, standardized by their real SDs, rather than treating score
estimates as the latent factors themselves.

## Appendix A

A complete R function `BASS(R, maxP, Print)` — eigen → varimax per level →
`t(T[[i-1]]) %*% S %*% T[[i]]` cross-level correlations. Historically the
first public R implementation of the method; sign convention is
column-sum-positive on the unrotated eigenvectors (contrast with our
primary-parent alignment, Invariant 4). Not used as an oracle (superseded by
Forbes's implementation, which is the test-backed contract — [[forbes2023]]).
