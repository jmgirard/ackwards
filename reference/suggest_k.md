# Suggest a maximum number of factors for bass-ackwards analysis

Runs five complementary selection criteria and reports their
recommendations. No single criterion is definitive; the goal is a
consensus range to inform your choice of `k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

## Usage

``` r
suggest_k(
  data,
  k_max = NULL,
  cor = "pearson",
  n_obs = NULL,
  n_iter = 20L,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  A data frame or numeric matrix (items in columns, observations in
  rows). Alternatively, a pre-computed **correlation matrix** may be
  supplied (a square, symmetric, numeric matrix with unit diagonal).
  When a correlation matrix is supplied, `n_obs` is required (PA and VSS
  need N), the `cor` argument is ignored, and the Comparison Data (CD)
  criterion is skipped (CD requires raw item distributions for
  resampling).

- k_max:

  Maximum number of components to test. Defaults to
  `min(ncol(data) - 1, 8)`. Increase if you expect a deeper hierarchy.

- cor:

  Correlation basis: `"pearson"` (default) or `"spearman"`. Should match
  the `cor` argument you plan to use in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  Ignored when `data` is a correlation matrix (the basis is already
  fixed).

- n_obs:

  Number of observations. Required when `data` is a pre-computed
  correlation matrix (PA and VSS need N). Ignored when raw data are
  supplied (N is determined from `nrow(data)`).

- n_iter:

  Number of Monte Carlo iterations for parallel analysis. Default `20`.
  Reduce to `5` for fast/exploratory runs; increase to `100+` for
  publication.

- seed:

  Integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  Comparison Data (CD) step. `NULL` (default) uses the current RNG
  state. **Note:** the parallel-analysis step uses
  [`psych::fa.parallel()`](https://rdrr.io/pkg/psych/man/fa.parallel.html),
  which does not respond reliably to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) – PA simulation
  results will vary across calls regardless of `seed`.

- ...:

  Reserved for future arguments.

## Value

An object of class `"suggest_k"`. Print it for a formatted summary; call
[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
on it for a diagnostic scree/criteria plot. The list contains:

- k_parallel_pc:

  Recommended k from PC-based parallel analysis.

- k_parallel_fa:

  Recommended k from FA-based parallel analysis (`NA_integer_` if no FA
  factor exceeded the random threshold).

- k_map:

  Recommended k from MAP.

- k_vss1:

  Recommended k from VSS complexity-1.

- k_vss2:

  Recommended k from VSS complexity-2.

- k_cd:

  Recommended k from Comparison Data (`NA_integer_` when EFAtools is not
  installed or CD fails).

- cd_available:

  Logical; `TRUE` when EFAtools was found and CD ran successfully.

- criteria:

  Data frame with one row per k: `k`, `ev_obs` (observed PC eigenvalue),
  `ev_obs_fa` (observed FA eigenvalue), `pa_pc_quant` / `pa_fa_quant`
  (95th-pct simulated eigenvalue for each basis), `pa_pc_suggested` /
  `pa_fa_suggested` (logical retention), `map`, `vss1`, `vss2`.

- k_max, n_obs, n_vars, cor:

  Metadata.

## Details

**Criteria computed:**

- **PA-PC** (Horn 1965, PC basis) – parallel analysis on
  principal-component eigenvalues. Compares observed eigenvalues to
  those from random correlation matrices; suggests retaining components
  whose eigenvalues exceed the 95th percentile of chance. Tends to
  overextract; treat as an upper bound.

- **PA-FA** (Horn 1965, FA basis) – parallel analysis using
  common-factor eigenvalues. More conservative than PA-PC and the better
  match for the EFA and ESEM engines in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

- **MAP** (Velicer 1976) – Minimum Average Partial criterion. Finds the
  k that minimises the average squared partial correlation remaining
  after extracting k components. Usually conservative.

- **VSS-1 / VSS-2** (Revelle & Rocklin 1979) – Very Simple Structure fit
  at complexities 1 and 2. Finds the k maximising the fit of a very
  simple loading structure.

- **CD** (Ruscio & Roche 2012; optional) – Comparison Data. Resamples
  from the observed item distributions to generate comparison eigenvalue
  profiles; retains factors until adding one no longer improves RMSE
  beyond chance. Requires the EFAtools package (install separately).

PA (both bases), MAP, and VSS share the same correlation matrix as
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
CD operates on the raw data matrix directly (required for resampling)
and is skipped gracefully when EFAtools is not installed.

## Interpreting the output

`k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
is a **maximum depth**, not a claim that exactly k factors exist. Users
commonly set k one or two levels above the consensus to watch
higher-level factors fragment – this is a feature of the method, not
overextraction.

## A note on overextraction

PA-PC in particular tends to recommend more factors than replicate
across independent samples, especially with correlated items (Forbes,
2023). PA-FA and CD are more conservative. Treat the full set of
criteria as a range: the true k is likely somewhere in the middle.

## References

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika*, 30, 179–185.

Revelle, W., & Rocklin, T. (1979). Very simple structure: An alternative
procedure for estimating the optimal number of interpretable factors.
*Multivariate Behavioral Research*, 14(4), 403–414.

Ruscio, J., & Roche, B. (2012). Determining the number of factors to
retain in an exploratory factor analysis using comparison data of a
known factorial structure. *Psychological Assessment*, 24(2), 282–292.

Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. *Psychometrika*, 41, 321–327.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)

## Examples

``` r
# \donttest{
sk <- suggest_k(bfi25)
#> ℹ Running parallel analysis (20 iterations, PC + FA)...
#> ✔ Running parallel analysis (20 iterations, PC + FA)... [190ms]
#> 
#> ℹ Running MAP and VSS...
#> CD: 125 rows with missing values removed (875 complete cases used).
#> ✔ Running MAP and VSS... [82ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [9.8s]
#> 
sk
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 25
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#> k = 1: PA-PC ✔ PA-FA ✔ MAP 0.0246 VSS-1 0.5121 VSS-2 0.0000 CD ✔
#> k = 2: PA-PC ✔ PA-FA ✔ MAP 0.0190 VSS-1 0.5761 VSS-2 0.6636 CD ✔
#> k = 3: PA-PC ✔ PA-FA ✔ MAP 0.0172 VSS-1 0.5969 VSS-2 0.7288 CD ✔
#> k = 4: PA-PC ✔ PA-FA ✔ MAP 0.0164 VSS-1 0.6198* VSS-2 0.7781 CD ✔
#> k = 5: PA-PC ✔ PA-FA ✔ MAP 0.0160* VSS-1 0.5730 VSS-2 0.7912* CD ✔
#> k = 6: PA-PC - PA-FA ✔ MAP 0.0170 VSS-1 0.5592 VSS-2 0.7530 CD ✔*
#> k = 7: PA-PC - PA-FA - MAP 0.0201 VSS-1 0.5698 VSS-2 0.7290 CD -
#> k = 8: PA-PC - PA-FA - MAP 0.0231 VSS-1 0.5615 VSS-2 0.7252 CD -
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 5
#> • PA-FA: k <= 6
#> • MAP: k = 5
#> • VSS-1: k = 4
#> • VSS-2: k = 5
#> • CD: k = 6
#> Consensus range: k = 4-6
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.
autoplot(sk)


# Faster exploratory run
suggest_k(bfi25, k_max = 6, n_iter = 5)
#> ℹ Running parallel analysis (5 iterations, PC + FA)...
#> ✔ Running parallel analysis (5 iterations, PC + FA)... [90ms]
#> 
#> ℹ Running MAP and VSS...
#> CD: 125 rows with missing values removed (875 complete cases used).
#> ✔ Running MAP and VSS... [74ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [8.9s]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 25
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-6
#> 
#> ── Criteria (k = 1-6) ──
#> 
#> k = 1: PA-PC ✔ PA-FA ✔ MAP 0.0246 VSS-1 0.5121 VSS-2 0.0000 CD ✔
#> k = 2: PA-PC ✔ PA-FA ✔ MAP 0.0190 VSS-1 0.5761 VSS-2 0.6636 CD ✔
#> k = 3: PA-PC ✔ PA-FA ✔ MAP 0.0172 VSS-1 0.5969 VSS-2 0.7288 CD ✔
#> k = 4: PA-PC ✔ PA-FA ✔ MAP 0.0164 VSS-1 0.6198* VSS-2 0.7781 CD ✔
#> k = 5: PA-PC ✔ PA-FA ✔ MAP 0.0160* VSS-1 0.5730 VSS-2 0.7912* CD ✔
#> k = 6: PA-PC - PA-FA ✔ MAP 0.0170 VSS-1 0.5592 VSS-2 0.7530 CD ✔*
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 5
#> • PA-FA: k <= 6
#> • MAP: k = 5
#> • VSS-1: k = 4
#> • VSS-2: k = 5
#> • CD: k = 6
#> Consensus range: k = 4-6
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.

# Correlation-matrix input (CD is skipped; n_obs required)
R <- cor(bfi25, use = "pairwise.complete.obs")
suggest_k(R, n_obs = 875L)
#> ℹ Comparison Data (CD) is skipped when a correlation matrix is supplied (CD
#>   requires raw item distributions for resampling).
#> ℹ Running parallel analysis (20 iterations, PC + FA)...
#> ✔ Running parallel analysis (20 iterations, PC + FA)... [190ms]
#> 
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [90ms]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 25
#> n: 875
#> Basis: (user-supplied matrix)
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#> k = 1: PA-PC ✔ PA-FA ✔ MAP 0.0246 VSS-1 0.5121 VSS-2 0.0000
#> k = 2: PA-PC ✔ PA-FA ✔ MAP 0.0190 VSS-1 0.5761 VSS-2 0.6636
#> k = 3: PA-PC ✔ PA-FA ✔ MAP 0.0172 VSS-1 0.5969 VSS-2 0.7288
#> k = 4: PA-PC ✔ PA-FA ✔ MAP 0.0164 VSS-1 0.6198* VSS-2 0.7781
#> k = 5: PA-PC ✔ PA-FA ✔ MAP 0.0160* VSS-1 0.5730 VSS-2 0.7912*
#> k = 6: PA-PC - PA-FA ✔ MAP 0.0170 VSS-1 0.5592 VSS-2 0.7530
#> k = 7: PA-PC - PA-FA - MAP 0.0201 VSS-1 0.5698 VSS-2 0.7290
#> k = 8: PA-PC - PA-FA - MAP 0.0231 VSS-1 0.5615 VSS-2 0.7252
#> + CD skipped (requires raw data; not available for matrix input).
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 5
#> • PA-FA: k <= 6
#> • MAP: k = 5
#> • VSS-1: k = 4
#> • VSS-2: k = 5
#> Consensus range: k = 4-6
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.
# }
```
