# Suggest a maximum number of factors for bass-ackwards analysis

Runs complementary factor-retention criteria and reports their
recommendations. No single criterion is definitive; the goal is a
consensus range to inform your choice of `k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

## Usage

``` r
suggest_k(
  data,
  k_max = NULL,
  criteria = c("pa_pc", "pa_fa", "map", "vss", "cd"),
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

  Maximum number of factors/components to *evaluate* when recommending a
  depth – not a depth itself. Defaults to `min(ncol(data) - 1, 8)`.
  Increase if you expect a deeper hierarchy. (Note: this is a distinct
  meaning from
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)'s
  `k_max`, which *is* the extraction depth; the two share a name because
  they're the same dial in the `suggest_k() -> ackwards()` workflow,
  just applied at different stages.)

- criteria:

  Character vector of criteria to compute. Any subset of
  `c("pa_pc", "pa_fa", "map", "vss", "cd")`; default is all five.
  `"pa_pc"` and `"pa_fa"` share one parallel-analysis call (both or
  neither are fast); `"map"` and `"vss"` share one VSS call. Criteria
  not requested are skipped entirely (no computation) and their `k_*`
  fields in the result are `NA`. `"vss"` selects both VSS-1 and VSS-2 as
  a unit.

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

  Recommended k from PC-based parallel analysis (`NA_integer_` when
  `"pa_pc"` not in `criteria`).

- k_parallel_fa:

  Recommended k from FA-based parallel analysis (`NA_integer_` if no FA
  factor exceeded the random threshold, or when `"pa_fa"` not in
  `criteria`).

- k_map:

  Recommended k from MAP (`NA_integer_` when `"map"` not in `criteria`).

- k_vss1:

  Recommended k from VSS complexity-1 (`NA_integer_` when `"vss"` not in
  `criteria`).

- k_vss2:

  Recommended k from VSS complexity-2 (`NA_integer_` when `"vss"` not in
  `criteria`).

- k_cd:

  Recommended k from Comparison Data (`NA_integer_` when `"cd"` not in
  `criteria`, EFAtools is not installed, or CD fails).

- cd_available:

  Logical; `TRUE` when `"cd"` was requested, EFAtools was found, and CD
  ran successfully.

- criteria:

  Data frame with one row per k: `k`, `ev_obs` (observed PC eigenvalue),
  `ev_obs_fa` (observed FA eigenvalue), `pa_pc_quant` / `pa_fa_quant`
  (95th-pct simulated eigenvalue for each basis), `pa_pc_suggested` /
  `pa_fa_suggested` (logical retention), `map`, `vss1`, `vss2`. Columns
  for non-requested criteria are `NA`.

- criteria_requested:

  Character vector of the criteria that were requested (and therefore
  computed).

- k_max, n_obs, n_vars, cor:

  Metadata.

## Details

**Criteria available (controlled by the `criteria` argument):**

- **PA-PC** (`"pa_pc"`) – Horn (1965) parallel analysis on PC
  eigenvalues. Compares observed eigenvalues to those from random
  correlation matrices; suggests retaining components whose eigenvalues
  exceed the 95th percentile of chance. Tends to overextract; treat as
  an upper bound.

- **PA-FA** (`"pa_fa"`) – Horn (1965) parallel analysis using
  common-factor eigenvalues. More conservative than PA-PC and the better
  match for the EFA and ESEM engines in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  Shares one
  [`psych::fa.parallel()`](https://rdrr.io/pkg/psych/man/fa.parallel.html)
  call with PA-PC.

- **MAP** (`"map"`) – Velicer (1976) Minimum Average Partial criterion.
  Finds the k that minimises the average squared partial correlation
  remaining after extracting k components. Usually conservative. Shares
  one [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html) call with
  VSS.

- **VSS** (`"vss"`) – Revelle & Rocklin (1979) Very Simple Structure fit
  at complexities 1 and 2 (VSS-1 and VSS-2). Finds the k maximising the
  fit of a very simple loading structure. Shares one
  [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html) call with
  MAP.

- **CD** (`"cd"`) – Ruscio & Roche (2012) Comparison Data. Resamples
  from the observed item distributions to generate comparison eigenvalue
  profiles; retains factors until adding one no longer improves RMSE
  beyond chance. Requires the EFAtools package (install separately).
  Skipped with an informational message when EFAtools is absent or when
  a correlation matrix is supplied.

PA (both bases), MAP, and VSS share the same correlation matrix as
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
CD operates on the raw data matrix directly (required for resampling).

## Interpreting the output

`k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
is a **maximum depth**, not a claim that exactly k factors exist. Users
commonly set k one or two levels above the consensus to watch
higher-level factors fragment – this is a feature of the method, not
overextraction.

Treating retention estimates as a range rather than a verdict has direct
support in the parallel-analysis literature: Lim and Jahng (2019)
recommend interpreting the PA estimate as a range of roughly plus or
minus one factor, resolved by interpretability, and Achim (2021) argues
that even that overstates PA's precision – the disagreement itself is
why `suggest_k()` reports several criteria and a consensus range, never
a single number.

## Ordinal (Likert) data

`suggest_k()` screens on the Pearson or Spearman basis by design and
never computes polychoric correlations itself, so when the raw data look
ordinal (at most 7 distinct integer values per column) it emits a
one-per-session warning – the same `detect_ordinal()` signal
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
and
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
use (Invariant 6: announce consequential defaults loudly). Crucially,
the advice points at the *final*
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
fit (`cor = "polychoric"`), **not** at `suggest_k()` itself: the
screening basis is intentional. The plausible-k *range* is robust to the
Pearson-vs-polychoric choice, so screening on Pearson and fitting the
final model polychoric is the recommended workflow – computing
polychoric correlations at the screening stage would only add cost and
NPD risk without changing the range. The warning is skipped for
correlation-matrix input (there are no items to inspect).

## A note on overextraction

PA-PC in particular tends to recommend more factors than replicate
across independent samples, especially with correlated items (Forbes,
2023). This is a long-standing observation in practice: Saucier (1997,
footnote 14) reported parallel analysis suggesting as many as 30 factors
in wide lexical item sets. PA-FA and CD are more conservative. Treat the
full set of criteria as a range: the true k is likely somewhere in the
middle.

## References

Achim, A. (2021). Determining the number of factors using parallel
analysis and its recent variants: Comment on Lim and Jahng (2019).
*Psychological Methods*, 26(1), 69–73.
[doi:10.1037/met0000269](https://doi.org/10.1037/met0000269)

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

Tong, L., Qu, W., & Zhang, Z. (2025). Comparison of the K1 rule,
parallel analysis, and the Bass-Ackward method on identifying the number
of factors in factor analysis. *Fudan Journal of the Humanities and
Social Sciences*, 18(1), 17–44.
[doi:10.1007/s40647-024-00423-2](https://doi.org/10.1007/s40647-024-00423-2)

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika*, 30, 179–185.

Lim, S., & Jahng, S. (2019). Determining the number of factors using
parallel analysis and its recent variants. *Psychological Methods*,
24(4), 452–467.
[doi:10.1037/met0000230](https://doi.org/10.1037/met0000230)

Revelle, W., & Rocklin, T. (1979). Very simple structure: An alternative
procedure for estimating the optimal number of interpretable factors.
*Multivariate Behavioral Research*, 14(4), 403–414.

Ruscio, J., & Roche, B. (2012). Determining the number of factors to
retain in an exploratory factor analysis using comparison data of a
known factorial structure. *Psychological Assessment*, 24(2), 282–292.

Saucier, G. (1997). Effects of variable selection on the factor
structure of person descriptors. *Journal of Personality and Social
Psychology*, 73(6), 1296–1312.
[doi:10.1037/0022-3514.73.6.1296](https://doi.org/10.1037/0022-3514.73.6.1296)

Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. *Psychometrika*, 41, 321–327.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md);
[`factorability()`](https://jmgirard.github.io/ackwards/reference/factorability.md),
[`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md),
and
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md),
the other pre-analysis diagnostics.

## Examples

``` r
# \donttest{
sk <- suggest_k(sim16)
#> ℹ Running parallel analysis (20 iterations, PC + FA)...
#> ✔ Running parallel analysis (20 iterations, PC + FA)... [182ms]
#> 
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [75ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [3.1s]
#> 
sk
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 16
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#>   k  PA-PC  PA-FA      MAP    VSS-1    VSS-2  CD
#>   1     ✔︎      ✔︎   0.0668   0.5505   0.0000   ✔︎ 
#>   2     ✔︎      ✔︎   0.0495   0.7439   0.7834   ✔︎ 
#>   3     ✔︎      ✔︎   0.0402   0.7685   0.8641   ✔︎ 
#>   4     ✔︎      ✔︎   0.0230*  0.7884*  0.9002*  ✔︎*
#>   5     -      -   0.0320   0.7881   0.8882   - 
#>   6     -      -   0.0425   0.7870   0.8633   - 
#>   7     -      -   0.0545   0.7867   0.8557   - 
#>   8     -      -   0.0714   0.7336   0.8372   - 
#>   ✔︎ retained   * optimal k   - not retained
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 4
#> • PA-FA: k <= 4
#> • MAP: k = 4
#> • VSS-1: k = 4
#> • VSS-2: k = 4
#> • CD: k = 4
#> Consensus: k = 4
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.
autoplot(sk)


# Run only MAP (fast; skips parallel analysis and CD)
suggest_k(sim16, criteria = "map")
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [40ms]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 16
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#>   k      MAP
#>   1  0.0668 
#>   2  0.0495 
#>   3  0.0402 
#>   4  0.0230*
#>   5  0.0320 
#>   6  0.0425 
#>   7  0.0545 
#>   8  0.0714 
#>   * optimal k
#> 
#> ── Recommendations ──
#> 
#> • MAP: k = 4
#> Consensus: k = 4
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.

# Run only the parallel-analysis criteria
suggest_k(sim16, criteria = c("pa_pc", "pa_fa"), n_iter = 5)
#> ℹ Running parallel analysis (5 iterations, PC + FA)...
#> ✔ Running parallel analysis (5 iterations, PC + FA)... [75ms]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 16
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#>   k  PA-PC  PA-FA
#>   1     ✔︎      ✔︎ 
#>   2     ✔︎      ✔︎ 
#>   3     ✔︎      ✔︎ 
#>   4     ✔︎      ✔︎ 
#>   5     -      - 
#>   6     -      - 
#>   7     -      - 
#>   8     -      - 
#>   ✔︎ retained   - not retained
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 4
#> • PA-FA: k <= 4
#> Consensus: k = 4
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.

# Faster exploratory run
suggest_k(sim16, k_max = 6, n_iter = 5)
#> ℹ Running parallel analysis (5 iterations, PC + FA)...
#> ✔ Running parallel analysis (5 iterations, PC + FA)... [80ms]
#> 
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [62ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [3.1s]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 16
#> n: 1,000
#> Basis: pearson
#> Tested k: 1-6
#> 
#> ── Criteria (k = 1-6) ──
#> 
#>   k  PA-PC  PA-FA      MAP    VSS-1    VSS-2  CD
#>   1     ✔︎      ✔︎   0.0668   0.5505   0.0000   ✔︎ 
#>   2     ✔︎      ✔︎   0.0495   0.7439   0.7834   ✔︎ 
#>   3     ✔︎      ✔︎   0.0402   0.7685   0.8641   ✔︎ 
#>   4     ✔︎      ✔︎   0.0230*  0.7884*  0.9002*  ✔︎*
#>   5     -      -   0.0320   0.7881   0.8882   - 
#>   6     -      -   0.0425   0.7870   0.8633   - 
#>   ✔︎ retained   * optimal k   - not retained
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 4
#> • PA-FA: k <= 4
#> • MAP: k = 4
#> • VSS-1: k = 4
#> • VSS-2: k = 4
#> • CD: k = 4
#> Consensus: k = 4
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.

# Correlation-matrix input (CD is skipped; n_obs required)
R <- cor(sim16)
suggest_k(R, n_obs = nrow(sim16))
#> ℹ CD is skipped when a correlation matrix is supplied (CD requires raw item
#>   distributions for resampling).
#> ℹ Running parallel analysis (20 iterations, PC + FA)...
#> ✔ Running parallel analysis (20 iterations, PC + FA)... [167ms]
#> 
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [77ms]
#> 
#> 
#> ── Factor / Component Count Suggestion (ackwards) ──────────────────────────────
#> Variables: 16
#> n: 1,000
#> Basis: (user-supplied matrix)
#> Tested k: 1-8
#> 
#> ── Criteria (k = 1-8) ──
#> 
#>   k  PA-PC  PA-FA      MAP    VSS-1    VSS-2
#>   1     ✔︎      ✔︎   0.0668   0.5505   0.0000 
#>   2     ✔︎      ✔︎   0.0495   0.7439   0.7834 
#>   3     ✔︎      ✔︎   0.0402   0.7685   0.8641 
#>   4     ✔︎      ✔︎   0.0230*  0.7884*  0.9002*
#>   5     -      -   0.0320   0.7881   0.8882 
#>   6     -      -   0.0425   0.7870   0.8633 
#>   7     -      -   0.0545   0.7867   0.8557 
#>   8     -      -   0.0714   0.7336   0.8372 
#>   ✔︎ retained   * optimal k   - not retained
#> + CD skipped (requires raw data; not available for matrix input).
#> 
#> ── Recommendations ──
#> 
#> • PA-PC: k <= 4
#> • PA-FA: k <= 4
#> • MAP: k = 4
#> • VSS-1: k = 4
#> • VSS-2: k = 4
#> Consensus: k = 4
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: k_max in ackwards() is a maximum depth. Setting k_max one or two levels
#> above the consensus to observe factor fragmentation is intentional.
#> Caution: PA-PC tends to overextract; structures may not replicate (Forbes,
#> 2023). PA-FA and CD are more conservative. Use the range.
# }
```
