# Tidy an ackwards object into a long data frame

Returns structured data from an `ackwards` object in tidy format. The
default (`what = "edges"`) returns the graph edge list that drives
diagrams.

## Usage

``` r
# S3 method for class 'ackwards'
tidy(
  x,
  what = c("edges", "loadings", "loadings_se", "variance", "fit", "nodes", "scores"),
  primary_only = FALSE,
  sort = c("none", "strength"),
  ...
)
```

## Arguments

- x:

  An `ackwards` object.

- what:

  What to extract:

  - `"edges"` *(default)* — one row per directed between-level edge:
    `from`, `to`, `level_from`, `level_to`, `r`, `is_primary`,
    `above_cut`.

  - `"loadings"` — one row per item × factor × level: `level`, `factor`,
    `item`, `loading`.

  - `"loadings_se"` — one row per item × factor × level with the
    rotation-aware standard error of each loading: `level`, `factor`,
    `item`, `se`. Only the ESEM engine (`engine = "esem"`) produces
    these; errors informatively for PCA/EFA objects, which carry no
    loading standard errors.

  - `"variance"` — one row per factor × level: `level`, `factor`,
    `variance_pct`, `cumulative_pct`.

  - `"fit"` — one row per fit index × level: `level`, `index`, `value`.
    For PCA objects the indices are eigenvalues; for EFA objects they
    are `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC`.

  - `"nodes"` — Forbes-extension pruning annotations (requires
    `prune != "none"` when the object was created). One row per factor
    across all levels: `id`, `level`, `pruned`, `prune_reason`. Returns
    an empty data frame with the same columns when no pruning was
    applied.

  - `"scores"` — long-format per-observation factor scores (requires
    `keep_scores = TRUE` at fit time or use
    [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
    for on-the-fly computation). Columns: `obs` (row index), `level`,
    `factor`, `score`.

- primary_only:

  For `what = "edges"` only. When `TRUE`, returns just each factor's
  primary-parent edge (`is_primary == TRUE`) — the lineage tree that the
  diagram draws as solid arrows. Default `FALSE` (all edges). Errors for
  any other value of `what`.

- sort:

  For `what = "edges"` only. One of `"none"` (default, natural order) or
  `"strength"` (descending `|r|`). Ignored for all other values of
  `what`.

- ...:

  Ignored.

## Value

A data frame (class `data.frame`).

## See also

[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)

## Examples

``` r
if (requireNamespace("psych", quietly = TRUE)) {
  x <- ackwards(psych::bfi[, 1:25], k_max = 5)
  tidy(x) # edges in natural order
  tidy(x, sort = "strength") # strongest edges first
  tidy(x, primary_only = TRUE) # just the primary-parent lineage
  tidy(x, what = "loadings")
  tidy(x, what = "variance")
}
#> Warning: ! 364 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#>    level factor variance_pct cumulative_pct
#> 1      1   m1f1        20.15          20.15
#> 2      2   m2f1        17.65          17.65
#> 3      2   m2f2        13.48          31.12
#> 4      3   m3f1        15.52          15.52
#> 5      3   m3f2        12.85          28.37
#> 6      3   m3f3        11.18          39.55
#> 7      4   m4f1        15.30          15.30
#> 8      4   m4f2        12.81          28.11
#> 9      4   m4f3        10.19          38.30
#> 10     4   m4f4         8.58          46.88
#> 11     5   m5f1        12.72          12.72
#> 12     5   m5f2        12.34          25.06
#> 13     5   m5f3        10.26          35.32
#> 14     5   m5f4         9.27          44.59
#> 15     5   m5f5         8.44          53.02
```
