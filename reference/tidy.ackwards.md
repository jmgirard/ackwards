# Tidy an ackwards object into a long data frame

Returns structured data from an `ackwards` object in tidy format. The
default (`what = "edges"`) returns the graph edge list that drives
diagrams.

## Usage

``` r
# S3 method for class 'ackwards'
tidy(
  x,
  what = c("edges", "loadings", "variance", "fit", "nodes", "scores"),
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

- ...:

  Ignored.

## Value

A data frame (class `data.frame`).

## See also

[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
