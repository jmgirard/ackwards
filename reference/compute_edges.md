# Compute between-level factor-score correlations

The centerpiece of the bass-ackwards algebra. For any engine whose
scoring is a **linear** map `S = Z W`, the cross-level correlation
matrix is:

## Usage

``` r
compute_edges(
  levels,
  R,
  edge_method = c("auto", "algebra", "scores"),
  pairs = c("adjacent", "all"),
  data = NULL,
  use = "pairwise.complete.obs",
  cut_show = 0.3,
  build_tidy = TRUE
)
```

## Arguments

- levels:

  Named list (indexed by k) of per-level objects produced by an engine.
  Each must contain a `scoring` sub-list with fields `linear`,
  `weights`, and `score_var`.

- R:

  Square correlation matrix (p x p). Required for the algebra path.

- edge_method:

  One of `"auto"` (algebra when possible, scores otherwise), `"algebra"`
  (force, and error if conditions are not met), or `"scores"` (always
  materialise).

- pairs:

  `"adjacent"` (classic Goldberg) or `"all"` (Forbes extension).

- data:

  Optional data frame / matrix of raw observations. Required only when
  `edge_method = "scores"` or the scores path is triggered.

- use:

  Passed to [`stats::cor()`](https://rdrr.io/r/stats/cor.html) when
  materialising scores.

- cut_show:

  Edges with `|r| >= cut_show` are flagged `above_cut` in the tidy
  tibble.

- build_tidy:

  Build the tidy edge data frame? `FALSE` returns `tidy = NULL` for
  matrices-only callers (lineage pass, `.cross_cor()`,
  `.boot_replicate()`), which would otherwise build and discard it
  (M60).

## Value

A list with:

- matrices:

  Named list of `(k_a x k_b)` edge matrices, keyed `"k_a:k_b"`.

- tidy:

  A data frame with one row per directed edge: `from`, `to`,
  `level_from`, `level_to`, `r`, `is_primary`, and `above_cut`. It is
  `NULL` when `build_tidy = FALSE`.

## Details

    E(a,b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}

where `R` is the correlation matrix passed in and
`D_x = diag(W_x' R W_x)` are the **actual** score variances (not assumed
to be 1). The algebra needs no scores. It is exact for that `R` under
the linear scoring of every engine: component score weights for PCA, and
ten Berge weights (regression weights as the fallback) for EFA and ESEM.
For the edges a fit reports, `R` is the matrix that
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
stores as `x$r`, so those edges are exact for `x$r`. For ESEM, lavaan
fits a different matrix in two settings. With `cor = "spearman"` and
`missing = "pairwise"` or `"listwise"`, `x$r` is a Spearman matrix, but
lavaan fits Pearson covariances of the raw data. With `cor = "pearson"`
and `missing = "pairwise"` on data with missing values, `x$r` holds
pairwise correlations. There, ML, MLR, and ULSMV fit the complete rows,
and WLSMV fits the pairwise covariances. The correlations of those
covariances differ slightly from `x$r`, because each item's variance
uses all of its observed rows. The split-half comparability check passes
a pooled `R` in place of `x$r`.

A pair goes to the scores branch when `edge_method = "scores"`, or when
`edge_method = "auto"` and either level's `scoring$linear` is not `TRUE`
or `R` is `NULL`. That branch correlates scores computed from `data`,
and it errors when `data` is `NULL`. Under `edge_method = "algebra"` the
same conditions raise an error instead. No shipped caller reaches the
scores branch. Each one passes `"auto"` or `"algebra"` with an `R`, and
every engine's scoring is linear. Only tests use it: the
algebra-vs-scores agreement tests and one error-path test.
