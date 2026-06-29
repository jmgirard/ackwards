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
  align = TRUE,
  use = "pairwise.complete.obs",
  cut_show = 0.3
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
  (force; errors if conditions not met), or `"scores"` (always
  materialise).

- pairs:

  `"adjacent"` (classic Goldberg) or `"all"` (Forbes extension).

- data:

  Optional data frame / matrix of raw observations. Required only when
  `edge_method = "scores"` or the scores path is triggered.

- align:

  Whether to sign-align edges to primary-parent lineage.

- use:

  Passed to [`stats::cor()`](https://rdrr.io/r/stats/cor.html) when
  materialising scores.

- cut_show:

  Edges with `|r| >= cut_show` are flagged `above_cut` in the tidy
  tibble.

## Value

A list with:

- matrices:

  Named list of `(k_a x k_b)` edge matrices, keyed `"k_a:k_b"`.

- tidy:

  A data frame with one row per directed edge: `from`, `to`,
  `level_from`, `level_to`, `r`, `is_primary`, `above_cut`.

## Details

    E(a,b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}

where `R` is the input correlation matrix and `D_x = diag(W_x' R W_x)`
are the **actual** score variances (not assumed to be 1). This avoids
materialising scores while remaining exact for PCA, EFA (regression /
Bartlett / tenBerge) – all of which produce linear score maps.

When the algebra cannot be used (nonlinear scoring, missing `R`, or the
user forces `edge_method = "scores"`), scores are materialised from
`data` instead.
