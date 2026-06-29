# Summarise an ackwards object

Returns a structured `summary_ackwards` object that, when printed, shows
per-level variance and fit indices, a readable lineage list, and (when
present) pruning annotations. More verbose than
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md);
designed for inspection and reporting.

## Usage

``` r
# S3 method for class 'ackwards'
summary(object, ...)
```

## Arguments

- object:

  An `ackwards` object.

- ...:

  Ignored.

## Value

An object of class `"summary_ackwards"`, printed via
[`print.summary_ackwards()`](https://jmgirard.github.io/ackwards/reference/print.summary_ackwards.md).

## See also

[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)

## Examples

``` r
if (requireNamespace("psych", quietly = TRUE)) {
  x <- ackwards(psych::bfi[, 1:25], k_max = 5)
  summary(x)
}
#> Warning: ! 364 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> 
#> ── Summary: Bass-Ackwards Analysis (ackwards) ──────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 2,800
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> k = 1: 1 factor (20.15% cumulative variance)
#> m1f1 20.15% eigenvalue 5.04
#> k = 2: 2 factors (31.12% cumulative variance)
#> m2f1 17.65% eigenvalue 5.04
#> m2f2 13.48% eigenvalue 2.74
#> k = 3: 3 factors (39.55% cumulative variance)
#> m3f1 15.52% eigenvalue 5.04
#> m3f2 12.85% eigenvalue 2.74
#> m3f3 11.18% eigenvalue 2.11
#> k = 4: 4 factors (46.88% cumulative variance)
#> m4f1 15.3% eigenvalue 5.04
#> m4f2 12.81% eigenvalue 2.74
#> m4f3 10.19% eigenvalue 2.11
#> m4f4 8.58% eigenvalue 1.83
#> k = 5: 5 factors (53.02% cumulative variance)
#> m5f1 12.72% eigenvalue 5.04
#> m5f2 12.34% eigenvalue 2.74
#> m5f3 10.26% eigenvalue 2.11
#> m5f4 9.27% eigenvalue 1.83
#> m5f5 8.44% eigenvalue 1.54
#> 
#> ── Lineage (primary parents) ──
#> 
#> m1f1 → m2f1, m2f2
#> m2f1 → m3f1, m3f3
#> m2f2 → m3f2
#> m3f1 → m4f1
#> m3f2 → m4f2
#> m3f3 → m4f3, m4f4
#> m4f1 → m5f2, m5f4
#> m4f2 → m5f1
#> m4f3 → m5f3
#> m4f4 → m5f5
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```
