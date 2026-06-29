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
x <- ackwards(bfi25, k_max = 5)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
summary(x)
#> 
#> ── Summary: Bass-Ackwards Analysis (ackwards) ──────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> k = 1: 1 factor (20.61% cumulative variance)
#> m1f1 20.61% eigenvalue 5.15
#> k = 2: 2 factors (31.82% cumulative variance)
#> m2f1 18.46% eigenvalue 5.15
#> m2f2 13.36% eigenvalue 2.8
#> k = 3: 3 factors (40.09% cumulative variance)
#> m3f1 16.39% eigenvalue 5.15
#> m3f2 12.42% eigenvalue 2.8
#> m3f3 11.28% eigenvalue 2.07
#> k = 4: 4 factors (47.29% cumulative variance)
#> m4f1 15.77% eigenvalue 5.15
#> m4f2 12.58% eigenvalue 2.8
#> m4f3 10.67% eigenvalue 2.07
#> m4f4 8.28% eigenvalue 1.8
#> k = 5: 5 factors (53.28% cumulative variance)
#> m5f1 12.57% eigenvalue 5.15
#> m5f2 12.51% eigenvalue 2.8
#> m5f3 10.75% eigenvalue 2.07
#> m5f4 9.47% eigenvalue 1.8
#> m5f5 7.98% eigenvalue 1.5
#> 
#> ── Lineage (primary parents) ──
#> 
#> m1f1 → m2f1, m2f2
#> m2f1 → m3f1, m3f3
#> m2f2 → m3f2
#> m3f1 → m4f1
#> m3f2 → m4f2
#> m3f3 → m4f3, m4f4
#> m4f1 → m5f1, m5f4
#> m4f2 → m5f2
#> m4f3 → m5f3
#> m4f4 → m5f5
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```
