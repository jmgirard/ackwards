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
x <- ackwards(sim16, k_max = 5)
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
#> k = 1: 1 factor (28.2% cumulative variance)
#> m1f1 28.2% eigenvalue 4.51
#> 
#> k = 2: 2 factors (46.5% cumulative variance)
#> m2f1 23.3% eigenvalue 4.51
#> m2f2 23.2% eigenvalue 2.93
#> 
#> k = 3: 3 factors (57.5% cumulative variance)
#> m3f1 23.0% eigenvalue 4.51
#> m3f2 17.5% eigenvalue 2.93
#> m3f3 16.9% eigenvalue 1.76
#> 
#> k = 4: 4 factors (67.7% cumulative variance)
#> m4f1 17.2% eigenvalue 4.51
#> m4f2 16.9% eigenvalue 2.93
#> m4f3 16.8% eigenvalue 1.76
#> m4f4 16.8% eigenvalue 1.63
#> 
#> k = 5: 5 factors (70.8% cumulative variance)
#> m5f1 17.1% eigenvalue 4.51
#> m5f2 16.9% eigenvalue 2.93
#> m5f3 16.8% eigenvalue 1.76
#> m5f4 16.7% eigenvalue 1.63
#> m5f5 3.3% eigenvalue 0.51
#> 
#> ── Lineage (primary parents) ──
#> 
#> m1f1 → m2f1, m2f2
#> m2f1 → m3f2, m3f3
#> m2f2 → m3f1
#> m3f1 → m4f3, m4f4
#> m3f2 → m4f1
#> m3f3 → m4f2
#> m4f1 → m5f1
#> m4f2 → m5f2
#> m4f3 → m5f3
#> m4f4 → m5f4, m5f5
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
```
