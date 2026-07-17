# Plot a comparability diagnostic

Renders a two-panel ggplot2 diagnostic for a
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
object: score comparability (r) and loading congruence (Tucker's phi)
for every factor at every level. Grey points are individual splits;
black points are the per-factor medians. Dashed and dotted reference
lines mark the conventional .90 / .95 benchmarks – visual guides, not
tests.

## Usage

``` r
# S3 method for class 'comparability'
autoplot(object, ...)
```

## Arguments

- object:

  A `comparability` object.

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

Requires the ggplot2 package.

## See also

[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cmp <- comparability(sim16, k_max = 5, n_splits = 5, seed = 1)
  autoplot(cmp)
}
#> ℹ Fitting 5 split-half replicates (pca, k = 1-5)...
#> ✔ Fitting 5 split-half replicates (pca, k = 1-5)... [276ms]
#> 

# }
```
