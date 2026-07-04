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
  cmp <- comparability(bfi25, k_max = 5, n_splits = 5, seed = 1)
  autoplot(cmp)
}
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> Warning: ! 25 columns look like ordinal/Likert items (<= 7 distinct integer values):
#>   "A1", "A2", "A3", "A4", "A5", "C1", …, "O4", and "O5".
#> ℹ Results use a "pearson" basis. Consider `cor = "polychoric"` for ordinal
#>   data.
#> This warning is displayed once per session.
#> ℹ Fitting 5 split-half replicates (pca, k = 1-5)...
#> ✔ Fitting 5 split-half replicates (pca, k = 1-5)... [774ms]
#> 

# }
```
