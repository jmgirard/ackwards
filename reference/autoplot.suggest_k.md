# Plot a suggest_k diagnostic

Renders a three-panel ggplot2 diagnostic for a `suggest_k` object: a
parallel-analysis/scree plot on top, a MAP panel in the middle, and a
VSS panel on the bottom. The recommended k for each criterion is marked
with a star-shaped point.

## Usage

``` r
# S3 method for class 'suggest_k'
autoplot(object, ...)
```

## Arguments

- object:

  A `suggest_k` object.

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

The scree panel shows both PC and FA observed eigenvalues alongside
their respective random-data thresholds. PA-PC compares the blue
"Observed (PC)" line to the dashed PA-PC threshold; PA-FA compares the
teal "Observed (FA)" line to the dotted PA-FA threshold. Reading the
lines on the same panel is informative but the two comparisons are
independent.

Requires the ggplot2 package.

## See also

[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)

## Examples

``` r
# \donttest{
if (requireNamespace("psych", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)) {
  sk <- suggest_k(psych::bfi[, 1:25], n_iter = 5)
  autoplot(sk)
}
#> ℹ Running parallel analysis (5 iterations, PC + FA)...
#> ✔ Running parallel analysis (5 iterations, PC + FA)... [100ms]
#> 
#> ℹ Running MAP and VSS...
#> CD: 364 rows with missing values removed (2436 complete cases used).
#> ✔ Running MAP and VSS... [120ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [23.8s]
#> 

# }
```
