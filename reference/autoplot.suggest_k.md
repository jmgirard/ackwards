# Plot a suggest_k diagnostic

Renders a ggplot2 diagnostic for a `suggest_k` object. The plot shows
one panel for each criterion group that was requested: Scree/PA (if
`"pa_pc"` or `"pa_fa"` were requested), MAP (if `"map"`), VSS (if
`"vss"`), and CD RMSE (if `"cd"` and EFAtools was available). When four
panels are shown the layout is a 2x2 grid; otherwise a single-column
layout is used. The recommended k for each criterion is marked with a
star-shaped point.

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
their respective random-data thresholds (for whichever PA bases were
requested). PA-PC compares the blue "Observed (PC)" line to the dashed
PA-PC threshold; PA-FA compares the teal "Observed (FA)" line to the
dotted PA-FA threshold.

The CD panel plots the mean RMSE between observed and comparison-data
eigenvalues at each k. CD uses a sequential one-sided Wilcoxon test
(Ruscio & Roche, 2012): a factor is retained while adding it
significantly reduces RMSE (default \\\alpha = 0.30\\); the starred k is
the last retained factor. The curve is shown only over the levels that
were actually computed; the starred k need not be the visible minimum of
the plotted curve.

Requires the ggplot2 package.

## See also

[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)

## Examples

``` r
# \donttest{
if (requireNamespace("ggplot2", quietly = TRUE)) {
  sk <- suggest_k(sim16, n_iter = 5)
  autoplot(sk)
}
#> ℹ Running parallel analysis (5 iterations, PC + FA)...
#> ✔ Running parallel analysis (5 iterations, PC + FA)... [97ms]
#> 
#> ℹ Running MAP and VSS...
#> ✔ Running MAP and VSS... [160ms]
#> 
#> ℹ Running Comparison Data (CD)...
#> ✔ Running Comparison Data (CD)... [4.1s]
#> 

# }
```
