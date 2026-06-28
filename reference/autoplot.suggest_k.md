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
if (FALSE) { # \dontrun{
sk <- suggest_k(psych::bfi[, 1:25], n_iter = 5)
autoplot(sk)
} # }
```
