# Suggest a maximum number of factors for bass-ackwards analysis

Runs two complementary criteria and reports their recommendations.
Neither alone is definitive; the goal is a consensus range to inform
your choice of `k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

## Usage

``` r
suggest_k(data, k_max = NULL, cor = "pearson", n_iter = 20L, ...)
```

## Arguments

- data:

  A data frame or numeric matrix (items in columns, observations in
  rows).

- k_max:

  Maximum number of components to test. Defaults to
  `min(ncol(data) - 1, 8)`. Increase if you expect a deeper hierarchy.

- cor:

  Correlation basis: `"pearson"` (default) or `"spearman"`. Should match
  the `cor` argument you plan to use in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

- n_iter:

  Number of Monte Carlo iterations for parallel analysis. Default `20`.
  Reduce to `5` for fast/exploratory runs; increase to `100+` for
  publication.

- ...:

  Reserved for future arguments.

## Value

An object of class `"suggest_k"`. Print it for a formatted summary. The
list contains:

- k_parallel:

  Recommended k from parallel analysis.

- k_map:

  Recommended k from MAP.

- criteria:

  Data frame with one row per k: `k`, `map` value, `pa_suggested`
  (logical; `TRUE` if k is within the parallel-analysis threshold).

- k_max, n_obs, n_vars, cor:

  Metadata.

## Details

**Criteria computed:**

- **Parallel analysis** (Horn 1965) – compares observed eigenvalues to
  those from random correlation matrices of the same size. Suggests the
  number of components whose eigenvalues exceed the 95th-percentile
  chance level.

- **MAP** (Velicer 1976) – the Minimum Average Partial criterion. Finds
  the number of components that minimises the average squared partial
  correlation after the components are partialled out.

Both use the same Pearson or Spearman correlation matrix as
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
Parallel analysis is implemented via
[`psych::fa.parallel()`](https://rdrr.io/pkg/psych/man/fa.parallel.html);
MAP via [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html).

## Interpreting the output

`k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
is a **maximum depth**, not a claim that exactly k factors exist. Users
commonly set k one or two levels above the consensus to watch
higher-level factors fragment – this is a feature of the method, not
overextraction.

## A note on overextraction

Parallel analysis in particular tends to recommend more factors than
replicate across independent samples, especially with correlated items
(Forbes, 2023). Treat these criteria as a starting range for
exploration, not a definitive stopping rule. Setting `k` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
one or two levels above the consensus is intentional – watching factors
fragment is part of the method.

## References

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000578](https://doi.org/10.1037/met0000578)

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika*, 30, 179–185.

Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. *Psychometrika*, 41, 321–327.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)

## Examples

``` r
if (FALSE) { # \dontrun{
suggest_k(psych::bfi[, 1:25])
suggest_k(psych::bfi[, 1:25], k_max = 6, n_iter = 5)
} # }
```
