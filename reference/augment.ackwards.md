# Augment data with factor scores from an ackwards object

Appends per-observation factor scores for every level to a data frame.
Score columns are named `.m{k}f{j}` (e.g., `.m1f1`, `.m2f1`, `.m2f2`),
matching the factor labels used throughout the object.

## Usage

``` r
# S3 method for class 'ackwards'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  An `ackwards` object.

- data:

  A data frame or numeric matrix with the same variables (columns) used
  to fit `x`. When `NULL` (default), uses pre-stored scores if available
  (requires `keep_scores = TRUE` at fit time).

- ...:

  Ignored.

## Value

A data frame. If `data` is supplied it is returned with score columns
appended. If `data` is `NULL` the return is a minimal data frame with a
`.obs` index column followed by score columns.

## Details

**Score computation.** Scores are `S = Z W / sqrt(score_var)`, where
`Z = scale(data)` (item z-scores), `W` is the per-level weight matrix
stored in the object, and `sqrt(score_var)` standardizes by the real
score standard deviations (Invariant 1: never assume unit variance). For
PCA the method is `"components"`; for EFA/ESEM it is `"tenBerge"`.

**Missing data.** Score projection applies weights row-wise and
propagates NAs listwise: any observation with at least one missing item
variable will produce `NA` scores at every level. This differs from
fitting, which uses pairwise-complete correlations. A warning is issued
if NA rows are detected. Use `na.omit(data)` before scoring if NA rows
are unwanted.

**Data source.** If `data` is supplied, scores are always recomputed
from it using the stored weights — this is how to score new
observations. If `data` is `NULL` and `keep_scores = TRUE` was set at
fit time, the stored scores are returned. If neither is available an
informative error is raised.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Score the training data on the fly (no keep_scores=TRUE needed)
x <- ackwards(psych::bfi[, 1:25], k_max = 5)
scores_df <- augment(x, data = psych::bfi[, 1:25])
head(scores_df[, grep("^\\.m", names(scores_df))])

# Or store at fit time and augment without re-supplying data
x2 <- ackwards(psych::bfi[, 1:25], k_max = 5, keep_scores = TRUE)
scores_df2 <- augment(x2)
} # }
```
