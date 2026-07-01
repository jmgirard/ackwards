# Augment data with factor scores from an ackwards object

Appends per-observation factor scores for every level to a data frame.
Score columns are named `.m{k}f{j}` (e.g., `.m1f1`, `.m2f1`, `.m2f2`),
matching the factor labels used throughout the object.

## Usage

``` r
# S3 method for class 'ackwards'
augment(x, data = NULL, append = TRUE, id_cols = NULL, ...)
```

## Arguments

- x:

  An `ackwards` object.

- data:

  A data frame or numeric matrix with the same variables (columns) used
  to fit `x`. When `NULL` (default), uses pre-stored scores if available
  (requires `keep_scores = TRUE` at fit time).

- append:

  Logical. When `TRUE` (default) the score columns are appended to
  `data` (or to a `.obs` index when `data` is `NULL`). When `FALSE` only
  the score columns are returned (plus any `id_cols`).

- id_cols:

  Optional character vector naming columns of `data` to carry through
  alongside the scores when `append = FALSE` (for example a subject
  identifier, so scores can be rejoined after filtering). Ignored – and
  an error – when `append = TRUE` (all columns are already kept) or when
  `data` is `NULL` (there are no source columns to carry). `NULL`
  (default) returns the bare score columns.

- ...:

  Ignored.

## Value

A data frame. With `append = TRUE`: the supplied `data` (or a `.obs`
index when `data` is `NULL`) with score columns appended. With
`append = FALSE`: only the score columns, optionally prefixed by the
`id_cols`. Row order and count always match the input.

## Details

**Score computation.** Scores are `S = Z W / sqrt(score_var)`, where
`Z = .standardize(data)` (item z-scores), `W` is the per-level weight
matrix stored in the object, and `sqrt(score_var)` standardizes by the
real score standard deviations (Invariant 1: never assume unit
variance). For PCA the method is `"components"`; for EFA/ESEM it is
`"tenBerge"`.

**Missing data.** Score projection applies weights row-wise and
propagates NAs listwise: any observation with at least one missing item
variable will produce `NA` scores at every level. This differs from
fitting, which uses pairwise-complete correlations. A warning is issued
if NA rows are detected. Use `na.omit(data)` before scoring if NA rows
are unwanted.

**Data source.** If `data` is supplied, scores are always recomputed
from it using the stored weights – this is how to score new
observations. If `data` is `NULL` and `keep_scores = TRUE` was set at
fit time, the stored scores are returned. If neither is available an
informative error is raised.

**Scores-only output (`append = FALSE`).** By default (`append = TRUE`)
the scores are appended to the supplied `data`, following the broom
convention. Set `append = FALSE` to return *only* the score columns –
convenient for feeding scores straight into
[`cor()`](https://rdrr.io/r/stats/cor.html),
[`lm()`](https://rdrr.io/r/stats/lm.html), or a clustering call without
dragging the item columns along. Because
[`augment()`](https://generics.r-lib.org/reference/augment.html) always
preserves row order and row count,
`cbind(data, augment(x, data, append = FALSE))` reproduces the appended
output exactly. For a rejoin that survives *filtering* the scores
afterwards, name identifier columns with `id_cols` so they travel with
the scores.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)

## Examples

``` r
# Score the training data on the fly (no keep_scores = TRUE needed)
x <- ackwards(bfi25, k_max = 5)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
scores_df <- augment(x, data = bfi25)
#> Warning: ! 125 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
head(scores_df[, startsWith(names(scores_df), ".m")])
#>         .m1f1      .m2f1      .m2f2       .m3f1      .m3f2       .m3f3
#> 1 -1.30030695 -1.4235382 -0.1052068 -0.85089105  0.2798278 -1.51077483
#> 2  0.01634436  0.5630258 -0.9989376  0.56828756 -0.9942347 -0.08567816
#> 3  0.59542725  0.5404323  0.2526270  0.67694066  0.3150119 -0.10907369
#> 4 -0.71828163 -0.5073496 -0.5700585  0.06274103 -0.2041248 -1.39504136
#> 5 -1.32474768 -1.2750967 -0.4286547 -1.42804149 -0.4773126 -0.05571538
#> 6  0.52117719  0.1647933  0.7867243  0.29379255  0.8458643 -0.07021866
#>        .m4f1      .m4f2       .m4f3       .m4f4       .m5f1        .m5f2
#> 1 -0.7568233  0.2637881 -1.60295498 -0.39893361  0.05797902  0.225048096
#> 2  0.5586824 -1.0007502 -0.03472584  0.09088737  0.43531863 -1.000292190
#> 3  0.6569910  0.4184690 -0.48684202  0.66308866  1.31520127  0.356916365
#> 4  0.1122025 -0.1138340 -1.83242813  0.39875175 -0.70662078 -0.005187861
#> 5 -1.3991091 -0.5747296  0.27674922 -0.71562603 -1.89416329 -0.511204049
#> 6  0.3599258  0.6416350  0.64476695 -1.28366300 -0.16261206  0.664345274
#>         .m5f3      .m5f4       .m5f5
#> 1 -1.41275114 -1.3997979 -0.65487376
#> 2 -0.04555192  0.3557413  0.06686548
#> 3 -0.26674609 -0.6561587  0.31185077
#> 4 -2.12703248  1.2472183  0.85014591
#> 5  0.05521847  0.1890915 -0.33729526
#> 6  0.56488072  0.7357341 -1.16642860

# Scores-only: just the .m{k}f{j} columns, ready for cor()/lm()
scores_only <- augment(x, data = bfi25, append = FALSE)
#> Warning: ! 125 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
round(cor(scores_only[, c(".m5f1", ".m5f2")]), 2)
#>       .m5f1 .m5f2
#> .m5f1     1    NA
#> .m5f2    NA     1

# Carry an identifier through for a safe post-filter rejoin
df <- data.frame(id = seq_len(nrow(bfi25)), bfi25)
scored <- augment(x, data = df, append = FALSE, id_cols = "id")
#> Warning: ! 125 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
head(scored[, c("id", ".m5f1")])
#>   id       .m5f1
#> 1  1  0.05797902
#> 2  2  0.43531863
#> 3  3  1.31520127
#> 4  4 -0.70662078
#> 5  5 -1.89416329
#> 6  6 -0.16261206

# Store at fit time and augment without re-supplying data
x2 <- ackwards(bfi25, k_max = 5, keep_scores = TRUE)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> Warning: ! 125 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
scores_df2 <- augment(x2)
```
