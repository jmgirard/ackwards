# Score new observations with a fitted ackwards model

Applies a fitted bass-ackwards model to new data — for example the
held-out test split of a cross-validation design — producing factor
scores for every level **without retraining**. This is the idiomatic
[`predict()`](https://rdrr.io/r/stats/predict.html) front door to the
same machinery as
[`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md):
the call `predict(object, newdata)` returns exactly
`augment(object, data = newdata, append = FALSE)`, a data frame holding
only the `.m{k}f{j}` score columns, one row per row of `newdata`.

## Usage

``` r
# S3 method for class 'ackwards'
predict(object, newdata, scaling = c("fit", "sample"), ...)
```

## Arguments

- object:

  An `ackwards` object.

- newdata:

  A data frame or numeric matrix with the same variables (columns) used
  to fit `object`. Required — to retrieve scores stored at fit time, use
  `augment(object)` instead.

- scaling:

  Which item means/SDs standardize `newdata`: `"fit"` (default, the
  training moments) or `"sample"` (`newdata`'s own moments). See
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md).

- ...:

  Ignored.

## Value

A data frame of factor scores (columns `.m{k}f{j}`, one per factor per
level), with one row per row of `newdata` in the original order.

## Details

Under the default `scaling = "fit"`, `newdata` is standardized by the
**fit-time** item means/SDs stored in the object before the stored
weight matrices are applied, so the new scores land on the same metric
as the training solution: an observation's score does not depend on
which other observations share its split, and train and test scores are
directly comparable. See
[`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
(section *Scoring new observations*) for the full semantics, the
`scaling = "sample"` alternative, and the non-Pearson-basis caveat.

`newdata` must contain the variables the model was fit on (matched by
column name, with extra columns ignored; a bare unnamed matrix is
matched positionally). Rows with missing items produce `NA` scores
(scoring does not impute).

## See also

[`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
for appending scores to the data (and the full scoring documentation),
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

## Examples

``` r
# Cross-validation: fit on a training split, score the test split
train <- bfi25[1:500, ]
test <- bfi25[501:1000, ]
x <- ackwards(train, k_max = 5)
#> Warning: ! 58 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
test_scores <- predict(x, test)
#> Warning: ! 67 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
head(test_scores)
#>         .m1f1      .m2f1      .m2f2      .m3f1      .m3f2     .m3f3      .m4f1
#> 1  1.31327917  1.2430095  0.4608252  0.6680509  0.1248161 1.4375648  0.7088025
#> 2  0.16529258  1.2321417 -1.9494198  1.0643443 -2.0400056 0.2062707  1.2228552
#> 3          NA         NA         NA         NA         NA        NA         NA
#> 4  0.66782843  0.3087139  0.8378206  0.3044256  0.8063030 0.2575057  0.4353962
#> 5 -0.02123811 -0.6503492  1.1686927 -0.8540634  1.0730075 0.4449248 -0.9429296
#> 6          NA         NA         NA         NA         NA        NA         NA
#>        .m4f2     .m4f3      .m4f4      .m5f1      .m5f2     .m5f3       .m5f4
#> 1  0.1179912 1.6867161 -0.3149617  0.1303632  1.6548929 1.5959910 -0.69755418
#> 2 -2.1590779 0.7967761 -1.0174953 -2.1501813  1.2787053 0.7224308  0.39711553
#> 3         NA        NA         NA         NA         NA        NA          NA
#> 4  0.7197966 0.5967368 -0.9248642  0.7224354  0.5193448 0.5700356  0.03266355
#> 5  1.1502964 0.1572142  0.5246276  1.1329384 -1.3643529 0.2357785  0.11707246
#> 6         NA        NA         NA         NA         NA        NA          NA
#>        .m5f5
#> 1 -0.8091157
#> 2 -1.2663365
#> 3         NA
#> 4 -1.0712061
#> 5  0.8709454
#> 6         NA

# Identical to the augment() spelling:
identical(test_scores, augment(x, data = test, append = FALSE))
#> Warning: ! 67 rows contain missing values and will produce NA scores.
#> ℹ Score projection applies weights row-wise and propagates NAs listwise; FIML
#>   estimation does not impute missing item responses.
#> ℹ Use `missing = "listwise"` when fitting (so the model and scores share the
#>   same complete rows), or call `na.omit(data)` before scoring.
#> [1] TRUE
```
