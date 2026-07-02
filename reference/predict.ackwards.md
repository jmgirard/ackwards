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
train <- sim16[1:500, ]
test <- sim16[501:1000, ]
x <- ackwards(train, k_max = 5)
test_scores <- predict(x, test)
head(test_scores)
#>        .m1f1       .m2f1      .m2f2       .m3f1       .m3f2      .m3f3
#> 1  0.9866825  1.10925975  0.2651224  0.22302480  0.83381091  0.7455003
#> 2 -1.3121549  0.41877448 -2.3452698 -2.31934440 -0.09306001  0.5863426
#> 3  1.0695420 -0.96685053  2.5697581  2.66347846 -0.91744168 -0.2973669
#> 4  0.4676516 -0.17389635  0.8617838  1.01558023 -0.99973200  0.8581230
#> 5  0.5951633  0.47949881  0.3594096  0.09180377  1.93864339 -1.3480924
#> 6  1.3390624  0.01471201  1.9282385  1.84818068  0.62412605 -0.5416109
#>         .m4f1      .m4f2       .m4f3       .m4f4      .m5f1      .m5f2
#> 1  0.70928336  0.8094370  0.45520440 -0.05611942  0.7250445  0.8034631
#> 2 -0.08628128  0.6373020 -1.40074580 -1.85689000 -0.1442605  0.6544344
#> 3 -0.78400610 -0.4193713  1.15211798  2.47154121 -0.7642160 -0.4206876
#> 4 -0.76118123  0.7194296 -0.07002827  1.34363823 -0.7490204  0.7191977
#> 5  1.29627226 -1.0086779  1.62903644 -1.17221196  1.2950204 -1.0113150
#> 6  0.25414070 -0.3714401  1.87301298  0.86141204  0.2659999 -0.3729043
#>        .m5f3       .m5f4      .m5f5
#> 1  0.4524271 -0.07636977  0.3129451
#> 2 -1.4086669 -1.75811718 -1.6661291
#> 3  1.1651161  2.42899090  0.7111726
#> 4 -0.0638014  1.32385148  0.3743907
#> 5  1.6229873 -1.17404983 -0.1157916
#> 6  1.8770026  0.83487823  0.3350342

# Identical to the augment() spelling:
identical(test_scores, augment(x, data = test, append = FALSE))
#> [1] TRUE
```
