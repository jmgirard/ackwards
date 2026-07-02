# Augment data with factor scores from an ackwards object

Appends per-observation factor scores for every level to a data frame.
Score columns are named `.m{k}f{j}` (e.g., `.m1f1`, `.m2f1`, `.m2f2`),
matching the factor labels used throughout the object.

## Usage

``` r
# S3 method for class 'ackwards'
augment(
  x,
  data = NULL,
  append = TRUE,
  id_cols = NULL,
  scaling = c("fit", "sample"),
  ...
)
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

- scaling:

  Which item means/SDs standardize `data` before the weights are
  applied. `"fit"` (default) uses the **fit-time** moments stored in the
  object – the correct choice for scoring new observations (e.g. a
  cross-validation test split) or subsets on the training metric.
  `"sample"` standardizes by the supplied data's own moments (the only
  option for objects fit from a correlation matrix, which carry no
  raw-data moments). Only used when `data` is supplied; passing it
  without `data` is an error (stored scores are returned exactly as
  computed at fit time).

- ...:

  Ignored.

## Value

A data frame. With `append = TRUE`: the supplied `data` (or a `.obs`
index when `data` is `NULL`) with score columns appended. With
`append = FALSE`: only the score columns, optionally prefixed by the
`id_cols`. Row order and count always match the input.

## Details

**Score computation.** Scores are `S = Z W / sqrt(score_var)`, where `Z`
is the item z-scores, `W` is the per-level weight matrix stored in the
object, and `sqrt(score_var)` standardizes by the real score standard
deviations (Invariant 1: never assume unit variance). For PCA the method
is `"components"`; for EFA/ESEM it is `"tenBerge"`. The `scaling`
argument controls which means/SDs build `Z`: by default the **fit-time**
moments stored in the object, so any data you score — the training data,
a subset of it, or entirely new observations — lands on the same metric
the model was estimated in.

**Scoring new observations (cross-validation).** Because scoring only
needs the stored weight matrices and the fit-time moments, you can fit
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
on a training split and score a held-out test split *without
retraining*: `augment(x, data = test_set)` (or, equivalently,
[`predict.ackwards()`](https://jmgirard.github.io/ackwards/reference/predict.ackwards.md)).
Under the default `scaling = "fit"` the test observations are
standardized by the *training* means/SDs, which is what "applying the
trained model" means: a test observation's score does not depend on
which other observations happen to share its split, and train and test
scores are directly comparable. `scaling = "sample"` instead
re-standardizes by the supplied data's own moments — a deliberate choice
when scoring a sample from a different population in its own metric, and
the only option for objects fit from a correlation matrix (which carry
no raw-data moments). For non-Pearson bases (polychoric, Spearman) the
usual caveat applies either way: the weights derive from the non-Pearson
`R` while `Z` is a linear standardization, so empirical score SDs are
close to but not exactly 1 (a one-time warning says so); train/test
comparability under `scaling = "fit"` is unaffected. For objects fit
with `missing = "fiml"` (PCA/EFA), the stored moments are the observed
(`na.rm`) means/SDs of the training items — the correlation matrix was
FIML-estimated, but scoring operates on observed responses, so the
observed moments are the consistent frame; incomplete rows still score
`NA` (scoring does not impute).

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
output exactly. A row that scores `NA` (because it is missing an item –
see *Missing data* above) is still returned in place, so the positional
alignment holds even with `NA` scores. For a rejoin that survives
*filtering* the scores afterwards, name identifier columns with
`id_cols` so they travel with the scores.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)

## Examples

``` r
# Score the training data on the fly (no keep_scores = TRUE needed)
x <- ackwards(sim16, k_max = 5)
scores_df <- augment(x, data = sim16)
head(scores_df[, startsWith(names(scores_df), ".m")])
#>        .m1f1       .m2f1      .m2f2      .m3f1      .m3f2       .m3f3
#> 1  0.7084102  2.08220676 -1.0812105 -1.0846753  1.2529408  1.66899528
#> 2 -0.5732057 -0.64530574 -0.1652010 -0.1361687 -0.6229163 -0.29658576
#> 3 -1.5222509  0.21929430 -2.3727764 -2.2976644 -0.5385154  0.79318695
#> 4  0.2786514  0.02349164  0.3706740  0.3590532  0.1184985 -0.07677641
#> 5 -1.0170680 -0.97317875 -0.4650366 -0.3494137 -1.4985488  0.09413691
#> 6 -0.2675764 -0.95846783  0.5804698  0.6811951 -1.3327340 -0.03027230
#>        .m4f1       .m4f2       .m4f3       .m4f4      .m5f1       .m5f2
#> 1  0.8367367  1.98120705  0.22462994 -1.61915885  0.8342794  1.98244560
#> 2 -0.4644816 -0.40243702 -0.47549374  0.22730576 -0.4843440 -0.40064269
#> 3  0.2058961  0.31402515 -3.08727273 -0.32514479  0.1577276  0.31322943
#> 4  0.2256783 -0.16054771  0.03139285  0.45632266  0.2093034 -0.15629258
#> 5 -1.4268829  0.06527855 -0.58548084  0.02156911 -1.3840396  0.05153121
#> 6 -1.0361522 -0.23837148 -0.35914294  1.20575073 -1.0246840 -0.24329695
#>        .m5f3      .m5f4      .m5f5
#> 1  0.2196692 -1.6038056 -0.2456481
#> 2 -0.4732311  0.3005452 -0.8525288
#> 3 -3.0886969 -0.1871828 -1.6902637
#> 4  0.0294733  0.5224934 -0.7365249
#> 5 -0.5787203 -0.1735409  2.2069843
#> 6 -0.3566370  1.1374441  0.8127603

# Scores-only: just the .m{k}f{j} columns, ready for cor()/lm()
scores_only <- augment(x, data = sim16, append = FALSE)
round(cor(scores_only[, c(".m5f1", ".m5f2")]), 2)
#>       .m5f1 .m5f2
#> .m5f1     1     0
#> .m5f2     0     1

# Carry an identifier through for a safe post-filter rejoin
df <- data.frame(id = seq_len(nrow(sim16)), sim16)
scored <- augment(x, data = df, append = FALSE, id_cols = "id")
head(scored[, c("id", ".m5f1")])
#>   id      .m5f1
#> 1  1  0.8342794
#> 2  2 -0.4843440
#> 3  3  0.1577276
#> 4  4  0.2093034
#> 5  5 -1.3840396
#> 6  6 -1.0246840

# Store at fit time and augment without re-supplying data
x2 <- ackwards(sim16, k_max = 5, keep_scores = TRUE)
scores_df2 <- augment(x2)

# Cross-validation: fit on a training split, score the test split on the
# training metric (no retraining; see also predict.ackwards())
train_idx <- seq_len(500)
x_train <- ackwards(sim16[train_idx, ], k_max = 5)
test_scores <- augment(x_train, data = sim16[-train_idx, ], append = FALSE)
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
```
