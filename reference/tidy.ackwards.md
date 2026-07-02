# Tidy an ackwards object into a long data frame

Returns structured data from an `ackwards` object in tidy format. The
default (`what = "edges"`) returns the graph edge list that drives
diagrams.

## Usage

``` r
# S3 method for class 'ackwards'
tidy(
  x,
  what = c("edges", "loadings", "variance", "fit", "nodes", "scores"),
  primary_only = FALSE,
  sort = c("none", "strength"),
  format = c("long", "wide"),
  conf_level = 0.95,
  ...
)
```

## Arguments

- x:

  An `ackwards` object.

- what:

  What to extract:

  - `"edges"` *(default)* – one row per directed between-level edge:
    `from`, `to`, `level_from`, `level_to`, `r`, `is_primary`,
    `above_cut`. If
    [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)
    has been run on the object, four bootstrap columns are appended:
    `se`, `lo`, `hi` (bootstrap standard error and percentile
    confidence-interval endpoints), and `n_boot_ok` (usable replicates).

  - `"loadings"` – one row per item x factor x level: `level`, `factor`,
    `item`, `loading`, `se`, `ci_lower`, `ci_upper`. `se`, `ci_lower`,
    and `ci_upper` are populated only for `engine = "esem"` (which
    produces rotation-aware loading SEs); they are `NA` for PCA and EFA.
    The confidence level is controlled by `conf_level`.

  - `"variance"` – one row per factor x level: `level`, `factor`,
    `proportion`, `cumulative`. Both are proportions of total item
    variance on a 0-1 scale (multiply by 100 for a percentage).

  - `"fit"` – one row per fit statistic x level: `level`, `statistic`,
    `value`. For PCA objects the statistics are eigenvalues; for EFA
    objects they are `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC` –
    where `chi` is the likelihood-ratio chi-square
    ([`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)'s
    `STATISTIC`), so `chi`, `dof`, `p_value`, `RMSEA`, and `TLI` all
    share one statistical framing (psych's residual-based *empirical*
    chi-square is a different statistic and is not reported); for ESEM
    they are `chi`, `dof`, `p_value`, `CFI`, `TLI`, `RMSEA`, `SRMR`,
    `BIC`. For ESEM under a scaled-test estimator (`"WLSMV"`/`"ULSMV"`
    for ordinal data, `"MLR"` for continuous), the whole row –
    `chi`/`dof`/ `p_value` **and** `CFI`/`TLI`/`RMSEA` – reports
    lavaan's mean-and-variance-adjusted ("scaled") variant, so every
    quantity shares one scaling. This matters most for WLSMV/ULSMV: the
    naive chi-square has no valid reference distribution (lavaan's own
    [`summary()`](https://rdrr.io/r/base/summary.html) labels its
    p-value "Unknown"), and the naive `CFI`/`TLI` are badly optimistic
    for ordinal data (Xia & Yang, 2019). `"ML"` has no scaled variant,
    so it reports the naive values (the correct ones for ML). `SRMR` has
    no scaled variant and is reported as-is. `BIC` is `NA` under
    WLSMV/ULSMV (no proper log-likelihood for a limited-information
    estimator) and populated under ML/MLR. Use `format = "wide"` for one
    row per **non-anchor** level (k \>= 2; the saturated 1-factor anchor
    is dropped, matching
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    `autoplot(what = "fit")`), one column per statistic. Conventional
    fit cutoffs (Hu & Bentler 1999) are shown as reference lines in
    `autoplot(what = "fit")` and inline in
    [`summary()`](https://rdrr.io/r/base/summary.html), but are not
    returned as a pass/fail column here – they are contested thresholds,
    report-only, and never gate anything (see those functions' docs).
    `format` is oriented to the EFA/ESEM model-fit statistics; for PCA
    the "statistics" are per-component eigenvalues.

  - `"nodes"` – Forbes-extension pruning annotations (requires
    `prune != "none"` when the object was created). One row per factor
    across all levels: `id`, `level`, `pruned`, `prune_reason`. Returns
    an empty data frame with the same columns when no pruning was
    applied.

  - `"scores"` – long-format per-observation factor scores (requires
    `keep_scores = TRUE` at fit time or use
    [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
    for on-the-fly computation). Columns: `obs` (row index), `level`,
    `factor`, `score`.

- primary_only:

  For `what = "edges"` only. When `TRUE`, returns just each factor's
  primary-parent edge (`is_primary == TRUE`) – the lineage tree that the
  diagram draws as solid arrows. Default `FALSE` (all edges). Errors for
  any other value of `what`.

- sort:

  For `what = "edges"` only. One of `"none"` (default, natural order) or
  `"strength"` (descending `|r|`). Ignored for all other values of
  `what`.

- format:

  For `what = "fit"` only. One of `"long"` (default, one row per
  statistic x level) or `"wide"` (one row per level, one column per
  statistic). Errors for all other values of `what`.

- conf_level:

  For `what = "loadings"` only. Confidence level for the loading
  intervals; default `0.95`. The intervals are computed as
  `loading ± qnorm((1 + conf_level) / 2) * se` and are `NA` for engines
  that carry no SEs (PCA, EFA). Errors for all other values of `what`.

- ...:

  Ignored.

## Value

A data frame (class `data.frame`).

## See also

[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)

## Examples

``` r
x <- ackwards(sim16, k_max = 5)
tidy(x) # edges in natural order
#>    from   to level_from level_to             r is_primary above_cut
#> 1  m1f1 m2f1          1        2  0.7072961126       TRUE      TRUE
#> 2  m1f1 m2f2          1        2  0.7069173991       TRUE      TRUE
#> 3  m2f1 m3f1          2        3 -0.0113423946      FALSE     FALSE
#> 4  m2f1 m3f2          2        3  0.6970929809       TRUE      TRUE
#> 5  m2f1 m3f3          2        3  0.7168910141       TRUE      TRUE
#> 6  m2f2 m3f1          2        3  0.9955932310       TRUE      TRUE
#> 7  m2f2 m3f2          2        3  0.0746157674      FALSE     FALSE
#> 8  m2f2 m3f3          2        3 -0.0568032178      FALSE     FALSE
#> 9  m3f1 m4f1          3        4 -0.0164580578      FALSE     FALSE
#> 10 m3f1 m4f2          3        4 -0.0048128349      FALSE     FALSE
#> 11 m3f1 m4f3          3        4  0.6639134885       TRUE      TRUE
#> 12 m3f1 m4f4          3        4  0.7476127666       TRUE      TRUE
#> 13 m3f2 m4f1          3        4  0.9406477359       TRUE      TRUE
#> 14 m3f2 m4f2          3        4  0.0317352150      FALSE     FALSE
#> 15 m3f2 m4f3          3        4  0.2627656201      FALSE     FALSE
#> 16 m3f2 m4f4          3        4 -0.2124357359      FALSE     FALSE
#> 17 m3f3 m4f1          3        4  0.0512578921      FALSE     FALSE
#> 18 m3f3 m4f2          3        4  0.9707784269       TRUE      TRUE
#> 19 m3f3 m4f3          3        4 -0.1715840980      FALSE     FALSE
#> 20 m3f3 m4f4          3        4  0.1597522194      FALSE     FALSE
#> 21 m4f1 m5f1          4        5  0.9997725525       TRUE      TRUE
#> 22 m4f1 m5f2          4        5  0.0020429352      FALSE     FALSE
#> 23 m4f1 m5f3          4        5 -0.0043948033      FALSE     FALSE
#> 24 m4f1 m5f4          4        5  0.0034107252      FALSE     FALSE
#> 25 m4f1 m5f5          4        5 -0.0204871256      FALSE     FALSE
#> 26 m4f2 m5f1          4        5 -0.0019639857      FALSE     FALSE
#> 27 m4f2 m5f2          4        5  0.9999850115       TRUE      TRUE
#> 28 m4f2 m5f3          4        5 -0.0025910717      FALSE     FALSE
#> 29 m4f2 m5f4          4        5 -0.0012689323      FALSE     FALSE
#> 30 m4f2 m5f5          4        5  0.0042184956      FALSE     FALSE
#> 31 m4f3 m5f1          4        5  0.0043781997      FALSE     FALSE
#> 32 m4f3 m5f2          4        5  0.0026035592      FALSE     FALSE
#> 33 m4f3 m5f3          4        5  0.9999840561       TRUE      TRUE
#> 34 m4f3 m5f4          4        5  0.0024297386      FALSE     FALSE
#> 35 m4f3 m5f5          4        5 -0.0001914960      FALSE     FALSE
#> 36 m4f4 m5f1          4        5 -0.0016399230      FALSE     FALSE
#> 37 m4f4 m5f2          4        5  0.0008901375      FALSE     FALSE
#> 38 m4f4 m5f3          4        5 -0.0023992589      FALSE     FALSE
#> 39 m4f4 m5f4          4        5  0.9962530528       TRUE      TRUE
#> 40 m4f4 m5f5          4        5  0.0864327285       TRUE     FALSE
tidy(x, sort = "strength") # strongest edges first
#>    from   to level_from level_to             r is_primary above_cut
#> 1  m4f2 m5f2          4        5  0.9999850115       TRUE      TRUE
#> 2  m4f3 m5f3          4        5  0.9999840561       TRUE      TRUE
#> 3  m4f1 m5f1          4        5  0.9997725525       TRUE      TRUE
#> 4  m4f4 m5f4          4        5  0.9962530528       TRUE      TRUE
#> 5  m2f2 m3f1          2        3  0.9955932310       TRUE      TRUE
#> 6  m3f3 m4f2          3        4  0.9707784269       TRUE      TRUE
#> 7  m3f2 m4f1          3        4  0.9406477359       TRUE      TRUE
#> 8  m3f1 m4f4          3        4  0.7476127666       TRUE      TRUE
#> 9  m2f1 m3f3          2        3  0.7168910141       TRUE      TRUE
#> 10 m1f1 m2f1          1        2  0.7072961126       TRUE      TRUE
#> 11 m1f1 m2f2          1        2  0.7069173991       TRUE      TRUE
#> 12 m2f1 m3f2          2        3  0.6970929809       TRUE      TRUE
#> 13 m3f1 m4f3          3        4  0.6639134885       TRUE      TRUE
#> 14 m3f2 m4f3          3        4  0.2627656201      FALSE     FALSE
#> 15 m3f2 m4f4          3        4 -0.2124357359      FALSE     FALSE
#> 16 m3f3 m4f3          3        4 -0.1715840980      FALSE     FALSE
#> 17 m3f3 m4f4          3        4  0.1597522194      FALSE     FALSE
#> 18 m4f4 m5f5          4        5  0.0864327285       TRUE     FALSE
#> 19 m2f2 m3f2          2        3  0.0746157674      FALSE     FALSE
#> 20 m2f2 m3f3          2        3 -0.0568032178      FALSE     FALSE
#> 21 m3f3 m4f1          3        4  0.0512578921      FALSE     FALSE
#> 22 m3f2 m4f2          3        4  0.0317352150      FALSE     FALSE
#> 23 m4f1 m5f5          4        5 -0.0204871256      FALSE     FALSE
#> 24 m3f1 m4f1          3        4 -0.0164580578      FALSE     FALSE
#> 25 m2f1 m3f1          2        3 -0.0113423946      FALSE     FALSE
#> 26 m3f1 m4f2          3        4 -0.0048128349      FALSE     FALSE
#> 27 m4f1 m5f3          4        5 -0.0043948033      FALSE     FALSE
#> 28 m4f3 m5f1          4        5  0.0043781997      FALSE     FALSE
#> 29 m4f2 m5f5          4        5  0.0042184956      FALSE     FALSE
#> 30 m4f1 m5f4          4        5  0.0034107252      FALSE     FALSE
#> 31 m4f3 m5f2          4        5  0.0026035592      FALSE     FALSE
#> 32 m4f2 m5f3          4        5 -0.0025910717      FALSE     FALSE
#> 33 m4f3 m5f4          4        5  0.0024297386      FALSE     FALSE
#> 34 m4f4 m5f3          4        5 -0.0023992589      FALSE     FALSE
#> 35 m4f1 m5f2          4        5  0.0020429352      FALSE     FALSE
#> 36 m4f2 m5f1          4        5 -0.0019639857      FALSE     FALSE
#> 37 m4f4 m5f1          4        5 -0.0016399230      FALSE     FALSE
#> 38 m4f2 m5f4          4        5 -0.0012689323      FALSE     FALSE
#> 39 m4f4 m5f2          4        5  0.0008901375      FALSE     FALSE
#> 40 m4f3 m5f5          4        5 -0.0001914960      FALSE     FALSE
tidy(x, primary_only = TRUE) # just the primary-parent lineage
#>    from   to level_from level_to          r is_primary above_cut
#> 1  m1f1 m2f1          1        2 0.70729611       TRUE      TRUE
#> 2  m1f1 m2f2          1        2 0.70691740       TRUE      TRUE
#> 3  m2f1 m3f2          2        3 0.69709298       TRUE      TRUE
#> 4  m2f1 m3f3          2        3 0.71689101       TRUE      TRUE
#> 5  m2f2 m3f1          2        3 0.99559323       TRUE      TRUE
#> 6  m3f1 m4f3          3        4 0.66391349       TRUE      TRUE
#> 7  m3f1 m4f4          3        4 0.74761277       TRUE      TRUE
#> 8  m3f2 m4f1          3        4 0.94064774       TRUE      TRUE
#> 9  m3f3 m4f2          3        4 0.97077843       TRUE      TRUE
#> 10 m4f1 m5f1          4        5 0.99977255       TRUE      TRUE
#> 11 m4f2 m5f2          4        5 0.99998501       TRUE      TRUE
#> 12 m4f3 m5f3          4        5 0.99998406       TRUE      TRUE
#> 13 m4f4 m5f4          4        5 0.99625305       TRUE      TRUE
#> 14 m4f4 m5f5          4        5 0.08643273       TRUE     FALSE
tidy(x, what = "loadings")
#>     level factor item      loading se ci_lower ci_upper
#> 1       1   m1f1   i1  0.461752431 NA       NA       NA
#> 2       1   m1f1   i2  0.485041500 NA       NA       NA
#> 3       1   m1f1   i3  0.514249968 NA       NA       NA
#> 4       1   m1f1   i4  0.515548605 NA       NA       NA
#> 5       1   m1f1   i5  0.575053347 NA       NA       NA
#> 6       1   m1f1   i6  0.568480393 NA       NA       NA
#> 7       1   m1f1   i7  0.589789932 NA       NA       NA
#> 8       1   m1f1   i8  0.508096441 NA       NA       NA
#> 9       1   m1f1   i9  0.563485394 NA       NA       NA
#> 10      1   m1f1  i10  0.551924658 NA       NA       NA
#> 11      1   m1f1  i11  0.546706411 NA       NA       NA
#> 12      1   m1f1  i12  0.563394632 NA       NA       NA
#> 13      1   m1f1  i13  0.519943165 NA       NA       NA
#> 14      1   m1f1  i14  0.525582317 NA       NA       NA
#> 15      1   m1f1  i15  0.502500732 NA       NA       NA
#> 16      1   m1f1  i16  0.481206189 NA       NA       NA
#> 17      2   m2f1   i1  0.670505125 NA       NA       NA
#> 18      2   m2f1   i2  0.688261361 NA       NA       NA
#> 19      2   m2f1   i3  0.666960388 NA       NA       NA
#> 20      2   m2f1   i4  0.695694541 NA       NA       NA
#> 21      2   m2f1   i5  0.677994613 NA       NA       NA
#> 22      2   m2f1   i6  0.671411581 NA       NA       NA
#> 23      2   m2f1   i7  0.688528565 NA       NA       NA
#> 24      2   m2f1   i8  0.646831586 NA       NA       NA
#> 25      2   m2f1   i9  0.142480573 NA       NA       NA
#> 26      2   m2f1  i10  0.106321144 NA       NA       NA
#> 27      2   m2f1  i11  0.099785669 NA       NA       NA
#> 28      2   m2f1  i12  0.132520254 NA       NA       NA
#> 29      2   m2f1  i13  0.045302307 NA       NA       NA
#> 30      2   m2f1  i14  0.056844031 NA       NA       NA
#> 31      2   m2f1  i15  0.027746191 NA       NA       NA
#> 32      2   m2f1  i16  0.002718709 NA       NA       NA
#> 33      2   m2f2   i1 -0.017672839 NA       NA       NA
#> 34      2   m2f2   i2 -0.002494047 NA       NA       NA
#> 35      2   m2f2   i3  0.060136416 NA       NA       NA
#> 36      2   m2f2   i4  0.033223911 NA       NA       NA
#> 37      2   m2f2   i5  0.135108279 NA       NA       NA
#> 38      2   m2f2   i6  0.132396786 NA       NA       NA
#> 39      2   m2f2   i7  0.145414946 NA       NA       NA
#> 40      2   m2f2   i8  0.071571268 NA       NA       NA
#> 41      2   m2f2   i9  0.654545268 NA       NA       NA
#> 42      2   m2f2  i10  0.674370338 NA       NA       NA
#> 43      2   m2f2  i11  0.673527623 NA       NA       NA
#> 44      2   m2f2  i12  0.664382532 NA       NA       NA
#> 45      2   m2f2  i13  0.690181088 NA       NA       NA
#> 46      2   m2f2  i14  0.686610283 NA       NA       NA
#> 47      2   m2f2  i15  0.683072674 NA       NA       NA
#> 48      2   m2f2  i16  0.677990466 NA       NA       NA
#> 49      3   m3f1   i1  0.018282326 NA       NA       NA
#> 50      3   m3f1   i2  0.030497878 NA       NA       NA
#> 51      3   m3f1   i3  0.088799850 NA       NA       NA
#> 52      3   m3f1   i4  0.065163580 NA       NA       NA
#> 53      3   m3f1   i5  0.088559024 NA       NA       NA
#> 54      3   m3f1   i6  0.085162586 NA       NA       NA
#> 55      3   m3f1   i7  0.100894156 NA       NA       NA
#> 56      3   m3f1   i8  0.021040215 NA       NA       NA
#> 57      3   m3f1   i9  0.635370326 NA       NA       NA
#> 58      3   m3f1  i10  0.654349429 NA       NA       NA
#> 59      3   m3f1  i11  0.656287051 NA       NA       NA
#> 60      3   m3f1  i12  0.643589738 NA       NA       NA
#> 61      3   m3f1  i13  0.702820782 NA       NA       NA
#> 62      3   m3f1  i14  0.704046359 NA       NA       NA
#> 63      3   m3f1  i15  0.701468675 NA       NA       NA
#> 64      3   m3f1  i16  0.697838670 NA       NA       NA
#> 65      3   m3f2   i1  0.132996558 NA       NA       NA
#> 66      3   m3f2   i2  0.167151079 NA       NA       NA
#> 67      3   m3f2   i3  0.189869691 NA       NA       NA
#> 68      3   m3f2   i4  0.181206856 NA       NA       NA
#> 69      3   m3f2   i5  0.775819255 NA       NA       NA
#> 70      3   m3f2   i6  0.776938349 NA       NA       NA
#> 71      3   m3f2   i7  0.767129561 NA       NA       NA
#> 72      3   m3f2   i8  0.784709189 NA       NA       NA
#> 73      3   m3f2   i9  0.260572503 NA       NA       NA
#> 74      3   m3f2  i10  0.245798186 NA       NA       NA
#> 75      3   m3f2  i11  0.220477487 NA       NA       NA
#> 76      3   m3f2  i12  0.267289841 NA       NA       NA
#> 77      3   m3f2  i13 -0.040980546 NA       NA       NA
#> 78      3   m3f2  i14 -0.070825420 NA       NA       NA
#> 79      3   m3f2  i15 -0.096078981 NA       NA       NA
#> 80      3   m3f2  i16 -0.122682914 NA       NA       NA
#> 81      3   m3f3   i1  0.806261359 NA       NA       NA
#> 82      3   m3f3   i2  0.798011727 NA       NA       NA
#> 83      3   m3f3   i3  0.747129970 NA       NA       NA
#> 84      3   m3f3   i4  0.795261223 NA       NA       NA
#> 85      3   m3f3   i5  0.192750257 NA       NA       NA
#> 86      3   m3f3   i6  0.182425580 NA       NA       NA
#> 87      3   m3f3   i7  0.216089072 NA       NA       NA
#> 88      3   m3f3   i8  0.139567888 NA       NA       NA
#> 89      3   m3f3   i9 -0.044575909 NA       NA       NA
#> 90      3   m3f3  i10 -0.080348554 NA       NA       NA
#> 91      3   m3f3  i11 -0.064812882 NA       NA       NA
#> 92      3   m3f3  i12 -0.064871464 NA       NA       NA
#> 93      3   m3f3  i13  0.114161326 NA       NA       NA
#> 94      3   m3f3  i14  0.159301070 NA       NA       NA
#> 95      3   m3f3  i15  0.143227502 NA       NA       NA
#> 96      3   m3f3  i16  0.134128155 NA       NA       NA
#> 97      4   m4f1   i1  0.128686176 NA       NA       NA
#> 98      4   m4f1   i2  0.149030299 NA       NA       NA
#> 99      4   m4f1   i3  0.143394315 NA       NA       NA
#> 100     4   m4f1   i4  0.164542205 NA       NA       NA
#> 101     4   m4f1   i5  0.802597931 NA       NA       NA
#> 102     4   m4f1   i6  0.810490299 NA       NA       NA
#> 103     4   m4f1   i7  0.791631093 NA       NA       NA
#> 104     4   m4f1   i8  0.832623843 NA       NA       NA
#> 105     4   m4f1   i9  0.082991593 NA       NA       NA
#> 106     4   m4f1  i10  0.069135718 NA       NA       NA
#> 107     4   m4f1  i11  0.049095872 NA       NA       NA
#> 108     4   m4f1  i12  0.101317088 NA       NA       NA
#> 109     4   m4f1  i13  0.090878291 NA       NA       NA
#> 110     4   m4f1  i14  0.058515420 NA       NA       NA
#> 111     4   m4f1  i15  0.037287536 NA       NA       NA
#> 112     4   m4f1  i16  0.005030840 NA       NA       NA
#> 113     4   m4f2   i1  0.813408272 NA       NA       NA
#> 114     4   m4f2   i2  0.814344929 NA       NA       NA
#> 115     4   m4f2   i3  0.782025429 NA       NA       NA
#> 116     4   m4f2   i4  0.809823496 NA       NA       NA
#> 117     4   m4f2   i5  0.165603797 NA       NA       NA
#> 118     4   m4f2   i6  0.150442425 NA       NA       NA
#> 119     4   m4f2   i7  0.190612643 NA       NA       NA
#> 120     4   m4f2   i8  0.098061225 NA       NA       NA
#> 121     4   m4f2   i9  0.067952773 NA       NA       NA
#> 122     4   m4f2  i10  0.031112689 NA       NA       NA
#> 123     4   m4f2  i11  0.043242680 NA       NA       NA
#> 124     4   m4f2  i12  0.039067952 NA       NA       NA
#> 125     4   m4f2  i13  0.010228110 NA       NA       NA
#> 126     4   m4f2  i14  0.057767882 NA       NA       NA
#> 127     4   m4f2  i15  0.039026995 NA       NA       NA
#> 128     4   m4f2  i16  0.034210703 NA       NA       NA
#> 129     4   m4f3   i1 -0.015409920 NA       NA       NA
#> 130     4   m4f3   i2  0.025695304 NA       NA       NA
#> 131     4   m4f3   i3  0.126581449 NA       NA       NA
#> 132     4   m4f3   i4  0.046794854 NA       NA       NA
#> 133     4   m4f3   i5  0.099123172 NA       NA       NA
#> 134     4   m4f3   i6  0.084120174 NA       NA       NA
#> 135     4   m4f3   i7  0.108694221 NA       NA       NA
#> 136     4   m4f3   i8  0.018603667 NA       NA       NA
#> 137     4   m4f3   i9  0.800528531 NA       NA       NA
#> 138     4   m4f3  i10  0.810953636 NA       NA       NA
#> 139     4   m4f3  i11  0.796816493 NA       NA       NA
#> 140     4   m4f3  i12  0.784529901 NA       NA       NA
#> 141     4   m4f3  i13  0.162504347 NA       NA       NA
#> 142     4   m4f3  i14  0.161065367 NA       NA       NA
#> 143     4   m4f3  i15  0.148774737 NA       NA       NA
#> 144     4   m4f3  i16  0.154760740 NA       NA       NA
#> 145     4   m4f4   i1  0.046208286 NA       NA       NA
#> 146     4   m4f4   i2  0.026498312 NA       NA       NA
#> 147     4   m4f4   i3  0.014558967 NA       NA       NA
#> 148     4   m4f4   i4  0.054441870 NA       NA       NA
#> 149     4   m4f4   i5  0.049164544 NA       NA       NA
#> 150     4   m4f4   i6  0.058020972 NA       NA       NA
#> 151     4   m4f4   i7  0.057083956 NA       NA       NA
#> 152     4   m4f4   i8  0.030583096 NA       NA       NA
#> 153     4   m4f4   i9  0.141224932 NA       NA       NA
#> 154     4   m4f4  i10  0.156811061 NA       NA       NA
#> 155     4   m4f4  i11  0.171594148 NA       NA       NA
#> 156     4   m4f4  i12  0.166644111 NA       NA       NA
#> 157     4   m4f4  i13  0.797841995 NA       NA       NA
#> 158     4   m4f4  i14  0.800352794 NA       NA       NA
#> 159     4   m4f4  i15  0.807231575 NA       NA       NA
#> 160     4   m4f4  i16  0.796319167 NA       NA       NA
#> 161     5   m5f1   i1  0.125845182 NA       NA       NA
#> 162     5   m5f1   i2  0.146281625 NA       NA       NA
#> 163     5   m5f1   i3  0.141457000 NA       NA       NA
#> 164     5   m5f1   i4  0.166373116 NA       NA       NA
#> 165     5   m5f1   i5  0.802636385 NA       NA       NA
#> 166     5   m5f1   i6  0.805601696 NA       NA       NA
#> 167     5   m5f1   i7  0.794174370 NA       NA       NA
#> 168     5   m5f1   i8  0.834242419 NA       NA       NA
#> 169     5   m5f1   i9  0.079330358 NA       NA       NA
#> 170     5   m5f1  i10  0.071440903 NA       NA       NA
#> 171     5   m5f1  i11  0.053199270 NA       NA       NA
#> 172     5   m5f1  i12  0.111309787 NA       NA       NA
#> 173     5   m5f1  i13  0.088459037 NA       NA       NA
#> 174     5   m5f1  i14  0.055792727 NA       NA       NA
#> 175     5   m5f1  i15  0.032789698 NA       NA       NA
#> 176     5   m5f1  i16  0.011575744 NA       NA       NA
#> 177     5   m5f2   i1  0.813880739 NA       NA       NA
#> 178     5   m5f2   i2  0.814971814 NA       NA       NA
#> 179     5   m5f2   i3  0.782834508 NA       NA       NA
#> 180     5   m5f2   i4  0.809628684 NA       NA       NA
#> 181     5   m5f2   i5  0.167503047 NA       NA       NA
#> 182     5   m5f2   i6  0.153331666 NA       NA       NA
#> 183     5   m5f2   i7  0.192001129 NA       NA       NA
#> 184     5   m5f2   i8  0.099430502 NA       NA       NA
#> 185     5   m5f2   i9  0.071729201 NA       NA       NA
#> 186     5   m5f2  i10  0.033692262 NA       NA       NA
#> 187     5   m5f2  i11  0.045365102 NA       NA       NA
#> 188     5   m5f2  i12  0.040036629 NA       NA       NA
#> 189     5   m5f2  i13  0.011914118 NA       NA       NA
#> 190     5   m5f2  i14  0.059428222 NA       NA       NA
#> 191     5   m5f2  i15  0.040979436 NA       NA       NA
#> 192     5   m5f2  i16  0.033839573 NA       NA       NA
#> 193     5   m5f3   i1 -0.018210287 NA       NA       NA
#> 194     5   m5f3   i2  0.022847977 NA       NA       NA
#> 195     5   m5f3   i3  0.123874092 NA       NA       NA
#> 196     5   m5f3   i4  0.043893857 NA       NA       NA
#> 197     5   m5f3   i5  0.095050271 NA       NA       NA
#> 198     5   m5f3   i6  0.079955315 NA       NA       NA
#> 199     5   m5f3   i7  0.104624657 NA       NA       NA
#> 200     5   m5f3   i8  0.014647210 NA       NA       NA
#> 201     5   m5f3   i9  0.799531025 NA       NA       NA
#> 202     5   m5f3  i10  0.810165899 NA       NA       NA
#> 203     5   m5f3  i11  0.796079656 NA       NA       NA
#> 204     5   m5f3  i12  0.783678485 NA       NA       NA
#> 205     5   m5f3  i13  0.160134023 NA       NA       NA
#> 206     5   m5f3  i14  0.158704883 NA       NA       NA
#> 207     5   m5f3  i15  0.146512650 NA       NA       NA
#> 208     5   m5f3  i16  0.152849170 NA       NA       NA
#> 209     5   m5f4   i1  0.049869940 NA       NA       NA
#> 210     5   m5f4   i2  0.030875158 NA       NA       NA
#> 211     5   m5f4   i3  0.018057780 NA       NA       NA
#> 212     5   m5f4   i4  0.039945216 NA       NA       NA
#> 213     5   m5f4   i5  0.050944187 NA       NA       NA
#> 214     5   m5f4   i6  0.080103321 NA       NA       NA
#> 215     5   m5f4   i7  0.048270051 NA       NA       NA
#> 216     5   m5f4   i8  0.025018104 NA       NA       NA
#> 217     5   m5f4   i9  0.171116784 NA       NA       NA
#> 218     5   m5f4  i10  0.162190188 NA       NA       NA
#> 219     5   m5f4  i11  0.168862429 NA       NA       NA
#> 220     5   m5f4  i12  0.139322561 NA       NA       NA
#> 221     5   m5f4  i13  0.802972834 NA       NA       NA
#> 222     5   m5f4  i14  0.806163196 NA       NA       NA
#> 223     5   m5f4  i15  0.820241598 NA       NA       NA
#> 224     5   m5f4  i16  0.763490591 NA       NA       NA
#> 225     5   m5f5   i1 -0.046701926 NA       NA       NA
#> 226     5   m5f5   i2 -0.054283683 NA       NA       NA
#> 227     5   m5f5   i3 -0.041636964 NA       NA       NA
#> 228     5   m5f5   i4  0.165490526 NA       NA       NA
#> 229     5   m5f5   i5 -0.002239278 NA       NA       NA
#> 230     5   m5f5   i6 -0.236088047 NA       NA       NA
#> 231     5   m5f5   i7  0.120061811 NA       NA       NA
#> 232     5   m5f5   i8  0.080680811 NA       NA       NA
#> 233     5   m5f5   i9 -0.315461030 NA       NA       NA
#> 234     5   m5f5  i10 -0.031706091 NA       NA       NA
#> 235     5   m5f5  i11  0.061565807 NA       NA       NA
#> 236     5   m5f5  i12  0.345595171 NA       NA       NA
#> 237     5   m5f5  i13 -0.018551750 NA       NA       NA
#> 238     5   m5f5  i14 -0.027424572 NA       NA       NA
#> 239     5   m5f5  i15 -0.110696540 NA       NA       NA
#> 240     5   m5f5  i16  0.417028618 NA       NA       NA
tidy(x, what = "variance")
#>    level factor proportion cumulative
#> 1      1   m1f1 0.28172978  0.2817298
#> 2      2   m2f1 0.23251408  0.2325141
#> 3      2   m2f2 0.23246133  0.4649754
#> 4      3   m3f1 0.23028614  0.2302861
#> 5      3   m3f2 0.17522787  0.4055140
#> 6      3   m3f3 0.16924345  0.5747575
#> 7      4   m4f1 0.17155300  0.1715530
#> 8      4   m4f2 0.16895577  0.3405088
#> 9      4   m4f3 0.16846890  0.5089777
#> 10     4   m4f4 0.16753770  0.6765154
#> 11     5   m5f1 0.17147388  0.1714739
#> 12     5   m5f2 0.16935421  0.3408281
#> 13     5   m5f3 0.16774286  0.5085710
#> 14     5   m5f4 0.16695127  0.6755222
#> 15     5   m5f5 0.03262035  0.7081426
tidy(x, what = "fit")
#>    level       statistic     value
#> 1      1 eigenvalue.m1f1 4.5076765
#> 2      2 eigenvalue.m2f1 4.5076765
#> 3      2 eigenvalue.m2f2 2.9319301
#> 4      3 eigenvalue.m3f1 4.5076765
#> 5      3 eigenvalue.m3f2 2.9319301
#> 6      3 eigenvalue.m3f3 1.7565127
#> 7      4 eigenvalue.m4f1 4.5076765
#> 8      4 eigenvalue.m4f2 2.9319301
#> 9      4 eigenvalue.m4f3 1.7565127
#> 10     4 eigenvalue.m4f4 1.6281266
#> 11     5 eigenvalue.m5f1 4.5076765
#> 12     5 eigenvalue.m5f2 2.9319301
#> 13     5 eigenvalue.m5f3 1.7565127
#> 14     5 eigenvalue.m5f4 1.6281266
#> 15     5 eigenvalue.m5f5 0.5060352
tidy(x, what = "fit", format = "wide")
#>   level eigenvalue.m2f1 eigenvalue.m2f2 eigenvalue.m3f1 eigenvalue.m3f2
#> 1     2        4.507677         2.93193              NA              NA
#> 2     3              NA              NA        4.507677         2.93193
#> 3     4              NA              NA              NA              NA
#> 4     5              NA              NA              NA              NA
#>   eigenvalue.m3f3 eigenvalue.m4f1 eigenvalue.m4f2 eigenvalue.m4f3
#> 1              NA              NA              NA              NA
#> 2        1.756513              NA              NA              NA
#> 3              NA        4.507677         2.93193        1.756513
#> 4              NA              NA              NA              NA
#>   eigenvalue.m4f4 eigenvalue.m5f1 eigenvalue.m5f2 eigenvalue.m5f3
#> 1              NA              NA              NA              NA
#> 2              NA              NA              NA              NA
#> 3        1.628127              NA              NA              NA
#> 4              NA        4.507677         2.93193        1.756513
#>   eigenvalue.m5f4 eigenvalue.m5f5
#> 1              NA              NA
#> 2              NA              NA
#> 3              NA              NA
#> 4        1.628127       0.5060352
```
