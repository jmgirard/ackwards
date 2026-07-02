# Bootstrap confidence intervals for between-level edges

Attaches nonparametric bootstrap standard errors and percentile
confidence intervals to every between-level correlation (edge) of a
fitted bass-ackwards hierarchy. Every edge
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
reports is a point estimate; `boot_edges()` quantifies its sampling
uncertainty, which matters most where a hard threshold consumes the
estimate –
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)'s
`|r| >= redundancy_r` redundancy rule, and the Forbes (2023) practice of
interpreting the strongest all-pairs edge.

## Usage

``` r
boot_edges(x, ...)

# S3 method for class 'ackwards'
boot_edges(x, data, n_boot = 1000L, conf = 0.95, seed = NULL, ...)
```

## Arguments

- x:

  An `ackwards` object fit with `engine = "pca"` or `"efa"` on a
  `"pearson"` or `"spearman"` basis. ESEM objects are not supported
  (refitting `n_boot` lavaan hierarchies needs its own performance
  treatment), nor are polychoric-basis objects (estimating a polychoric
  matrix in every resample is slow and unstable) or objects fit from a
  correlation matrix (resampling needs rows).

- ...:

  Reserved for future arguments.

- data:

  The raw item data the model was fit on. Required – the `ackwards`
  object deliberately does not store raw data (light core), so it must
  be re-supplied here. Columns are matched by name against the fit; a
  warning is issued if the data do not look like the fit data.

- n_boot:

  Number of bootstrap replicates. Default `1000L`; larger values (2000+)
  are advisable for published interval endpoints (Efron & Tibshirani,
  1993). Each replicate refits the full hierarchy, so cost scales
  linearly.

- conf:

  Confidence level for the percentile intervals. Default `0.95`.

- seed:

  Integer seed for reproducible resampling. `NULL` (default) uses the
  current RNG state.

## Value

`x`, invisibly modified: the `$boot` element is populated with

- edges:

  Data frame with one row per edge (aligned with
  `tidy(x, what = "edges")`): `from`, `to`, `level_from`, `level_to`,
  `r` (the full-sample point estimate), `se` (bootstrap standard error),
  `lo`, `hi` (percentile interval endpoints), and `n_boot_ok` (usable
  replicates for that edge).

- n_boot, conf, seed:

  The request.

After calling `boot_edges()`, `tidy(x, what = "edges")` gains `se`,
`lo`, and `hi` columns, and `print(x)`/`summary(x)` note the interval
coverage.

## Details

For each of `n_boot` replicates, `n` rows are resampled with
replacement, the correlation matrix is recomputed with the same basis
and missing-data routine used at fit time, and the full hierarchy is
refit. Each replicate level is then **anchored to the full-sample
solution** – its factors matched (greedy max-\|r\| with removal) and
sign-oriented against the full-sample factors on the full-sample
correlation matrix – before its edges are computed. Without anchoring,
factor label switching and sign flipping across replicates would corrupt
the pooled edge distributions; this is the same matching machinery
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
uses.

All resample indices are drawn upfront from `seed`, so results are
reproducible and identical whether replicates run serially or in
parallel. Replicate fits are dispatched through future.apply when it is
installed and the user has set a
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
(serial otherwise, as in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)'s
ESEM engine).

## What the intervals do and do not fix

Per-edge intervals make sampling uncertainty **visible**: an edge whose
interval straddles
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)'s
`redundancy_r` threshold should not be treated as decisively above or
below it. They do **not** correct the selection bias of scanning many
edges for the strongest one – the maximum of hundreds of correlations
capitalizes on chance even when every individual interval is honest.
Treat the intervals as per-edge error bars, not a familywise inference.

## Failed replicates

A replicate whose hierarchy fails to converge (in full or at some
levels) contributes `NA` to the affected edges and is dropped from their
distributions – convergence is data, not an error. The usable replicate
count is reported per edge in `n_boot_ok`, and a message summarises any
shortfall.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. Chapman & Hall.

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

## See also

[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) for
the thresholded rules the intervals contextualise,
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
for split-half replicability of the factors themselves,
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
for the augmented edge table.

## Examples

``` r
# \donttest{
x <- ackwards(bfi25, k_max = 3)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
x <- boot_edges(x, bfi25, n_boot = 100, seed = 1)
#> ℹ Fitting 100 bootstrap replicates (PCA, k = 1-3)...
#> ✔ Fitting 100 bootstrap replicates (PCA, k = 1-3)... [3s]
#> 
x$boot$edges
#>   from   to level_from level_to           r         se          lo          hi
#> 1 m1f1 m2f1          1        2  0.87806671 0.02742806  0.80759107  0.91632427
#> 2 m1f1 m2f2          1        2  0.47853824 0.04807772  0.40043672  0.58969785
#> 3 m2f1 m3f1          2        3  0.90501945 0.03590910  0.81103593  0.95162512
#> 4 m2f1 m3f2          2        3 -0.08491543 0.02438128 -0.12484971 -0.03685198
#> 5 m2f1 m3f3          2        3  0.41680832 0.07425865  0.29177975  0.57911195
#> 6 m2f2 m3f1          2        3 -0.03221772 0.03653585 -0.07908454  0.07245676
#> 7 m2f2 m3f2          2        3  0.96337357 0.02841986  0.89415596  0.99601862
#> 8 m2f2 m3f3          2        3  0.26622056 0.11056458  0.02974731  0.44542779
#>   n_boot_ok
#> 1       100
#> 2       100
#> 3       100
#> 4       100
#> 5       100
#> 6       100
#> 7       100
#> 8       100
head(tidy(x)) # now carries se / lo / hi
#>   from   to level_from level_to           r is_primary above_cut         se
#> 1 m1f1 m2f1          1        2  0.87806671       TRUE      TRUE 0.02742806
#> 2 m1f1 m2f2          1        2  0.47853824       TRUE      TRUE 0.04807772
#> 3 m2f1 m3f1          2        3  0.90501945       TRUE      TRUE 0.03590910
#> 4 m2f1 m3f2          2        3 -0.08491543      FALSE     FALSE 0.02438128
#> 5 m2f1 m3f3          2        3  0.41680832       TRUE      TRUE 0.07425865
#> 6 m2f2 m3f1          2        3 -0.03221772      FALSE     FALSE 0.03653585
#>            lo          hi n_boot_ok
#> 1  0.80759107  0.91632427       100
#> 2  0.40043672  0.58969785       100
#> 3  0.81103593  0.95162512       100
#> 4 -0.12484971 -0.03685198       100
#> 5  0.29177975  0.57911195       100
#> 6 -0.07908454  0.07245676       100
# }
```
