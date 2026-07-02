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
    `above_cut`.

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
x <- ackwards(bfi25, k_max = 5)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
tidy(x) # edges in natural order
#>    from   to level_from level_to            r is_primary above_cut
#> 1  m1f1 m2f1          1        2  0.878066713       TRUE      TRUE
#> 2  m1f1 m2f2          1        2  0.478538241       TRUE      TRUE
#> 3  m2f1 m3f1          2        3  0.905019446       TRUE      TRUE
#> 4  m2f1 m3f2          2        3 -0.084915431      FALSE     FALSE
#> 5  m2f1 m3f3          2        3  0.416808316       TRUE      TRUE
#> 6  m2f2 m3f1          2        3 -0.032217716      FALSE     FALSE
#> 7  m2f2 m3f2          2        3  0.963373570       TRUE      TRUE
#> 8  m2f2 m3f3          2        3  0.266220556      FALSE     FALSE
#> 9  m3f1 m4f1          3        4  0.996745750       TRUE      TRUE
#> 10 m3f1 m4f2          3        4 -0.006661471      FALSE     FALSE
#> 11 m3f1 m4f3          3        4  0.042095931      FALSE     FALSE
#> 12 m3f1 m4f4          3        4  0.068421245      FALSE     FALSE
#> 13 m3f2 m4f1          3        4  0.016233327      FALSE     FALSE
#> 14 m3f2 m4f2          3        4  0.988176397       TRUE      TRUE
#> 15 m3f2 m4f3          3        4  0.018059508      FALSE     FALSE
#> 16 m3f2 m4f4          3        4 -0.151386067      FALSE     FALSE
#> 17 m3f3 m4f1          3        4 -0.070214243      FALSE     FALSE
#> 18 m3f3 m4f2          3        4  0.061809300      FALSE     FALSE
#> 19 m3f3 m4f3          3        4  0.861695922       TRUE      TRUE
#> 20 m3f3 m4f4          3        4  0.498728091       TRUE      TRUE
#> 21 m4f1 m5f1          4        5  0.806394128       TRUE      TRUE
#> 22 m4f1 m5f2          4        5 -0.001363543      FALSE     FALSE
#> 23 m4f1 m5f3          4        5  0.006259133      FALSE     FALSE
#> 24 m4f1 m5f4          4        5  0.589479673       TRUE      TRUE
#> 25 m4f1 m5f5          4        5 -0.046916846      FALSE     FALSE
#> 26 m4f2 m5f1          4        5  0.030991419      FALSE     FALSE
#> 27 m4f2 m5f2          4        5  0.998569303       TRUE      TRUE
#> 28 m4f2 m5f3          4        5  0.015523208      FALSE     FALSE
#> 29 m4f2 m5f4          4        5 -0.040546852      FALSE     FALSE
#> 30 m4f2 m5f5          4        5 -0.003723118      FALSE     FALSE
#> 31 m4f3 m5f1          4        5 -0.105831031      FALSE     FALSE
#> 32 m4f3 m5f2          4        5 -0.006424971      FALSE     FALSE
#> 33 m4f3 m5f3          4        5  0.984798330       TRUE      TRUE
#> 34 m4f3 m5f4          4        5  0.135975178      FALSE     FALSE
#> 35 m4f3 m5f5          4        5  0.021012193      FALSE     FALSE
#> 36 m4f4 m5f1          4        5  0.191253503      FALSE     FALSE
#> 37 m4f4 m5f2          4        5 -0.010261396      FALSE     FALSE
#> 38 m4f4 m5f3          4        5  0.025504732      FALSE     FALSE
#> 39 m4f4 m5f4          4        5 -0.185238698      FALSE     FALSE
#> 40 m4f4 m5f5          4        5  0.963510734       TRUE      TRUE
tidy(x, sort = "strength") # strongest edges first
#>    from   to level_from level_to            r is_primary above_cut
#> 1  m4f2 m5f2          4        5  0.998569303       TRUE      TRUE
#> 2  m3f1 m4f1          3        4  0.996745750       TRUE      TRUE
#> 3  m3f2 m4f2          3        4  0.988176397       TRUE      TRUE
#> 4  m4f3 m5f3          4        5  0.984798330       TRUE      TRUE
#> 5  m4f4 m5f5          4        5  0.963510734       TRUE      TRUE
#> 6  m2f2 m3f2          2        3  0.963373570       TRUE      TRUE
#> 7  m2f1 m3f1          2        3  0.905019446       TRUE      TRUE
#> 8  m1f1 m2f1          1        2  0.878066713       TRUE      TRUE
#> 9  m3f3 m4f3          3        4  0.861695922       TRUE      TRUE
#> 10 m4f1 m5f1          4        5  0.806394128       TRUE      TRUE
#> 11 m4f1 m5f4          4        5  0.589479673       TRUE      TRUE
#> 12 m3f3 m4f4          3        4  0.498728091       TRUE      TRUE
#> 13 m1f1 m2f2          1        2  0.478538241       TRUE      TRUE
#> 14 m2f1 m3f3          2        3  0.416808316       TRUE      TRUE
#> 15 m2f2 m3f3          2        3  0.266220556      FALSE     FALSE
#> 16 m4f4 m5f1          4        5  0.191253503      FALSE     FALSE
#> 17 m4f4 m5f4          4        5 -0.185238698      FALSE     FALSE
#> 18 m3f2 m4f4          3        4 -0.151386067      FALSE     FALSE
#> 19 m4f3 m5f4          4        5  0.135975178      FALSE     FALSE
#> 20 m4f3 m5f1          4        5 -0.105831031      FALSE     FALSE
#> 21 m2f1 m3f2          2        3 -0.084915431      FALSE     FALSE
#> 22 m3f3 m4f1          3        4 -0.070214243      FALSE     FALSE
#> 23 m3f1 m4f4          3        4  0.068421245      FALSE     FALSE
#> 24 m3f3 m4f2          3        4  0.061809300      FALSE     FALSE
#> 25 m4f1 m5f5          4        5 -0.046916846      FALSE     FALSE
#> 26 m3f1 m4f3          3        4  0.042095931      FALSE     FALSE
#> 27 m4f2 m5f4          4        5 -0.040546852      FALSE     FALSE
#> 28 m2f2 m3f1          2        3 -0.032217716      FALSE     FALSE
#> 29 m4f2 m5f1          4        5  0.030991419      FALSE     FALSE
#> 30 m4f4 m5f3          4        5  0.025504732      FALSE     FALSE
#> 31 m4f3 m5f5          4        5  0.021012193      FALSE     FALSE
#> 32 m3f2 m4f3          3        4  0.018059508      FALSE     FALSE
#> 33 m3f2 m4f1          3        4  0.016233327      FALSE     FALSE
#> 34 m4f2 m5f3          4        5  0.015523208      FALSE     FALSE
#> 35 m4f4 m5f2          4        5 -0.010261396      FALSE     FALSE
#> 36 m3f1 m4f2          3        4 -0.006661471      FALSE     FALSE
#> 37 m4f3 m5f2          4        5 -0.006424971      FALSE     FALSE
#> 38 m4f1 m5f3          4        5  0.006259133      FALSE     FALSE
#> 39 m4f2 m5f5          4        5 -0.003723118      FALSE     FALSE
#> 40 m4f1 m5f2          4        5 -0.001363543      FALSE     FALSE
tidy(x, primary_only = TRUE) # just the primary-parent lineage
#>    from   to level_from level_to         r is_primary above_cut
#> 1  m1f1 m2f1          1        2 0.8780667       TRUE      TRUE
#> 2  m1f1 m2f2          1        2 0.4785382       TRUE      TRUE
#> 3  m2f1 m3f1          2        3 0.9050194       TRUE      TRUE
#> 4  m2f1 m3f3          2        3 0.4168083       TRUE      TRUE
#> 5  m2f2 m3f2          2        3 0.9633736       TRUE      TRUE
#> 6  m3f1 m4f1          3        4 0.9967458       TRUE      TRUE
#> 7  m3f2 m4f2          3        4 0.9881764       TRUE      TRUE
#> 8  m3f3 m4f3          3        4 0.8616959       TRUE      TRUE
#> 9  m3f3 m4f4          3        4 0.4987281       TRUE      TRUE
#> 10 m4f1 m5f1          4        5 0.8063941       TRUE      TRUE
#> 11 m4f1 m5f4          4        5 0.5894797       TRUE      TRUE
#> 12 m4f2 m5f2          4        5 0.9985693       TRUE      TRUE
#> 13 m4f3 m5f3          4        5 0.9847983       TRUE      TRUE
#> 14 m4f4 m5f5          4        5 0.9635107       TRUE      TRUE
tidy(x, what = "loadings")
#>     level factor item       loading se ci_lower ci_upper
#> 1       1   m1f1   A1 -0.3061961995 NA       NA       NA
#> 2       1   m1f1   A2  0.5402566804 NA       NA       NA
#> 3       1   m1f1   A3  0.6024065155 NA       NA       NA
#> 4       1   m1f1   A4  0.4483209246 NA       NA       NA
#> 5       1   m1f1   A5  0.6459661082 NA       NA       NA
#> 6       1   m1f1   C1  0.3809485870 NA       NA       NA
#> 7       1   m1f1   C2  0.4073284603 NA       NA       NA
#> 8       1   m1f1   C3  0.3662536411 NA       NA       NA
#> 9       1   m1f1   C4 -0.4996871502 NA       NA       NA
#> 10      1   m1f1   C5 -0.5604936869 NA       NA       NA
#> 11      1   m1f1   E1 -0.4424410074 NA       NA       NA
#> 12      1   m1f1   E2 -0.6361801372 NA       NA       NA
#> 13      1   m1f1   E3  0.5552669387 NA       NA       NA
#> 14      1   m1f1   E4  0.5873197321 NA       NA       NA
#> 15      1   m1f1   E5  0.5711174529 NA       NA       NA
#> 16      1   m1f1   N1 -0.4328876727 NA       NA       NA
#> 17      1   m1f1   N2 -0.3711237506 NA       NA       NA
#> 18      1   m1f1   N3 -0.3678223061 NA       NA       NA
#> 19      1   m1f1   N4 -0.4898668897 NA       NA       NA
#> 20      1   m1f1   N5 -0.3149134890 NA       NA       NA
#> 21      1   m1f1   O1  0.3536240085 NA       NA       NA
#> 22      1   m1f1   O2 -0.1696017175 NA       NA       NA
#> 23      1   m1f1   O3  0.4269569343 NA       NA       NA
#> 24      1   m1f1   O4 -0.0450614665 NA       NA       NA
#> 25      1   m1f1   O5 -0.2428184266 NA       NA       NA
#> 26      2   m2f1   A1 -0.2974353300 NA       NA       NA
#> 27      2   m2f1   A2  0.6192390413 NA       NA       NA
#> 28      2   m2f1   A3  0.6903578930 NA       NA       NA
#> 29      2   m2f1   A4  0.4356493140 NA       NA       NA
#> 30      2   m2f1   A5  0.6611472978 NA       NA       NA
#> 31      2   m2f1   C1  0.3410425462 NA       NA       NA
#> 32      2   m2f1   C2  0.4253111751 NA       NA       NA
#> 33      2   m2f1   C3  0.3332485909 NA       NA       NA
#> 34      2   m2f1   C4 -0.3372972262 NA       NA       NA
#> 35      2   m2f1   C5 -0.3773528609 NA       NA       NA
#> 36      2   m2f1   E1 -0.5092361988 NA       NA       NA
#> 37      2   m2f1   E2 -0.5951785530 NA       NA       NA
#> 38      2   m2f1   E3  0.6624432525 NA       NA       NA
#> 39      2   m2f1   E4  0.6300976720 NA       NA       NA
#> 40      2   m2f1   E5  0.6122939866 NA       NA       NA
#> 41      2   m2f1   N1 -0.0810478503 NA       NA       NA
#> 42      2   m2f1   N2 -0.0073626233 NA       NA       NA
#> 43      2   m2f1   N3  0.0061178407 NA       NA       NA
#> 44      2   m2f1   N4 -0.2068825134 NA       NA       NA
#> 45      2   m2f1   N5 -0.0186525984 NA       NA       NA
#> 46      2   m2f1   O1  0.4207181492 NA       NA       NA
#> 47      2   m2f1   O2 -0.0552477704 NA       NA       NA
#> 48      2   m2f1   O3  0.5201614551 NA       NA       NA
#> 49      2   m2f1   O4  0.1184274579 NA       NA       NA
#> 50      2   m2f1   O5 -0.1884883764 NA       NA       NA
#> 51      2   m2f2   A1 -0.0940951697 NA       NA       NA
#> 52      2   m2f2   A2 -0.0072648510 NA       NA       NA
#> 53      2   m2f2   A3 -0.0078860370 NA       NA       NA
#> 54      2   m2f2   A4  0.1374848610 NA       NA       NA
#> 55      2   m2f2   A5  0.1367386516 NA       NA       NA
#> 56      2   m2f2   C1  0.1702904231 NA       NA       NA
#> 57      2   m2f2   C2  0.0707924089 NA       NA       NA
#> 58      2   m2f2   C3  0.1538835145 NA       NA       NA
#> 59      2   m2f2   C4 -0.4252903243 NA       NA       NA
#> 60      2   m2f2   C5 -0.4788597467 NA       NA       NA
#> 61      2   m2f2   E1  0.0098264827 NA       NA       NA
#> 62      2   m2f2   E2 -0.2373345574 NA       NA       NA
#> 63      2   m2f2   E3 -0.0551730838 NA       NA       NA
#> 64      2   m2f2   E4  0.0711582430 NA       NA       NA
#> 65      2   m2f2   E5  0.0699682532 NA       NA       NA
#> 66      2   m2f2   N1 -0.7558899633 NA       NA       NA
#> 67      2   m2f2   N2 -0.7620266161 NA       NA       NA
#> 68      2   m2f2   N3 -0.7798628127 NA       NA       NA
#> 69      2   m2f2   N4 -0.6440660642 NA       NA       NA
#> 70      2   m2f2   N5 -0.6238482902 NA       NA       NA
#> 71      2   m2f2   O1 -0.0330059179 NA       NA       NA
#> 72      2   m2f2   O2 -0.2530424511 NA       NA       NA
#> 73      2   m2f2   O3 -0.0622301881 NA       NA       NA
#> 74      2   m2f2   O4 -0.3114665925 NA       NA       NA
#> 75      2   m2f2   O5 -0.1615608763 NA       NA       NA
#> 76      3   m3f1   A1 -0.3303815400 NA       NA       NA
#> 77      3   m3f1   A2  0.6115763242 NA       NA       NA
#> 78      3   m3f1   A3  0.7183590386 NA       NA       NA
#> 79      3   m3f1   A4  0.4230101783 NA       NA       NA
#> 80      3   m3f1   A5  0.7055300255 NA       NA       NA
#> 81      3   m3f1   C1  0.0812624767 NA       NA       NA
#> 82      3   m3f1   C2  0.1603434517 NA       NA       NA
#> 83      3   m3f1   C3  0.1206567726 NA       NA       NA
#> 84      3   m3f1   C4 -0.0851470599 NA       NA       NA
#> 85      3   m3f1   C5 -0.2032909896 NA       NA       NA
#> 86      3   m3f1   E1 -0.5777625170 NA       NA       NA
#> 87      3   m3f1   E2 -0.6574242246 NA       NA       NA
#> 88      3   m3f1   E3  0.6728364031 NA       NA       NA
#> 89      3   m3f1   E4  0.7349322733 NA       NA       NA
#> 90      3   m3f1   E5  0.4897306828 NA       NA       NA
#> 91      3   m3f1   N1 -0.0964706076 NA       NA       NA
#> 92      3   m3f1   N2 -0.0493075848 NA       NA       NA
#> 93      3   m3f1   N3 -0.0351947964 NA       NA       NA
#> 94      3   m3f1   N4 -0.2229928522 NA       NA       NA
#> 95      3   m3f1   N5 -0.0334413067 NA       NA       NA
#> 96      3   m3f1   O1  0.2989168015 NA       NA       NA
#> 97      3   m3f1   O2  0.1271938761 NA       NA       NA
#> 98      3   m3f1   O3  0.3963740430 NA       NA       NA
#> 99      3   m3f1   O4  0.0159519776 NA       NA       NA
#> 100     3   m3f1   O5 -0.0169511178 NA       NA       NA
#> 101     3   m3f2   A1 -0.1039099166 NA       NA       NA
#> 102     3   m3f2   A2 -0.0290454598 NA       NA       NA
#> 103     3   m3f2   A3 -0.0102562844 NA       NA       NA
#> 104     3   m3f2   A4  0.1153472243 NA       NA       NA
#> 105     3   m3f2   A5  0.1425063110 NA       NA       NA
#> 106     3   m3f2   C1  0.0020184624 NA       NA       NA
#> 107     3   m3f2   C2 -0.1012245093 NA       NA       NA
#> 108     3   m3f2   C3  0.0144123380 NA       NA       NA
#> 109     3   m3f2   C4 -0.2572858865 NA       NA       NA
#> 110     3   m3f2   C5 -0.3556383574 NA       NA       NA
#> 111     3   m3f2   E1 -0.0172032498 NA       NA       NA
#> 112     3   m3f2   E2 -0.2539180175 NA       NA       NA
#> 113     3   m3f2   E3 -0.0665042347 NA       NA       NA
#> 114     3   m3f2   E4  0.1151820137 NA       NA       NA
#> 115     3   m3f2   E5 -0.0218610730 NA       NA       NA
#> 116     3   m3f2   N1 -0.7497924042 NA       NA       NA
#> 117     3   m3f2   N2 -0.7737882151 NA       NA       NA
#> 118     3   m3f2   N3 -0.7913133924 NA       NA       NA
#> 119     3   m3f2   N4 -0.6367982646 NA       NA       NA
#> 120     3   m3f2   N5 -0.6213997894 NA       NA       NA
#> 121     3   m3f2   O1 -0.1172407466 NA       NA       NA
#> 122     3   m3f2   O2 -0.1377080454 NA       NA       NA
#> 123     3   m3f2   O3 -0.1499302915 NA       NA       NA
#> 124     3   m3f2   O4 -0.3708421618 NA       NA       NA
#> 125     3   m3f2   O5 -0.0506242222 NA       NA       NA
#> 126     3   m3f3   A1 -0.0174112818 NA       NA       NA
#> 127     3   m3f3   A2  0.1518303860 NA       NA       NA
#> 128     3   m3f3   A3  0.0944272838 NA       NA       NA
#> 129     3   m3f3   A4  0.1502168591 NA       NA       NA
#> 130     3   m3f3   A5  0.0833234089 NA       NA       NA
#> 131     3   m3f3   C1  0.6421892585 NA       NA       NA
#> 132     3   m3f3   C2  0.6516225810 NA       NA       NA
#> 133     3   m3f3   C3  0.5404779292 NA       NA       NA
#> 134     3   m3f3   C4 -0.6767739810 NA       NA       NA
#> 135     3   m3f3   C5 -0.5363850435 NA       NA       NA
#> 136     3   m3f3   E1  0.0292443610 NA       NA       NA
#> 137     3   m3f3   E2 -0.0522024215 NA       NA       NA
#> 138     3   m3f3   E3  0.1148393306 NA       NA       NA
#> 139     3   m3f3   E4 -0.0605784374 NA       NA       NA
#> 140     3   m3f3   E5  0.4011960574 NA       NA       NA
#> 141     3   m3f3   N1 -0.1377348228 NA       NA       NA
#> 142     3   m3f3   N2 -0.0682444634 NA       NA       NA
#> 143     3   m3f3   N3 -0.0701159278 NA       NA       NA
#> 144     3   m3f3   N4 -0.1418965086 NA       NA       NA
#> 145     3   m3f3   N5 -0.0987360250 NA       NA       NA
#> 146     3   m3f3   O1  0.3364546181 NA       NA       NA
#> 147     3   m3f3   O2 -0.4367816880 NA       NA       NA
#> 148     3   m3f3   O3  0.3567679371 NA       NA       NA
#> 149     3   m3f3   O4  0.1739417937 NA       NA       NA
#> 150     3   m3f3   O5 -0.4257258220 NA       NA       NA
#> 151     4   m4f1   A1 -0.3313259612 NA       NA       NA
#> 152     4   m4f1   A2  0.6039367656 NA       NA       NA
#> 153     4   m4f1   A3  0.7147046980 NA       NA       NA
#> 154     4   m4f1   A4  0.4273539113 NA       NA       NA
#> 155     4   m4f1   A5  0.7021785549 NA       NA       NA
#> 156     4   m4f1   C1  0.0422291268 NA       NA       NA
#> 157     4   m4f1   C2  0.1231215140 NA       NA       NA
#> 158     4   m4f1   C3  0.0958302985 NA       NA       NA
#> 159     4   m4f1   C4 -0.0487939384 NA       NA       NA
#> 160     4   m4f1   C5 -0.1797876827 NA       NA       NA
#> 161     4   m4f1   E1 -0.5697606325 NA       NA       NA
#> 162     4   m4f1   E2 -0.6516925126 NA       NA       NA
#> 163     4   m4f1   E3  0.6532615369 NA       NA       NA
#> 164     4   m4f1   E4  0.7414915202 NA       NA       NA
#> 165     4   m4f1   E5  0.4594144547 NA       NA       NA
#> 166     4   m4f1   N1 -0.0941837756 NA       NA       NA
#> 167     4   m4f1   N2 -0.0553300908 NA       NA       NA
#> 168     4   m4f1   N3 -0.0406377287 NA       NA       NA
#> 169     4   m4f1   N4 -0.2240453960 NA       NA       NA
#> 170     4   m4f1   N5 -0.0274511387 NA       NA       NA
#> 171     4   m4f1   O1  0.2567842820 NA       NA       NA
#> 172     4   m4f1   O2  0.1711555605 NA       NA       NA
#> 173     4   m4f1   O3  0.3497320811 NA       NA       NA
#> 174     4   m4f1   O4 -0.0131031326 NA       NA       NA
#> 175     4   m4f1   O5  0.0321115985 NA       NA       NA
#> 176     4   m4f2   A1 -0.0955211947 NA       NA       NA
#> 177     4   m4f2   A2 -0.0446683588 NA       NA       NA
#> 178     4   m4f2   A3 -0.0303498078 NA       NA       NA
#> 179     4   m4f2   A4  0.0645877449 NA       NA       NA
#> 180     4   m4f2   A5  0.1316413255 NA       NA       NA
#> 181     4   m4f2   C1  0.0167405670 NA       NA       NA
#> 182     4   m4f2   C2 -0.1023277576 NA       NA       NA
#> 183     4   m4f2   C3 -0.0046958003 NA       NA       NA
#> 184     4   m4f2   C4 -0.2673087704 NA       NA       NA
#> 185     4   m4f2   C5 -0.3481249429 NA       NA       NA
#> 186     4   m4f2   E1 -0.0441517470 NA       NA       NA
#> 187     4   m4f2   E2 -0.2654751733 NA       NA       NA
#> 188     4   m4f2   E3 -0.0311162926 NA       NA       NA
#> 189     4   m4f2   E4  0.0942070452 NA       NA       NA
#> 190     4   m4f2   E5  0.0007009622 NA       NA       NA
#> 191     4   m4f2   N1 -0.7661581249 NA       NA       NA
#> 192     4   m4f2   N2 -0.7746854036 NA       NA       NA
#> 193     4   m4f2   N3 -0.7952344769 NA       NA       NA
#> 194     4   m4f2   N4 -0.6311055866 NA       NA       NA
#> 195     4   m4f2   N5 -0.6549982895 NA       NA       NA
#> 196     4   m4f2   O1 -0.0363850621 NA       NA       NA
#> 197     4   m4f2   O2 -0.2257924242 NA       NA       NA
#> 198     4   m4f2   O3 -0.0594071797 NA       NA       NA
#> 199     4   m4f2   O4 -0.3140180240 NA       NA       NA
#> 200     4   m4f2   O5 -0.1535961414 NA       NA       NA
#> 201     4   m4f3   A1 -0.0525499338 NA       NA       NA
#> 202     4   m4f3   A2  0.2327708767 NA       NA       NA
#> 203     4   m4f3   A3  0.1881019714 NA       NA       NA
#> 204     4   m4f3   A4  0.3507576734 NA       NA       NA
#> 205     4   m4f3   A5  0.1387975421 NA       NA       NA
#> 206     4   m4f3   C1  0.6448308647 NA       NA       NA
#> 207     4   m4f3   C2  0.7160906394 NA       NA       NA
#> 208     4   m4f3   C3  0.6569091667 NA       NA       NA
#> 209     4   m4f3   C4 -0.6930818547 NA       NA       NA
#> 210     4   m4f3   C5 -0.6037711831 NA       NA       NA
#> 211     4   m4f3   E1  0.1188657289 NA       NA       NA
#> 212     4   m4f3   E2 -0.0205889366 NA       NA       NA
#> 213     4   m4f3   E3  0.0107474091 NA       NA       NA
#> 214     4   m4f3   E4  0.0203836671 NA       NA       NA
#> 215     4   m4f3   E5  0.3631598043 NA       NA       NA
#> 216     4   m4f3   N1 -0.0736905655 NA       NA       NA
#> 217     4   m4f3   N2 -0.0526571866 NA       NA       NA
#> 218     4   m4f3   N3 -0.0430977492 NA       NA       NA
#> 219     4   m4f3   N4 -0.1628048739 NA       NA       NA
#> 220     4   m4f3   N5  0.0287274682 NA       NA       NA
#> 221     4   m4f3   O1  0.0816460319 NA       NA       NA
#> 222     4   m4f3   O2 -0.1504248925 NA       NA       NA
#> 223     4   m4f3   O3  0.0713838170 NA       NA       NA
#> 224     4   m4f3   O4 -0.0068410908 NA       NA       NA
#> 225     4   m4f3   O5 -0.0894989674 NA       NA       NA
#> 226     4   m4f4   A1  0.0210757701 NA       NA       NA
#> 227     4   m4f4   A2 -0.0071811608 NA       NA       NA
#> 228     4   m4f4   A3 -0.0312817116 NA       NA       NA
#> 229     4   m4f4   A4 -0.2526735344 NA       NA       NA
#> 230     4   m4f4   A5  0.0098017544 NA       NA       NA
#> 231     4   m4f4   C1  0.1773942500 NA       NA       NA
#> 232     4   m4f4   C2  0.0993324605 NA       NA       NA
#> 233     4   m4f4   C3 -0.0372129122 NA       NA       NA
#> 234     4   m4f4   C4 -0.1332430159 NA       NA       NA
#> 235     4   m4f4   C5 -0.0144851966 NA       NA       NA
#> 236     4   m4f4   E1 -0.2214795546 NA       NA       NA
#> 237     4   m4f4   E2 -0.1279460322 NA       NA       NA
#> 238     4   m4f4   E3  0.3075220245 NA       NA       NA
#> 239     4   m4f4   E4 -0.0639678547 NA       NA       NA
#> 240     4   m4f4   E5  0.2415681995 NA       NA       NA
#> 241     4   m4f4   N1 -0.0671574532 NA       NA       NA
#> 242     4   m4f4   N2  0.0423634074 NA       NA       NA
#> 243     4   m4f4   N3  0.0267094752 NA       NA       NA
#> 244     4   m4f4   N4  0.0434481315 NA       NA       NA
#> 245     4   m4f4   N5 -0.1702988944 NA       NA       NA
#> 246     4   m4f4   O1  0.5742195385 NA       NA       NA
#> 247     4   m4f4   O2 -0.5638093117 NA       NA       NA
#> 248     4   m4f4   O3  0.6486197349 NA       NA       NA
#> 249     4   m4f4   O4  0.3976634657 NA       NA       NA
#> 250     4   m4f4   O5 -0.6754313044 NA       NA       NA
#> 251     5   m5f1   A1  0.0646226924 NA       NA       NA
#> 252     5   m5f1   A2  0.2672782045 NA       NA       NA
#> 253     5   m5f1   A3  0.3820508243 NA       NA       NA
#> 254     5   m5f1   A4  0.1363548060 NA       NA       NA
#> 255     5   m5f1   A5  0.4646610567 NA       NA       NA
#> 256     5   m5f1   C1  0.0583256948 NA       NA       NA
#> 257     5   m5f1   C2  0.0431850468 NA       NA       NA
#> 258     5   m5f1   C3  0.0025303501 NA       NA       NA
#> 259     5   m5f1   C4 -0.0160723150 NA       NA       NA
#> 260     5   m5f1   C5 -0.1858360861 NA       NA       NA
#> 261     5   m5f1   E1 -0.6841842145 NA       NA       NA
#> 262     5   m5f1   E2 -0.7215885907 NA       NA       NA
#> 263     5   m5f1   E3  0.6535774970 NA       NA       NA
#> 264     5   m5f1   E4  0.7051005673 NA       NA       NA
#> 265     5   m5f1   E5  0.5747070478 NA       NA       NA
#> 266     5   m5f1   N1  0.0536435676 NA       NA       NA
#> 267     5   m5f1   N2  0.0840702661 NA       NA       NA
#> 268     5   m5f1   N3 -0.0250658759 NA       NA       NA
#> 269     5   m5f1   N4 -0.3018914765 NA       NA       NA
#> 270     5   m5f1   N5 -0.1617027467 NA       NA       NA
#> 271     5   m5f1   O1  0.3013971766 NA       NA       NA
#> 272     5   m5f1   O2  0.0929484238 NA       NA       NA
#> 273     5   m5f1   O3  0.3788252574 NA       NA       NA
#> 274     5   m5f1   O4 -0.1774937648 NA       NA       NA
#> 275     5   m5f1   O5 -0.0328080987 NA       NA       NA
#> 276     5   m5f2   A1 -0.1256711315 NA       NA       NA
#> 277     5   m5f2   A2 -0.0285960200 NA       NA       NA
#> 278     5   m5f2   A3 -0.0162766653 NA       NA       NA
#> 279     5   m5f2   A4  0.0760973347 NA       NA       NA
#> 280     5   m5f2   A5  0.1383137204 NA       NA       NA
#> 281     5   m5f2   C1  0.0051847203 NA       NA       NA
#> 282     5   m5f2   C2 -0.1083356030 NA       NA       NA
#> 283     5   m5f2   C3 -0.0088518694 NA       NA       NA
#> 284     5   m5f2   C4 -0.2594918952 NA       NA       NA
#> 285     5   m5f2   C5 -0.3346991836 NA       NA       NA
#> 286     5   m5f2   E1 -0.0258185253 NA       NA       NA
#> 287     5   m5f2   E2 -0.2470504064 NA       NA       NA
#> 288     5   m5f2   E3 -0.0418378502 NA       NA       NA
#> 289     5   m5f2   E4  0.0823275105 NA       NA       NA
#> 290     5   m5f2   E5 -0.0233822649 NA       NA       NA
#> 291     5   m5f2   N1 -0.7788020049 NA       NA       NA
#> 292     5   m5f2   N2 -0.7867916446 NA       NA       NA
#> 293     5   m5f2   N3 -0.7961909316 NA       NA       NA
#> 294     5   m5f2   N4 -0.6172249847 NA       NA       NA
#> 295     5   m5f2   N5 -0.6445215425 NA       NA       NA
#> 296     5   m5f2   O1 -0.0425566548 NA       NA       NA
#> 297     5   m5f2   O2 -0.2240602326 NA       NA       NA
#> 298     5   m5f2   O3 -0.0652192570 NA       NA       NA
#> 299     5   m5f2   O4 -0.2953825855 NA       NA       NA
#> 300     5   m5f2   O5 -0.1521555346 NA       NA       NA
#> 301     5   m5f3   A1  0.0465825257 NA       NA       NA
#> 302     5   m5f3   A2  0.1721850162 NA       NA       NA
#> 303     5   m5f3   A3  0.1362534957 NA       NA       NA
#> 304     5   m5f3   A4  0.3037543655 NA       NA       NA
#> 305     5   m5f3   A5  0.1144396535 NA       NA       NA
#> 306     5   m5f3   C1  0.6581766639 NA       NA       NA
#> 307     5   m5f3   C2  0.7081232447 NA       NA       NA
#> 308     5   m5f3   C3  0.6471358944 NA       NA       NA
#> 309     5   m5f3   C4 -0.6954819104 NA       NA       NA
#> 310     5   m5f3   C5 -0.6299162649 NA       NA       NA
#> 311     5   m5f3   E1  0.0546615449 NA       NA       NA
#> 312     5   m5f3   E2 -0.0833393990 NA       NA       NA
#> 313     5   m5f3   E3  0.0438759226 NA       NA       NA
#> 314     5   m5f3   E4  0.0615241916 NA       NA       NA
#> 315     5   m5f3   E5  0.4279181649 NA       NA       NA
#> 316     5   m5f3   N1 -0.0373994780 NA       NA       NA
#> 317     5   m5f3   N2 -0.0198162758 NA       NA       NA
#> 318     5   m5f3   N3 -0.0472910889 NA       NA       NA
#> 319     5   m5f3   N4 -0.2100690587 NA       NA       NA
#> 320     5   m5f3   N5 -0.0124671510 NA       NA       NA
#> 321     5   m5f3   O1  0.0943080836 NA       NA       NA
#> 322     5   m5f3   O2 -0.1481708602 NA       NA       NA
#> 323     5   m5f3   O3  0.0825441466 NA       NA       NA
#> 324     5   m5f3   O4 -0.0744790766 NA       NA       NA
#> 325     5   m5f3   O5 -0.0870502848 NA       NA       NA
#> 326     5   m5f4   A1 -0.6608493873 NA       NA       NA
#> 327     5   m5f4   A2  0.6619262832 NA       NA       NA
#> 328     5   m5f4   A3  0.6899489715 NA       NA       NA
#> 329     5   m5f4   A4  0.5197398678 NA       NA       NA
#> 330     5   m5f4   A5  0.5565035029 NA       NA       NA
#> 331     5   m5f4   C1 -0.0028200627 NA       NA       NA
#> 332     5   m5f4   C2  0.1502576422 NA       NA       NA
#> 333     5   m5f4   C3  0.1500251371 NA       NA       NA
#> 334     5   m5f4   C4 -0.0644977888 NA       NA       NA
#> 335     5   m5f4   C5 -0.0427313221 NA       NA       NA
#> 336     5   m5f4   E1 -0.0394679906 NA       NA       NA
#> 337     5   m5f4   E2 -0.1191348980 NA       NA       NA
#> 338     5   m5f4   E3  0.2320593807 NA       NA       NA
#> 339     5   m5f4   E4  0.2806622823 NA       NA       NA
#> 340     5   m5f4   E5 -0.0014964851 NA       NA       NA
#> 341     5   m5f4   N1 -0.2452915782 NA       NA       NA
#> 342     5   m5f4   N2 -0.2121788489 NA       NA       NA
#> 343     5   m5f4   N3 -0.0344889973 NA       NA       NA
#> 344     5   m5f4   N4  0.0426402534 NA       NA       NA
#> 345     5   m5f4   N5  0.1637513829 NA       NA       NA
#> 346     5   m5f4   O1  0.0656490791 NA       NA       NA
#> 347     5   m5f4   O2  0.1181424307 NA       NA       NA
#> 348     5   m5f4   O3  0.1232896124 NA       NA       NA
#> 349     5   m5f4   O4  0.2602278049 NA       NA       NA
#> 350     5   m5f4   O5  0.0454013494 NA       NA       NA
#> 351     5   m5f5   A1 -0.1205757842 NA       NA       NA
#> 352     5   m5f5   A2  0.0618886113 NA       NA       NA
#> 353     5   m5f5   A3  0.0205631842 NA       NA       NA
#> 354     5   m5f5   A4 -0.1966166381 NA       NA       NA
#> 355     5   m5f5   A5  0.0233730897 NA       NA       NA
#> 356     5   m5f5   C1  0.1546256305 NA       NA       NA
#> 357     5   m5f5   C2  0.1035115928 NA       NA       NA
#> 358     5   m5f5   C3 -0.0275059206 NA       NA       NA
#> 359     5   m5f5   C4 -0.1318524854 NA       NA       NA
#> 360     5   m5f5   C5  0.0267485097 NA       NA       NA
#> 361     5   m5f5   E1 -0.1033688412 NA       NA       NA
#> 362     5   m5f5   E2 -0.0128878820 NA       NA       NA
#> 363     5   m5f5   E3  0.2324427219 NA       NA       NA
#> 364     5   m5f5   E4 -0.1531437588 NA       NA       NA
#> 365     5   m5f5   E5  0.1247753442 NA       NA       NA
#> 366     5   m5f5   N1 -0.1348113658 NA       NA       NA
#> 367     5   m5f5   N2 -0.0213668932 NA       NA       NA
#> 368     5   m5f5   N3  0.0188382240 NA       NA       NA
#> 369     5   m5f5   N4  0.1122029266 NA       NA       NA
#> 370     5   m5f5   N5 -0.1197031854 NA       NA       NA
#> 371     5   m5f5   O1  0.5458112825 NA       NA       NA
#> 372     5   m5f5   O2 -0.5793620831 NA       NA       NA
#> 373     5   m5f5   O3  0.6188115571 NA       NA       NA
#> 374     5   m5f5   O4  0.4968108232 NA       NA       NA
#> 375     5   m5f5   O5 -0.6850859719 NA       NA       NA
tidy(x, what = "variance")
#>    level factor proportion cumulative
#> 1      1   m1f1 0.20609551  0.2060955
#> 2      2   m2f1 0.18456586  0.1845659
#> 3      2   m2f2 0.13360875  0.3181746
#> 4      3   m3f1 0.16388269  0.1638827
#> 5      3   m3f2 0.12421872  0.2881014
#> 6      3   m3f3 0.11276933  0.4008707
#> 7      4   m4f1 0.15768276  0.1576828
#> 8      4   m4f2 0.12576460  0.2834474
#> 9      4   m4f3 0.10665811  0.3901055
#> 10     4   m4f4 0.08276426  0.4728697
#> 11     5   m5f1 0.12571964  0.1257196
#> 12     5   m5f2 0.12505250  0.2507721
#> 13     5   m5f3 0.10750049  0.3582726
#> 14     5   m5f4 0.09473744  0.4530101
#> 15     5   m5f5 0.07975781  0.5327679
tidy(x, what = "fit")
#>    level       statistic    value
#> 1      1 eigenvalue.m1f1 5.152388
#> 2      2 eigenvalue.m2f1 5.152388
#> 3      2 eigenvalue.m2f2 2.801977
#> 4      3 eigenvalue.m3f1 5.152388
#> 5      3 eigenvalue.m3f2 2.801977
#> 6      3 eigenvalue.m3f3 2.067404
#> 7      4 eigenvalue.m4f1 5.152388
#> 8      4 eigenvalue.m4f2 2.801977
#> 9      4 eigenvalue.m4f3 2.067404
#> 10     4 eigenvalue.m4f4 1.799974
#> 11     5 eigenvalue.m5f1 5.152388
#> 12     5 eigenvalue.m5f2 2.801977
#> 13     5 eigenvalue.m5f3 2.067404
#> 14     5 eigenvalue.m5f4 1.799974
#> 15     5 eigenvalue.m5f5 1.497454
tidy(x, what = "fit", format = "wide")
#>   level eigenvalue.m2f1 eigenvalue.m2f2 eigenvalue.m3f1 eigenvalue.m3f2
#> 1     2        5.152388        2.801977              NA              NA
#> 2     3              NA              NA        5.152388        2.801977
#> 3     4              NA              NA              NA              NA
#> 4     5              NA              NA              NA              NA
#>   eigenvalue.m3f3 eigenvalue.m4f1 eigenvalue.m4f2 eigenvalue.m4f3
#> 1              NA              NA              NA              NA
#> 2        2.067404              NA              NA              NA
#> 3              NA        5.152388        2.801977        2.067404
#> 4              NA              NA              NA              NA
#>   eigenvalue.m4f4 eigenvalue.m5f1 eigenvalue.m5f2 eigenvalue.m5f3
#> 1              NA              NA              NA              NA
#> 2              NA              NA              NA              NA
#> 3        1.799974              NA              NA              NA
#> 4              NA        5.152388        2.801977        2.067404
#>   eigenvalue.m5f4 eigenvalue.m5f5
#> 1              NA              NA
#> 2              NA              NA
#> 3              NA              NA
#> 4        1.799974        1.497454
```
