# Display the salient items for each factor

Returns, per level and factor, the items whose absolute loading meets or
exceeds `cut`, sorted by descending absolute loading. This gives a
concise reading of "what each factor is about" without printing a full
item-by-factor matrix, which does not scale well to large `k` or many
items.

## Usage

``` r
top_items(x, level = NULL, cut = 0.3, n = NULL, sort = TRUE)
```

## Arguments

- x:

  An `ackwards` object.

- level:

  Integer vector selecting which level(s) to include. `NULL` (default)
  returns all levels.

- cut:

  Absolute-loading threshold. Items with `|loading| >= cut` are shown.
  Default `0.3`, matching the `cut_show` plotting default.

- n:

  Maximum number of items to show per factor. `NULL` (default) shows all
  items meeting the cut. When set, the top-`n` items by `|loading|` are
  kept after applying the cut.

- sort:

  Logical. When `TRUE` (default), items within each factor are ordered
  by descending `|loading|`. Set to `FALSE` to keep the original item
  order (useful when items have a meaningful sequence).

## Value

An object of class `"top_items"`. Print it for a grouped level-by-factor
cli listing. The underlying data frame (one row per selected item) is
accessible via `$data` and contains columns `level`, `factor`, `item`,
and `loading`. The values equal the corresponding
`tidy(x, what = "loadings")` rows (after filtering and optional
sorting).

## Details

Loadings are signed and reflect the object's primary-parent sign
alignment (see
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)).
Items near the cut threshold may appear for one sign orientation but not
the other; this is expected and informative.

## See also

[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md),
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)

## Examples

``` r
if (requireNamespace("psych", quietly = TRUE)) {
  x <- ackwards(psych::bfi[, 1:25], k_max = 5)
  top_items(x)
  top_items(x, level = 5, cut = 0.4, n = 5)
}
#> Warning: ! 364 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.4
#> Top-n: 5
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> N1 [0.793]
#> N3 [0.793]
#> N2 [0.788]
#> N4 [0.646]
#> N5 [0.633]
#> m5f2
#> E2 [-0.722]
#> E4 [0.704]
#> E1 [-0.677]
#> E3 [0.632]
#> E5 [0.582]
#> m5f3
#> C2 [0.738]
#> C3 [0.674]
#> C4 [-0.667]
#> C1 [0.651]
#> C5 [-0.617]
#> m5f4
#> A2 [0.711]
#> A3 [0.677]
#> A1 [-0.627]
#> A5 [0.564]
#> A4 [0.525]
#> m5f5
#> O5 [-0.681]
#> O3 [0.636]
#> O2 [-0.594]
#> O1 [0.588]
#> O4 [0.499]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```
