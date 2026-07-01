# Display the salient items for each factor

Returns, per level and factor, the items whose absolute loading meets or
exceeds `cut`, sorted by descending absolute loading. This gives a
concise reading of "what each factor is about" without printing a full
item-by-factor matrix, which does not scale well to large `k` or many
items.

## Usage

``` r
top_items(
  x,
  level = NULL,
  cut = 0.3,
  n = NULL,
  sort = TRUE,
  by = c("factor", "item"),
  show_labels = TRUE
)
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

  Logical. When `TRUE` (default), items within each group are ordered by
  descending `|loading|`. Set to `FALSE` to keep the original order
  (useful when items have a meaningful sequence).

- by:

  One of `"factor"` (default) or `"item"`. `"factor"` groups the listing
  by factor (the salient items *of* each factor – "what is this factor
  about?"). `"item"` inverts the grouping to list, for each item, the
  factors it loads on – which makes cross-loadings legible ("where does
  this item go?"). `n` and `sort` apply within whichever unit `by`
  selects.

- show_labels:

  Logical. When `TRUE` (default) and the data carried a variable-label
  attribute at fit time (see
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)),
  items are shown as `label (id)`; items without a label fall back to
  the bare id. Set to `FALSE` to always show the bare `m{k}f{j}`-style
  item ids.

## Value

An object of class `"top_items"`. Print it for a grouped cli listing.
The underlying data frame (one row per selected item) is accessible via
`$data` and contains columns `level`, `factor`, `item`, and `loading`
(plus `label` when labels are available). The values equal the
corresponding `tidy(x, what = "loadings")` rows (after filtering and
optional sorting).

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
x <- ackwards(bfi25, k_max = 5)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
top_items(x)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.3
#> Top-n: all
#> 
#> ── Level 1 (1 factor) ──
#> 
#> m1f1
#> A5 [0.646]
#> E2 [-0.636]
#> A3 [0.602]
#> E4 [0.587]
#> E5 [0.571]
#> C5 [-0.560]
#> E3 [0.555]
#> A2 [0.540]
#> C4 [-0.500]
#> N4 [-0.490]
#> A4 [0.448]
#> E1 [-0.442]
#> N1 [-0.433]
#> O3 [0.427]
#> C2 [0.407]
#> C1 [0.381]
#> N2 [-0.371]
#> N3 [-0.368]
#> C3 [0.366]
#> O1 [0.354]
#> N5 [-0.315]
#> A1 [-0.306]
#> 
#> ── Level 2 (2 factors) ──
#> 
#> m2f1
#> A3 [0.690]
#> E3 [0.662]
#> A5 [0.661]
#> E4 [0.630]
#> A2 [0.619]
#> E5 [0.612]
#> E2 [-0.595]
#> O3 [0.520]
#> E1 [-0.509]
#> A4 [0.436]
#> C2 [0.425]
#> O1 [0.421]
#> C5 [-0.377]
#> C1 [0.341]
#> C4 [-0.337]
#> C3 [0.333]
#> m2f2
#> N3 [-0.780]
#> N2 [-0.762]
#> N1 [-0.756]
#> N4 [-0.644]
#> N5 [-0.624]
#> C5 [-0.479]
#> C4 [-0.425]
#> O4 [-0.311]
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> E4 [0.735]
#> A3 [0.718]
#> A5 [0.706]
#> E3 [0.673]
#> E2 [-0.657]
#> A2 [0.612]
#> E1 [-0.578]
#> E5 [0.490]
#> A4 [0.423]
#> O3 [0.396]
#> A1 [-0.330]
#> m3f2
#> N3 [-0.791]
#> N2 [-0.774]
#> N1 [-0.750]
#> N4 [-0.637]
#> N5 [-0.621]
#> O4 [-0.371]
#> C5 [-0.356]
#> m3f3
#> C4 [-0.677]
#> C2 [0.652]
#> C1 [0.642]
#> C3 [0.540]
#> C5 [-0.536]
#> O2 [-0.437]
#> O5 [-0.426]
#> E5 [0.401]
#> O3 [0.357]
#> O1 [0.336]
#> 
#> ── Level 4 (4 factors) ──
#> 
#> m4f1
#> E4 [0.741]
#> A3 [0.715]
#> A5 [0.702]
#> E3 [0.653]
#> E2 [-0.652]
#> A2 [0.604]
#> E1 [-0.570]
#> E5 [0.459]
#> A4 [0.427]
#> O3 [0.350]
#> A1 [-0.331]
#> m4f2
#> N3 [-0.795]
#> N2 [-0.775]
#> N1 [-0.766]
#> N5 [-0.655]
#> N4 [-0.631]
#> C5 [-0.348]
#> O4 [-0.314]
#> m4f3
#> C2 [0.716]
#> C4 [-0.693]
#> C3 [0.657]
#> C1 [0.645]
#> C5 [-0.604]
#> E5 [0.363]
#> A4 [0.351]
#> m4f4
#> O5 [-0.675]
#> O3 [0.649]
#> O1 [0.574]
#> O2 [-0.564]
#> O4 [0.398]
#> E3 [0.308]
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> E2 [-0.722]
#> E4 [0.705]
#> E1 [-0.684]
#> E3 [0.654]
#> E5 [0.575]
#> A5 [0.465]
#> A3 [0.382]
#> O3 [0.379]
#> N4 [-0.302]
#> O1 [0.301]
#> m5f2
#> N3 [-0.796]
#> N2 [-0.787]
#> N1 [-0.779]
#> N5 [-0.645]
#> N4 [-0.617]
#> C5 [-0.335]
#> m5f3
#> C2 [0.708]
#> C4 [-0.695]
#> C1 [0.658]
#> C3 [0.647]
#> C5 [-0.630]
#> E5 [0.428]
#> A4 [0.304]
#> m5f4
#> A3 [0.690]
#> A2 [0.662]
#> A1 [-0.661]
#> A5 [0.557]
#> A4 [0.520]
#> m5f5
#> O5 [-0.685]
#> O3 [0.619]
#> O2 [-0.579]
#> O1 [0.546]
#> O4 [0.497]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
top_items(x, level = 5, cut = 0.4, n = 5)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.4
#> Top-n: 5
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> E2 [-0.722]
#> E4 [0.705]
#> E1 [-0.684]
#> E3 [0.654]
#> E5 [0.575]
#> m5f2
#> N3 [-0.796]
#> N2 [-0.787]
#> N1 [-0.779]
#> N5 [-0.645]
#> N4 [-0.617]
#> m5f3
#> C2 [0.708]
#> C4 [-0.695]
#> C1 [0.658]
#> C3 [0.647]
#> C5 [-0.630]
#> m5f4
#> A3 [0.690]
#> A2 [0.662]
#> A1 [-0.661]
#> A5 [0.557]
#> A4 [0.520]
#> m5f5
#> O5 [-0.685]
#> O3 [0.619]
#> O2 [-0.579]
#> O1 [0.546]
#> O4 [0.497]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.

# Invert the grouping to read cross-loadings item-by-item
top_items(x, level = 3, cut = 0.25, by = "item")
#> 
#> ── Salient factors by item (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.25
#> Top-n: all
#> 
#> ── Level 3 (3 factors) ──
#> 
#> A1
#> m3f1 [-0.330]
#> A2
#> m3f1 [0.612]
#> A3
#> m3f1 [0.718]
#> A4
#> m3f1 [0.423]
#> A5
#> m3f1 [0.706]
#> C1
#> m3f3 [0.642]
#> C2
#> m3f3 [0.652]
#> C3
#> m3f3 [0.540]
#> C4
#> m3f3 [-0.677]
#> m3f2 [-0.257]
#> C5
#> m3f3 [-0.536]
#> m3f2 [-0.356]
#> E1
#> m3f1 [-0.578]
#> E2
#> m3f1 [-0.657]
#> m3f2 [-0.254]
#> E3
#> m3f1 [0.673]
#> E4
#> m3f1 [0.735]
#> E5
#> m3f1 [0.490]
#> m3f3 [0.401]
#> N1
#> m3f2 [-0.750]
#> N2
#> m3f2 [-0.774]
#> N3
#> m3f2 [-0.791]
#> N4
#> m3f2 [-0.637]
#> N5
#> m3f2 [-0.621]
#> O1
#> m3f3 [0.336]
#> m3f1 [0.299]
#> O2
#> m3f3 [-0.437]
#> O3
#> m3f1 [0.396]
#> m3f3 [0.357]
#> O4
#> m3f2 [-0.371]
#> O5
#> m3f3 [-0.426]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```
