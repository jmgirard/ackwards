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
  items are shown as `id: label`; items without a label fall back to the
  bare id. Set to `FALSE` to always show the bare `m{k}f{j}`-style item
  ids.

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

If [factor
labels](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md)
have been attached, the factor dimension is shown as `label (id)`
wherever it appears – the group headers under `by = "factor"` and the
body entries under `by = "item"`.

## See also

[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md),
[`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md),
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)

## Examples

``` r
# Fit the raw dataset (not na.omit(), which would drop the column
# attributes): bfi25's IPIP item labels are then captured and printed as
# `code: label`. `missing = "listwise"` handles the NAs cleanly. A 10-item
# subset keeps the example fast; use all items in practice.
x <- ackwards(bfi25[, 1:10], k_max = 3, cor = "polychoric", missing = "listwise")
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
#> A3: Know how to comfort others [0.683]
#> C4: Do things in a half-way manner [-0.652]
#> A2: Inquire about others' well-being [0.648]
#> A5: Make people feel at ease [0.634]
#> C5: Waste my time [-0.628]
#> C2: Continue until everything is perfect [0.614]
#> A4: Love children [0.593]
#> C3: Do things according to a plan [0.572]
#> C1: Am exacting in my work [0.554]
#> A1: Am indifferent to the feelings of others [-0.407]
#> 
#> ── Level 2 (2 factors) ──
#> 
#> m2f1
#> C4: Do things in a half-way manner [-0.773]
#> C2: Continue until everything is perfect [0.738]
#> C1: Am exacting in my work [0.734]
#> C5: Waste my time [-0.699]
#> C3: Do things according to a plan [0.657]
#> 
#> m2f2
#> A3: Know how to comfort others [0.846]
#> A2: Inquire about others' well-being [0.786]
#> A5: Make people feel at ease [0.737]
#> A1: Am indifferent to the feelings of others [-0.626]
#> A4: Love children [0.567]
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> C4: Do things in a half-way manner [-0.801]
#> C1: Am exacting in my work [0.740]
#> C5: Waste my time [-0.711]
#> C2: Continue until everything is perfect [0.701]
#> C3: Do things according to a plan [0.636]
#> 
#> m3f2
#> A3: Know how to comfort others [0.812]
#> A5: Make people feel at ease [0.762]
#> A4: Love children [0.728]
#> A2: Inquire about others' well-being [0.644]
#> 
#> m3f3
#> A1: Am indifferent to the feelings of others [-0.876]
#> A2: Inquire about others' well-being [0.464]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
top_items(x, level = 3, cut = 0.4, n = 5)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.4
#> Top-n: 5
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> C4: Do things in a half-way manner [-0.801]
#> C1: Am exacting in my work [0.740]
#> C5: Waste my time [-0.711]
#> C2: Continue until everything is perfect [0.701]
#> C3: Do things according to a plan [0.636]
#> 
#> m3f2
#> A3: Know how to comfort others [0.812]
#> A5: Make people feel at ease [0.762]
#> A4: Love children [0.728]
#> A2: Inquire about others' well-being [0.644]
#> 
#> m3f3
#> A1: Am indifferent to the feelings of others [-0.876]
#> A2: Inquire about others' well-being [0.464]
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
#> A1: Am indifferent to the feelings of others
#> m3f3 [-0.876]
#> 
#> A2: Inquire about others' well-being
#> m3f2 [0.644]
#> m3f3 [0.464]
#> 
#> A3: Know how to comfort others
#> m3f2 [0.812]
#> m3f3 [0.287]
#> 
#> A4: Love children
#> m3f2 [0.728]
#> 
#> A5: Make people feel at ease
#> m3f2 [0.762]
#> 
#> C1: Am exacting in my work
#> m3f1 [0.740]
#> 
#> C2: Continue until everything is perfect
#> m3f1 [0.701]
#> m3f2 [0.269]
#> 
#> C3: Do things according to a plan
#> m3f1 [0.636]
#> 
#> C4: Do things in a half-way manner
#> m3f1 [-0.801]
#> 
#> C5: Waste my time
#> m3f1 [-0.711]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```
