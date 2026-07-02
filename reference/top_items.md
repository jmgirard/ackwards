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
# Fit the raw dataset (not na.omit(), which would drop the column
# attributes): bfi25's IPIP item labels are then captured and printed as
# `label (code)`. `missing = "listwise"` handles the NAs cleanly.
x <- ackwards(bfi25, k_max = 5, cor = "polychoric", missing = "listwise")
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
#> Make people feel at ease (A5) [0.685]
#> Find it difficult to approach others (E2) [-0.656]
#> Know how to comfort others (A3) [0.651]
#> Make friends easily (E4) [0.607]
#> Inquire about others' well-being (A2) [0.598]
#> Take charge (E5) [0.598]
#> Waste my time (C5) [-0.590]
#> Know how to captivate people (E3) [0.582]
#> Do things in a half-way manner (C4) [-0.546]
#> Often feel blue (N4) [-0.489]
#> Love children (A4) [0.484]
#> Don't talk a lot (E1) [-0.460]
#> Carry the conversation to a higher level (O3) [0.457]
#> Am exacting in my work (C1) [0.450]
#> Get angry easily (N1) [-0.444]
#> Continue until everything is perfect (C2) [0.443]
#> Do things according to a plan (C3) [0.411]
#> Am full of ideas (O1) [0.391]
#> Get irritated easily (N2) [-0.366]
#> Have frequent mood swings (N3) [-0.364]
#> Am indifferent to the feelings of others (A1) [-0.344]
#> Panic easily (N5) [-0.319]
#> 
#> ── Level 2 (2 factors) ──
#> 
#> m2f1
#> Know how to comfort others (A3) [0.730]
#> Make people feel at ease (A5) [0.694]
#> Know how to captivate people (E3) [0.684]
#> Inquire about others' well-being (A2) [0.661]
#> Make friends easily (E4) [0.638]
#> Take charge (E5) [0.630]
#> Find it difficult to approach others (E2) [-0.614]
#> Carry the conversation to a higher level (O3) [0.565]
#> Don't talk a lot (E1) [-0.525]
#> Love children (A4) [0.471]
#> Continue until everything is perfect (C2) [0.468]
#> Am full of ideas (O1) [0.467]
#> Am exacting in my work (C1) [0.421]
#> Waste my time (C5) [-0.411]
#> Do things in a half-way manner (C4) [-0.389]
#> Do things according to a plan (C3) [0.365]
#> Am indifferent to the feelings of others (A1) [-0.335]
#> 
#> m2f2
#> Have frequent mood swings (N3) [-0.806]
#> Get angry easily (N1) [-0.790]
#> Get irritated easily (N2) [-0.789]
#> Often feel blue (N4) [-0.677]
#> Panic easily (N5) [-0.659]
#> Waste my time (C5) [-0.491]
#> Do things in a half-way manner (C4) [-0.438]
#> Spend time reflecting on things (O4) [-0.358]
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> Make friends easily (E4) [0.784]
#> Know how to comfort others (A3) [0.750]
#> Make people feel at ease (A5) [0.735]
#> Find it difficult to approach others (E2) [-0.698]
#> Know how to captivate people (E3) [0.664]
#> Inquire about others' well-being (A2) [0.650]
#> Don't talk a lot (E1) [-0.611]
#> Take charge (E5) [0.507]
#> Love children (A4) [0.490]
#> Carry the conversation to a higher level (O3) [0.355]
#> Am indifferent to the feelings of others (A1) [-0.354]
#> 
#> m3f2
#> Have frequent mood swings (N3) [-0.812]
#> Get irritated easily (N2) [-0.791]
#> Get angry easily (N1) [-0.777]
#> Often feel blue (N4) [-0.680]
#> Panic easily (N5) [-0.650]
#> Waste my time (C5) [-0.428]
#> Spend time reflecting on things (O4) [-0.414]
#> Do things in a half-way manner (C4) [-0.348]
#> 
#> m3f3
#> Am exacting in my work (C1) [0.658]
#> Continue until everything is perfect (C2) [0.642]
#> Do things in a half-way manner (C4) [-0.609]
#> Will not probe deeply into a subject (O5) [-0.553]
#> Avoid difficult reading material (O2) [-0.520]
#> Carry the conversation to a higher level (O3) [0.504]
#> Am full of ideas (O1) [0.493]
#> Do things according to a plan (C3) [0.476]
#> Waste my time (C5) [-0.453]
#> Take charge (E5) [0.392]
#> Spend time reflecting on things (O4) [0.319]
#> 
#> ── Level 4 (4 factors) ──
#> 
#> m4f1
#> Make friends easily (E4) [0.783]
#> Know how to comfort others (A3) [0.732]
#> Make people feel at ease (A5) [0.725]
#> Find it difficult to approach others (E2) [-0.709]
#> Know how to captivate people (E3) [0.681]
#> Don't talk a lot (E1) [-0.633]
#> Inquire about others' well-being (A2) [0.626]
#> Take charge (E5) [0.485]
#> Love children (A4) [0.442]
#> Carry the conversation to a higher level (O3) [0.384]
#> Am indifferent to the feelings of others (A1) [-0.354]
#> 
#> m4f2
#> Have frequent mood swings (N3) [-0.824]
#> Get angry easily (N1) [-0.798]
#> Get irritated easily (N2) [-0.793]
#> Panic easily (N5) [-0.700]
#> Often feel blue (N4) [-0.659]
#> Waste my time (C5) [-0.357]
#> Spend time reflecting on things (O4) [-0.332]
#> 
#> m4f3
#> Continue until everything is perfect (C2) [0.735]
#> Do things in a half-way manner (C4) [-0.714]
#> Do things according to a plan (C3) [0.688]
#> Am exacting in my work (C1) [0.680]
#> Waste my time (C5) [-0.626]
#> Love children (A4) [0.397]
#> Take charge (E5) [0.372]
#> 
#> m4f4
#> Will not probe deeply into a subject (O5) [-0.704]
#> Carry the conversation to a higher level (O3) [0.674]
#> Am full of ideas (O1) [0.621]
#> Avoid difficult reading material (O2) [-0.586]
#> Spend time reflecting on things (O4) [0.499]
#> Know how to captivate people (E3) [0.301]
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> Find it difficult to approach others (E2) [-0.752]
#> Make friends easily (E4) [0.747]
#> Don't talk a lot (E1) [-0.701]
#> Know how to captivate people (E3) [0.677]
#> Take charge (E5) [0.597]
#> Make people feel at ease (A5) [0.487]
#> Know how to comfort others (A3) [0.415]
#> Carry the conversation to a higher level (O3) [0.394]
#> Often feel blue (N4) [-0.314]
#> Am full of ideas (O1) [0.302]
#> 
#> m5f2
#> Have frequent mood swings (N3) [-0.825]
#> Get angry easily (N1) [-0.810]
#> Get irritated easily (N2) [-0.805]
#> Panic easily (N5) [-0.688]
#> Often feel blue (N4) [-0.646]
#> Waste my time (C5) [-0.344]
#> Spend time reflecting on things (O4) [-0.317]
#> 
#> m5f3
#> Continue until everything is perfect (C2) [0.735]
#> Do things in a half-way manner (C4) [-0.716]
#> Am exacting in my work (C1) [0.690]
#> Do things according to a plan (C3) [0.679]
#> Waste my time (C5) [-0.652]
#> Take charge (E5) [0.448]
#> Love children (A4) [0.345]
#> 
#> m5f4
#> Am indifferent to the feelings of others (A1) [-0.704]
#> Know how to comfort others (A3) [0.703]
#> Inquire about others' well-being (A2) [0.692]
#> Make people feel at ease (A5) [0.580]
#> Love children (A4) [0.522]
#> 
#> m5f5
#> Will not probe deeply into a subject (O5) [-0.705]
#> Carry the conversation to a higher level (O3) [0.655]
#> Am full of ideas (O1) [0.604]
#> Avoid difficult reading material (O2) [-0.595]
#> Spend time reflecting on things (O4) [0.551]
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
#> Find it difficult to approach others (E2) [-0.752]
#> Make friends easily (E4) [0.747]
#> Don't talk a lot (E1) [-0.701]
#> Know how to captivate people (E3) [0.677]
#> Take charge (E5) [0.597]
#> 
#> m5f2
#> Have frequent mood swings (N3) [-0.825]
#> Get angry easily (N1) [-0.810]
#> Get irritated easily (N2) [-0.805]
#> Panic easily (N5) [-0.688]
#> Often feel blue (N4) [-0.646]
#> 
#> m5f3
#> Continue until everything is perfect (C2) [0.735]
#> Do things in a half-way manner (C4) [-0.716]
#> Am exacting in my work (C1) [0.690]
#> Do things according to a plan (C3) [0.679]
#> Waste my time (C5) [-0.652]
#> 
#> m5f4
#> Am indifferent to the feelings of others (A1) [-0.704]
#> Know how to comfort others (A3) [0.703]
#> Inquire about others' well-being (A2) [0.692]
#> Make people feel at ease (A5) [0.580]
#> Love children (A4) [0.522]
#> 
#> m5f5
#> Will not probe deeply into a subject (O5) [-0.705]
#> Carry the conversation to a higher level (O3) [0.655]
#> Am full of ideas (O1) [0.604]
#> Avoid difficult reading material (O2) [-0.595]
#> Spend time reflecting on things (O4) [0.551]
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
#> Am indifferent to the feelings of others (A1)
#> m3f1 [-0.354]
#> 
#> Inquire about others' well-being (A2)
#> m3f1 [0.650]
#> 
#> Know how to comfort others (A3)
#> m3f1 [0.750]
#> 
#> Love children (A4)
#> m3f1 [0.490]
#> 
#> Make people feel at ease (A5)
#> m3f1 [0.735]
#> 
#> Am exacting in my work (C1)
#> m3f3 [0.658]
#> 
#> Continue until everything is perfect (C2)
#> m3f3 [0.642]
#> 
#> Do things according to a plan (C3)
#> m3f3 [0.476]
#> 
#> Do things in a half-way manner (C4)
#> m3f3 [-0.609]
#> m3f2 [-0.348]
#> 
#> Waste my time (C5)
#> m3f3 [-0.453]
#> m3f2 [-0.428]
#> m3f1 [-0.252]
#> 
#> Don't talk a lot (E1)
#> m3f1 [-0.611]
#> 
#> Find it difficult to approach others (E2)
#> m3f1 [-0.698]
#> 
#> Know how to captivate people (E3)
#> m3f1 [0.664]
#> 
#> Make friends easily (E4)
#> m3f1 [0.784]
#> 
#> Take charge (E5)
#> m3f1 [0.507]
#> m3f3 [0.392]
#> 
#> Get angry easily (N1)
#> m3f2 [-0.777]
#> 
#> Get irritated easily (N2)
#> m3f2 [-0.791]
#> 
#> Have frequent mood swings (N3)
#> m3f2 [-0.812]
#> 
#> Often feel blue (N4)
#> m3f2 [-0.680]
#> 
#> Panic easily (N5)
#> m3f2 [-0.650]
#> 
#> Am full of ideas (O1)
#> m3f3 [0.493]
#> m3f1 [0.253]
#> 
#> Avoid difficult reading material (O2)
#> m3f3 [-0.520]
#> 
#> Carry the conversation to a higher level (O3)
#> m3f3 [0.504]
#> m3f1 [0.355]
#> 
#> Spend time reflecting on things (O4)
#> m3f2 [-0.414]
#> m3f3 [0.319]
#> 
#> Will not probe deeply into a subject (O5)
#> m3f3 [-0.553]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```
