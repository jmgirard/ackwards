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
# `code: label`. `missing = "listwise"` handles the NAs cleanly.
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
#> A5: Make people feel at ease [0.685]
#> E2: Find it difficult to approach others [-0.656]
#> A3: Know how to comfort others [0.651]
#> E4: Make friends easily [0.607]
#> A2: Inquire about others' well-being [0.598]
#> E5: Take charge [0.598]
#> C5: Waste my time [-0.590]
#> E3: Know how to captivate people [0.582]
#> C4: Do things in a half-way manner [-0.546]
#> N4: Often feel blue [-0.489]
#> A4: Love children [0.484]
#> E1: Don't talk a lot [-0.460]
#> O3: Carry the conversation to a higher level [0.457]
#> C1: Am exacting in my work [0.450]
#> N1: Get angry easily [-0.444]
#> C2: Continue until everything is perfect [0.443]
#> C3: Do things according to a plan [0.411]
#> O1: Am full of ideas [0.391]
#> N2: Get irritated easily [-0.366]
#> N3: Have frequent mood swings [-0.364]
#> A1: Am indifferent to the feelings of others [-0.344]
#> N5: Panic easily [-0.319]
#> 
#> ── Level 2 (2 factors) ──
#> 
#> m2f1
#> A3: Know how to comfort others [0.730]
#> A5: Make people feel at ease [0.694]
#> E3: Know how to captivate people [0.684]
#> A2: Inquire about others' well-being [0.661]
#> E4: Make friends easily [0.638]
#> E5: Take charge [0.630]
#> E2: Find it difficult to approach others [-0.614]
#> O3: Carry the conversation to a higher level [0.565]
#> E1: Don't talk a lot [-0.525]
#> A4: Love children [0.471]
#> C2: Continue until everything is perfect [0.468]
#> O1: Am full of ideas [0.467]
#> C1: Am exacting in my work [0.421]
#> C5: Waste my time [-0.411]
#> C4: Do things in a half-way manner [-0.389]
#> C3: Do things according to a plan [0.365]
#> A1: Am indifferent to the feelings of others [-0.335]
#> 
#> m2f2
#> N3: Have frequent mood swings [-0.806]
#> N1: Get angry easily [-0.790]
#> N2: Get irritated easily [-0.789]
#> N4: Often feel blue [-0.677]
#> N5: Panic easily [-0.659]
#> C5: Waste my time [-0.491]
#> C4: Do things in a half-way manner [-0.438]
#> O4: Spend time reflecting on things [-0.358]
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> E4: Make friends easily [0.784]
#> A3: Know how to comfort others [0.750]
#> A5: Make people feel at ease [0.735]
#> E2: Find it difficult to approach others [-0.698]
#> E3: Know how to captivate people [0.664]
#> A2: Inquire about others' well-being [0.650]
#> E1: Don't talk a lot [-0.611]
#> E5: Take charge [0.507]
#> A4: Love children [0.490]
#> O3: Carry the conversation to a higher level [0.355]
#> A1: Am indifferent to the feelings of others [-0.354]
#> 
#> m3f2
#> N3: Have frequent mood swings [-0.812]
#> N2: Get irritated easily [-0.791]
#> N1: Get angry easily [-0.777]
#> N4: Often feel blue [-0.680]
#> N5: Panic easily [-0.650]
#> C5: Waste my time [-0.428]
#> O4: Spend time reflecting on things [-0.414]
#> C4: Do things in a half-way manner [-0.348]
#> 
#> m3f3
#> C1: Am exacting in my work [0.658]
#> C2: Continue until everything is perfect [0.642]
#> C4: Do things in a half-way manner [-0.609]
#> O5: Will not probe deeply into a subject [-0.553]
#> O2: Avoid difficult reading material [-0.520]
#> O3: Carry the conversation to a higher level [0.504]
#> O1: Am full of ideas [0.493]
#> C3: Do things according to a plan [0.476]
#> C5: Waste my time [-0.453]
#> E5: Take charge [0.392]
#> O4: Spend time reflecting on things [0.319]
#> 
#> ── Level 4 (4 factors) ──
#> 
#> m4f1
#> E4: Make friends easily [0.783]
#> A3: Know how to comfort others [0.732]
#> A5: Make people feel at ease [0.725]
#> E2: Find it difficult to approach others [-0.709]
#> E3: Know how to captivate people [0.681]
#> E1: Don't talk a lot [-0.633]
#> A2: Inquire about others' well-being [0.626]
#> E5: Take charge [0.485]
#> A4: Love children [0.442]
#> O3: Carry the conversation to a higher level [0.384]
#> A1: Am indifferent to the feelings of others [-0.354]
#> 
#> m4f2
#> N3: Have frequent mood swings [-0.824]
#> N1: Get angry easily [-0.798]
#> N2: Get irritated easily [-0.793]
#> N5: Panic easily [-0.700]
#> N4: Often feel blue [-0.659]
#> C5: Waste my time [-0.357]
#> O4: Spend time reflecting on things [-0.332]
#> 
#> m4f3
#> C2: Continue until everything is perfect [0.735]
#> C4: Do things in a half-way manner [-0.714]
#> C3: Do things according to a plan [0.688]
#> C1: Am exacting in my work [0.680]
#> C5: Waste my time [-0.626]
#> A4: Love children [0.397]
#> E5: Take charge [0.372]
#> 
#> m4f4
#> O5: Will not probe deeply into a subject [-0.704]
#> O3: Carry the conversation to a higher level [0.674]
#> O1: Am full of ideas [0.621]
#> O2: Avoid difficult reading material [-0.586]
#> O4: Spend time reflecting on things [0.499]
#> E3: Know how to captivate people [0.301]
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> E2: Find it difficult to approach others [-0.752]
#> E4: Make friends easily [0.747]
#> E1: Don't talk a lot [-0.701]
#> E3: Know how to captivate people [0.677]
#> E5: Take charge [0.597]
#> A5: Make people feel at ease [0.487]
#> A3: Know how to comfort others [0.415]
#> O3: Carry the conversation to a higher level [0.394]
#> N4: Often feel blue [-0.314]
#> O1: Am full of ideas [0.302]
#> 
#> m5f2
#> N3: Have frequent mood swings [-0.825]
#> N1: Get angry easily [-0.810]
#> N2: Get irritated easily [-0.805]
#> N5: Panic easily [-0.688]
#> N4: Often feel blue [-0.646]
#> C5: Waste my time [-0.344]
#> O4: Spend time reflecting on things [-0.317]
#> 
#> m5f3
#> C2: Continue until everything is perfect [0.735]
#> C4: Do things in a half-way manner [-0.716]
#> C1: Am exacting in my work [0.690]
#> C3: Do things according to a plan [0.679]
#> C5: Waste my time [-0.652]
#> E5: Take charge [0.448]
#> A4: Love children [0.345]
#> 
#> m5f4
#> A1: Am indifferent to the feelings of others [-0.704]
#> A3: Know how to comfort others [0.703]
#> A2: Inquire about others' well-being [0.692]
#> A5: Make people feel at ease [0.580]
#> A4: Love children [0.522]
#> 
#> m5f5
#> O5: Will not probe deeply into a subject [-0.705]
#> O3: Carry the conversation to a higher level [0.655]
#> O1: Am full of ideas [0.604]
#> O2: Avoid difficult reading material [-0.595]
#> O4: Spend time reflecting on things [0.551]
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
#> E2: Find it difficult to approach others [-0.752]
#> E4: Make friends easily [0.747]
#> E1: Don't talk a lot [-0.701]
#> E3: Know how to captivate people [0.677]
#> E5: Take charge [0.597]
#> 
#> m5f2
#> N3: Have frequent mood swings [-0.825]
#> N1: Get angry easily [-0.810]
#> N2: Get irritated easily [-0.805]
#> N5: Panic easily [-0.688]
#> N4: Often feel blue [-0.646]
#> 
#> m5f3
#> C2: Continue until everything is perfect [0.735]
#> C4: Do things in a half-way manner [-0.716]
#> C1: Am exacting in my work [0.690]
#> C3: Do things according to a plan [0.679]
#> C5: Waste my time [-0.652]
#> 
#> m5f4
#> A1: Am indifferent to the feelings of others [-0.704]
#> A3: Know how to comfort others [0.703]
#> A2: Inquire about others' well-being [0.692]
#> A5: Make people feel at ease [0.580]
#> A4: Love children [0.522]
#> 
#> m5f5
#> O5: Will not probe deeply into a subject [-0.705]
#> O3: Carry the conversation to a higher level [0.655]
#> O1: Am full of ideas [0.604]
#> O2: Avoid difficult reading material [-0.595]
#> O4: Spend time reflecting on things [0.551]
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
#> m3f1 [-0.354]
#> 
#> A2: Inquire about others' well-being
#> m3f1 [0.650]
#> 
#> A3: Know how to comfort others
#> m3f1 [0.750]
#> 
#> A4: Love children
#> m3f1 [0.490]
#> 
#> A5: Make people feel at ease
#> m3f1 [0.735]
#> 
#> C1: Am exacting in my work
#> m3f3 [0.658]
#> 
#> C2: Continue until everything is perfect
#> m3f3 [0.642]
#> 
#> C3: Do things according to a plan
#> m3f3 [0.476]
#> 
#> C4: Do things in a half-way manner
#> m3f3 [-0.609]
#> m3f2 [-0.348]
#> 
#> C5: Waste my time
#> m3f3 [-0.453]
#> m3f2 [-0.428]
#> m3f1 [-0.252]
#> 
#> E1: Don't talk a lot
#> m3f1 [-0.611]
#> 
#> E2: Find it difficult to approach others
#> m3f1 [-0.698]
#> 
#> E3: Know how to captivate people
#> m3f1 [0.664]
#> 
#> E4: Make friends easily
#> m3f1 [0.784]
#> 
#> E5: Take charge
#> m3f1 [0.507]
#> m3f3 [0.392]
#> 
#> N1: Get angry easily
#> m3f2 [-0.777]
#> 
#> N2: Get irritated easily
#> m3f2 [-0.791]
#> 
#> N3: Have frequent mood swings
#> m3f2 [-0.812]
#> 
#> N4: Often feel blue
#> m3f2 [-0.680]
#> 
#> N5: Panic easily
#> m3f2 [-0.650]
#> 
#> O1: Am full of ideas
#> m3f3 [0.493]
#> m3f1 [0.253]
#> 
#> O2: Avoid difficult reading material
#> m3f3 [-0.520]
#> 
#> O3: Carry the conversation to a higher level
#> m3f3 [0.504]
#> m3f1 [0.355]
#> 
#> O4: Spend time reflecting on things
#> m3f2 [-0.414]
#> m3f3 [0.319]
#> 
#> O5: Will not probe deeply into a subject
#> m3f3 [-0.553]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```
