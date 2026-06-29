# Interpreting and Labeling Factors

``` r

library(ackwards)
bfi <- na.omit(psych::bfi[, 1:25])
# BFI items are ordinal, so use polychoric correlations (matches the intro and
# visualization vignettes). This also avoids the ordinal-detection warning.
x <- ackwards(bfi, k_max = 5, cor = "polychoric")
```

Fitting an `ackwards` model gives you a hierarchy of factors with stable
but opaque IDs — `m1f1`, `m2f1`, `m2f2`, and so on. Turning those IDs
into something you can reason about is *interpretive* work: you read
each factor’s loadings, decide what construct it represents, and give it
a name. This is the part of the analysis that requires judgment, and it
is harder in a bass-ackwards hierarchy than in a single flat factor
solution, because the same construct can appear, split, and merge across
levels.

This article covers the full workflow:

1.  **Read** each factor with
    [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md).
2.  **Understand** the sign convention so you don’t misread a factor.
3.  **Name** factors in a way that respects the hierarchy.
4.  **Apply** your names to the diagram with
    [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
    and `autoplot(node_labels = ...)`.

## Reading a factor with `top_items()`

A factor is defined by the items that load strongly on it.
[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
lists, for each factor, the items whose absolute loading meets a
threshold, sorted from strongest to weakest. This is far easier to read
than a full item-by-factor matrix, especially at deeper levels.

``` r

top_items(x, level = 5)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.3
#> Top-n: all
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> N1 [0.829]
#> N2 [0.818]
#> N3 [0.817]
#> N4 [0.671]
#> N5 [0.657]
#> C5 [0.332]
#> m5f2
#> E2 [-0.752]
#> E4 [0.734]
#> E1 [-0.707]
#> E3 [0.636]
#> E5 [0.603]
#> A5 [0.455]
#> A3 [0.375]
#> N4 [-0.365]
#> O3 [0.363]
#> m5f3
#> C2 [0.759]
#> C4 [-0.719]
#> C3 [0.701]
#> C1 [0.673]
#> C5 [-0.651]
#> E5 [0.366]
#> m5f4
#> A2 [0.743]
#> A3 [0.705]
#> A1 [-0.677]
#> A5 [0.594]
#> A4 [0.540]
#> m5f5
#> O5 [-0.708]
#> O3 [0.676]
#> O1 [0.641]
#> O2 [-0.621]
#> O4 [0.548]
#> E3 [0.311]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```

The default `cut = 0.3` shows every item with `|loading| >= 0.3`. Raise
it to isolate the items that most define each factor, or lower it to
surface weaker cross-loadings.

``` r

top_items(x, level = 5, cut = 0.45)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.45
#> Top-n: all
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> N1 [0.829]
#> N2 [0.818]
#> N3 [0.817]
#> N4 [0.671]
#> N5 [0.657]
#> m5f2
#> E2 [-0.752]
#> E4 [0.734]
#> E1 [-0.707]
#> E3 [0.636]
#> E5 [0.603]
#> A5 [0.455]
#> m5f3
#> C2 [0.759]
#> C4 [-0.719]
#> C3 [0.701]
#> C1 [0.673]
#> C5 [-0.651]
#> m5f4
#> A2 [0.743]
#> A3 [0.705]
#> A1 [-0.677]
#> A5 [0.594]
#> A4 [0.540]
#> m5f5
#> O5 [-0.708]
#> O3 [0.676]
#> O1 [0.641]
#> O2 [-0.621]
#> O4 [0.548]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```

When a factor has many salient items, `n` caps the list at the strongest
few so you can see the core of each factor at a glance:

``` r

top_items(x, level = 5, cut = 0.3, n = 4)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.3
#> Top-n: 4
#> 
#> ── Level 5 (5 factors) ──
#> 
#> m5f1
#> N1 [0.829]
#> N2 [0.818]
#> N3 [0.817]
#> N4 [0.671]
#> m5f2
#> E2 [-0.752]
#> E4 [0.734]
#> E1 [-0.707]
#> E3 [0.636]
#> m5f3
#> C2 [0.759]
#> C4 [-0.719]
#> C3 [0.701]
#> C1 [0.673]
#> m5f4
#> A2 [0.743]
#> A3 [0.705]
#> A1 [-0.677]
#> A5 [0.594]
#> m5f5
#> O5 [-0.708]
#> O3 [0.676]
#> O1 [0.641]
#> O2 [-0.621]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```

By default items are sorted by descending `|loading|`. If your items
have a meaningful order (for example, a numbered scale you want to keep
in sequence), set `sort = FALSE`.

[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
returns an object whose `$data` field is a plain data frame — a
filtered, sorted slice of `tidy(x, what = "loadings")` — so you can
compute on it if you need to:

``` r

ti <- top_items(x, level = 5, cut = 0.3)
head(ti$data)
#>   level factor item   loading
#> 1     5   m5f1   N1 0.8292885
#> 2     5   m5f1   N2 0.8181280
#> 3     5   m5f1   N3 0.8170949
#> 4     5   m5f1   N4 0.6711163
#> 5     5   m5f1   N5 0.6566396
#> 6     5   m5f1   C5 0.3321598
```

### Cross-loadings are signal

Items that load on more than one factor are not noise to be suppressed —
they often tell you how two factors relate. Lowering the cut reveals
them:

``` r

top_items(x, level = 3, cut = 0.25)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.25
#> Top-n: all
#> 
#> ── Level 3 (3 factors) ──
#> 
#> m3f1
#> E4 [0.760]
#> A3 [0.733]
#> A5 [0.724]
#> E3 [0.653]
#> A2 [0.652]
#> E2 [-0.649]
#> E1 [-0.604]
#> A4 [0.536]
#> E5 [0.506]
#> O3 [0.345]
#> A1 [-0.320]
#> N4 [-0.271]
#> m3f2
#> N3 [0.807]
#> N1 [0.794]
#> N2 [0.794]
#> N4 [0.711]
#> N5 [0.632]
#> C5 [0.436]
#> O4 [0.389]
#> C4 [0.364]
#> E2 [0.298]
#> m3f3
#> C1 [0.657]
#> C2 [0.634]
#> C4 [-0.593]
#> O1 [0.540]
#> O5 [-0.540]
#> O3 [0.525]
#> O2 [-0.510]
#> C3 [0.479]
#> E5 [0.430]
#> C5 [-0.425]
#> O4 [0.351]
#> E3 [0.276]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```

An item that appears under two factors at the same level marks a point
where the constructs overlap. Whether that overlap is substantively
meaningful or a sign of overextraction is a judgment call — see
[`vignette("ackwards-suggest-k")`](https://jmgirard.github.io/ackwards/articles/ackwards-suggest-k.md)
for the overextraction discussion.

## The sign convention: negative does not mean “low”

Before naming anything, understand how `ackwards` orients factors.
Loadings are **sign-aligned to the primary parent** (see
[`?ackwards`](https://jmgirard.github.io/ackwards/reference/ackwards.md)):
the level-1 factor is anchored so its loadings sum positive, and every
deeper factor is flipped so its correlation with its primary parent is
positive. This makes the diagram readable, but it has a consequence for
interpretation.

A factor’s *sign* is arbitrary in the sense that flipping every loading
and the factor’s orientation describes the same dimension. So a column
of negative loadings does **not** mean “low” on that construct — it
means the construct’s positive pole was oriented the other way by the
alignment step. Read the *pattern* of items, not the bare sign:

``` r

top_items(x, level = 2)
#> 
#> ── Salient items by factor (ackwards) ──────────────────────────────────────────
#> Engine: pca
#> Cut: |loading| >= 0.3
#> Top-n: all
#> 
#> ── Level 2 (2 factors) ──
#> 
#> m2f1
#> E3 [0.694]
#> A3 [0.667]
#> E5 [0.662]
#> A5 [0.644]
#> E4 [0.614]
#> A2 [0.611]
#> O3 [0.589]
#> E2 [-0.572]
#> O1 [0.500]
#> E1 [-0.478]
#> A4 [0.473]
#> C2 [0.466]
#> C1 [0.450]
#> C4 [-0.414]
#> C5 [-0.384]
#> C3 [0.359]
#> O5 [-0.300]
#> m2f2
#> N3 [-0.808]
#> N2 [-0.797]
#> N1 [-0.795]
#> N4 [-0.726]
#> N5 [-0.629]
#> C5 [-0.439]
#> O4 [-0.400]
#> C4 [-0.359]
#> E2 [-0.339]
#> ────────────────────────────────────────────────────────────────────────────────
#> Loadings reflect primary-parent sign alignment. Use tidy(x, what = "loadings")
#> for the full matrix.
```

If a Neuroticism factor shows up with negative loadings on the anxiety
items, that is an orientation artifact of the alignment, not a “low
anxiety” factor. Name it for the construct (Neuroticism), and if the
sign matters for your downstream use, flip the scores yourself. The
values shown by
[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
and `tidy(what = "loadings")` are the aligned values, so they are
consistent with the diagram and the edge table.

## Naming in a hierarchy

In a flat factor solution you name each factor once. In a bass-ackwards
hierarchy you are naming factors at *every* level, and the levels are
related — so the names should be related too. The lineage tells you how.

``` r

summary(x)
#> 
#> ── Summary: Bass-Ackwards Analysis (ackwards) ──────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: polychoric
#> n: 2,436
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> k = 1: 1 factor (22.9% cumulative variance)
#> m1f1 22.9% eigenvalue 5.73
#> k = 2: 2 factors (34.74% cumulative variance)
#> m2f1 20.28% eigenvalue 5.73
#> m2f2 14.46% eigenvalue 2.96
#> k = 3: 3 factors (43.92% cumulative variance)
#> m3f1 17.1% eigenvalue 5.73
#> m3f2 14.05% eigenvalue 2.96
#> m3f3 12.76% eigenvalue 2.29
#> k = 4: 4 factors (51.77% cumulative variance)
#> m4f1 16.98% eigenvalue 5.73
#> m4f2 13.72% eigenvalue 2.96
#> m4f3 11.29% eigenvalue 2.29
#> m4f4 9.78% eigenvalue 1.96
#> k = 5: 5 factors (58.33% cumulative variance)
#> m5f1 13.6% eigenvalue 5.73
#> m5f2 13.38% eigenvalue 2.96
#> m5f3 11.36% eigenvalue 2.29
#> m5f4 10.27% eigenvalue 1.96
#> m5f5 9.72% eigenvalue 1.64
#> 
#> ── Lineage (primary parents) ──
#> 
#> m1f1 → m2f1, m2f2
#> m2f1 → m3f1, m3f3
#> m2f2 → m3f2
#> m3f1 → m4f1
#> m3f2 → m4f2
#> m3f3 → m4f3, m4f4
#> m4f1 → m5f2, m5f4
#> m4f2 → m5f1
#> m4f3 → m5f3
#> m4f4 → m5f5
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```

The lineage list (`m1f1 → m2f1, m2f2 → ...`) shows each factor’s primary
children. Use it to name **top-down**:

- **Upper-level factors are broader.** A level-2 factor that is the
  primary parent of two level-3 factors is the construct they share.
  Name it for the blend, not for one child. In the `bfi` data above, the
  level-1 factor is a single broad dimension; at level 2 it separates
  into a Neuroticism factor (`m2f2`, defined by the N items) and a broad
  factor blending the remaining four trait families (`m2f1`). The Big
  Five themselves do not appear cleanly until level 5 — so a substantive
  name for `m2f1` (“broad well-adjustment”, say) is necessarily coarser
  than the names you give its descendants.

- **A split is a refinement, not a contradiction.** When a parent factor
  splits into two children, the children carve up the parent’s content.
  Their names should read as specializations of the parent’s name, so
  that reading down a branch tells a coherent story.

- **Watch for factors that reorganize.** A child whose strongest parent
  is at the *other* side of the level above (a crossing edge), or a
  factor whose primary edge is weak (`|r|` well below the near-1.0
  values of stable dimensions), is a place where the structure is
  genuinely rearranging. The edge table makes these visible:

``` r

# Primary-parent edges, weakest last: the bottom rows are where structure shifts
tidy(x, what = "edges", primary_only = TRUE, sort = "strength") |> tail()
#>    from   to level_from level_to         r is_primary above_cut
#> 9  m4f1 m5f2          4        5 0.7879035       TRUE      TRUE
#> 10 m3f3 m4f3          3        4 0.7150725       TRUE      TRUE
#> 11 m3f3 m4f4          3        4 0.6989590       TRUE      TRUE
#> 12 m4f1 m5f4          4        5 0.6157559       TRUE      TRUE
#> 13 m2f1 m3f3          2        3 0.5767812       TRUE      TRUE
#> 14 m1f1 m2f2          1        2 0.4864530       TRUE      TRUE
```

Edges with `|r|` near 1.0 are factors that pass through nearly unchanged
— name the child the same as the parent. The smaller `|r|` values at the
bottom flag where a new, distinct construct is emerging and deserves its
own name.

## Applying your names

Once you have decided on names, attach them to the diagram.
[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
takes a `node_labels` argument: a named character vector mapping factor
IDs to display strings.

Typing that vector out by hand is tedious and error-prone, so
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
generates it for you, in the same order the diagram uses, and prints an
editable `c(...)` literal you can paste straight into your script:

``` r

label_template(x)
#> `label_template()` scaffold (id style):
#> c(
#>   "m1f1" = "m1f1",
#>   "m2f1" = "m2f1",
#>   "m2f2" = "m2f2",
#>   "m3f1" = "m3f1",
#>   "m3f2" = "m3f2",
#>   "m3f3" = "m3f3",
#>   "m4f1" = "m4f1",
#>   "m4f2" = "m4f2",
#>   "m4f3" = "m4f3",
#>   "m4f4" = "m4f4",
#>   "m5f1" = "m5f1",
#>   "m5f2" = "m5f2",
#>   "m5f3" = "m5f3",
#>   "m5f4" = "m5f4",
#>   "m5f5" = "m5f5"
#> )
```

Copy that literal, fill in your names, and pass it to
[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md).
Unspecified IDs keep their default `m{k}f{j}` label, so you can label
just the level you care about:

``` r

autoplot(x, node_labels = c(
  m5f1 = "Neuroticism",
  m5f2 = "Extraversion",
  m5f3 = "Conscientiousness",
  m5f4 = "Agreeableness",
  m5f5 = "Openness"
))
```

![](ackwards-interpret_files/figure-html/node-labels-1.png)

### The Forbes letter convention

Forbes (2023) labels nodes by level-letter and within-level index — `A1`
for the single level-1 factor, `B1`/`B2` at level 2, and so on.
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
produces this convention directly with `style = "forbes"`:

``` r

autoplot(x, node_labels = label_template(x, style = "forbes"))
#> `label_template()` scaffold (forbes style):
#> c(
#>   "m1f1" = "A1",
#>   "m2f1" = "B1",
#>   "m2f2" = "B2",
#>   "m3f1" = "C1",
#>   "m3f2" = "C2",
#>   "m3f3" = "C3",
#>   "m4f1" = "D1",
#>   "m4f2" = "D2",
#>   "m4f3" = "D3",
#>   "m4f4" = "D4",
#>   "m5f1" = "E1",
#>   "m5f2" = "E2",
#>   "m5f3" = "E3",
#>   "m5f4" = "E4",
#>   "m5f5" = "E5"
#> )
```

![](ackwards-interpret_files/figure-html/label-forbes-1.png)

This is useful when you want to refer to nodes by position rather than
by substantive name — for example, in a methods section that walks
through the hierarchy before interpreting it.

### Starting from a blank slate

If you would rather supply every label yourself with no defaults showing
through, `style = "blank"` gives you an all-empty scaffold to fill in:

``` r

labs <- label_template(x, style = "blank")
#> `label_template()` scaffold (blank style):
#> c(
#>   "m1f1" = "",
#>   "m2f1" = "",
#>   "m2f2" = "",
#>   "m3f1" = "",
#>   "m3f2" = "",
#>   "m3f3" = "",
#>   "m4f1" = "",
#>   "m4f2" = "",
#>   "m4f3" = "",
#>   "m4f4" = "",
#>   "m5f1" = "",
#>   "m5f2" = "",
#>   "m5f3" = "",
#>   "m5f4" = "",
#>   "m5f5" = ""
#> )
labs["m5f1"] <- "N"
labs["m5f2"] <- "E"
labs["m5f3"] <- "C"
labs["m5f4"] <- "A"
labs["m5f5"] <- "O"
autoplot(x, node_labels = labs)
```

![](ackwards-interpret_files/figure-html/label-blank-1.png)

## Where to go next

This article is about *what* to put on the diagram. For *how* the
diagram looks — colours, edge thresholds, monochrome and publication
styling, level labels, and the Forbes pruned-diagram mode — see
[`vignette("ackwards-visualization")`](https://jmgirard.github.io/ackwards/articles/ackwards-visualization.md).
