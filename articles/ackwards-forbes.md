# The Forbes Extension: Skip-Level Connections and Pruning

The classic Goldberg (2006) method computes between-level factor-score
correlations only for *adjacent* levels: 1↔︎2, 2↔︎3, 3↔︎4, and so on.
Forbes (2023) extended this in two ways: computing correlations across
*all* level pairs (not just adjacent), and using those extra connections
to identify and flag **redundant** or **artefactual** factors in the
hierarchy.

This vignette covers both extensions: `pairs = "all"` and `prune`.

## The limitation of adjacent-only edges

An adjacent-only hierarchy shows you the immediate parent–child
relationships. What it cannot tell you is whether a factor at level k is
essentially the same construct as a factor several levels up — a sign
that the intermediate levels are adding noise rather than resolution.

Consider a factor that appears at k = 2, k = 3, and k = 4 and correlates
\> 0.97 with its counterpart at every adjacent level. The adjacent-only
diagram shows three consecutive arrows, each nearly perfect. But nothing
in the diagram directly flags that the k = 3 factor is redundant: you
could skip straight from k = 2 to k = 4 without losing any information.

**Skip-level correlations** make this visible by computing the
correlation between every pair of levels, not just neighbors.

## Setup

``` r

library(ackwards)
bfi <- na.omit(bfi25)
```

## `pairs = "all"`: computing every between-level correlation

Adding `pairs = "all"` extends the edge table from adjacent pairs only
to every combination of levels.

``` r

# Classic adjacent-only
x_adj <- ackwards(bfi, k_max = 5, cor = "polychoric")

# All pairs
x_all <- ackwards(bfi, k_max = 5, cor = "polychoric", pairs = "all")

# How many edges?
nrow(tidy(x_adj, what = "edges")) # adjacent only
#> [1] 40
nrow(tidy(x_all, what = "edges")) # all pairs
#> [1] 85
```

With k = 5, the adjacent-only model has 40 edges (1×2 + 2×3 + 3×4 +
4×5). The all-pairs model adds every non-adjacent pair — 1↔︎3, 1↔︎4, 1↔︎5,
2↔︎4, 2↔︎5, 3↔︎5 — for 85 edges total.

### Reading the skip-level edge table

The table below keeps the non-adjacent edges (levels more than one
apart) with `|r| >= 0.5`, strongest first — drawn from
`tidy(x_all, what = "edges")`:

| from | to   | level_from | level_to |     r |
|:-----|:-----|-----------:|---------:|------:|
| m3f2 | m5f2 |          3 |        5 |  0.98 |
| m2f2 | m4f2 |          2 |        4 | -0.97 |
| m2f2 | m5f2 |          2 |        5 | -0.97 |
| m2f1 | m4f1 |          2 |        4 |  0.85 |
| m3f1 | m5f1 |          3 |        5 |  0.82 |
| m1f1 | m3f1 |          1 |        3 |  0.77 |
| m1f1 | m4f1 |          1 |        4 |  0.75 |
| m3f3 | m5f3 |          3 |        5 |  0.73 |
| m2f1 | m5f1 |          2 |        5 |  0.69 |
| m3f3 | m5f5 |          3 |        5 |  0.68 |
| m1f1 | m5f1 |          1 |        5 |  0.61 |
| m3f1 | m5f4 |          3 |        5 |  0.56 |

Strongest skip-level edges (\|r\| \>= 0.5, non-adjacent levels) {.table}

Several factors connect across two or more levels with correlations
above 0.90. `m3f2` (level 3, factor 2) correlates 0.98 with `m5f1`
(level 5, factor 1), jumping *two* levels. This tells you that m3f2 and
m5f1 are essentially the same construct — the intermediate levels are
just refinements within a stable dimension.

## Pruning: identifying redundant factors

The `prune` argument uses the skip-level correlations to automatically
flag factors that may not be adding genuine information to the
hierarchy.

### `prune = "redundant"`: chains of near-identical factors

A **redundant chain** is a sequence of factors connected by near-perfect
correlations (\|r\| ≥ 0.9 by default) across levels. If m2f2 → m3f2 →
m4f2 → m5f1 all share r \> 0.97, the intermediate nodes m3f2 and m4f2
are flagged as redundant: they repeat rather than refine the same
dimension.

``` r

x_prune <- ackwards(bfi,
  k_max = 5, cor = "polychoric",
  pairs = "all", prune = "redundant"
)
#> ℹ Redundancy pruning (|r| ≥ 0.9) flagged 6 nodes.
#> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
#>   `x$prune$chains`.
```

``` r

tidy(x_prune, what = "nodes")
#>      id level pruned prune_reason
#> 1  m1f1     1  FALSE         <NA>
#> 2  m2f1     2  FALSE         <NA>
#> 3  m2f2     2   TRUE    redundant
#> 4  m3f1     3  FALSE         <NA>
#> 5  m3f2     3   TRUE    redundant
#> 6  m3f3     3  FALSE         <NA>
#> 7  m4f1     4   TRUE    redundant
#> 8  m4f2     4   TRUE    redundant
#> 9  m4f3     4   TRUE    redundant
#> 10 m4f4     4   TRUE    redundant
#> 11 m5f1     5  FALSE         <NA>
#> 12 m5f2     5  FALSE         <NA>
#> 13 m5f3     5  FALSE         <NA>
#> 14 m5f4     5  FALSE         <NA>
#> 15 m5f5     5  FALSE         <NA>
```

Six factors are flagged as redundant. The entire k = 4 level (m4f1
through m4f4) is flagged, as are m2f2 and m3f2. This is a striking
finding: for this dataset and this k, the four-factor level adds little
beyond what you already know from k = 3 and k = 5.

The flagged factors are not removed from the object — `prune` is purely
a diagnostic annotation, not a deletion. You can still inspect their
loadings, use their scores, and include them in the diagram. Pruning
flags guide interpretation; they do not alter the model.

### The pruned-factor diagram

For presentations and publications it is cleaner to omit the flagged
factors entirely and draw direct connections from each retained factor
to its single strongest kept ancestor — even when that ancestor is
several levels away. This is the Forbes (2023) pruned-factor diagram,
activated by `drop_pruned = TRUE`.

Forbes (2023) presents two variants: one with correlation labels on each
arrow and one without. The first is useful when the strength of each
spanning connection matters to the interpretation; the second is cleaner
for presentations. Both are reproduced below using the same publication
style — black lines, uniform width, plain line ends, and no legend.

**With correlation labels** (`show_r = TRUE`):

``` r

autoplot(x_prune,
  drop_pruned = TRUE, show_r = TRUE,
  color_pos = "black", color_neg = "black",
  edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE
)
```

![](ackwards-forbes_files/figure-html/drop-pruned-1.png)

**Without labels** (cleaner for slides or when the exact values are not
the focus):

``` r

autoplot(x_prune,
  drop_pruned = TRUE,
  color_pos = "black", color_neg = "black",
  edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE
)
```

![](ackwards-forbes_files/figure-html/drop-pruned-unlabeled-1.png)

Level 4 is entirely pruned, leaving a visible gap in the y-axis. The gap
is intentional: it shows *which* level was removed. Spanning arrows
bridge directly from level 3 factors to level 5 factors (and from level
2 to level 5 where intermediate levels were flagged). Solid lines
indicate strong connections (\|r\| ≥ 0.5 by default); dashed lines
indicate weaker ones — the same `cut_strong` threshold used throughout.

To close the gaps and compact the layout while retaining the original
level numbers on the axis:

``` r

autoplot(x_prune,
  drop_pruned = TRUE, compress_levels = TRUE,
  color_pos = "black", color_neg = "black",
  edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE
)
```

![](ackwards-forbes_files/figure-html/drop-pruned-compressed-1.png)

The level labels still read “1 factor”, “2 factors”, “3 factors”, “5
factors” so readers know which levels were retained; the uniform
vertical spacing makes the diagram easier to read in constrained page
layouts.

For further cosmetic customization — colours, node labels, arrowheads,
and more — see
[`vignette("ackwards-visualization")`](https://jmgirard.github.io/ackwards/articles/ackwards-visualization.md).

### `prune = "artefact"`: factors defined by structural similarity

An **artefact** factor is one whose loading pattern is more similar to a
factor at a *non-adjacent* level than to its own-level neighbors.
Similarity is measured by Tucker’s congruence coefficient (φ):

``` math
\phi(F_a, F_b) = \frac{\sum_i \lambda_{ia}\lambda_{ib}}
{\sqrt{\sum_i \lambda_{ia}^2 \cdot \sum_i \lambda_{ib}^2}}
```

φ ranges from −1 to +1, with values \> 0.95 indicating near-identical
loading patterns regardless of sign. A factor at k = 3 that has φ \>
0.95 with a factor at k = 1 is suspected to be an artefact of the
rotation rather than a genuine new dimension.

``` r

x_art <- ackwards(bfi,
  k_max = 5, cor = "polychoric",
  pairs = "all", prune = "artefact"
)
#> ℹ Artefact mode: Tucker's computed for all cross-level factor pairs.
#> ℹ Inspect `x$prune$phi` to identify potential artefacts; removal is a
#>   researcher judgment (Forbes, 2023).
tidy(x_art, what = "nodes")
#>      id level pruned prune_reason
#> 1  m1f1     1  FALSE         <NA>
#> 2  m2f1     2  FALSE         <NA>
#> 3  m2f2     2  FALSE         <NA>
#> 4  m3f1     3  FALSE         <NA>
#> 5  m3f2     3  FALSE         <NA>
#> 6  m3f3     3  FALSE         <NA>
#> 7  m4f1     4  FALSE         <NA>
#> 8  m4f2     4  FALSE         <NA>
#> 9  m4f3     4  FALSE         <NA>
#> 10 m4f4     4  FALSE         <NA>
#> 11 m5f1     5  FALSE         <NA>
#> 12 m5f2     5  FALSE         <NA>
#> 13 m5f3     5  FALSE         <NA>
#> 14 m5f4     5  FALSE         <NA>
#> 15 m5f5     5  FALSE         <NA>
```

For the BFI, the artefact criterion flags nothing — no factor at any
level has a loading pattern more similar to a factor from a non-adjacent
level than to its own-level solution. This is a good result for a
well-validated instrument: it suggests the rotation is producing
genuinely distinct factors at each level, not recycling old ones.

On data with weaker or noisier factor structure, or with rotations that
struggle to separate highly correlated factors, the artefact criterion
will often flag some nodes. A node can be flagged by one criterion but
not the other; the two criteria capture different flavors of redundancy:

- `"redundant"`: this factor appears at multiple levels with
  near-identical *score correlations* — it persists unchanged as k
  increases.
- `"artefact"`: this factor’s *loading pattern* closely resembles a
  factor from a different, non-adjacent level — it looks like a copy
  rather than a refinement.

## Tuning the thresholds

The redundancy criterion has an adjustable `redundancy_r` threshold
(default `0.90`) matching Forbes (2023). The artefact criterion has no
auto-flag threshold — `prune = "artefact"` computes Tucker’s φ for
researcher inspection; no factors are auto-flagged.

For the BFI, the result is the same across a wide range of thresholds
because the redundant chains all have correlations \> 0.97 — the
flagging is unambiguous. With your own data you may find borderline
cases where the threshold matters:

Refitting at a few `redundancy_r` thresholds and counting flagged
factors at each shows how sensitive the pruning is:

    #> ℹ Redundancy pruning (|r| ≥ 0.8) flagged 9 nodes.
    #> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
    #>   `x$prune$chains`.
    #> ℹ Redundancy pruning (|r| ≥ 0.85) flagged 8 nodes.
    #> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
    #>   `x$prune$chains`.
    #> ℹ Redundancy pruning (|r| ≥ 0.9) flagged 6 nodes.
    #> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
    #>   `x$prune$chains`.
    #> ℹ Redundancy pruning (|r| ≥ 0.95) flagged 6 nodes.
    #> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
    #>   `x$prune$chains`.

| redundancy_r | n_flagged |
|-------------:|----------:|
|         0.80 |         9 |
|         0.85 |         8 |
|         0.90 |         6 |
|         0.95 |         6 |

Factors flagged redundant at each redundancy_r threshold {.table}

For the BFI all thresholds agree: the flagged factors are robustly
redundant, not borderline cases. In noisier datasets or smaller samples
you will typically see the count increase as you lower the threshold.

## Practical interpretation

The Forbes extension does not change the core bass-ackwards analysis. It
enriches it with two questions:

1.  **Do any factors persist unchanged across multiple levels?**
    (`pairs = "all"`) Skip-level correlations near 1.0 indicate stable
    dimensions that survive changes in k — exactly the kind of robust
    construct you want to report.

2.  **Are there levels where the factor structure is just reorganizing
    rather than genuinely differentiating?** (`prune = "redundant"`)
    Flagged levels can often be removed from the k range without losing
    interpretive content.

A common workflow: fit with `pairs = "all"` first to examine the full
picture, then apply `prune = "redundant"` to identify which levels add
the most new information, and use that to guide your focus in reporting.

## References

Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research
in Personality*, *40*(4), 347–358.

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg’s bass-ackward method.
*Psychological Methods*. <https://doi.org/10.1037/met0000546>

Tucker, L. R. (1951). *A method for synthesis of factor analysis
studies* (Personnel Research Section Report No. 984). Department of the
Army.
