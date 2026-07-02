# Flag redundant or artifactual factors (Forbes 2023 extension)

`prune()` never removes anything from an `ackwards` object – it only
annotates factors with `pruned`/`prune_reason` flags in `x$prune$nodes`
(flag-only, never removes; the object keeps every level). Because it is
a separate, cheap step from extraction, you can re-prune with new
thresholds without re-running
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md):


      x <- ackwards(bfi25, k_max = 6, engine = "esem")  # expensive
      x |> prune("redundant")                            # cheap, repeatable
      x |> prune("redundant", redundancy_r = 0.95)        # no re-extraction

`prune()` is an S3 generic (rather than a plain function) so it coexists
with the `prune` generics already defined by recursive-partitioning
packages (e.g.
[`rpart::prune`](https://rdrr.io/pkg/rpart/man/prune.rpart.html))
regardless of package load order.

## Usage

``` r
prune(x, ...)

# S3 method for class 'ackwards'
prune(
  x,
  rules = "none",
  manual = NULL,
  redundancy_r = 0.9,
  redundancy_phi = NULL,
  min_items = 3L,
  orphan_r = 0.5,
  ...
)
```

## Arguments

- x:

  An `ackwards` object.

- ...:

  Reserved for future methods/arguments.

- rules:

  Character vector controlling which auto-rules run. Default `"none"`
  (no auto rule; combine with `manual` for pure manual pruning, or call
  `prune(x)` with no arguments to clear any existing pruning). Options:

  - `"redundant"` – identify chains of factors connected by
    primary-parent links with `|r| >= redundancy_r` (and optionally
    `phi > redundancy_phi`). Applies Forbes's (2023) retention rule:
    keep the bottom node when the chain reaches level `k_max` (most
    specific); keep the top node otherwise. Flagged nodes get
    `pruned = TRUE` and `prune_reason = "redundant"` in `x$prune$nodes`.

  - `"artifact"` (or the alias `"artefact"`, normalized to `"artifact"`)
    – compute Tucker's congruence coefficient (phi) for all cross-level
    factor pairs and store in `x$prune$phi`, plus structural signals
    (`few_items`/`orphan`/`split_merge`) in `x$prune$structural`. No
    factors are auto-flagged; artifact identification requires judgment
    (Forbes, 2023; Wicherts et al., 2016).

- manual:

  Character vector of factor labels (e.g. `c("m4f3", "m4f4")`) to flag
  directly, in addition to or instead of an auto rule. Standalone manual
  pruning is supported: `prune(x, manual = c("m4f3"))`. Unknown labels
  error. A node already flagged by an auto rule keeps that
  `prune_reason`; only otherwise-unflagged manual nodes get
  `prune_reason = "manual"`.

- redundancy_r:

  Scalar in `(0, 1]`. Adjacent primary-parent `|r|` threshold for
  redundancy chains. Default `0.9` (Forbes, 2023).

- redundancy_phi:

  Scalar in `(0, 1]`, `NULL` (default, auto), or `NA` (explicit
  opt-out). When `NULL`:

  - `x$engine == "pca"` – no phi filter. Component scores are
    *determinate* (exact linear functions of the data, with no
    factor-score indeterminacy), so the score correlation `|r|` is the
    true correlation between the components themselves; phi adds nothing
    that `|r|` does not already capture.

  - `x$engine` is `"efa"` or `"esem"` – automatically set to `0.95`
    (Lorenzo-Seva & ten Berge, 2006). Factor-score indeterminacy off-PCA
    means `|r|`-alone is liberal; the conjunctive phi criterion is the
    conservative default. A cli message announces the resolved value.
    Pass `NA` to disable phi filtering regardless of engine. Pass a
    numeric value to override on any engine.

- min_items:

  Minimum number of items for which a factor must be the primary loader
  (highest `|loading|`). Factors with fewer than `min_items` primary
  items are flagged `few_items = TRUE` in `x$prune$structural`. Only
  used when `rules` includes `"artifact"`. Default `3L` – a factor
  defined by one or two items is under-identified and frequently an
  extraction artifact rather than a replicable construct (the classic
  "three-indicator rule"; Forbes, 2023, Fig. 2).

- orphan_r:

  Threshold for the `orphan` structural signal. A factor whose maximum
  **adjacent-level** `|r|` (to the immediately shallower and deeper
  levels) falls below `orphan_r` is flagged `orphan = TRUE` in
  `x$prune$structural` – it does not connect to the neighbouring
  solutions and so does not replicate across the hierarchy. Only used
  when `rules` includes `"artifact"`. Default `0.5` – a moderate
  correlation; a factor that shares less than a quarter of its variance
  with every neighbour is a structural outlier worth inspecting.

## Value

`x`, with `$prune` populated (replacing any prior pruning).

## References

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
coefficient as a meaningful index of factor similarity. *Methodology*,
2(2), 57–64.
[doi:10.1027/1614-2241.2.2.57](https://doi.org/10.1027/1614-2241.2.2.57)

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
(`what = "nodes"`),
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
(`drop_pruned`)

## Examples

``` r
# sim16 has a planted redundant chain + overextraction artifact at k = 5,
# so the prune rules always have a finding to show (and no ordinal warning).
x <- ackwards(sim16, k_max = 5)

xp <- prune(x, "redundant")
#> ℹ Redundancy pruning (|r| ≥ 0.9) flagged 7 nodes.
#> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
#>   `x$prune$chains`.
xp$prune$nodes
#>      id level pruned prune_reason
#> 1  m1f1     1  FALSE         <NA>
#> 2  m2f1     2  FALSE         <NA>
#> 3  m2f2     2  FALSE         <NA>
#> 4  m3f1     3   TRUE    redundant
#> 5  m3f2     3   TRUE    redundant
#> 6  m3f3     3   TRUE    redundant
#> 7  m4f1     4   TRUE    redundant
#> 8  m4f2     4   TRUE    redundant
#> 9  m4f3     4   TRUE    redundant
#> 10 m4f4     4   TRUE    redundant
#> 11 m5f1     5  FALSE         <NA>
#> 12 m5f2     5  FALSE         <NA>
#> 13 m5f3     5  FALSE         <NA>
#> 14 m5f4     5  FALSE         <NA>
#> 15 m5f5     5  FALSE         <NA>

# Re-prune with a new threshold -- no re-extraction needed
prune(x, "redundant", redundancy_r = 0.95)
#> ℹ Redundancy pruning (|r| ≥ 0.95) flagged 6 nodes.
#> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
#>   `x$prune$chains`.
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 28.2% variance
#> ✔ k = 2: 2 factors, 46.5% variance
#> ✔ k = 3: 3 factors, 57.5% variance
#> ✔ k = 4: 4 factors, 67.7% variance
#> ✔ k = 5: 5 factors, 70.8% variance
#> 
#> ── Edges ──
#> 
#> 13 of 40 edges have |r| ≥ 0.3
#> 
#> ── Pruning ──
#> 
#> Redundancy (|r| ≥ 0.95): 6 nodes flagged
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: Pruning is interpretive relabeling, not re-estimation. Flagged nodes
#> remain in the object; all edges are preserved. Inspect with `x$prune$nodes` and
#> `tidy(x, what = "nodes")`.
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.

# Manual pruning: standalone, or mixed with an auto rule
prune(x, manual = "m4f2")
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 28.2% variance
#> ✔ k = 2: 2 factors, 46.5% variance
#> ✔ k = 3: 3 factors, 57.5% variance
#> ✔ k = 4: 4 factors, 67.7% variance
#> ✔ k = 5: 5 factors, 70.8% variance
#> 
#> ── Edges ──
#> 
#> 13 of 40 edges have |r| ≥ 0.3
#> 
#> ── Pruning ──
#> 
#> Manual: 1 node explicitly flagged (m4f2)
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: Pruning is interpretive relabeling, not re-estimation. Flagged nodes
#> remain in the object; all edges are preserved. Inspect with `x$prune$nodes` and
#> `tidy(x, what = "nodes")`.
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
prune(x, "redundant", manual = "m4f2")
#> ℹ Redundancy pruning (|r| ≥ 0.9) flagged 7 nodes.
#> ℹ Nodes are retained in the object; inspect with `x$prune$nodes` and
#>   `x$prune$chains`.
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 28.2% variance
#> ✔ k = 2: 2 factors, 46.5% variance
#> ✔ k = 3: 3 factors, 57.5% variance
#> ✔ k = 4: 4 factors, 67.7% variance
#> ✔ k = 5: 5 factors, 70.8% variance
#> 
#> ── Edges ──
#> 
#> 13 of 40 edges have |r| ≥ 0.3
#> 
#> ── Pruning ──
#> 
#> Redundancy (|r| ≥ 0.9): 7 nodes flagged
#> Manual: 1 node explicitly flagged (m4f2)
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: Pruning is interpretive relabeling, not re-estimation. Flagged nodes
#> remain in the object; all edges are preserved. Inspect with `x$prune$nodes` and
#> `tidy(x, what = "nodes")`.
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
```
