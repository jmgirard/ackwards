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
  redundancy_criterion = c("direct", "adjacent"),
  min_items = 3L,
  orphan_r = 0.5,
  near_margin = 0.1,
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

  - `"redundant"` – identify chains of factors connected by score
    correlations `|r| >= redundancy_r` (and optionally
    `phi > redundancy_phi`), using `redundancy_criterion` (default
    `"direct"`, faithful to Forbes). Applies Forbes's (2023) retention
    rule: keep the bottom node when the chain reaches level `k_max`
    (most specific); keep the top node otherwise. Flagged nodes get
    `pruned = TRUE` and `prune_reason = "redundant"` in `x$prune$nodes`.

  - `"artifact"` (or the alias `"artefact"`, normalized to `"artifact"`)
    – compute Tucker's congruence coefficient (phi) for all cross-level
    factor pairs and store in `x$prune$phi`, structural signals
    (`few_items`/`orphan`/`split_merge`) in `x$prune$structural`, and
    the **near-redundant band** (`near_margin`, see Details) in
    `x$prune$near_redundant`. No factors are auto-flagged; artifact
    identification requires judgment (Forbes, 2023; Wicherts et al.,
    2016).

- manual:

  Character vector of factor labels (e.g. `c("m4f3", "m4f4")`) to flag
  directly, in addition to or instead of an auto rule. Standalone manual
  pruning is supported: `prune(x, manual = c("m4f3"))`. Unknown labels
  error. A node already flagged by an auto rule keeps that
  `prune_reason`; only otherwise-unflagged manual nodes get
  `prune_reason = "manual"`.

- redundancy_r:

  Scalar in `(0, 1]`. Score-correlation `|r|` threshold for redundancy
  chains. Default `0.9` (Forbes, 2023).

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

- redundancy_criterion:

  How redundancy chains are traced. One of:

  - `"direct"` (default) – chase upward via the **direct (skip-level)**
    correlation between a factor and each ancestor level, continuing
    while `|r| >= redundancy_r` contiguously. This is Forbes's (2023)
    published `ChaseCorrPaths` rule and the honest operationalization of
    "the same construct" (two factor scores share `>= redundancy_r^2`
    variance *directly*). It reproduces her AMH applied example exactly.

  - `"adjacent"` – trace **adjacent primary-parent** links only (each
    consecutive level `|r| >= redundancy_r`). This was the pre-M53
    default; because correlation is non-transitive it can both over- and
    under-flag versus `"direct"` in deep (many-level) hierarchies, so it
    is retained only as an opt-in. On shallow/transitive hierarchies the
    two agree.

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

- near_margin:

  Scalar in `(0, 1]`. Width of the **near-redundant band** reported by
  `"artifact"` mode (see Details). A cross-level pair is flagged
  near-redundant when its direct `|r|` sits in
  `[redundancy_r - near_margin, redundancy_r)` **or** its Tucker `phi`
  sits in `[redundancy_phi - near_margin, redundancy_phi)` – i.e. within
  `near_margin` *below* a redundancy threshold – while the pair is not
  itself fully redundant. Only used when `rules` includes `"artifact"`.
  Default `0.1`.

## Value

`x`, with `$prune` populated (replacing any prior pruning).

## Details

**The `"direct"` criterion is a star anchored on the leaf, not a walk.**
Take a three-level chain candidate with deepest leaf `m3f1` and
shallower factors `m2f1` (level 2) and `m1f1` (level 1). Under the
default `"direct"` criterion the chain `m1f1 -> m2f1 -> m3f1` forms when
**both** direct-to-leaf correlations `|r(m1f1, m3f1)|` and
`|r(m2f1, m3f1)|` meet `redundancy_r` – every member is judged by its
own direct correlation to the *same* deepest factor (a star centred on
the leaf). It does **not** require the adjacent hop `|r(m1f1, m2f1)|` to
meet the threshold (that is the `"adjacent"` criterion), and it does
**not** screen every ancestor pair against every other (there is no
all-pairs test). So an ancestor can join on a strong direct link to the
leaf even where its adjacent hop to the next chain member is weak –
which is exactly why `r_to_prev` (below) can sit under `redundancy_r`.

**Reading `x$prune$chains` under `redundancy_criterion = "direct"`.**
The `r_to_prev` and `phi_to_prev` columns report the **adjacent-level**
correlation and congruence between consecutive chain members (for
continuity of display), but chain *membership* is decided by the
**direct (skip-level)** correlation to the chain's deepest factor. A
direct chain can therefore legitimately contain a link whose `r_to_prev`
is *below* `redundancy_r` – the stronger direct link is what justified
it. The `endpoint_r` column gives the direct root-to-leaf correlation as
an at-a-glance cross-check. Under `redundancy_criterion = "adjacent"`,
`r_to_prev` *is* the criterion and always meets `redundancy_r`.

**The near-redundant band (`"artifact"` mode).** `prune("redundant")`
drops *full* redundancy – pairs at or above the thresholds. Forbes
(2023) uses the artifact flags mainly for the messier band *just below*
them: a pair that correlates, say, `r = 0.89` and shares a loading
pattern `phi = 0.93` is not quite redundant but is a candidate
re-rotation worth a second look. Artifact mode surfaces this as
`x$prune$near_redundant` – a data frame of every cross-level pair that
is **not** itself fully redundant yet has its direct (skip-level) `|r|`
**or** its Tucker `phi` within `near_margin` *below* the corresponding
threshold (`redundancy_r` / `redundancy_phi`). Columns: `from`, `to`,
`level_from`, `level_to`, `r` (direct score correlation), `phi` (loading
congruence), and the logical band flags `near_r` / `near_phi`. Under the
PCA engine `redundancy_phi` is `NULL` (component scores are
determinate), so only the `|r|` band applies; under EFA/ESEM
`redundancy_phi` auto-resolves to `0.95` (announced via cli) and the
`phi` band is active too. Like phi and the structural signals, the band
is **report-only** – it never drops or mutates the kept node set.

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
#> ℹ Redundancy pruning (direct criterion, |r| ≥ 0.9) flagged 7 nodes.
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
#> ℹ Redundancy pruning (direct criterion, |r| ≥ 0.95) flagged 6 nodes.
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
#> Redundancy (direct, |r| ≥ 0.95): 6 nodes flagged
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
#> ℹ Redundancy pruning (direct criterion, |r| ≥ 0.9) flagged 7 nodes.
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
#> Redundancy (direct, |r| ≥ 0.9): 7 nodes flagged
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
