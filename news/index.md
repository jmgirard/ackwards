# Changelog

## ackwards 0.0.0.9000 (dev)

### Milestone 8 — Plot customization (Waves 1 & 2)

- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `show_r = FALSE` / `r_digits = 2L`: when `show_r = TRUE`, draws
  `round(r, r_digits)` as a text label at each edge midpoint. Pairs
  naturally with `mono = TRUE` for greyscale journal figures.
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `mono = FALSE`: monochrome mode draws all edges in black with
  `linewidth` encoding `|r|` and `linetype` encoding sign (`solid` =
  positive, `dashed` = negative). Drops the `cut_strong` strong/weak
  linetype distinction (linewidth already conveys magnitude, so a second
  strength encoding on linetype is redundant).
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `show_level_labels = TRUE` / `level_label_size = 3`: draws “1
  factor”, “2 factors”, … to the left of the diagram at each level’s `y`
  position (using `coord_cartesian(clip = "off")` for overflow-safe
  placement, which was already set). Default is `TRUE` (was missing
  before).
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `node_labels = NULL`: a named character vector mapping factor
  IDs (e.g. `"m5f1"`) to custom display strings. Unspecified factors
  keep their `m{k}f{j}` label. Warns if a supplied name matches no
  factor ID in the object.
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `primary_only = FALSE`: when `TRUE`, filters the edge table to
  `is_primary == TRUE` before drawing, producing a clean tree of primary
  parent links. Because skip-level edges are never primary, this also
  suppresses curved skip arcs.
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `drop_pruned = FALSE`: activates the Forbes
  2023. pruned-view rendering path. Pruned nodes are removed from the
        diagram entirely; each retained node is connected to its single
        strongest kept ancestor by a straight arrow, even when the arrow
        spans multiple level gaps (no curved arcs in this mode).
        Requires `prune != "none"` at fit time; errors with a clear
        message if pruning annotations are absent. `show_r`
        auto-defaults to `TRUE` in this mode (Forbes convention).
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  gains `compress_levels = FALSE`: used together with
  `drop_pruned = TRUE`, closes the vertical gaps left by pruned levels
  so retained nodes are evenly spaced; level axis labels still show the
  original level numbers so each row remains identifiable.
- Added private helper `.drop_pruned_nodes(x, nodes, compress_levels)`
  in `layout.R`: computes the kept-only node set and a reduced
  primary-edge table for the drop-pruned rendering path. For each kept
  node it selects the single strongest `|r|` edge to any kept shallower
  node, recomputing primary parentage on the reduced graph (the original
  `is_primary` column is not reused). When `compress_levels = TRUE`,
  re-indexes `y = -rank(level)` while preserving the original `level`
  column for display labels.
- `show_r` default changed from `FALSE` to `NULL` to support the
  `drop_pruned` auto-default; calling `autoplot(x)` without `show_r`
  produces the same behaviour as before (`FALSE`).

### M6/M7 conformance fixes (second-pass review)

- `augment(x, data=)` and `scores = TRUE` now warn when `data` contains
  rows with missing values: score projection applies weights row-wise
  and propagates NAs listwise, even though fitting used
  pairwise-complete correlations. The warning advises `na.omit(data)` to
  remove incomplete rows before scoring.
  [`?augment.ackwards`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
  documents the listwise-vs-pairwise asymmetry.
- `augment(x, data=)` now errors with a clear message when `data`
  contains non-numeric columns (previously would error inside
  [`scale()`](https://rdrr.io/r/base/scale.html) with a cryptic
  message).
- `detect_ordinal()` is now called once in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  and its result reused for both the Invariant-6 console warning and
  `meta$ordinal_warned` (was called twice, scanning all columns each
  time).
- The non-Pearson basis warning in `.compute_scores()` now uses
  `.frequency = "once"` so repeated
  [`augment()`](https://generics.r-lib.org/reference/augment.html) calls
  on the same object do not produce duplicate warnings in a session.
- Tests added: non-Pearson basis warns on scoring (Fix 3);
  `keep_fits = TRUE` only stores converged levels under truncation (Fix
  4, labelled B4); NA-data warns and produces NA scores for
  [`augment()`](https://generics.r-lib.org/reference/augment.html) and
  `scores = TRUE` (Fix 1); non-numeric
  [`augment()`](https://generics.r-lib.org/reference/augment.html) data
  errors cleanly (Fix 5).
- DESIGN.md §12 and CLAUDE.md dependency section updated to match
  DESCRIPTION: `Imports` is `cli`, `generics`, `rlang`, `stats`, `utils`
  (`methods` is not imported; no `methods::` usage exists). `Suggests`
  reflects the actual set (`psych`, `GPArotation`, `lavaan`, `ggplot2`,
  `testthat`, `knitr`, `rmarkdown`, `covr`);
  `ggraph`/`igraph`/`tidygraph`, `EGAnet`/`paran`, `clue`, `future` are
  not in DESCRIPTION and the §12 table no longer lists them.
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  uses
  [`psych::fa.parallel`](https://rdrr.io/pkg/psych/man/fa.parallel.html) +
  [`psych::vss`](https://rdrr.io/pkg/psych/man/VSS.html) (not separate
  `EGAnet`/`paran` deps).

### Milestone 6 post-review fixes

- `augment(x, data=)` now validates column names/count before
  projecting. Named data: errors if any model variable is missing;
  silently subsets extra columns (so passing the full dataset works when
  the model used a subset). Unnamed data (bare matrix): errors if column
  count doesn’t match.
- `.compute_scores()` now warns when the scoring basis is non-Pearson
  (polychoric/Spearman): empirical score SDs will differ from 1.0
  because the raw projection uses Pearson z-scores while `score_var`
  comes from the model’s non-Pearson R (Invariant 6).
- [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`glance()`](https://generics.r-lib.org/reference/glance.html), and
  [`augment()`](https://generics.r-lib.org/reference/augment.html) are
  now re-exported from `generics` so they work after
  [`library(ackwards)`](https://jmgirard.github.io/ackwards/) without
  also loading `broom` or `generics` (consistent with
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)).
- Tests added: EFA/ESEM scores+fits coverage (B1), scores truncation for
  k_eff \< k (B3), three augment column-validation tests (B5). 526
  total.

### Milestone 6 — Storage materialization + cfQ cleanup

- `scores = TRUE` in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  now stores factor scores in `x$scores` as a named list of `n x k_j`
  matrices (one per level), standardized by real score SDs via
  `sqrt(diag(W'RW))` (Invariant 1 — never assume unit variance). Column
  names match the `m{k}f{j}` factor labels used throughout the object.
- Added
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
  (broom generic): appends score columns (`.m{k}f{j}`) to a supplied
  data frame by applying the stored weight matrices to new or training
  data. When `data = NULL`, uses pre-stored scores (requires
  `scores = TRUE` at fit time). Gives an informative error if neither is
  available.
- `keep_fits = TRUE` in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  now stores the raw per-level engine fit objects (psych or lavaan) in
  `x$fits` as a named list indexed by level.
- Added `tidy(x, what = "scores")`: long-format per-observation scores
  with columns `obs`, `level`, `factor`, `score` (requires
  `scores = TRUE` at fit time; for on-the-fly computation use
  [`augment()`](https://generics.r-lib.org/reference/augment.html)).
- `rotation = "cfQ"` (oblique) is now documented as **unsupported** (not
  “not yet implemented”) across all engines. Oblique rotation confounds
  the between-level score correlations that are the method’s core output
  (Goldberg, 2006; Kim & Eaton, 2015). Only `cfT` (orthogonal CF,
  approx. varimax) is available.
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  now errors immediately with a clear message rather than forwarding to
  the engine.
- All three engines (`pca_levels`, `efa_levels`, `esem_levels`) now
  return `list(levels = ..., fits = ...)` consistently, matching the
  pre-existing ESEM convention. `keep_fits = FALSE` (default) keeps
  `fits = NULL`.

### Milestone 5 — Forbes (2023) extension

- Added `pairs = "all"` to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md):
  computes between-level factor-score correlations for *every* pair of
  levels, not just adjacent ones. Reveals associations that span
  multiple levels and is required for redundancy pruning. Default
  remains `"adjacent"` (classic Goldberg). Adjacent-level edges are
  identical between the two modes.
- Added `prune = "redundant"` to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md):
  identifies chains of factors connected by adjacent primary-parent
  links with `|r| >= redundancy_r` (default 0.9) and applies
  Forbes’s (2023) retention rule — keep the bottom node when the chain
  reaches level `k` (most specific, best-defined), else keep the top
  node (broadest manifestation). Pruning is *flag-only*: flagged nodes
  are retained in the object with `pruned = TRUE` /
  `prune_reason = "redundant"` in `x$prune$nodes`; nothing is removed.
- Added `prune = "artefact"` to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md):
  computes Tucker’s congruence coefficient (φ) for all cross-level
  factor pairs and stores results in `x$prune$phi` for researcher
  inspection. No factors are auto-flagged; artefact identification
  requires researcher judgment (Forbes, 2023).
- Added `redundancy_phi` argument (default `NULL`): when non-`NULL`,
  Tucker’s φ must also exceed the threshold for a primary-parent link to
  qualify as a redundancy chain link (conjunctive with `redundancy_r`).
- Additive enrichments over the published method: redundancy chains
  report both `r` and `phi` per link (report-first / flag-second so
  borderline cases stay visible), plus a direct *endpoint* `r` from the
  all-levels edge matrix with a flag when it disagrees with the chain
  criterion (correlation is non-transitive; a clean adjacent chain is
  neither necessary nor sufficient for endpoint identity).
- Setting `prune` to anything other than `"none"` automatically upgrades
  `pairs` to `"all"` with a `cli_inform()` message (Invariant 6).
- [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  now renders skip-level edges as curved arcs
  ([`ggplot2::geom_curve()`](https://ggplot2.tidyverse.org/reference/geom_segment.html);
  curvature tunable via `curvature` arg) when `pairs = "all"`. Pruned
  (redundant) nodes are displayed with a distinct fill (`color_pruned`,
  default `"grey80"`). Both features degrade gracefully when not
  applicable.
- New meta fields: `meta$pairs`, `meta$prune`, `meta$redundancy_r`,
  `meta$redundancy_phi`.
- [`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
  now shows a **Pruning** section when `prune != "none"`: flagged-node
  count per rule, and an explicit caveat that pruning is interpretive
  relabeling, not re-estimation.
- `tidy(x, what = "nodes")` exposes the Forbes pruning node-annotation
  table (`id`, `level`, `pruned`, `prune_reason`). Returns an empty
  frame with the same columns when no pruning was applied (safe to call
  unconditionally).
- `x$prune$phi` now includes an `abs_phi` column alongside `phi` to aid
  artefact screening when sign alignment is ambiguous for non-primary
  cross-level factor pairs.
- **Bug fix:** `match_parents()` previously used `clue::solve_LSAP`
  (bipartite bijection), which is ill-posed for bass-ackwards
  hierarchies — adjacent levels always have more children than parents,
  so the padding row returned indices beyond `nrow(E)`. Replaced with
  greedy per-column argmax (`apply(abs(E), 2, which.max)`), which is
  correct for many-to-one parent assignment. `clue` removed from
  `Suggests`.

### Previous milestones

- Added ESEM engine (`method = "esem"`) using
  [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) (requires
  lavaan \>= 0.6-13) with WLSMV ordinal estimation, rotation-aware SEs
  in `loadings_se`, and per-level fit indices (CFI, TLI, RMSEA, SRMR).
  Self-computed tenBerge weights keep
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  on the algebra path for ESEM as for EFA.

- Added `cor = "polychoric"` as a general basis option for all three
  engines. PCA/EFA use
  [`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  with NPD smoothing; ESEM uses lavaan’s own polychoric matrix
  (extracted from WLSMV sampstat) to ensure consistency with the fitted
  model.

- Added `estimator` argument for the ESEM engine (`NULL` auto-selects
  `"WLSMV"` for `cor = "polychoric"`, `"ML"` otherwise;
  user-overridable).

- Added `loadings_se` field to the §4 level contract (p × k matrix;
  `NULL` for PCA and EFA engines that do not produce rotation-aware
  SEs).

- Updated ordinal-detection warning: no longer fires when
  `cor = "polychoric"` is already set; message updated to remove “future
  release” language.

- Added EFA engine (`method = "efa"`) using
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) with tenBerge
  scoring weights, keeping the W’RW algebra path in
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md).
  Convergence failures truncate the hierarchy at the deepest converged
  level (warn + valid object). Heywood cases warn but do not truncate.
  Orthogonal rotation only for now; `fm` controls extraction method
  (`"minres"` default, `"ml"`, `"pa"`).

- Added
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  implementing Goldberg’s (2006) bass-ackwards method with a PCA engine.

- Added
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  computing between-level factor-score correlations via W’RW algebra.

- Added
  [`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
  [`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
  and
  [`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)
  for output.

- Added
  [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)
  computing a Sugiyama-style barycenter layout for bass-ackwards
  diagrams.

- Added
  [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  and
  [`plot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  for ggplot2-based hierarchy diagrams (requires ggplot2; edges coloured
  by sign, scaled by \|r\|, solid/dashed by strength).

- Added
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  with parallel analysis and MAP (Velicer) criteria to guide choice of
  `k`.

- Fixed: `scoring$basis` in both PCA and EFA engines now reflects the
  actual `cor=` argument passed by the user rather than always reporting
  `"pearson"`.

- Fixed: EFA tenBerge weight fallback now labels
  `scoring$method = "regression"` when `.tenBerge_weights()` fails and
  psych regression weights are used instead (Invariant 6).

- Added `tidy(x, what = "fit")` returning one row per fit index × level
  (`chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC` for EFA; eigenvalues
  for PCA).

- Fixed: `rotation = "cfQ"` now errors immediately rather than silently
  running oblimin.

- Fixed: `scores = TRUE` / `keep_fits = TRUE` now warn that storage is
  not yet implemented.

- Fixed:
  [`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
  cut threshold now reads from the stored object rather than a param.

- Fixed: `k = 1` now errors; at least two levels are required for a
  hierarchy.

- Fixed: `kappa` default corrected from `1 / (2 * p)` to `1 / p`, which
  is the varimax-equivalent Crawford-Ferguson value (Crawford &
  Ferguson, 1970; Browne, 2001; Kim & Eaton, 2015). The previous value
  (equamax-adjacent) was silently wrong.
