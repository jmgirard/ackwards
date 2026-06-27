# ackwards 0.0.0.9000 (dev)

## Milestone 6 post-review fixes

* `augment(x, data=)` now validates column names/count before projecting.
  Named data: errors if any model variable is missing; silently subsets extra
  columns (so passing the full dataset works when the model used a subset).
  Unnamed data (bare matrix): errors if column count doesn't match.
* `.compute_scores()` now warns when the scoring basis is non-Pearson
  (polychoric/Spearman): empirical score SDs will differ from 1.0 because
  the raw projection uses Pearson z-scores while `score_var` comes from the
  model's non-Pearson R (Invariant 6).
* `tidy()`, `glance()`, and `augment()` are now re-exported from `generics`
  so they work after `library(ackwards)` without also loading `broom` or
  `generics` (consistent with `autoplot()`).
* Tests added: EFA/ESEM scores+fits coverage (B1), scores truncation for
  k_eff < k (B3), three augment column-validation tests (B5). 526 total.

## Milestone 6 — Storage materialization + cfQ cleanup

* `scores = TRUE` in `ackwards()` now stores factor scores in `x$scores` as a
  named list of `n x k_j` matrices (one per level), standardized by real score
  SDs via `sqrt(diag(W'RW))` (Invariant 1 — never assume unit variance).
  Column names match the `m{k}f{j}` factor labels used throughout the object.
* Added `augment.ackwards()` (broom generic): appends score columns (`.m{k}f{j}`)
  to a supplied data frame by applying the stored weight matrices to new or
  training data. When `data = NULL`, uses pre-stored scores (requires
  `scores = TRUE` at fit time). Gives an informative error if neither is
  available.
* `keep_fits = TRUE` in `ackwards()` now stores the raw per-level engine fit
  objects (psych or lavaan) in `x$fits` as a named list indexed by level.
* Added `tidy(x, what = "scores")`: long-format per-observation scores with
  columns `obs`, `level`, `factor`, `score` (requires `scores = TRUE` at fit
  time; for on-the-fly computation use `augment()`).
* `rotation = "cfQ"` (oblique) is now documented as **unsupported** (not
  "not yet implemented") across all engines. Oblique rotation confounds the
  between-level score correlations that are the method's core output (Goldberg,
  2006; Kim & Eaton, 2015). Only `cfT` (orthogonal CF, approx. varimax) is
  available. `ackwards()` now errors immediately with a clear message rather
  than forwarding to the engine.
* All three engines (`pca_levels`, `efa_levels`, `esem_levels`) now return
  `list(levels = ..., fits = ...)` consistently, matching the pre-existing ESEM
  convention. `keep_fits = FALSE` (default) keeps `fits = NULL`.

## Milestone 5 — Forbes (2023) extension

* Added `pairs = "all"` to `ackwards()`: computes between-level factor-score
  correlations for *every* pair of levels, not just adjacent ones. Reveals
  associations that span multiple levels and is required for redundancy pruning.
  Default remains `"adjacent"` (classic Goldberg). Adjacent-level edges are
  identical between the two modes.
* Added `prune = "redundant"` to `ackwards()`: identifies chains of factors
  connected by adjacent primary-parent links with `|r| >= redundancy_r` (default
  0.9) and applies Forbes's (2023) retention rule — keep the bottom node when the
  chain reaches level `k` (most specific, best-defined), else keep the top node
  (broadest manifestation). Pruning is *flag-only*: flagged nodes are retained in
  the object with `pruned = TRUE` / `prune_reason = "redundant"` in
  `x$prune$nodes`; nothing is removed.
* Added `prune = "artefact"` to `ackwards()`: computes Tucker's congruence
  coefficient (φ) for all cross-level factor pairs and stores results in
  `x$prune$phi` for researcher inspection. No factors are auto-flagged; artefact
  identification requires researcher judgment (Forbes, 2023).
* Added `redundancy_phi` argument (default `NULL`): when non-`NULL`, Tucker's φ
  must also exceed the threshold for a primary-parent link to qualify as a
  redundancy chain link (conjunctive with `redundancy_r`).
* Additive enrichments over the published method: redundancy chains report both
  `r` and `phi` per link (report-first / flag-second so borderline cases stay
  visible), plus a direct *endpoint* `r` from the all-levels edge matrix with a
  flag when it disagrees with the chain criterion (correlation is non-transitive;
  a clean adjacent chain is neither necessary nor sufficient for endpoint identity).
* Setting `prune` to anything other than `"none"` automatically upgrades `pairs`
  to `"all"` with a `cli_inform()` message (Invariant 6).
* `autoplot.ackwards()` now renders skip-level edges as curved arcs
  (`ggplot2::geom_curve()`; curvature tunable via `curvature` arg) when
  `pairs = "all"`. Pruned (redundant) nodes are displayed with a distinct fill
  (`color_pruned`, default `"grey80"`). Both features degrade gracefully when
  not applicable.
* New meta fields: `meta$pairs`, `meta$prune`, `meta$redundancy_r`,
  `meta$redundancy_phi`.
* `print.ackwards()` now shows a **Pruning** section when `prune != "none"`:
  flagged-node count per rule, and an explicit caveat that pruning is
  interpretive relabeling, not re-estimation.
* `tidy(x, what = "nodes")` exposes the Forbes pruning node-annotation table
  (`id`, `level`, `pruned`, `prune_reason`). Returns an empty frame with the
  same columns when no pruning was applied (safe to call unconditionally).
* `x$prune$phi` now includes an `abs_phi` column alongside `phi` to aid
  artefact screening when sign alignment is ambiguous for non-primary
  cross-level factor pairs.
* **Bug fix:** `match_parents()` previously used `clue::solve_LSAP` (bipartite
  bijection), which is ill-posed for bass-ackwards hierarchies — adjacent levels
  always have more children than parents, so the padding row returned indices
  beyond `nrow(E)`. Replaced with greedy per-column argmax
  (`apply(abs(E), 2, which.max)`), which is correct for many-to-one
  parent assignment. `clue` removed from `Suggests`.

## Previous milestones

* Added ESEM engine (`method = "esem"`) using `lavaan::efa()` (requires lavaan
  >= 0.6-13) with WLSMV ordinal estimation, rotation-aware SEs in `loadings_se`,
  and per-level fit indices (CFI, TLI, RMSEA, SRMR). Self-computed tenBerge
  weights keep `compute_edges()` on the algebra path for ESEM as for EFA.
* Added `cor = "polychoric"` as a general basis option for all three engines.
  PCA/EFA use `psych::polychoric()` with NPD smoothing; ESEM uses lavaan's own
  polychoric matrix (extracted from WLSMV sampstat) to ensure consistency with
  the fitted model.
* Added `estimator` argument for the ESEM engine (`NULL` auto-selects `"WLSMV"`
  for `cor = "polychoric"`, `"ML"` otherwise; user-overridable).
* Added `loadings_se` field to the §4 level contract (p × k matrix; `NULL` for
  PCA and EFA engines that do not produce rotation-aware SEs).
* Updated ordinal-detection warning: no longer fires when `cor = "polychoric"`
  is already set; message updated to remove "future release" language.

* Added EFA engine (`method = "efa"`) using `psych::fa()` with tenBerge scoring weights,
  keeping the W'RW algebra path in `compute_edges()`. Convergence failures truncate the
  hierarchy at the deepest converged level (warn + valid object). Heywood cases warn but
  do not truncate. Orthogonal rotation only for now; `fm` controls extraction method
  (`"minres"` default, `"ml"`, `"pa"`).
* Added `ackwards()` implementing Goldberg's (2006) bass-ackwards method with a PCA engine.
* Added `compute_edges()` computing between-level factor-score correlations via W'RW algebra.
* Added `print.ackwards()`, `tidy.ackwards()`, and `glance.ackwards()` for output.
* Added `ba_layout()` computing a Sugiyama-style barycenter layout for bass-ackwards diagrams.
* Added `autoplot.ackwards()` and `plot.ackwards()` for ggplot2-based hierarchy diagrams
  (requires ggplot2; edges coloured by sign, scaled by |r|, solid/dashed by strength).
* Added `suggest_k()` with parallel analysis and MAP (Velicer) criteria to guide choice of `k`.
* Fixed: `scoring$basis` in both PCA and EFA engines now reflects the actual `cor=` argument
  passed by the user rather than always reporting `"pearson"`.
* Fixed: EFA tenBerge weight fallback now labels `scoring$method = "regression"` when
  `.tenBerge_weights()` fails and psych regression weights are used instead (Invariant 6).
* Added `tidy(x, what = "fit")` returning one row per fit index × level (`chi`, `dof`,
  `p_value`, `RMSEA`, `TLI`, `BIC` for EFA; eigenvalues for PCA).
* Fixed: `rotation = "cfQ"` now errors immediately rather than silently running oblimin.
* Fixed: `scores = TRUE` / `keep_fits = TRUE` now warn that storage is not yet implemented.
* Fixed: `print.ackwards()` cut threshold now reads from the stored object rather than a param.
* Fixed: `k = 1` now errors; at least two levels are required for a hierarchy.
* Fixed: `kappa` default corrected from `1 / (2 * p)` to `1 / p`, which is the
  varimax-equivalent Crawford-Ferguson value (Crawford & Ferguson, 1970; Browne, 2001;
  Kim & Eaton, 2015). The previous value (equamax-adjacent) was silently wrong.
