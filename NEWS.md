# ackwards 0.0.0.9000 (dev)

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
