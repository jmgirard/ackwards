# ackwards 0.1.0

First public release. Licensed under MIT.

## Extraction engines

Three engines are available via the `engine` argument to `ackwards()`:

* **PCA** (`engine = "pca"`) — the original Goldberg (2006) bass-ackwards method.
  Fastest; never fails to converge; exact W′RW algebra for between-level edges.
* **EFA** (`engine = "efa"`) — exploratory factor analysis via `psych::fa()` with
  tenBerge scoring weights and the same closed-form algebra path.
* **ESEM** (`engine = "esem"`) — EFA estimated in a SEM framework via `lavaan::efa()`
  with WLSMV ordinal estimation, per-level fit indices (CFI, TLI, RMSEA, SRMR), and
  rotation-aware standard errors. Required for the clinical/HiTOP workflow.

All engines use **varimax** (orthogonal) rotation — the only rotation compatible with the
W′RW algebra. Oblique rotation is unsupported by design.

## Choosing k — `suggest_k()`

`suggest_k()` runs five complementary retention criteria and returns a consensus range:

* **PA-PC** and **PA-FA** — parallel analysis on the PC and FA eigenvalue bases
  (`psych::fa.parallel(fa = "both")`).
* **MAP** (Velicer 1976) and **VSS-1/VSS-2** (Revelle & Rocklin 1979) via `psych::vss()`.
* **Comparison Data** (Ruscio & Roche 2012) via `EFAtools::CD()` — skipped gracefully when
  `EFAtools` is not installed.

`autoplot.suggest_k()` produces a three-panel ggplot2 diagnostic (scree/PA, MAP, VSS).

## Between-level edges — exact W′RW algebra

Between-level factor-score correlations are computed via closed-form `W′RW` algebra for all
linear engines (PCA, EFA, ESEM with tenBerge weights) — no scores need to be materialized.
Standardization uses real score SDs from `sqrt(diag(W′RW))` rather than assuming unit
variance. A materialized-scores fallback is available and a cross-check test verifies the
two routes agree within tolerance.

## Ordinal / Likert data

Pass `cor = "polychoric"` to use polychoric correlations as the basis for PCA or EFA.
For ESEM, lavaan's own polychoric matrix is used for edge computation (consistent with the
fitted model). The WLSMV estimator is the default for `cor = "polychoric"` with ESEM.
`ackwards()` detects likely-ordinal columns (≤ 7 distinct integer values) and emits a
suppressible warning when `cor = "pearson"` is set.

## Missing data — `missing` argument

Three strategies via `missing` on `ackwards()`:

* `"pairwise"` (default) — pairwise-complete correlations; warns when NAs are present.
* `"listwise"` — complete cases only, consistent N across fit and edges.
* `"fiml"` — Full Information ML via lavaan (ESEM ML/MLR only); edges derived from the
  FIML saturated model.

## Forbes (2023) extension

* `pairs = "all"` computes between-level correlations for every level pair, revealing
  skip-level associations.
* `prune = "redundant"` flags chains where adjacent `|r| ≥ 0.9` (Forbes retention rule).
* `prune = "artefact"` computes Tucker's φ for researcher inspection; no auto-flagging.
* `autoplot()` renders pruned views (`drop_pruned = TRUE`) and skip-level arcs.

## Visualization — `autoplot.ackwards()`

A ggplot2 layered diagram with extensive customization: edge color by sign, width/alpha by
`|r|`, solid/dashed by strength; monochrome mode (`mono`); APA-formatted correlation labels
(`show_r`, `r_digits`, `r_label_size`); custom node labels (`node_labels`); primary-only
mode (`primary_only`); level axis labels (`show_level_labels`); arrowhead toggle
(`show_arrows`); fixed edge width (`edge_linewidth`); legend toggle (`legend`); Forbes-style
pruned view (`drop_pruned`, `compress_levels`).

## Tidy interface and scoring

* `print.ackwards()` — compact cli summary card.
* `summary.ackwards()` — per-level variance, fit indices, lineage list, pruning notes.
* `tidy(x, what = "edges" / "loadings" / "variance" / "nodes" / "scores")` — long-format
  tidy tibbles.
* `glance(x)` — one-row model-level summary.
* `augment(x, data = ...)` — appends factor scores to a data frame; recomputes from stored
  weights when scores were not materialized at fit time.

---

## Development history

## Milestone 16 — Estimator-aware missing-data handling

## Milestone 16 — Estimator-aware missing-data handling

* New `missing` argument on `ackwards()`. Controls how incomplete rows are
  handled across all three engines. `"pairwise"` (default) preserves existing
  behaviour — pairwise-complete correlations — and emits a warning when missing
  rows are detected. `"listwise"` reduces data to complete cases before all
  downstream steps, so the correlation matrix, the engine fit, and the edges
  are all consistent; `n_obs` in the result reflects the reduced N. `"fiml"`
  (ESEM ML/MLR only) passes `missing = "fiml"` to `lavaan::efa()` and derives
  edge correlations from lavaan's FIML-estimated saturated model, ensuring
  fit-vs-edges consistency under missingness. `"fiml"` errors clearly for PCA,
  EFA, and WLSMV/ULSMV estimators. (#M16)

* Fixed a latent inconsistency in the ESEM continuous (ML/MLR) path: the
  engine fit used lavaan's default listwise deletion while edge correlations
  were computed from a separately computed pairwise matrix. With
  `missing = "listwise"`, the data are pre-filtered so both use identical rows.
  With `missing = "fiml"`, the edge R is extracted from lavaan's FIML saturated
  model rather than a pairwise `stats::cor()`. (#M16)

* For ESEM with `cor = "polychoric"` (WLSMV/ULSMV), `missing = "pairwise"` now
  passes `missing = "available.cases"` to `lavaan::efa()` rather than lavaan's
  listwise default. This uses every row that contributes to each pair of
  polychoric thresholds (MCAR-valid) and makes `n_obs` in the result reflect the
  full sample size, not just complete cases. (#M16)

* `ackwards()` result objects now carry `$meta$missing` (the `missing` argument
  value used) and `$meta$n_complete` (the number of complete cases in the
  original data, regardless of deletion mode). (#M16)

## Milestone 15 — Naming clarity & consistency pass

* `ackwards()` argument `k` renamed to `k_max`. Aligns with the pre-existing
  `suggest_k(k_max=)` parameter and `$k_max` result-object field, making the
  maximum-depth semantic explicit and consistent across the full API. (#M15)

* `ackwards()` argument `method` renamed to `engine`. Resolves the overload with
  `compute_edges()`'s internal `edge_method` argument and matches DESIGN.md's
  "three engines" prose (`"pca"`, `"efa"`, `"esem"`). (#M15)

* `ackwards()` argument `scores` renamed to `keep_scores`. Mirrors the existing
  `keep_fits` argument and removes ambiguity between "store scores" (the storage
  option) and "factor scores" (the output concept). (#M15)

* `ackwards()` argument `align` renamed to `align_signs`. Clarifies that the
  argument controls sign orientation of loadings, not structural alignment or
  rotation. (#M15)

* Result-object field `$method` renamed to `$engine`. Result-object field
  `$cor_type` renamed to `$cor`. All S3 methods (`print`, `summary`, `glance`,
  `tidy`, `augment`, `autoplot`) updated to read from the new field names. (#M15)

* `compute_edges()` argument `method` renamed to `edge_method`. Disambiguates
  from the (now-removed) `method` engine argument on `ackwards()`. (#M15)

* No behaviour changes. All existing functionality is preserved; only surface
  names changed.

## Milestone 14 — Dedicated `suggest_k()` vignette

* New vignette `vignette("ackwards-suggest-k")` — "Choosing k: How Many
  Factors?" — a narrative/educational companion to the `suggest_k()` reference
  page. Covers all five criteria with pros/cons and bias direction, engine-to-
  criterion pairing (PA-FA for efa/esem; PA-PC for pca), all four arguments
  including the ordinal→Pearson caveat and the PA non-reproducibility note, and
  a worked k-selection for the BFI. Listed first in the pkgdown "Deep dives"
  section. (#M14)

* Intro vignette (`vignette("ackwards-intro")`) Step 1 trimmed to the default
  `suggest_k()` call plus a pointer to the new vignette; criteria prose moved
  there. (#M14)

* README corrected: `suggest_k()` description updated from "parallel analysis
  and Velicer's MAP criterion" to the accurate five-criteria description;
  "Both criteria agree" prose updated to "The criteria converge"; new vignette
  added to the "Learn more" table. (#M14)

* Vignette worked-example prose corrected: the initial version mis-stated
  PA-PC and PA-FA recommendations and assumed CD was absent. The section now
  reads from the actual rendered output (PA-PC=5, PA-FA=6, MAP=5, VSS-1=4,
  VSS-2=5, CD=8 on `bfi`), explains the unusual PA-FA > PA-PC relationship,
  and documents why CD=8 is a high outlier on large correlated samples rather
  than a directive to extract 8 factors. (#M14)

* CD summary table: "Conservative, accurate" changed to "Accurate in
  simulation; can over-retain on large, correlated samples". (#M14)

* Tests added for three previously-uncovered branches in `autoplot.suggest_k()`
  and `print.suggest_k()`: `k_parallel_fa = NA` (no FA star rendered),
  `cd_available = FALSE` (no CD vline; "requires EFAtools" note in print). (#M14)

## Milestone 13 — Rotation honesty: remove dead `kappa` argument

* The `kappa` argument to `ackwards()` has been removed. It was accepted and
  stored in `x$meta` but was never passed to any of the three engines — all
  engines already hardcoded `cfT → "varimax"`. Crawford-Ferguson with κ = 1/p
  is numerically identical to varimax (Crawford & Ferguson, 1970; Browne, 2001),
  and no reference paper (Goldberg 2006; Kim & Eaton 2015; Forbush et al. 2024)
  varies kappa. Exposing it implied quartimax/equamax were reasonable alternatives
  for this method; they are not. (#M13)

* The `rotation` argument to `ackwards()` has been removed. Varimax is the only
  valid rotation for this method, so exposing it as a choice was misleading.
  "cfT" (Crawford-Ferguson family name) has been renamed to "varimax" everywhere:
  the result object's `$rotation` field, `print()` output, and all engine
  internals. The `@section Defaults` explains why orthogonal rotation is
  load-bearing (T'T = I so T' = T^-1, enabling closed-form W'RW algebra) and
  that varimax is what Goldberg (2006), Kim & Eaton (2015), and Forbush et al.
  (2024) all used. (#M13)

## Milestone 12 — Best-practice `suggest_k()` expansion

* `suggest_k()` now runs five complementary criteria instead of two. New
  criteria: **PA-FA** (factor-eigenvalue parallel analysis, the model-consistent
  companion to PA-PC), **VSS-1/VSS-2** (Very Simple Structure; already computed
  from `psych::vss()` but now surfaced), and **CD** (Comparison Data, Ruscio &
  Roche 2012, via `EFAtools::CD()` — skipped gracefully when `EFAtools` is not
  installed). (#M12)

* The `suggest_k` object has been enriched with new fields: `k_parallel_pc`,
  `k_parallel_fa` (replaces `k_parallel`), `k_vss1`, `k_vss2`, `k_cd`,
  `cd_available`; and the `criteria` table now includes `ev_obs`, `pa_pc_quant`,
  `pa_pc_suggested`, `pa_fa_quant`, `pa_fa_suggested`, `vss1`, `vss2`. (#M12)

* New `seed` argument to `suggest_k()` for reproducibility of the stochastic
  parallel-analysis and Comparison Data steps. (#M12)

* New `autoplot.suggest_k()` method produces a three-panel ggplot2 diagnostic:
  a parallel-analysis/scree plot, a MAP panel, and a VSS panel. Optimal k for
  each criterion is marked with a star-shaped point. (#M12)

* `print.suggest_k()` redesigned as a multi-criterion table showing PA-PC,
  PA-FA, MAP, VSS-1, VSS-2, and CD (when available) side by side, with a
  consensus range summary and overextraction caution. (#M12)

* `EFAtools` added to `Suggests` in DESCRIPTION. Never imported directly;
  gated by `rlang::is_installed()`. (#M12)

## Milestone 11 — Edge-label polish + `show_r` decoupling

* `show_r` now defaults to `FALSE` in all `autoplot.ackwards()` views.
  Previously, `show_r = NULL` resolved to `TRUE` when `drop_pruned = TRUE`
  (coupling two orthogonal concerns). The default is now explicit: pass
  `show_r = TRUE` to label edges. The Forbes vignette demonstrates both the
  labeled and unlabeled variants, matching the two-figure treatment in Forbes
  (2023). (#M11)

* Edge correlation labels now use APA-style formatting via the new internal
  helper `.format_r()`: the leading zero is stripped (`.23` not `0.23`),
  trailing zeros are padded to `r_digits` decimal places (`.30` not `.3`),
  and `±1` formats as `1.00` / `-1.00`. (#M11)

* Edge labels are now rendered with `geom_label` (white background, no border)
  offset perpendicular to each edge's direction, so they clear near-vertical
  arrows and the arrowhead regardless of edge angle. The previous flat
  horizontal nudge was ineffective for near-vertical segments. (#M11)

* New `r_label_size` argument to `autoplot.ackwards()` controls the font size
  of edge correlation labels when `show_r = TRUE`. Default `2.5` preserves
  existing behaviour. (#M11)

* `.format_r()` no longer produces `"-.00"` when a small negative correlation
  rounds to zero at the requested precision; sign is suppressed in that case.
  (#M11)

## Milestone 10 — Conformance + robustness

### Wave 2 follow-up: conformance fixes

* Fixed `print.summary_ackwards()` spuriously printing "Redundant: 0 nodes
  flagged" when only `prune = "artefact"` was requested. Root cause:
  `.summary_prune()` always returned a `redundant` field (a length-0 character
  vector), and `print` gated on `!is.null(redundant)` rather than on the
  `rules` slot. Fix: `.summary_prune()` now carries `rules` through from the
  prune slot; both pruning display blocks in `print` gate on `"redundant" %in%
  p$rules` / `"artefact" %in% p$rules` respectively — consistent with how
  `print.ackwards()` handles the same slot.
* `.summary_lineage()` now uses `which(e$is_primary & ...)` instead of direct
  logical indexing to guard against `NA` values in `is_primary` (skip-level
  edges from `pairs = "all"` carry `NA` there). Parents are now sorted
  explicitly by level then label rather than relying on tidy-table insertion
  order, giving stable top-down ordering regardless of pairs mode.
* 7 new tests added to `test-print.R` covering: artefact-only prune (no
  spurious redundant line), redundant+artefact combined, truncated ESEM
  hierarchy, polychoric basis, PCA eigenvalue presence in fit table, lineage
  content correctness (specific parent→children check), and the print method.
* Heywood fixture in `test-esem.R` now starts with a `skip_if` pre-check
  that verifies the fixture still produces `theta <= 0` under the installed
  lavaan version, making the test skip gracefully rather than silently pass
  if a future lavaan changes its residual-variance bounding behaviour.

### Wave 2: summary method

* Added `summary.ackwards()` and `print.summary_ackwards()` — the previously
  documented but unimplemented summary method (DESIGN.md §6/§10). Returns a
  structured `summary_ackwards` object that, when printed, shows:
  - **Per-level block**: per-factor variance % and cumulative variance for every
    level. PCA levels additionally show eigenvalues; EFA/ESEM levels (k ≥ 2)
    append a fit-index line (RMSEA/TLI/chi/df for EFA; CFI/TLI/RMSEA/SRMR for
    ESEM).
  - **Lineage list**: `m1f1 → m2f1, m2f2` readable chains walking the
    adjacent primary-parent edges top-down.
  - **Pruning section**: when `prune != "none"`, lists flagged-node IDs per
    rule and notes that pruning is interpretive, not re-estimation.
  - Honesty caveat (series of linked solutions, not a fitted hierarchy).
  The returned object carries `variance`, `fit`, `lineage`, and `prune` fields
  for programmatic access, matching the `print.suggest_k` pattern.
* Added `summary.ackwards` and `print.summary_ackwards` to `_pkgdown.yml`
  reference and mentioned `summary()` in the intro vignette Step 3 section.

### Wave 1: engine robustness

* **ESEM improper-solution warning.** The ESEM engine now warns when any residual
  variance in the lavaan theta matrix is at or below zero (Heywood case). lavaan
  clamps negative residual variances to 0 by default; the check catches both
  clamped-to-zero and unconstrained-negative values. The object is still returned
  at the full requested depth — the warning is diagnostic, not truncating
  (Invariant 7). Parity with the EFA engine, which has warned on Heywood cases
  since M3.
* **`cor = "spearman"` + `method = "esem"` inconsistency warning.** `ackwards()`
  now warns when this combination is requested: lavaan fits a Pearson-ML model on
  raw data while `compute_edges()` uses the Spearman correlation matrix for
  scoring, mixing bases silently. The warning suggests `cor = "polychoric"` for
  ordinal data or `cor = "pearson"` for a fully consistent continuous path.
  Emitted once per session (`.frequency = "once"`).
* **DESIGN.md §8 reconciled.** `suggest_k()` documentation in §8 now lists only
  the two implemented criteria (parallel analysis and MAP/Velicer) and explicitly
  marks Empirical Kaiser Criterion (EKC) and EGA (`{EGAnet}`) as out of scope,
  consistent with the §12/§14 no-extra-deps decision.

## Milestone 9 — Visualization round 2 + vignette restructure

### Wave 2: vignette restructure

* New `vignette("ackwards-visualization")`: a dedicated visual reference for
  all `autoplot.ackwards()` presentation options, with rendered before/after
  figures for each argument. Covers thresholds, colours, monochrome mode,
  correlation labels, custom node labels, primary-only mode, level labels,
  arrowheads, edge width, and legend; closes with a publication-ready worked
  example.
* `vignette("ackwards-forbes")`: removed the curved-arc `autoplot(x_all)`
  figure (skip-level edges are better understood via the edge table); the
  annotated grey-box view is replaced by the Forbes-styled `drop_pruned`
  diagram as the primary figure, using `color_pos = color_neg = "black"`,
  `edge_linewidth = 0.6`, `show_arrows = FALSE`, `legend = FALSE`. Added a
  cross-reference to the visualization vignette for cosmetic options.
* `vignette("ackwards-intro")`: trimmed the "Adjusting the diagram" section to
  prose pointing to `vignette("ackwards-visualization")`; fixed a stale code
  comment that promised to "relabel the k = 5 factors" when the code did not.
  Added the visualization vignette to the "Next steps" table.

### Wave 1: new autoplot args

* `autoplot.ackwards()` gains `show_arrows = TRUE`: when `FALSE`, edges are
  drawn with plain line ends (`arrow = NULL`) instead of closed arrowheads.
  Applies to both straight and curved edge layers. Forbes (2023) figures use
  plain line ends throughout.
* `autoplot.ackwards()` gains `edge_linewidth = NULL`: `NULL` (default) keeps
  the existing behaviour of mapping `|r|` to `linewidth` via a continuous scale.
  A numeric value draws every edge at that constant width, removes the
  `linewidth` aesthetic mapping, and drops the `|r|` linewidth legend. Forbes
  figures use uniform thin lines (≈ 0.5–0.6).
* `autoplot.ackwards()` gains `legend = TRUE`: when `FALSE`, suppresses all
  plot legends (`legend.position = "none"`). Useful with `color_pos = color_neg
  = "black"` to remove the otherwise redundant "Direction" key.
* **Forbes (2023) publication style** is now reproducible without a new mode.
  Combine existing colour-mode args with the three new args:
  `autoplot(xp, drop_pruned = TRUE, color_pos = "black", color_neg = "black",`
  `         edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE)`.
  The colour mode's `cut_strong`-based solid/dashed distinction matches
  Forbes's weak/secondary dashing; `mono` is not needed and uses different
  linetype semantics (sign, not strength).

## Milestone 8 — Plot customization (Waves 1 & 2)

* `autoplot.ackwards()` gains `show_r = FALSE` / `r_digits = 2L`: when
  `show_r = TRUE`, draws `round(r, r_digits)` as a text label at each edge
  midpoint. Pairs naturally with `mono = TRUE` for greyscale journal figures.
* `autoplot.ackwards()` gains `mono = FALSE`: monochrome mode draws all edges
  in black with `linewidth` encoding `|r|` and `linetype` encoding sign
  (`solid` = positive, `dashed` = negative). Drops the `cut_strong`
  strong/weak linetype distinction (linewidth already conveys magnitude, so a
  second strength encoding on linetype is redundant).
* `autoplot.ackwards()` gains `show_level_labels = TRUE` / `level_label_size = 3`:
  draws "1 factor", "2 factors", … to the left of the diagram at each level's
  `y` position (using `coord_cartesian(clip = "off")` for overflow-safe placement,
  which was already set). Default is `TRUE` (was missing before).
* `autoplot.ackwards()` gains `node_labels = NULL`: a named character vector
  mapping factor IDs (e.g. `"m5f1"`) to custom display strings. Unspecified
  factors keep their `m{k}f{j}` label. Warns if a supplied name matches no
  factor ID in the object.
* `autoplot.ackwards()` gains `primary_only = FALSE`: when `TRUE`, filters
  the edge table to `is_primary == TRUE` before drawing, producing a clean
  tree of primary parent links. Because skip-level edges are never primary,
  this also suppresses curved skip arcs.
* `autoplot.ackwards()` gains `drop_pruned = FALSE`: activates the Forbes
  (2023) pruned-view rendering path. Pruned nodes are removed from the
  diagram entirely; each retained node is connected to its single strongest
  kept ancestor by a straight arrow, even when the arrow spans multiple
  level gaps (no curved arcs in this mode). Requires `prune != "none"` at
  fit time; errors with a clear message if pruning annotations are absent.
  `show_r` auto-defaults to `TRUE` in this mode (Forbes convention).
* `autoplot.ackwards()` gains `compress_levels = FALSE`: used together with
  `drop_pruned = TRUE`, closes the vertical gaps left by pruned levels so
  retained nodes are evenly spaced; level axis labels still show the original
  level numbers so each row remains identifiable.
* Added private helper `.drop_pruned_nodes(x, nodes, compress_levels)` in
  `layout.R`: computes the kept-only node set and a reduced primary-edge
  table for the drop-pruned rendering path. For each kept node it selects the
  single strongest `|r|` edge to any kept shallower node, recomputing primary
  parentage on the reduced graph (the original `is_primary` column is not
  reused). When `compress_levels = TRUE`, re-indexes `y = -rank(level)` while
  preserving the original `level` column for display labels.
* `show_r` default changed from `FALSE` to `NULL` to support the
  `drop_pruned` auto-default; calling `autoplot(x)` without `show_r`
  produces the same behaviour as before (`FALSE`).

## M6/M7 conformance fixes (second-pass review)

* `augment(x, data=)` and `scores = TRUE` now warn when `data` contains rows
  with missing values: score projection applies weights row-wise and propagates
  NAs listwise, even though fitting used pairwise-complete correlations. The
  warning advises `na.omit(data)` to remove incomplete rows before scoring.
  `?augment.ackwards` documents the listwise-vs-pairwise asymmetry.
* `augment(x, data=)` now errors with a clear message when `data` contains
  non-numeric columns (previously would error inside `scale()` with a
  cryptic message).
* `detect_ordinal()` is now called once in `ackwards()` and its result reused
  for both the Invariant-6 console warning and `meta$ordinal_warned` (was
  called twice, scanning all columns each time).
* The non-Pearson basis warning in `.compute_scores()` now uses
  `.frequency = "once"` so repeated `augment()` calls on the same object
  do not produce duplicate warnings in a session.
* Tests added: non-Pearson basis warns on scoring (Fix 3); `keep_fits = TRUE`
  only stores converged levels under truncation (Fix 4, labelled B4);
  NA-data warns and produces NA scores for `augment()` and `scores = TRUE`
  (Fix 1); non-numeric `augment()` data errors cleanly (Fix 5).
* DESIGN.md §12 and CLAUDE.md dependency section updated to match DESCRIPTION:
  `Imports` is `cli`, `generics`, `rlang`, `stats`, `utils` (`methods` is not
  imported; no `methods::` usage exists). `Suggests` reflects the actual set
  (`psych`, `GPArotation`, `lavaan`, `ggplot2`, `testthat`, `knitr`,
  `rmarkdown`, `covr`); `ggraph`/`igraph`/`tidygraph`, `EGAnet`/`paran`,
  `clue`, `future` are not in DESCRIPTION and the §12 table no longer lists
  them. `suggest_k()` uses `psych::fa.parallel` + `psych::vss` (not separate
  `EGAnet`/`paran` deps).

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
