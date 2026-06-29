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

## Factor interpretation

* **`top_items(x, level, cut, n, sort)`** — lists the salient items for each factor,
  filtered to `|loading| >= cut` (default `0.3`) and sorted by descending absolute loading.
  Arguments `level`, `n`, and `sort` allow subsetting levels, capping items per factor, and
  controlling order. Returns a `top_items` S3 object with a grouped cli print method.

* **`label_template(x, style)`** — generates the named-character-vector scaffold for
  `autoplot(x, node_labels = ...)`. Styles: `"id"` (default), `"forbes"` (level-letter +
  within-level index: `"A1"`, `"B1"`, `"B2"`, …), `"blank"`. Factor IDs are in canonical
  layout order; an editable `c(...)` literal is also printed for copy-paste.

* See `vignette("ackwards-interpret")` for the end-to-end naming workflow: reading a factor
  with `top_items()`, hierarchy-aware naming (parent vs. child, blends), sign-alignment
  caveat, and the `top_items()` → `label_template()` → `autoplot(node_labels = ...)` round-trip.

## Tidy interface and scoring

* `print.ackwards()` — compact cli summary card.
* `summary.ackwards()` — per-level variance, fit indices, lineage list, pruning notes.
* `tidy(x, what = "edges" / "loadings" / "variance" / "nodes" / "scores")` — long-format
  tidy tibbles.
* `glance(x)` — one-row model-level summary.
* `augment(x, data = ...)` — appends factor scores to a data frame; recomputes from stored
  weights when scores were not materialized at fit time.

