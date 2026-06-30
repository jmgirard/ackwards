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

## Dependencies

* **`psych` moved to `Imports`** — the default PCA and EFA engines require psych; placing it in
  Suggests meant a mandatory install prompt for core functionality. The SEM (`lavaan`) and
  plotting (`ggplot2`) stacks remain in Suggests.
* **`GPArotation` removed** — varimax rotation routes through base `stats::varimax`; GPArotation
  was never loaded on any supported path.
* **`bfi25` example dataset** bundled — a 1 000-row, 25-item subset of the SAPA/IPIP Big Five
  data (sampled from `psych::bfi`). Used throughout examples and vignettes so they run without
  reaching into psych's namespace. See `?bfi25` for provenance and `@source`.

## `autoplot.suggest_k()` — CD panel

`suggest_k()` now stores `cd_rmse` (column means of `EFAtools::CD()`'s RMSE eigenvalue matrix)
in the returned object. When CD is available, `autoplot(suggest_k(...))` renders a four-panel
2×2 grid including a dedicated **"CD (RMSE, minimize)"** panel with the retention threshold
marked — giving CD its own diagnostic space rather than sharing the MAP panel. The three-panel
single-column layout is unchanged when EFAtools is absent.

## Correlation-matrix input

`ackwards()` and `suggest_k()` now accept a pre-computed **correlation matrix** in place of
raw item data, detected automatically from the matrix shape (square, symmetric, unit diagonal).

* **Engine gating** — only `"pca"` and `"efa"` are supported; `"esem"` errors clearly (lavaan
  requires raw data for WLSMV estimation and per-level fit indices).
* **`n_obs` argument** — new optional argument on both functions. Required for `engine = "efa"`
  with a matrix (psych needs N for chi-square/RMSEA/TLI); optional for `engine = "pca"` (stored
  as `NA` when omitted). Ignored (with a warning) when raw data are supplied.
* **Edge correctness** — edges from `ackwards(R, ...)` match `ackwards(data, ...)` exactly for
  the same correlation matrix: both routes use the closed-form `W'RW` algebra.
* **`cor` and `missing` arguments** — ignored (with a warning if set) for matrix input. `$cor`
  is stored as `NA_character_` and printed as `"(user-supplied matrix)"`.
* **Score paths blocked** — `keep_scores = TRUE`, `augment()`, and `tidy(what = "scores")`
  error clearly: individual-level scoring requires row-level item responses.
* **CD skipped in `suggest_k()`** — Comparison Data resamples raw item distributions and is
  skipped with an info note when a matrix is supplied.
* **Validation** — square, numeric, finite, symmetric, unit diagonal, `|r| <= 1`, no NA;
  synthesises `V1..Vp` dimnames when absent; warns (does not auto-smooth) when not
  positive-definite; errors with a clear message when a covariance matrix is detected.

## Tidy interface and scoring

* `print.ackwards()` — compact cli summary card.
* `summary.ackwards()` — per-level variance, fit indices, lineage list, pruning notes.
* `tidy(x, what = "edges" / "loadings" / "loadings_se" / "variance" / "nodes" / "scores")` —
  long-format tidy tibbles. For `what = "edges"`, `primary_only = TRUE` returns just each
  factor's primary-parent edge (the lineage tree) and `sort = "strength"` orders by `|r|`.
  `what = "loadings_se"` returns the rotation-aware loading standard errors (ESEM only).
* `glance(x)` — one-row model-level summary.
* `augment(x, data = ...)` — appends factor scores to a data frame; recomputes from stored
  weights when scores were not materialized at fit time.

## Test coverage

Test suite reaches **100% line coverage**. Genuinely unreachable defensive branches
(engine fallbacks, lavaan error/warning handlers) are excluded via `# nocov` rather than
brittle environment-dependent tests.

