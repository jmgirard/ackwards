# ackwards (development)

## Vignette corrections: `ackwards-intro` and `ackwards-suggest-k`

`vignette("ackwards-intro")` hardcoded a cumulative-variance jump ("22.9% →
34.7%") that had drifted from the code's actual output (23.2% → 35.5%); it is
now computed inline from the fitted object so it cannot drift again. The
lineage walkthrough had `m4f1`'s primary children backwards (claimed
`m5f2`/`m5f4`; the edge table shows `m5f1`/`m5f4`), and the diagram narrative
misattributed which traits differentiate at which level (Conscientiousness/
Openness split at k = 4, not Agreeableness/Extraversion, which split at
k = 5) — both corrected to match a live run. `vignette("ackwards-suggest-k")`
now clarifies why the printed "Recommendations" block shows six lines for
five criteria (`"vss"` reports both VSS-1 and VSS-2).

## Guard against `cor = "polychoric"` with an incompatible ESEM estimator

`ackwards(engine = "esem", cor = "polychoric", estimator = "ML")` (or
`"MLR"`) now errors immediately with a clear explanation instead of failing
many calls deep inside the per-level ESEM fit and surfacing a misleading
"failed to build at least 2 converged levels... check your data for
multicollinearity" abort. Polychoric correlations mark every item `ordered`
for lavaan, and lavaan itself does not support ML/MLR estimation on ordered
indicators. `"WLSMV"`/`"ULSMV"` with a continuous `cor` remains allowed (a
valid, if atypical, continuous WLS/ADF estimator).

## ESEM fit indices: honest p-value and BIC

Under `estimator = "WLSMV"`/`"ULSMV"`, lavaan's naive chi-square test has no
valid reference distribution for these limited-information estimators —
lavaan's own `summary()` reports that p-value as "Unknown" (`NA`). `chi`,
`dof`, and `p_value` now fall back to lavaan's mean-and-variance-adjusted
("scaled") test whenever the naive p-value is undefined, which does have a
genuine null distribution; ML/MLR (whose naive p-value is already valid) are
unaffected. A genuinely saturated level (`dof = 0`) still reports `NA` — there
is no test to perform on a model that fits perfectly by construction.

`BIC` is now a first-class ESEM fit index (previously silently absent):
populated under `"ML"`/`"MLR"` (a proper log-likelihood exists), and `NA`
under `"WLSMV"`/`"ULSMV"` (no proper log-likelihood for these estimators —
genuinely inapplicable, not a bug). `tidy(x, what = "fit")` and `glance(x)`
both surface it consistently across engines.

## `tidy(what = "fit", cutoffs = TRUE)`: no more always-NA `_meets` columns

`format = "wide"` no longer emits a `{index}_meets` column for indices with no
defined threshold (`chi`, `dof`, `p_value`, `BIC`, PCA eigenvalues) — the pivot
previously generated one for every index, and those columns were always `NA`.
`format = "long"` is unaffected (`meets` is still `NA` for those rows, as
documented).

## Citation hygiene

`citation("ackwards")` now returns a single software entry (Girard) instead
of also listing Goldberg (2006) as a package author — the method paper is
cited in `DESCRIPTION` and roxygen `@references`, not as a software
co-citation. `ackwards()`'s `@references` now also lists Forbes (2023)
alongside Goldberg (2006) and Waller (2007), since `ackwards()` implements
all three (original method, exact edge algebra, and the extended method).
The README's Citation section is updated to guide users to the right
paper(s) for what they used.

## Documentation: no internal milestone numbers

User-facing docs (`NEWS.md`, `README.md`, vignettes) no longer reference internal
development-milestone tags — they are meaningless outside this repo's own process.
A regression test now guards against reintroducing them.

## CD criterion: bug fix and honest framing

`suggest_k()` now correctly handles `EFAtools::CD`'s output. `EFAtools::CD`
fills its `RMSE_eigenvalues` matrix only up to the level it actually tested;
trailing columns were literal zeros, causing the plotted RMSE curve to dive
spuriously to zero at higher k and `which.min()` to land on a fake value.
The `cd_rmse` vector now has those unfilled positions set to `NA`, so the curve
terminates at the last genuinely computed level.

The CD plot panel is relabeled from `"CD (RMSE, minimize)"` to
`"CD (RMSE; sequential test)"` to reflect how `EFAtools::CD` actually works: a
sequential one-sided Wilcoxon test (default α = 0.30) that stops at the first
non-significant RMSE improvement. The starred k is the last retained factor and
need not be the visible minimum of the curve. The `autoplot.suggest_k()` and
`suggest_k()` documentation has been updated accordingly.

## Per-level fit indices and loading SEs as first-class output

**`glance()` now carries fit indices** from the deepest converged level. A
consistent five-column set (`CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`) is present
for all engines; indices unavailable for a given engine are `NA` (e.g. `CFI`
and `SRMR` for EFA; all five for PCA).

**`tidy(what = "fit")` gains two new arguments:**
- `format = "wide"`: returns one row per non-anchor level (k >= 2) with one
  column per index — the natural shape for reporting tables. The saturated
  1-factor anchor is dropped, matching `summary()` and `autoplot(what = "fit")`.
- `cutoffs = TRUE`: appends a `meets` column flagging each index against Hu &
  Bentler (1999) conventional thresholds (CFI/TLI >= .95, RMSEA <= .06,
  SRMR <= .08). Indices without a defined threshold (chi, BIC, eigenvalues)
  get `NA`. Thresholds are report-only and never gate any behaviour.

**`tidy(what = "loadings")` now includes `se`, `ci_lower`, and `ci_upper`
columns** for all engines. For ESEM these are populated from the rotation-aware
loading SEs; for PCA and EFA they are `NA`. A `conf_level` argument (default
`0.95`) controls the interval width. This replaces the need to compute CIs by
hand, and supersedes the previous `tidy(what = "loadings_se")` accessor, which
has been **removed** (the package is unreleased, so no deprecation cycle is
needed).

**`autoplot(x, what = "fit")`** produces a two-panel ggplot2 chart of per-level
fit indices (CFI/TLI in the top panel; RMSEA/SRMR in the bottom panel) with
Hu & Bentler (1999) reference lines. PCA objects return an informative empty
plot. `what = "hierarchy"` (default) is unchanged.

**`summary()` fit lines** now append a check (✔) or cross (✘) to each
index compared against its Hu & Bentler (1999) threshold.

**`print()` and `summary()` footer** extended: per-level fit indices describe
how well a k-factor model fits the items at that level — they do not validate
the edges or the hierarchy itself.

**New vignette section** in `vignette("ackwards-engines")`: "Per-level fit:
what it tells you (and what it doesn't)" — covers the level-vs-hierarchy
distinction, reporting with the new wide table and fit plot, the ESEM
cost/benefit tradeoff, and when to care about fit in exploratory vs publication
workflows.

## Faster ESEM on large item sets

The ESEM engine no longer recomputes lavaan's sample statistics (thresholds, the
polychoric correlation matrix, and the asymptotic weight matrix) at every level.
These depend only on the data, so they are now computed once at the first level and
reused for every deeper level via lavaan's `slotSampleStats=`. Solutions are
identical; the redundant recomputation — which dominated run time with hundreds of
items — is removed.

The independent per-level model fits can now also run in parallel. When the
optional `future.apply` package is installed, `ackwards()` dispatches the ESEM
levels through the `future` framework. The default plan is sequential (no behaviour
change); call `future::plan(future::multisession, workers = N)` before `ackwards()`
to parallelize. Results are reproducible across plans when `seed` is supplied.
`future.apply` added to `Suggests`. See the "Performance with many items" section of
`vignette("ackwards-engines")`.

## Tucker's φ default for non-PCA redundancy (⚠ resolved-default change)

`redundancy_phi = NULL` (the default) now **auto-resolves** based on the extraction engine:

* `engine = "pca"` — no φ filter (the closed-form W′RW algebra is exact; `|r|` alone is
  sufficient, as in the original Waller 2007 method).
* `engine = "efa"` or `"esem"` — automatically applies `redundancy_phi = 0.95` (Lorenzo-Seva
  & ten Berge, 2006). Factor-score indeterminacy off-PCA makes `|r|`-only redundancy too
  liberal; the conjunctive φ criterion is the conservative default.

A cli message is emitted whenever auto-resolution applies (Invariant 6: loud defaults). The
resolved value is stored in `x$prune$redundancy_phi`.

**Existing non-PCA calls using `prune = "redundant"` will become more conservative** — some
previously flagged redundant chains may no longer meet the φ criterion. To restore the old
behaviour, pass `redundancy_phi = NA` (explicit opt-out; no φ filter regardless of engine).
Explicit numeric values (e.g., `redundancy_phi = 0.8`) override on any engine.

## Structural artefact signals (Forbes extension)

`ackwards(..., prune = "artefact")` now computes structural signals alongside Tucker's φ
table. Two new arguments on `ackwards()`:

* `min_items = 3L` — minimum primary-item count per factor; factors with fewer items are
  flagged `few_items = TRUE`.
* `orphan_r = 0.5` — minimum adjacent-level `|r|` required for a factor to be considered
  "connected"; factors below this threshold across all adjacent levels are flagged
  `orphan = TRUE`.

A third signal, `split_merge`, is `TRUE` when a factor at level k draws its primary items
from multiple distinct primary factors at level k−1 (items that were in separate groups at
the shallower level merged into a single factor at the deeper level).

All signals are **flag-only** (artefact identification requires researcher judgment;
Forbes, 2023). Results are stored in `x$prune$structural` (one row per factor × level);
`print()` and `summary()` report the count of flagged factors. Both `min_items` (positive
integer) and `orphan_r` (number in `[0, 1]`) are validated with a clear error. See
`vignette("ackwards-forbes")` for a worked example.

## `suggest_k()` criterion selection

`suggest_k()` gains a `criteria` argument that controls which retention criteria are
computed. Any subset of `c("pa_pc", "pa_fa", "map", "vss", "cd")` may be requested;
the default runs all five (identical to prior behaviour). Criteria not requested are
skipped entirely (no computation), delivering a real speed win:

* `"pa_pc"` / `"pa_fa"` share one `psych::fa.parallel()` call — both skipped if neither
  is requested.
* `"map"` / `"vss"` share one `psych::vss()` call — both skipped if neither is
  requested.
* `"cd"` runs `EFAtools::CD()` only when explicitly requested **and** `EFAtools` is
  installed; if requested but unavailable an informational message is emitted.

`k_*` fields for non-requested criteria are `NA_integer_`; the `criteria` data frame
uses `NA` for non-run columns (stable schema). The returned object gains a
`criteria_requested` field. `print()` and `autoplot()` render only the requested
criteria; the consensus range is computed from requested criteria only.

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

`autoplot.suggest_k()` produces a ggplot2 diagnostic: a four-panel 2×2 grid (scree/PA,
MAP, VSS, CD RMSE) when `EFAtools` is installed, or a three-panel single-column layout
(scree/PA, MAP, VSS) otherwise. The returned object also carries `cd_rmse` (column means
of CD's RMSE eigenvalue matrix) for use in the CD panel.

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

## Documentation

Vignette comparison tables reworked for legibility:

* `ackwards-engines` and `ackwards-ordinal` — stacked long-format `kable` tables replaced
  by **wide gt tables**: one row per item/edge, one column per engine/basis, plus an explicit
  Δ column that is the teaching point (EFA/polychoric attenuation). Factor and sign
  correspondence is asserted before differencing; edge tables surface primary-parent
  disagreements as NA rather than hiding them.
* `ackwards-forbes` — `prune-nodes` raw `tidy()` print replaced by a styled gt table;
  narrated counts (`n_adj`, `n_all`, `n_redundant`, top skip-edge) are now computed via
  inline R expressions so they can never go stale.
* `gt` added to Suggests (vignette table formatting only; never loaded by core).

## Test coverage

Test suite reaches **100% line coverage**. Genuinely unreachable defensive branches
(engine fallbacks, lavaan error/warning handlers) are excluded via `# nocov` rather than
brittle environment-dependent tests.

