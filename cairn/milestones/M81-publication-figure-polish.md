# M81: Publication-figure polish — item lists, per-node box sizes, manual factor ordering

- **Status:** planned
- **Priority:** normal
- **Depends on:** M80
- **Driving RR:** —
- **Principles touched:** GP5 (lean install — ggplot2 stays Suggests); works under D-015 (light `ba_layout`, layered DAG not a tree) and D-029 (persistent factor labels)

## Goal

Add three opt-in `autoplot()`/`ba_layout()` controls for publication figures: item lists under the deepest-level factors, per-node box sizes, and manual left-to-right factor ordering.

## Scope

**In:**
- (a) `autoplot(x, show_items = TRUE)` draws the top-`n_items` observed items (by `|loading|`) beneath each **deepest-level (`k_max`)** factor box, reusing `top_items()` extraction (including variable labels when present); default `show_items = FALSE` leaves output byte-identical.
- (b) `node_width`/`node_height` accept a **named numeric vector** keyed by factor id (per-box override); a bare scalar reproduces current output.
- (c) `ba_layout(x, order = <deepest-level permutation>)` fixes the deepest level's left-to-right order; upper levels re-center above their primary children (D-015 intact); `autoplot(..., order =)` passes through.

**Out:**
- Independently pinning an **upper-level** factor's position (breaks centered-above-children) → candidate row if Forbes later needs it; deepest-level order already reaches any primary-forest arrangement.
- Auto-size-to-text box sizing → rejected at plan gate (nchar is not reliable text extent); per-node override is the mechanism.
- Item lists under levels other than `k_max` → out; deepest-only per Forbes's "lowest-level factors".

## Acceptance criteria

- [ ] AC1: `autoplot(x, show_items = TRUE)` adds a text layer listing the top-`n_items` items (by `|loading|`, `top_items()` semantics; variable labels when the fit carried them) under each `k_max` box, in both `direction` values; `show_items = FALSE` (default) produces a plot with no such layer. Verified by inspecting the built plot's layer data.
- [ ] AC2: named-vector `node_width`/`node_height` render each named box at its override size and unnamed boxes at the scalar default (per-row `width`/`height` in the tile layer); a bare scalar is unchanged from current output. Verified by tile-layer data assertions.
- [ ] AC3: `ba_layout(x, order = p)` yields deepest-level x-coordinates whose left-to-right rank equals `p`, with upper-level nodes at the mean-x of their primary children; a non-permutation or wrong-id `order` errors; an `order` entry for a non-deepest level warns and is ignored. Verified by coordinate assertions + error/warning tests.
- [ ] AC4: `autoplot(x, order = p)` passes `order` through so the rendered node x-coordinates reflect the manual order. Verified by node-coordinate assertion on the built plot.
- [ ] AC5: `ackwards-visualization` vignette + roxygen `@examples` demonstrate all three features; NEWS.md entry added; `precompute.R` regenerated and vignette-freshness clean; no new exports (pkgdown reference unchanged).
- [ ] AC6: `Rscript tools/dod-gate.R` clean (check 0/0/0, coverage, style, lint, pkgdown).

## Coverage

- AC1 → T3
- AC2 → T2
- AC3 → T1
- AC4 → T1
- AC5 → T4
- AC6 → T1, T2, T3, T4

## Tasks

- [ ] T1 — Manual ordering (`ba_layout`, `R/layout.R:45`): add `order =` (character vector or list keyed by level; only the `k_max` entry honored). Validate it is a permutation of the deepest level's ids; error otherwise; warn+ignore non-deepest entries. When supplied, set the deepest ordinal to the manual order and run `.assign_x` once (skip the two-candidate keep-best, which collapses when the deepest order is fixed). Add `order =` to `autoplot.ackwards()` and pass through. Tests + roxygen (`@param order`, `@examples`).
- [ ] T2 — Per-node box sizes (`R/autoplot.R`): let `node_width`/`node_height` be a scalar **or** a named numeric vector; resolve to a per-node `width`/`height` column consumed by `geom_tile` in `.ba_add_nodes()`; make the `min_sep < node_width` warning vector-safe (compare against `max`). Tests + roxygen update.
- [ ] T3 — Item lists (`R/autoplot.R`): add `show_items = FALSE`, `n_items = 5`, `item_cut = 0.3`. When on, extract deepest-level items via the shared `top_items()` path and draw a stacked text layer beyond the `k_max` boxes (below when vertical, right when horizontal); extend the plot margin. Tests + roxygen + one `@example`.
- [ ] T4 — Docs: NEWS.md entry; update `vignettes/ackwards-visualization.Rmd.orig` to show all three; `Rscript vignettes/precompute.R`; commit regenerated `.Rmd` + `assets/` (revert unrelated timing churn per M61/M75); confirm `tools/check-vignette-freshness.R` clean. Run `Rscript tools/dod-gate.R` at the gate.

## Work log

- 2026-07-24: created by /milestone-plan.
