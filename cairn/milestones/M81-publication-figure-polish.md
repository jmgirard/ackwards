# M81: Publication-figure polish — item lists, per-node box sizes, manual factor ordering

- **Status:** review
- **Priority:** normal
- **Depends on:** M80
- **Driving RR:** —
- **Principles touched:** GP5 (lean install — ggplot2 stays Suggests); works under D-015 (light `ba_layout`, layered DAG not a tree) and D-029 (persistent factor labels)
- **Branch/PR:** m81-publication-figure-polish / https://github.com/jmgirard/ackwards/pull/86

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

- [x] AC1: `autoplot(x, show_items = TRUE)` adds a text layer listing the top-`n_items` items (by `|loading|`, `top_items()` semantics; variable labels when the fit carried them) under each `k_max` box, in both `direction` values; `show_items = FALSE` (default) produces a plot with no such layer. Verified by inspecting the built plot's layer data.
- [x] AC2: named-vector `node_width`/`node_height` render each named box at its override size and unnamed boxes at the scalar default (per-row `width`/`height` in the tile layer); a bare scalar is unchanged from current output. Verified by tile-layer data assertions.
- [x] AC3: `ba_layout(x, order = p)` yields deepest-level x-coordinates whose left-to-right rank equals `p`, with upper-level nodes at the mean-x of their primary children; a non-permutation or wrong-id `order` errors; an `order` entry for a non-deepest level warns and is ignored. Verified by coordinate assertions + error/warning tests.
- [x] AC4: `autoplot(x, order = p)` passes `order` through so the rendered node x-coordinates reflect the manual order. Verified by node-coordinate assertion on the built plot.
- [x] AC5: `ackwards-visualization` vignette + roxygen `@examples` demonstrate all three features; NEWS.md entry added; `precompute.R` regenerated and vignette-freshness clean; no new exports (pkgdown reference unchanged).
- [x] AC6: `Rscript tools/dod-gate.R` clean (check 0/0/0, coverage, style, lint, pkgdown).

## Coverage

- AC1 → T3
- AC2 → T2
- AC3 → T1
- AC4 → T1
- AC5 → T4
- AC6 → T1, T2, T3, T4

## Tasks

- [x] T1 — Manual ordering (`ba_layout`, `R/layout.R:45`): add `order =` (character vector or list keyed by level; only the `k_max` entry honored). Validate it is a permutation of the deepest level's ids; error otherwise; warn+ignore non-deepest entries. When supplied, set the deepest ordinal to the manual order and run `.assign_x` once (skip the two-candidate keep-best, which collapses when the deepest order is fixed). Add `order =` to `autoplot.ackwards()` and pass through. Tests + roxygen (`@param order`, `@examples`).
- [x] T2 — Per-node box sizes (`R/autoplot.R`): let `node_width`/`node_height` be a scalar **or** a named numeric vector; resolve to a per-node `width`/`height` column consumed by `geom_tile` in `.ba_add_nodes()`; make the `min_sep < node_width` warning vector-safe (compare against `max`). Tests + roxygen update.
- [x] T3 — Item lists (`R/autoplot.R`): add `show_items = FALSE`, `n_items = 5`, `item_cut = 0.3`. When on, extract deepest-level items via the shared `top_items()` path and draw a stacked text layer beyond the `k_max` boxes (below when vertical, right when horizontal); extend the plot margin. Tests + roxygen + one `@example`.
- [x] T4 — Docs: NEWS.md entry; update `vignettes/ackwards-visualization.Rmd.orig` to show all three; `Rscript vignettes/precompute.R`; commit regenerated `.Rmd` + `assets/` (revert unrelated timing churn per M61/M75); confirm `tools/check-vignette-freshness.R` clean. Run `Rscript tools/dod-gate.R` at the gate.

## Work log

- 2026-07-24: created by /milestone-plan.
- 2026-07-24: T1 done — `order=` on `ba_layout()` (deepest-level permutation, list-keyed or bare vector; validates, warns+ignores non-deepest, errors on non-permutation) + `autoplot()` pass-through. `.resolve_manual_order()` helper; 6 tests. Layout suite 330 pass.
- 2026-07-24: T2 done — `node_width`/`node_height` accept a scalar or named-vector (per-box override, unnamed fall back to 0.8/0.4). `.resolve_node_size()` helper; per-node `nw`/`nh` mapped in `geom_tile`; edge faces + overlap warning vector-safe. 4 tests; 340 pass.
- 2026-07-24: T3 done — `show_items`/`n_items`/`item_cut` list the top deepest-level items (via `top_items()`, labels when present) below boxes (vertical) or right (horizontal); `.ba_item_text()` helper, margin widened. 4 tests; 346 pass. Visually verified item + box figures.
- 2026-07-24: T4 done — NEWS entry; visualization vignette gains order/box-size/item sections (`.orig` edited, precompute regenerated, unrelated churn reverted, freshness OK). `.ba_item_text` vectorised for 100% coverage. DoD gate clean: check 0/0/0, coverage 100%, lint clean, pkgdown complete. All tasks done → status review.
- 2026-07-24: review — fresh evidence gathered, consistency gate + full DoD gate re-run clean, three-lens fan-out returned zero findings. PR #86.

## Review

**Acceptance criteria — fresh evidence (2026-07-24):**

- AC1 ✓ — `test-layout.R` cases "draws top-n items under the deepest level", "draws no item layer (default unchanged)", "works in horizontal direction", "too-high item_cut draws no items", ".ba_item_text uses variable labels" all pass; item IDs match the `top_items(level=k_max)` oracle.
- AC2 ✓ — cases ".resolve_node_size resolves scalars/named/errors", "per-node node_width overrides", "bare-scalar node sizes", "overlap warning vector-safe" pass; named box renders at 1.9, unnamed at 0.8 default.
- AC3 ✓ — six `order=` cases pass: deepest x-rank equals supplied permutation, list form ignores non-deepest (warn), non-permutation/dup/wrong-type error, list omitting deepest errors.
- AC4 ✓ — "autoplot() forwards order= to ba_layout()" passes (built-plot node x reflects manual order).
- AC5 ✓ — `check-vignette-freshness.R` clean; `document()` produces no diff; `NAMESPACE` + `_pkgdown.yml` unchanged in the branch diff (no new exports); NEWS entry present; visualization vignette + roxygen `@examples` cover all three.
- AC6 ✓ — `Rscript tools/dod-gate.R` → GATE PASSED (check 0/0/0, coverage 100%, style/lint clean, pkgdown index complete). Full layout suite 349 pass.

**Consistency gate:** `cairn_validate` exit 0 (all CHECKs pass incl. coverage-complete; 91 pre-existing "dangling id tokens" advisories re archived M10–M53, not introduced here). No principle changed → `cairn_impact` skipped. Toolchain slot: `document()` no-diff, full `check()` 0/0/0, pkgdown, NEWS — all clean.

**Independent three-lens review — zero actioned findings:**
- [O] diff-bug (Opus): no defect; traced rank math (`.resolve_manual_order`), `top_items` sort dependency (`.ba_item_text`), per-node edge faces (`.attach_coords`), `.resolve_node_size` validation, `drop_pruned`/`compress_levels` interaction, defaults-unchanged. One non-finding dropped per taxonomy: long item lists can visually crowd siblings in horizontal layout (accepted aesthetic tradeoff, mirrors vertical downward extension).
- [S] blame-history: no regression; `order=NULL` branch byte-identical to master (M80 keep-best intact), manual bypass faithful to M80's own "deepest level is the only free variable" design; vignette churn scoped correctly (M61/M75).
- [S] prior-review: prior evidence exists (M80 F1 keep-best, M79 F1 secondary-linewidth) — neither regressed; GitHub PR-comment probe empty.
- Scorer: no findings to score (no-op).

Noted (not actioned): horizontal item-list crowding for long lists — a known aesthetic limitation of hand-tuned figures; becomes a candidate only if it bites in practice.
