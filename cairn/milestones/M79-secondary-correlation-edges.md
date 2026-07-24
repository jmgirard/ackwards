# M79: Secondary-correlation edges in the pruned/publication view

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, GP2
- **Branch/PR:** m79-secondary-correlation-edges / #84

## Goal

In the pruned (`drop_pruned`) view, optionally draw the between-level
correlation edges the single-strongest-ancestor primary view hides: every kept
cross-level factor pair with `|r| >= cut_show` that is not already a primary
edge. Render them in a channel distinct from the sign encoding (dimmed +
thinner) so the sign dash/colour stays intact (Forbes website-review item E).

## Scope

**In:**
- Extend `.drop_pruned_nodes` (`layout.R:173`) to also return a **secondary**
  edge set: from the all-pairs table it already builds
  (`compute_edges(pairs = "all")` — IP1, one edge path; weights are sign-aligned),
  take every kept cross-level pair (`from`/`to` both kept, `level_from` shallower
  than the child) and subtract the per-node primary edges it already selects. The
  set therefore includes *both* cross-branch second parents (e.g. `f1 -> d2`) and
  same-lineage skip arcs (e.g. `f1 -> b1` when the primary path is
  `f1 -> d1 -> b1`) — consistent with D-032/D-017 that a direct skip-level `r` is
  a distinct, non-transitive fact, not implied by the path.
- In `autoplot`'s `drop_pruned` branch (`autoplot.R:398-429`), gate the secondary
  layer behind a new `show_secondary` argument (default `FALSE`). Filter the
  secondary set by the resolved `cut_show`, attach coordinates + sign annotation,
  and draw it as a **separate** `geom_segment` layer at reduced alpha (~0.4) and a
  constant thin linewidth *below* the primary magnitude range; the layer inherits
  the existing sign colour/linetype scales (via `sign_group`) so sign still reads,
  and never maps the magnitude channel. Reuse `cut_show` (no separate threshold);
  by construction each node's secondaries are weaker than its primary (= strongest
  kept ancestor), so no secondary ever shows while that node's primary is hidden.
- Tests via `ggplot_build()`/`layer_data()` introspection (the repo has **no
  vdiffr**; plots are asserted this way, e.g. `test-layout.R`). roxygen `@param` +
  runnable `@example`; showcase in the precomputed `ackwards-visualization`
  vignette; NEWS; `_pkgdown.yml` only if the export surface changes.

**Out:** layout crossing reduction / edge-label dodging -> M80; manual factor
ordering / box sizing / item lists -> F (candidate, depends on M80); a separate
secondary-only threshold -> candidate (reuse `cut_show` for now). Shares the
`autoplot`/`layout` surface with M80 but is logically independent — whichever
lands first, the second rebases.

## Acceptance criteria

- [x] AC1: With `show_secondary = TRUE` on a pruned object, `autoplot` draws a
  secondary edge for every kept cross-level pair `|r| >= cut_show` that is not the
  primary edge — including at least one cross-branch pair and one same-lineage
  skip pair (verified on the returned secondary set + `layer_data`).
- [x] AC2: Secondary edges render dimmed (reduced alpha) + thinner and inherit
  the sign encoding; the sign channel (colour/dash) is byte-identical to the
  `show_secondary = FALSE` render (verified via `layer_data` colour/linetype).
- [x] AC3: The default (`show_secondary = FALSE`) reproduces the current pruned
  view exactly — no extra layer, existing edge `layer_data` unchanged.
- [x] AC4: `devtools::test()` clean; `devtools::document()` no diff; the
  `ackwards-visualization` vignette re-precomputed with its freshness stamp
  current; `_pkgdown.yml` updated iff the surface changed; `devtools::check()`
  clean.

## Coverage

- AC1 -> T1, T2
- AC2 -> T2
- AC3 -> T2
- AC4 -> T3

## Tasks

- [x] T1: (test-first) extend `.drop_pruned_nodes` (`layout.R:173`) to return the
  secondary-edge set (all non-primary kept cross-level pairs, unfiltered by cut);
  unit-test the selection at the helper level — a planted cross-branch pair and a
  planted same-lineage skip both present, every primary edge absent.
- [x] T2: add the `show_secondary` arg + the dimmed/thinner secondary
  `geom_segment` layer in the `drop_pruned` branch (`autoplot.R:398-429`);
  `layer_data` tests for the new render (secondary layer alpha/linewidth present;
  sign colour/linetype intact) and for the `FALSE` default (no extra layer,
  existing edge layer unchanged).
- [x] T3: roxygen `@param` + `@example`; showcase in
  `ackwards-visualization.Rmd.orig` then re-run `Rscript vignettes/precompute.R`
  (revert timing-only churn line-by-line per M75/M61, confirm freshness stamp);
  NEWS.md entry; `_pkgdown.yml` iff exports change; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-24: re-planned (/milestone-plan M79) — resolved the secondary-edge
  definition (all non-primary above-cut, incl. same-lineage skips), the visual
  channel (dimmed + thinner, sign preserved), and the vignette showcase; corrected
  the test method from the non-existent "vdiffr snapshot" to the repo's
  `ggplot_build()`/`layer_data()` introspection pattern.
- 2026-07-24: T1 — `.drop_pruned_nodes` now returns `secondary` (all non-primary
  kept cross-level pairs); helper-level selection test green (bfi25 k_max=4 →
  9 secondary edges, both adjacent + skip kinds present).
- 2026-07-24: T2 — `show_secondary` arg (default FALSE) + a distinct secondary
  `geom_segment` layer (alpha 0.4, linewidth 0.3, no arrowheads) drawn beneath
  the primary edges, inheriting the sign colour/linetype. Render + default +
  linetype-sign tests green. Minor discovered fix (M77 lesson): reworded
  prune.R roxygen `` `r = 0.89` `` → `` `|r| = 0.89` `` to clear a non-fatal
  document() markdown warning.
- 2026-07-24: T3 — roxygen `@param`/`@example`; `ackwards-visualization` vignette
  subsection (re-precomputed, timing churn on the 5 untouched vignettes reverted
  per M61/M75, freshness stamp current); NEWS entry; no `_pkgdown.yml` change
  (new arg, not a new export). Gate: check 0/0/0, coverage 100% (layout +
  autoplot 100%), lint 0, style applied, pkgdown OK. (dod-gate.R exit 144 after
  a clean check was resource pressure from covr running back-to-back with check;
  each piece passes standalone.)

## Decisions

## Review

**AC evidence (fresh, 2026-07-24):**
- AC1 — `test-layout.R` M79 tests green: helper test asserts primary/secondary
  disjoint, union = all kept cross-level pairs, both an adjacent second-parent
  (gap 1) and a same-lineage skip (gap ≥ 2) present; render test asserts the
  drawn secondary count == `sum(|r| ≥ cut_show)` on `$secondary`. FAIL 0.
- AC2 — render test: secondary `layer_data` colours ⊆ sign palette
  (`{#2166AC,#D6604D}`), primary layer colours byte-identical with/without
  secondary; `sign_by="linetype"` test: secondary linetype set == primary's.
- AC3 — default test: `drop_pruned` alone has exactly one `GeomSegment` layer,
  no dimmed (alpha 0.4) layer; matches the prior pruned view.
- AC4 — layout suite FAIL 0 (249 assertions); `document()` no diff; vignette
  freshness OK; NEWS entry present, no milestone numbers; no `_pkgdown.yml`
  change (new arg, not a new export); `check()` 0 err/0 warn/0 note.

**Consistency gate:** `cairn_validate` all checks passed (108 advisory warnings
= pre-existing archive cross-refs, not gate failures); coverage-completeness map
clean. `document()` no diff (generated files not hand-edited). `pkgdown::check_pkgdown()`
"No problems found." coverage 100% (layout + autoplot 100%; diagnostic only).
lint 0, styled. No `IPn`/`GPn` added or changed (IP1/GP2 worked-under) → `cairn_impact` skipped.

**Independent three-lens review** (fresh-context, parallel, distinct evidence bases):
- blame-history (Sonnet): **0 findings** — secondary set reuses M42's single
  `compute_edges(pairs="all")` table (IP1/D-004 intact), doesn't disturb
  primary-edge selection, preserves "no unmapped aesthetic" + default-unchanged.
- prior-review (Sonnet): **0 findings** — vignette churn correctly scoped to the
  one touched vignette + its 2 assets (M61/M75), no M59 one-line branch, the
  prune.R edit is the M77 fix; PR-comment probe returned `[]` (no thread walk).
- diff-bug (Opus): 2 raised. Scorer (Sonnet, independent):
  - **F1 (85, actioned → fixed):** secondary `linewidth` was hardcoded `0.3`,
    so a user `edge_linewidth < 0.3` (e.g. `0.2`) rendered the *primary* thinner
    than the secondary — inverting the documented "thinner" contract. Fixed:
    `sec_linewidth <- min(0.3, 0.6 * width_val)`, strictly below the primary in
    every mode (mapped range floor 0.4, or any constant primary width). Regression
    test added (`edge_linewidth = 0.2` → secondary 0.12 < primary 0.2).
  - **F2 (40, logged, not actioned):** `show_r` labels only primary edges, not
    secondary. Reviewer + scorer agree this is an intentional de-emphasis; no
    roxygen promises secondary labels. Below the 80 action bar.

**Post-fix gate:** layout suite 252 assertions FAIL 0; `check()` 0/0/0;
coverage 100% (layout + autoplot); lint 0; styled.
