# M79: Secondary-correlation edges in the pruned/publication view

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, GP2
- **Branch/PR:** m79-secondary-correlation-edges

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

- [ ] AC1: With `show_secondary = TRUE` on a pruned object, `autoplot` draws a
  secondary edge for every kept cross-level pair `|r| >= cut_show` that is not the
  primary edge — including at least one cross-branch pair and one same-lineage
  skip pair (verified on the returned secondary set + `layer_data`).
- [ ] AC2: Secondary edges render dimmed (reduced alpha) + thinner and inherit
  the sign encoding; the sign channel (colour/dash) is byte-identical to the
  `show_secondary = FALSE` render (verified via `layer_data` colour/linetype).
- [ ] AC3: The default (`show_secondary = FALSE`) reproduces the current pruned
  view exactly — no extra layer, existing edge `layer_data` unchanged.
- [ ] AC4: `devtools::test()` clean; `devtools::document()` no diff; the
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
- [ ] T2: add the `show_secondary` arg + the dimmed/thinner secondary
  `geom_segment` layer in the `drop_pruned` branch (`autoplot.R:398-429`);
  `layer_data` tests for the new render (secondary layer alpha/linewidth present;
  sign colour/linetype intact) and for the `FALSE` default (no extra layer,
  existing edge layer unchanged).
- [ ] T3: roxygen `@param` + `@example`; showcase in
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

## Decisions

## Review
