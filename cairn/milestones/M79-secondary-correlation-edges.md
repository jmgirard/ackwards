# M79: Secondary-correlation edges in the pruned/publication view

- **Status:** planned
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, GP2
- **Branch/PR:** —

## Goal

In the pruned (`drop_pruned`) view, draw secondary between-level edges (`|r| >= cut_show`) that the primary hierarchy path does not already imply, using a visual channel distinct from the sign-encoding dash (candidate E).

## Scope

**In:**
- In the `drop_pruned` render path (`.drop_pruned_nodes`, `layout.R:157-173`; `autoplot.R:398-400`), after the primary-parent edges are drawn, add **secondary edges**: cross-level factor pairs with `|r| >= cut_show` (default `0.3`) that are **not** already the primary path between those two nodes (e.g. `f1 -> d2` when the primary path is `f1 -> d1 -> b1`).
- Compute the secondary set via `compute_edges()` (IP1 — one edge path; the weights are already sign-aligned), not a second correlation route.
- Render them in a channel **distinct from the dash** — dashed already encodes *sign* (`autoplot`), so secondary edges need a different visual channel (lighter/thinner stroke, or a separate colour/alpha), leaving the sign dash intact.
- Gate behind an argument (e.g. `show_secondary`, default `FALSE`) with roxygen; vdiffr snapshot; NEWS; `_pkgdown.yml` if the surface changes.

**Out:** layout crossing reduction → M80; manual factor ordering / box sizing / item lists → F (candidate). This milestone shares the `autoplot`/`layout` surface with M80 but is logically independent — whichever lands first, the second rebases.

## Acceptance criteria

- [ ] AC1: With `show_secondary = TRUE` on a pruned object, `autoplot` draws secondary edges for cross-level pairs `|r| >= cut_show` not implied by the primary path, in a channel distinct from the sign dash (verified on the resolved edge set + vdiffr).
- [ ] AC2: Sign encoding (dashed = negative) is preserved and not conflated with the secondary channel.
- [ ] AC3: The default (`show_secondary = FALSE`) reproduces the current pruned view unchanged (existing vdiffr snapshot unchanged).
- [ ] AC4: `devtools::test()` clean; `devtools::document()` no diff; `_pkgdown.yml` updated if surface changed; `devtools::check()` clean.

## Coverage

- AC1 → T1, T2
- AC2 → T2
- AC3 → T2
- AC4 → T3

## Tasks

- [ ] T1: (test-first) compute the secondary-edge set — cross-level `|r| >= cut_show` minus the primary path — in `.drop_pruned_nodes`/`autoplot` (`layout.R:157-173`, `autoplot.R:398-400`); unit-test the selection logic (planted secondary edge present, primary-path edge excluded).
- [ ] T2: Add the `show_secondary` arg + distinct aesthetic (not dashed); vdiffr snapshot for the new render; assert sign dash intact.
- [ ] T3: roxygen `@param`; `_pkgdown.yml`; NEWS.md entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.

## Decisions

## Review
