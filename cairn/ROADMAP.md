# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-11 (migrated from Lineage B via /cairn-init)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M54 | Export `amh_cor` as a bundled dataset | planned | — | normal | milestones/M54-amh-cor-dataset.md |
<!-- M01–M53 are done/dropped (entombed in cairn/legacy/MILESTONES.md); new work from M54. -->

## Candidates

- Wire `amh_cor` into a vignette (worked Forbes AMH example) once the dataset ships (M54) — added 2026-07-12 — M54 Out
- Owner-only 0.1.0 CRAN release tail: win-builder / R-hub remote checks, then interactive `devtools::submit_cran()`; update README "on CRAN" phrasing when it lands; if CRAN bounces, patch branches from the `v0.1.0` tag — added 2026-07-11 — legacy/ROADMAP.md
- `comparability()` ESEM engine/basis extension: split-half factor comparability per level per factor; feasible but demand-gated (2·n_splits lavaan hierarchies per call; per-half convergence handling) — added 2026-07-11 — DESIGN §14.35 / M46
- `boot_edges()` engine/basis extension: WLSMV/polychoric bootstrap edges; expensive (n_boot × (k_max−1) fits) and statistically wrinkly (resample can drop a response category) — keep off schedule until asked — added 2026-07-11 — DESIGN §14.36 / M47
- Full DESIGN §14 → DECISIONS.md extraction (Compromise B): lift the entire embedded decision log out of DESIGN.md into DECISIONS.md and repoint every inline `§14.x` reference; run as its own focused pass — added 2026-07-11 — M20 migration gate
- Formalize IP/GP principles via /design-interview: ackwards' hard constraints live as CLAUDE.md "Invariants"; a design-interview run would elicit and number them as DESIGN IP<n>/GP<n> — added 2026-07-11 — M20 migration gate
