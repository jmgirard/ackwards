# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-24 (CRAN accepted 0.1.1 — release-tail candidate row updated; validate green, all inboxes empty)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M77 | Near-redundant band flag in artifact mode + fix artifact vignette example | done | M76 | normal | milestones/archive/M77-near-redundant-band.md |
| M78 | Gap-tolerant (skip-level) direct redundancy chains — gated on ChaseCorrPaths semantics | done | — | normal | milestones/archive/M78-skip-level-redundancy-chains.md |
| M79 | Secondary-correlation edges in the pruned/publication view | done | — | normal | milestones/archive/M79-secondary-correlation-edges.md |
| M80 | Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging | done | — | normal | milestones/archive/M80-deep-hierarchy-layout-quality.md |
| M81 | Publication-figure polish — item lists, per-node box sizes, manual factor ordering | done | M80 | normal | milestones/archive/M81-publication-figure-polish.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent done rows. -->

## Candidates

- Owner-only post-M55 release tail: 0.1.1 **ACCEPTED on CRAN 2026-07-24** (submitted 2026-07-17, tarball from master `0a5da58`). Remaining acceptance follow-up: **tag `v0.1.1`** at `0a5da58` and push it, and refresh the README to confirm the on-CRAN state (line 39 already reads "released version from CRAN"; no CRAN status badge yet — add via `README.Rmd` + re-render if wanted). — added 2026-07-12, updated 2026-07-24
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23
