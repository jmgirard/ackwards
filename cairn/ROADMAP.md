# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-24 (M79 merged + archived, PR #84; terminal-row retention keeps the 5 most recent done rows — M74 row pruned)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M75 | Beautify suggest_k() print output as an aligned criteria table | done | — | normal | milestones/archive/M75-suggest-k-print-table.md |
| M76 | Clarify redundancy-criterion + oblique-rotation prose and reframe artifact-mode DoF wording | done | — | normal | milestones/archive/M76-redundancy-oblique-dof-prose.md |
| M77 | Near-redundant band flag in artifact mode + fix artifact vignette example | done | M76 | normal | milestones/archive/M77-near-redundant-band.md |
| M78 | Gap-tolerant (skip-level) direct redundancy chains — gated on ChaseCorrPaths semantics | done | — | normal | milestones/archive/M78-skip-level-redundancy-chains.md |
| M79 | Secondary-correlation edges in the pruned/publication view | done | — | normal | milestones/archive/M79-secondary-correlation-edges.md |
| M80 | Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging | in-progress | — | normal | milestones/M80-deep-hierarchy-layout-quality.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent done rows. -->

## Candidates

- Owner-only post-M55 release tail: 0.1.1 **resubmitted 2026-07-17** (tarball from master `0a5da58`, after the 2026-07-16 withdrawal; CRAN-SUBMISSION committed) — **owner: confirm the CRAN submission email if pending; on acceptance tag `v0.1.1` + update README "on CRAN" phrasing; if CRAN bounces again, plan the next resubmission milestone** — added 2026-07-12, updated 2026-07-17
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80 (all planned 2026-07-23).** F and H remain below.

- **[F] feature — publication-figure polish (plan as ONE milestone, 3 tasks — owner-confirmed 2026-07-23).** Three asks: (a) list the top items under the lowest-level factors; (b) adjustable box sizes to fit manual node labels; (c) manual factor ordering per level to untangle/arrange the plot (override `ba_layout`'s ordinal pass). **Depends on M80** — (c) overrides the same ordinal pass M80 rewrites, so plan F after M80 lands. Ready to plan. — added 2026-07-23
- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23
