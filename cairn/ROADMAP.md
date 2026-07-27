# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-27 (release-tail row corrected — `v0.1.1` tag was already pushed; row repointed at the 0.2.0 submission and the NEWS misfiling recorded)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M79 | Secondary-correlation edges in the pruned/publication view | done | — | normal | milestones/archive/M79-secondary-correlation-edges.md |
| M80 | Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging | done | — | normal | milestones/archive/M80-deep-hierarchy-layout-quality.md |
| M81 | Publication-figure polish — item lists, per-node box sizes, manual factor ordering | done | M80 | normal | milestones/archive/M81-publication-figure-polish.md |
| M82 | macOS oldrel (arm64) CI visibility + skip-filter repair | done | — | normal | milestones/archive/M82-macos-oldrel-ci-visibility.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent done rows. -->

## Candidates

- Owner-only post-M55 release tail: 0.1.1 **ACCEPTED on CRAN 2026-07-24** (submitted 2026-07-17, tarball from master `0a5da58`). **Acceptance follow-up is complete** — annotated tag `v0.1.1` exists at `0a5da58` and is on origin, README refreshed and CRAN status badge added at `cc0bdbc` (the earlier "tag it and push it" item was already done when written; corrected 2026-07-27). **0.2.0 release prep done 2026-07-27** at `f03ea2f` (version bump, NEWS 0.2.0 section, cran-comments leading with the `r-oldrel-macos-arm64` ERROR it clears; the misfiled M75/M77/M79/M80/M81 bullets moved out of the 0.1.1 section, which is now byte-identical to the accepted tarball's). Local gate 0/0/0 + 100% coverage; win-builder R-devel submitted. **Submission itself is the owner's action** and had not been taken as of that date. — added 2026-07-12, updated 2026-07-27
- Further CRAN-flavour coverage beyond M82 (grouped): a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its own push-to-master job) and a macOS **x86_64** row (`macos-13`) — the 0.1.1 failure was arm64-specific and `r-oldrel-macos-x86_64` passed. Promote either on any CRAN macOS failure M82's job could not have caught. — added 2026-07-27
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23
