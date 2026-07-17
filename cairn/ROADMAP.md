# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-17 (M66 merged + archived; terminal-row retention pruned M61)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M66 | Upstream re-oracle watch (weekly CRAN-current CI) | done | — | normal | milestones/archive/M66-upstream-reoracle-watch.md |
| M65 | Mechanize the precomputed-vignette staleness guard | done | — | normal | milestones/archive/M65-vignette-staleness-guard.md |
| M64 | DESIGN §14 → DECISIONS.md extraction (hybrid + entomb) | done | — | normal | milestones/archive/M64-design-s14-extraction.md |
| M63 | User-facing text accuracy pass (Goldberg 2006 citations + PCA n_obs message) | done | — | normal | milestones/archive/M63-text-accuracy-pass.md |
| M62 | Worked Forbes (2023) AMH example vignette | done | — | normal | milestones/archive/M62-forbes2023-vignette.md |
<!-- M01–M59 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent. -->

## Candidates

- Draft the author-owned Intro + Discussion prose for the BRM manuscript (scholarly argument/framing) — M56 scaffold shipped 2026-07-12, now actionable — M56 Out
- Owner-only post-M55 release tail: 0.1.1 **resubmitted 2026-07-17** (tarball from master `0a5da58`, after the 2026-07-16 withdrawal; CRAN-SUBMISSION committed) — **owner: confirm the CRAN submission email if pending; on acceptance tag `v0.1.1` + update README "on CRAN" phrasing; if CRAN bounces again, plan the next resubmission milestone** — added 2026-07-12, updated 2026-07-17
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16
