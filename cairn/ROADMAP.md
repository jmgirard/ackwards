# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-16 (M64 done)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M64 | DESIGN §14 → DECISIONS.md extraction (hybrid + entomb) | done | — | normal | milestones/archive/M64-design-s14-extraction.md |
| M63 | User-facing text accuracy pass (Goldberg 2006 citations + PCA n_obs message) | done | — | normal | milestones/archive/M63-text-accuracy-pass.md |
| M62 | Worked Forbes (2023) AMH example vignette | done | — | normal | milestones/archive/M62-forbes2023-vignette.md |
| M61 | Enrich suggest_k() docs with verified k-selection citations | done | — | normal | milestones/archive/M61-suggest-k-citations.md |
| M60 | De-duplicate the setup path (audit bucket 3) | done | — | normal | milestones/archive/M60-bucket3-dedup.md |
<!-- M01–M59 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent. -->

## Candidates

- Draft the author-owned Intro + Discussion prose for the BRM manuscript (scholarly argument/framing) — M56 scaffold shipped 2026-07-12, now actionable — M56 Out
- Owner-only post-M55 release tail: 0.1.1 **submitted to CRAN, awaiting response (as of 2026-07-16)**; NB the submitted tarball predates the lavaan 0.7 compatibility fix (PR #66, 2026-07-16) — if CRAN checks against lavaan 0.7-2 the ESEM tests fail and it will bounce; a resubmission includes that fix; on acceptance tag `v0.1.1` + update README "on CRAN" phrasing; if CRAN bounces again, plan the next resubmission milestone — added 2026-07-12 — supersedes the 0.1.0 tail row (0.1.0 submitted; 2026-07-12 reviewer feedback became M55, which absorbed the remote-check steps and superseded the patch-branch-from-tag guidance at its plan gate)
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16
- Formalize IP/GP principles via /design-interview: ackwards' hard constraints live as CLAUDE.md "Invariants"; a design-interview run would elicit and number them as DESIGN IP<n>/GP<n> — added 2026-07-11 — M20 migration gate — NB: M57 adds Invariant #8 (oracle-backed numerics) as interim home; fold it in at IP/GP time
