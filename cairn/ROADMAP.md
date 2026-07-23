# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-23 (M72 merged + archived; M66 row pruned under terminal-row retention)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M70 | Author + verify the 5 default-rationale backer notes and wire them into DESIGN §9 + roxygen | done | — | normal | milestones/archive/M70-reference-notes-backers.md |
| M71 | Author + verify the 5 application source notes as citation precedents + Forbes drift-watch | in-progress | M70 | normal | milestones/M71-reference-notes-applications.md |
| M72 | Correct the stale Forbes contract line + author the Goldberg/Forbes departures ledger | done | — | normal | milestones/archive/M72-source-departures-ledger.md |
| M69 | Author + verify the 8 secondary-methods source notes against their shelf PDFs | done | — | normal | milestones/archive/M69-reference-verification-secondary-sources.md |
| M67 | Re-verify the 9 single-source reference pages against their shelf PDFs | done | — | normal | milestones/archive/M67-reference-verification-method-pages.md |
| M68 | Re-verify the 3 collapsed synthesis pages against their member sources | done | M67 | normal | milestones/archive/M68-reference-verification-collapsed-pages.md |
<!-- M01–M59 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent. -->

## Candidates

- Draft the author-owned Intro + Discussion prose for the BRM manuscript (scholarly argument/framing) — M56 scaffold shipped 2026-07-12, now actionable — M56 Out
- Owner-only post-M55 release tail: 0.1.1 **resubmitted 2026-07-17** (tarball from master `0a5da58`, after the 2026-07-16 withdrawal; CRAN-SUBMISSION committed) — **owner: confirm the CRAN submission email if pending; on acceptance tag `v0.1.1` + update README "on CRAN" phrasing; if CRAN bounces again, plan the next resubmission milestone** — added 2026-07-12, updated 2026-07-17
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16
