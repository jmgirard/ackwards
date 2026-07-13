# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-13 (M58 shipped)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M59 | De-duplicate console output & plot builders | planned | M58 | normal | milestones/M59-dedup-console-output.md |
| M58 | Consolidate input-validation helpers & fix two drift bugs | done | — | normal | milestones/archive/M58-consolidate-validation-helpers.md |
| M57 | Ossify oracles — reproducible, catalogued oracle discipline | done | — | normal | milestones/archive/M57-ossify-oracles.md |
| M56 | BRM manuscript: reproducible Quarto scaffold + worked example | done | — | normal | milestones/archive/M56-brm-manuscript.md |
| M55 | Address CRAN 0.1.0 feedback and resubmit as 0.1.1 | done | — | high | milestones/archive/M55-cran-resubmission.md |
| M54 | Export `forbes2023` as a bundled dataset | done | — | normal | milestones/archive/M54-amh-cor-dataset.md |
<!-- M01–M53 are done/dropped (entombed in cairn/legacy/MILESTONES.md); new work from M55. -->

## Candidates

- Remaining dedup pass (readability-only; no correctness/perf gain — DESIGN §3 "measure before optimizing", the setup-path recomputes aren't hotspots): shared per-engine variance / score-var / positive-manifold-flip helpers; drop `ackwards()`'s doubly-built edge tidy tibble + discarded `.align_signs()$edges`; thread `eigen(R)`'s smallest eigenvalue instead of recomputing 2–3×; remove the dead `k_eff>=1L` guard + vestigial all-TRUE convergence meta; `.fit_levels_muffled()` + convergence-shortfall reporter shared by comparability/boot_edges — added 2026-07-12 — bucket 3 of the 2026-07-12 codebase audit, deferred at the M58/M59 plan gate
- Draft the author-owned Intro + Discussion prose for the BRM manuscript (scholarly argument/framing) — M56 scaffold shipped 2026-07-12, now actionable — M56 Out
- Wire `forbes2023` into a vignette (worked Forbes AMH example) once the dataset ships (M54) — added 2026-07-12 — M54 Out
- Owner-only post-M55 release tail: interactive `devtools::submit_cran()` of 0.1.1, tag `v0.1.1`, update README "on CRAN" phrasing when accepted; if CRAN bounces again, plan the next resubmission milestone — added 2026-07-12 — supersedes the 0.1.0 tail row (0.1.0 submitted; 2026-07-12 reviewer feedback became M55, which absorbed the remote-check steps and superseded the patch-branch-from-tag guidance at its plan gate)
- `comparability()` ESEM engine/basis extension: split-half factor comparability per level per factor; feasible but demand-gated (2·n_splits lavaan hierarchies per call; per-half convergence handling) — added 2026-07-11 — DESIGN §14.35 / M46
- `boot_edges()` engine/basis extension: WLSMV/polychoric bootstrap edges; expensive (n_boot × (k_max−1) fits) and statistically wrinkly (resample can drop a response category) — keep off schedule until asked — added 2026-07-11 — DESIGN §14.36 / M47
- Full DESIGN §14 → DECISIONS.md extraction (Compromise B): lift the entire embedded decision log out of DESIGN.md into DECISIONS.md and repoint every inline `§14.x` reference; run as its own focused pass — added 2026-07-11 — M20 migration gate
- Formalize IP/GP principles via /design-interview: ackwards' hard constraints live as CLAUDE.md "Invariants"; a design-interview run would elicit and number them as DESIGN IP<n>/GP<n> — added 2026-07-11 — M20 migration gate — NB: M57 adds Invariant #8 (oracle-backed numerics) as interim home; fold it in at IP/GP time
