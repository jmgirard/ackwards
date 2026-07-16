# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-16 (M60 done)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M60 | De-duplicate the setup path (audit bucket 3) | done | — | normal | milestones/archive/M60-bucket3-dedup.md |
| M59 | De-duplicate console output & plot builders | done | M58 | normal | milestones/archive/M59-dedup-console-output.md |
| M58 | Consolidate input-validation helpers & fix two drift bugs | done | — | normal | milestones/archive/M58-consolidate-validation-helpers.md |
| M57 | Ossify oracles — reproducible, catalogued oracle discipline | done | — | normal | milestones/archive/M57-ossify-oracles.md |
| M56 | BRM manuscript: reproducible Quarto scaffold + worked example | done | — | normal | milestones/archive/M56-brm-manuscript.md |
<!-- M01–M55 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent. -->

## Candidates

- Positive-manifold flip helper across engines: the k=1 sign-flip test is identical but each engine's follow-up diverges (PCA negates `fit$weights`, EFA reuses a `flip` flag in its fallback, ESEM leaves `L_se` alone) and it overlaps `.align_signs()` — a shared helper saves ~2 lines/engine while hiding the coupling — dropped from M60 at the 2026-07-13 plan gate as low-value/awkward
- Shortfall-reporter helper shared by comparability/boot_edges: same two-bullet cli template but divergent wording, nouns, and computed stats — would need every string parameterized for little gain — dropped from M60 at the 2026-07-13 plan gate
- Draft the author-owned Intro + Discussion prose for the BRM manuscript (scholarly argument/framing) — M56 scaffold shipped 2026-07-12, now actionable — M56 Out
- Wire `forbes2023` into a vignette (worked Forbes AMH example) once the dataset ships (M54) — added 2026-07-12 — M54 Out
- Owner-only post-M55 release tail: interactive `devtools::submit_cran()` of 0.1.1, tag `v0.1.1`, update README "on CRAN" phrasing when accepted; if CRAN bounces again, plan the next resubmission milestone — added 2026-07-12 — supersedes the 0.1.0 tail row (0.1.0 submitted; 2026-07-12 reviewer feedback became M55, which absorbed the remote-check steps and superseded the patch-branch-from-tag guidance at its plan gate)
- `comparability()` ESEM engine/basis extension: split-half factor comparability per level per factor; feasible but demand-gated (2·n_splits lavaan hierarchies per call; per-half convergence handling) — added 2026-07-11 — DESIGN §14.35 / M46
- `boot_edges()` engine/basis extension: WLSMV/polychoric bootstrap edges; expensive (n_boot × (k_max−1) fits) and statistically wrinkly (resample can drop a response category) — keep off schedule until asked — added 2026-07-11 — DESIGN §14.36 / M47
- Full DESIGN §14 → DECISIONS.md extraction (Compromise B): lift the entire embedded decision log out of DESIGN.md into DECISIONS.md and repoint every inline `§14.x` reference; run as its own focused pass — added 2026-07-11 — M20 migration gate
- Formalize IP/GP principles via /design-interview: ackwards' hard constraints live as CLAUDE.md "Invariants"; a design-interview run would elicit and number them as DESIGN IP<n>/GP<n> — added 2026-07-11 — M20 migration gate — NB: M57 adds Invariant #8 (oracle-backed numerics) as interim home; fold it in at IP/GP time
- Enrich `suggest_k()` roxygen with two verified-at-ingest citations: Lim & Jahng (2019) "PA estimate as a ±1 range" (matches the advisory-range stance; pair with Achim 2021 counterpoint) and Saucier (1997) fn 14 (PA suggested ~30 factors in wide lexical sets — in-the-wild "PA-PC overextracts" case); reference notes in cairn/references/ — added 2026-07-16 — comparability-citations hotfix (#61) follow-up
