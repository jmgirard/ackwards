# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-27 (release tail reconciled with git — `v0.2.0` tagged + GitHub release published + `use_dev_version()` run; the 91 dangling-id-token advisories are all pre-migration M≤53 ids entombed in `cairn/legacy/`)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M80 | Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging | done | — | normal | milestones/archive/M80-deep-hierarchy-layout-quality.md |
| M81 | Publication-figure polish — item lists, per-node box sizes, manual factor ordering | done | M80 | normal | milestones/archive/M81-publication-figure-polish.md |
| M82 | macOS oldrel (arm64) CI visibility + skip-filter repair | done | — | normal | milestones/archive/M82-macos-oldrel-ci-visibility.md |
| M83 | Windows parallel-testthat crash — measure, diagnose, mitigate | done | — | high | milestones/archive/M83-windows-parallel-testthat-crash.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent done rows. -->

## Candidates

- Owner-only release tail: **0.2.0 ACCEPTED on CRAN 2026-07-27** (released content = master `98b7366`) — tagged `v0.2.0` at `98b7366` and pushed, GitHub release published, `usethis::use_dev_version()` run (`9181435`; DESCRIPTION now 0.2.0.9000). 0.1.1 accepted 2026-07-24, tagged `v0.1.1`. The 0.2.0 road: submitted `07b6f41` → auto-rejected the same day (noLD ERROR), fixed by degenerate-level truncation (`6ad4879`) and a bit-identity test tolerance (`da8a1eb`), verified on R-hub `atlas` + `nold` and re-checked on win-builder, resubmitted and accepted. This also cleared the ATLAS/noLD issue-database entry (deadline was 2026-08-21) and the `r-oldrel-macos-arm64` ERROR. **Remaining: push `9181435` to origin (master is ahead 1) and commit the `CRAN-SUBMISSION` deletion sitting unstaged in the worktree; then this row retires.** — updated 2026-07-27
- Further CRAN-flavour coverage beyond M82 (grouped): a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its own push-to-master job) and a macOS **x86_64** row (`macos-13`) — the 0.1.1 failure was arm64-specific and `r-oldrel-macos-x86_64` passed. Promote either on any CRAN macOS failure M82's job could not have caught. **Promoted-by-evidence 2026-07-27: the alternative-numerics gap fired** — 0.2.0 was auto-rejected on a noLD ERROR no green macOS/Windows/Ubuntu matrix could have caught, because the defect needed degenerate rather than merely ill-conditioned arithmetic. Add R-hub `atlas` + `nold` to `PROFILE.md`'s release-walk as a pre-submission step (the `rhub.yaml` workflow already exists; `rhub::rhub_check(platforms = c("atlas", "nold"))` is the invocation). — added 2026-07-27, updated 2026-07-27
- Untested axis from M83: **disabling testthat parallelism outright** (`Config/testthat/parallel: false` / `TESTTHAT_PARALLEL=false`), as distinct from the worker counts M83 measured. Evidence it matters: at `TESTTHAT_CPUS=1` the crashes still read `testthat subprocess exited`, so the parallel machinery still spawned a subprocess — the crashing component was never removed in any sweep. Plausibly the actual fix, and the only candidate that eliminates rather than reduces. Costs M48's speedup (27s parallel vs 81s serial locally) and, unlike M83's mitigation, changes **tarball content** (DESCRIPTION), so it needs release re-verification. Promote if a CRAN Windows flavour fails with the -1073741819 signature. — added 2026-07-27
- Upstream report for the M83 Windows crash, if it isolates to a dependency's compiled code (`EFAtools` or `mnormt`; `psych`/`lavaan`/`GPArotation` are pure R): file with a minimal reproducer. Held out of M83 at the 2026-07-27 plan gate because a maintainer's timeline must not gate our CI health. Promote when M83's bisection names a component. — added 2026-07-27
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23
