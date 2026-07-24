# M77: Near-redundant band flag in artifact mode + fix the artifact vignette example

- **Status:** planned
- **Priority:** normal
- **Depends on:** M76
- **Driving RR:** —
- **Principles touched:** GP2, IP6
- **Branch/PR:** —

## Goal

Add a report-only "near-redundant" band signal to `prune("artifact")` flagging factor pairs just below the redundancy thresholds, and fix the artifact-mode vignette so its illustrated pair is genuinely near-redundant (candidate D).

## Scope

**In:**
- In `prune(x, "artifact")`, compute a **near-redundant band** signal: cross-level factor pairs whose direct `|r|` and/or Tucker `phi` fall within a margin *below* `redundancy_r` / `redundancy_phi` but are **not** already flagged by `prune("redundant")` (Forbes uses the flags mainly for this just-below band — e.g. `r=.89`/`phi=.93`, a messy re-rotation — since full redundancy is already dropped in the redundancy stage).
- Store it report-only in `x$prune` (GP2 — never drops or mutates the kept set), alongside the existing `x$prune$phi` evidence.
- Expose the band width as an argument (e.g. `near_margin`) with a documented default; if the default is engine-dependent, announce the auto-resolved value via cli (IP6).
- Fix the artifact-mode worked example in `ackwards-forbes.Rmd.orig` — the current `phi ~= .99` pair would already be flagged by `prune("redundant")`, so it is a poor illustration; use a genuinely near-redundant pair.
- roxygen + `_pkgdown.yml` reference if the surface changes; NEWS; re-precompute the vignette.

**Out:** the prose rewrite of `redundancy_criterion`/oblique → M76; the skip-level chain change → M78.

## Acceptance criteria

- [ ] AC1: On a fixture with a planted near-redundant pair, `prune(x, "artifact")` returns a band signal identifying pairs with direct `|r|` or `phi` within `near_margin` below `redundancy_r`/`redundancy_phi` and not already flagged redundant.
- [ ] AC2: The band is report-only — `prune("artifact")` never drops or mutates the kept node set because of it (GP2).
- [ ] AC3: `near_margin` has a documented default; an engine-dependent auto-resolution is announced via cli (IP6).
- [ ] AC4: The artifact vignette example illustrates a genuinely near-redundant pair (one *not* already dropped by `prune("redundant")`).
- [ ] AC5: `devtools::test()` clean; `devtools::document()` no diff; `_pkgdown.yml` updated if any export changed (`pkgdown::check_pkgdown()`); vignette re-precomputed + freshness passes; `devtools::check()` clean.

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T2, T3
- AC4 → T4
- AC5 → T3, T4, T5

## Tasks

- [ ] T1: (test-first) tests: fixture (planted near-redundant pair, e.g. a `data-raw/` matrix or `sim16` at a chosen `k`) → `prune("artifact")` flags it; report-only invariant (kept set unchanged); boundary edges (`NA` phi, exactly-at-threshold, disabled phi).
- [ ] T2: Implement the near-redundant band in artifact mode (`prune.R`); add `near_margin` arg + default resolution (mirror `redundancy_phi`'s auto-resolve/announce pattern where engine-dependent).
- [ ] T3: roxygen `@param`/`@details`; `devtools::document()`; update `_pkgdown.yml` if the surface changed.
- [ ] T4: Fix the artifact vignette example + prose; re-run `Rscript vignettes/precompute.R`; `git checkout --` untouched assets.
- [ ] T5: NEWS.md entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.

## Decisions

## Review
