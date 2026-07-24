# M77: Near-redundant band flag in artifact mode + fix the artifact vignette example

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M76
- **Driving RR:** —
- **Principles touched:** GP2, IP6
- **Branch/PR:** m77-near-redundant-band

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

- [x] T1: (test-first) tests: fixture (planted near-redundant pair, e.g. a `data-raw/` matrix or `sim16` at a chosen `k`) → `prune("artifact")` flags it; report-only invariant (kept set unchanged); boundary edges (`NA` phi, exactly-at-threshold, disabled phi).
- [x] T2: Implement the near-redundant band in artifact mode (`prune.R`); add `near_margin` arg + default resolution (mirror `redundancy_phi`'s auto-resolve/announce pattern where engine-dependent).
- [x] T3: roxygen `@param`/`@details`; `devtools::document()`; update `_pkgdown.yml` if the surface changed.
- [ ] T4: Fix the artifact vignette example + prose; re-run `Rscript vignettes/precompute.R`; `git checkout --` untouched assets.
- [ ] T5: NEWS.md entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: question gate — confirmed OR band semantics + `near_margin = 0.1` default (both recommended).
- 2026-07-23: T1+T2 — added `.near_redundant_pairs()` + `near_margin` arg to `prune()`; artifact mode now stores `x$prune$near_redundant`; `redundancy_phi` auto-resolve extended to artifact mode; cli announces the band. Unit + integration tests pass (bfi25 k=5 polychoric band includes the planted m1f1↔m2f1 r=.89/φ=.94 pair; sim16 empty-band case).

## Decisions

- **M77-D1: Near-redundant band = OR of just-below r/φ windows, pair-level, report-only.** `prune(x, "artifact")` computes `x$prune$near_redundant`: cross-level pairs that are **not** fully redundant (`|r| ≥ redundancy_r` and, where a φ threshold is active, `φ > redundancy_phi`) yet have **at least one** signal in its just-below window — `|r| ∈ [redundancy_r − near_margin, redundancy_r)` **or** `φ ∈ [redundancy_phi − near_margin, redundancy_phi)`. `|r|` is the direct/skip-level score correlation (same quantity `redundancy_criterion="direct"` chases); φ is signed (mirrors the redundant rule's `φ > threshold` test). `near_margin` default `0.1` (absolute, fixed — not engine-dependent). "Not already flagged" is operationalized at the **pair** level (does the pair itself meet full redundancy) rather than against the chain/retention node flags, since redundancy flags nodes and correlation is non-transitive. PCA (`redundancy_phi` NULL) → `|r|` band only. `redundancy_phi` auto-resolution (0.95 for EFA/ESEM) now also fires in artifact mode, since the band consumes it (IP6-announced). Report-only (GP2): never drops/mutates the kept set.

## Review
