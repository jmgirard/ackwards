# M91: Near-singular guard for the within-level score correlation

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP6, GP1, GP2
- **Resolves:** —
- **Surface tier:** user-facing, because it adds a warning to `tidy()` output and changes the user docs
- **Branch/PR:** m091-near-singular-score-cor-guard

## Goal

If a level's within-level score correlation is nearly singular, the partialled `beta` values that level feeds always come with a warning.

## Scope

**In:** a smallest-eigenvalue check on Φ_s inside `.partialled_edges()` (R/utils.R:216), at 1e-2. The check warns and keeps the values. Only `tidy(what = "edges")` raises it, because `beta` is the value that collinearity makes unstable. `.tidy_edges()` dedups both the new warning and the existing singular warning per level. Roxygen, NEWS, and a correction to the `beta` entry's scores-path sentence (R/tidy.R:35).

GP1 trade: 1e-2 is a numerical guard that the package chose, not a published cutoff. At 1e-2, two scores at one level correlate about .99. Every real oblique fit probed at plan time stayed above 0.14.

**Out:** a stored `meta` flag for Φ_s is not queued, because `tidy()` recomputes Φ_s on every call and the warning recurs. The manuscript's oblique wording is M92. Intervals on `beta` and `r2` from `boot_edges()` stay in the existing `[low]` candidate row. The R-level check `.near_singular_check()` and its 1e-4 threshold keep their behavior.

## Acceptance criteria

- [x] AC1: Take a level above the deepest level, so that it starts at least one stored pair. Its within-level score correlation Φ_s has a smallest eigenvalue below 1e-2, and `solve()` still inverts Φ_s. Then `tidy(what = "edges")` keeps `beta` finite on that level's edges. It emits exactly one cli warning for that level, and the message names `k = <level>`, the smallest eigenvalue, and `beta`. This holds for a level-2 plant (a 2 by 2 Φ_s) and a level-3 plant (a 3 by 3 Φ_s) of `ackwards(sim16, k_max = 4)`. Each plant rebuilds every stored edge matrix that involves the planted level from the planted weights. The level-2 plant's smallest eigenvalue lies between 5e-3 and 1e-2. The level-3 plant's lies between 1e-8 and 1e-2. On both, `beta` equals `solve(Phi_s, E)` within a relative tolerance of 1e-8.
- [x] AC2: The near-singular warning stays silent in these cases. A planted level-2 control whose smallest eigenvalue lies between 1e-2 and 1.5e-2 raises no warning from `tidy(what = "edges")`. On the AC1 level-2 plant, `tidy(what = "variance")` and `summary()` raise no warning and keep `r2` finite. The healthy fits are `ackwards(sim16, k_max = 4)` and `ackwards(bfi25, k_max = 8, rotation = "promax")` with `engine = "pca"` and `"efa"`. On each, `tidy(what = "edges")`, `tidy(what = "variance")`, and `summary()` raise no warning.
- [x] AC3: The existing singular path takes precedence, and each path warns once per level per call. If `solve()` rejects Φ_s, `beta` and `r2` stay `NA` with the existing "cannot be inverted" warning and no near-singular warning. Under `pairs = "all"`, a near-singular level that starts two or more stored pairs raises exactly one near-singular warning per `tidy(what = "edges")` call.
- [x] AC4: Two help sections each state five facts. They are the `beta` entry of `?tidy.ackwards` and the Caution tier of "When to trust the result" in `?ackwards`. The facts are the 1e-2 smallest-eigenvalue cutoff, that the package chose it rather than a published rule, and that `beta` is still reported. The other two are that a warning names the level, and that `tidy()` raises it, not `ackwards()`. The opening sentence of that section names `tidy()` as well as `ackwards()`. The `beta` entry no longer says that the `r` of an `ackwards()` object can come from materialised scores. `NEWS.md` carries an entry under the development version.
- [x] AC5: On the committed tree, `devtools::check()` reports 0 errors, 0 warnings, and 0 notes.

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T1, T2
- AC4 → T3
- AC5 → T3

## Tasks

- [x] T1: Write the tests first in `tests/testthat/test-partialled-edges.R`, beside the singular block (line 106). Rebuild each plant's `score_var` and edge matrices from the planted weights, as the `lm()` oracle test does (lines 60 to 72). Add the level-2 and level-3 plants and the level-2 control. Add the three healthy fits, with `skip_if_not_installed("GPArotation")` on the EFA promax fit. Add the precedence case and the `pairs = "all"` case. In every case, assert the measured eigenvalue, the warning count, and which message. Update the direct helper test (line 146) to the new return value.
- [x] T2: In `.partialled_edges()` (R/utils.R:216), compute the smallest eigenvalue of Φ_s after a successful `solve()`. Below the 1e-2 cutoff, warn, unless the caller turns the warning off. `.tidy_variance()` (R/tidy.R:443) turns it off. Return a flag, so that `.tidy_edges()` dedups both warnings per level. Its dedup keys on `anyNA(B)` today (R/tidy.R:248). Define the cutoff once, as an internal constant beside the helper.
- [x] T3: Update the roxygen in R/tidy.R (`beta`) and R/ackwards.R (trust section opening at line 300, Caution tier). Before you correct the scores-path sentence (R/tidy.R:35), check that `ackwards()` always passes `edge_method = "auto"` with no data (R/ackwards.R:940 and 1023). Correct the matching sentence in DESIGN's Known limitations, marked `corrected M91`. Run `devtools::document()` and add the NEWS entry. Run `Rscript tools/check-prose.R` on each edited doc file, then `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-30: created by /milestone-plan, from the candidate row added at the M89 review (finding 14) and re-rated at the M90 review (finding F10).
- 2026-09-30: criteria audit ran in full mode twice (fresh-context Opus readers). The first pass returned 12 findings on an earlier draft and the second returned 7 on this one. All were fixed at the plan gate, and no criterion carries an open finding. AC1 now excludes the deepest level, asserts the warning count, and uses a relative tolerance. The plants rebuild their edge matrices, a level-3 plant joins, and the scores-path plant is gone (no `ackwards()` object takes that path). Plants and controls sit at the cutoff's edge. Test wording moved from the criteria to T1. AC5 narrowed from the gate script to `devtools::check()`, and the GP1 trade is stated.
- 2026-09-30: plan gate chose a 1e-2 cutoff over the 1e-4 item-matrix cutoff. At 4.5e-4, a planted level reached a `beta` of about 33 with no warning. It also rejected a variance-inflation cutoff of 10, which needs an outside source and is disputed. Falsified by a real fit that warns at 1e-2 while its `beta` stays stable across resamples.
- 2026-09-30: plan gate chose to warn from the edge table only, not also from the variance table and `summary()`. On consistent plants, `r2` stayed between 0.49 and 0.99 down to an eigenvalue of 5e-9. Falsified by a fit where `r2` from a near-singular level departs from the `lm()` R-squared oracle.
- 2026-09-30: plan chose to warn and keep the values over returning `NA`. The package reports values beside its cautions, and `r2` stays valid.
- 2026-09-30: implement started on branch `m091-near-singular-score-cor-guard`. The question gate was skipped, because the plan gate fixed the cutoff, the warning site, and the keep-values behavior.
- 2026-09-30: T1 tests written in `test-partialled-edges.R`. Before the code change, 4 of 11 tests fail as intended: the level-2, level-3, and `pairs = "all"` plants, and the helper's new `status` field. The silence and precedence tests pass on the old code. Plants sit at a smallest eigenvalue of 0.0083 (level 2), 2.5e-7 (level 3), and 0.0111 (control).
- 2026-09-30: T2 done. `.partialled_edges()` gains `warn_near` and a `status` field ("ok", "near_singular", "singular"), with the cutoff in `.phi_s_near_singular` (R/utils.R). `.tidy_edges()` dedups on `status`, and `.tidy_variance()` passes `warn_near = FALSE`. `test-partialled-edges.R` passes 11 of 11. Planted cutoffs of 0.02 and 0.005 each turned tests red (the control, then the level-2 and `pairs = "all"` plants), and the cutoff is restored. Full suite with `TESTTHAT_CPUS=8`: 805 tests, 0 failed, 2 skipped (both "On Mac").
- 2026-09-30: T3 done. The `beta` entry of `?tidy.ackwards` and the trust section of `?ackwards` (opening sentence and a new Caution item) state the cutoff, that the package chose it, that `beta` is still reported, that the warning names the level, and that `tidy()` raises it. The scores-path sentence is gone from `?tidy.ackwards`, after a read of R/ackwards.R:937 and :1020 and R/compute_edges.R:78. DESIGN's Known-limitations entry is corrected in place. NEWS entry added. With `warn_near = FALSE` removed from `.tidy_variance()`, the silence test failed twice (variance table and `summary()`), so the NEWS claim is test-backed. `tools/check-prose.R` is clean on the edited files. `Rscript tools/dod-gate.R` on commit `247a948`: GATE PASSED (check 0 err/0 warn/0 note, coverage 100%, styler, lintr, prose, pkgdown clean).
- 2026-09-30: claim audit: 50 claims read, 4 corrected — NEWS.md, R/ackwards.R, R/tidy.R, R/utils.R, man/ackwards.Rd, man/tidy.ackwards.Rd, tests/testthat/test-partialled-edges.R. The corrections: two test comments (the rcond reasoning, and "untouched levels"), the `plant_weights()` header, and the trust section's opening, which now also covers the singular warning that `tidy()` raises. The same reader re-read the four once and refined two wordings, which were applied. Only roxygen and comments changed after the gate. Full suite afterwards: 805 tests, 0 failed. Status set to review. Falsified by a near-singular level whose returned `beta` departs from an independent regression fit beyond tolerance.
- 2026-09-30: review started (checkpoint, half-done). AC1 to AC4 have evidence and ticks. The gate script for AC5 and the three reviewers are still running.

## Decisions

## Review

Review run 2026-09-30 on branch head `3a32040`. `origin/master` had not moved since the branch was cut (`9f161ec`), and no PR existed, so the resume route was (d).

- **AC1 evidence.** `devtools::test(filter = "partialled-edges")`: 11 tests, 128 expectations, 0 failed, 0 skipped. The level-2 and level-3 plant tests assert one warning whose text names `k = <level>`, the eigenvalue, and `beta`. They also assert a finite `beta` equal to `solve(Phi_s, E)` at tolerance 1e-8. An independent probe measured the planted smallest eigenvalues. Level 2 is 0.00834 (inside 5e-3 to 1e-2), and level 3 is 2.5e-7 (inside 1e-8 to 1e-2). In a scratch copy of the package, a planted cutoff of 5e-3 turned the level-2 plant test and the `pairs = "all"` test red.
- **AC2 evidence.** The same run passes the silence test and both healthy-fit tests. The probe measured the level-2 control at 0.01106 (inside 1e-2 to 1.5e-2). On the level-2 plant, the test asserts zero warnings and finite `r2` from `tidy(what = "variance")` and `summary()`. The probe measured the smallest eigenvalue over levels 2 and up for each healthy fit. It is 1 for `sim16` varimax, 0.391 for `bfi25` PCA promax, and 0.288 for `bfi25` EFA promax. In the scratch copy, a planted cutoff of 2e-2 turned the silence test red. Passing `warn_near = TRUE` from `.tidy_variance()` also turned it red.
- **AC3 evidence.** The same run passes the precedence test and the `pairs = "all"` test. The precedence test plants an exact copy of a level-2 weight column, so `solve()` rejects Φ_s. Each of `tidy(what = "edges")` and `tidy(what = "variance")` then raises one "cannot be inverted" warning and no near-singular warning, with `beta` and `r2` all `NA`. Under `pairs = "all"`, level 2 starts the stored pairs `2:3` and `2:4`. Two calls each raise exactly one near-singular warning. In the scratch copy, the old dedup on `anyNA(B)` turned the `pairs = "all"` test red.
- **AC4 evidence.** `tools::Rd2txt()` renders both help sections from the committed `man/` files. The `beta` entry of `?tidy.ackwards` states all five facts. It gives the cutoff "below `1e-2`" and calls it "a numerical guard", then says "It is not a published rule". It says "still reports `beta`" and "a warning that names the level". It also says "Fitting with `ackwards()` does not raise it". The Caution item in "When to trust the result" of `?ackwards` states the same five facts in the same words. The opening sentence of that section names both `ackwards()` and `tidy()` as the sources of its diagnostics. A grep for "materialis" in `R/tidy.R` and `man/tidy.ackwards.Rd` finds nothing. `NEWS.md` carries the entry under "# ackwards (development version)".
- **AC5 evidence.** `Rscript tools/dod-gate.R` on the committed code tree (`3a32040`, same code as checkpoint `35ec60d`) printed "GATE PASSED". Its `devtools::check()` line reads 0 errors, 0 warnings, and 0 notes in 62 s. The same run reports coverage 100%, styler and lintr clean, prose clean, vignette freshness clean, and the pkgdown reference index complete.
- **Consistency gate.** `cairn_validate.py` exits 0. Its 16 advisory warnings are all work-log format lines in M84. No DESIGN principle changed, so `cairn_impact` does not apply. `devtools::document()` leaves no diff. The branch does not touch README.Rmd or README.md, and README.md is newer. `pkgdown::check_pkgdown()` passes inside the gate. NEWS has the entry and names no milestone. The branch adds no new file, so no `.Rbuildignore` entry is owed.
