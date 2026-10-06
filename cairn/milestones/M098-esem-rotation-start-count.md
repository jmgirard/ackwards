# M098: Counted lavaan rotation failures on every ESEM rotation

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP6, IP7
- **Resolves:** —
- **Surface tier:** user-facing — it changes the warnings and levels that `ackwards(engine = "esem")` returns
- **Branch/PR:** `m098-esem-rotation-start-count`

## Goal

ESEM counts lavaan's failed rotation starts under every rotation, varimax included, and ends the
hierarchy at a level where every start failed.

## Scope

**In:** These changes are in scope.

- `.esem_rotation_args()` (`R/engine_esem.R:64`) turns lavaan's rotation warnings on for varimax.
  The k = 1 level (`"none"`) does not change.
- `.esem_fit_one()` counts the "rotation algorithm did not converge" warnings at each level. It
  reads the start total from the fit's lavaan options (`rotation.args$rstarts`).
- If some starts failed, the level is kept with a counted warning. If every start failed, the
  hierarchy ends at the level before. One start counts as every start, as for PCA in M90.
- The same rule applies to `"oblimin"` and `"geomin"`.
- Tests, the `rotation` paragraph of `?ackwards`, NEWS, DESIGN §4, and the code comments at
  `R/engine_esem.R:57` and `:366`.

**Out:** These items are not in this milestone.

- EFA and PCA varimax. psych calls `stats::varimax()`, which reports no convergence, so there is
  nothing to show.
- An `ackwards()` argument for lavaan's `rstarts` or `max_iter`. No user asked for one.
- The deprecated-formal detection in `.esem_rotation_args()`. Its `[low]` candidate row stays.

## Acceptance criteria

- [x] AC1: If `f` contains `"rotation_args"`, `.esem_rotation_args("varimax", f)` returns
      `list(rotation = list("varimax", warn = TRUE))`. Else it returns
      `list(rotation = "varimax", rotation.args = list(warn = TRUE))`.
      `.esem_rotation_args("none", f)` returns `list(rotation = "none")` for `f` with and
      without `"rotation_args"`. A unit test in `tests/testthat/test-oblique-branches.R` asserts
      these four results with `expect_identical()`, and its oblimin and geomin assertions stay.
- [x] AC2: Some starts fail at a level, but not all. Then `ackwards()` keeps that level, and
      exactly one of its warnings contains both `k = 3` and `"<n> of <N> random starts"`. Here
      `<n>` is the count of lavaan non-convergence warnings at that level, and `<N>` is
      `rotation.args$rstarts` from the lavaan options of that fit. A test plants this case at
      k = 3 of `k_max = 4` for `rotation = "varimax"` and for `rotation = "oblimin"` (ML, sim16).
      The plant forwards the real `lavaan::efa` fit and then raises a fixed count `n` of
      warnings in lavaan's exact text, with `0 < n < N`. The test asserts that level 3 is
      present and that the warning shows the planted `n`.
- [x] AC3: Every start fails at level k. Then `ackwards()` returns levels `1..k-1` only. One
      warning contains `k = 3`, contains `"<N> of <N>"` with `<N>` read as in AC2, and says
      that the kept rotation did not converge. A test plants `max_iter = 2` at k = 3 of
      `k_max = 4` in the real lavaan call for four fits: varimax, oblimin, and geomin on ML
      (sim16), and varimax on WLSMV (`bfi25`, `cor = "polychoric"`). For each fit it asserts
      that `names(x$levels)` is `c("1", "2")` and the warning text above. These plants replace
      the M90 test "lavaan's rotation warning is shown under an oblique rotation only".
- [x] AC4: Two named default fits raise no rotation warning. A test asserts that no warning
      that contains `"random starts"` or `"rotation algorithm"` comes from
      `ackwards(sim16, k_max = 4, engine = "esem", seed = 1)` or from
      `ackwards(bfi25, k_max = 4, engine = "esem", cor = "polychoric", seed = 1)`. Separately,
      `tests/testthat/test-baseline-m89.R` passes with its fixture unchanged. That file
      suppresses warnings, so it shows unchanged values and level names only.
- [x] AC5: The rule is documented in three places. The `rotation` paragraph of `?ackwards`
      says that ESEM reports failed lavaan rotation starts under every rotation, varimax
      included, keeps a level where some starts failed, and ends the hierarchy at the level
      before one where every start failed. A NEWS.md bullet under the development version says
      the same and names varimax as newly covered. DESIGN §4 states the keep and truncate rule
      and that `<N>` is lavaan's `rstarts`.
- [x] AC6: `Rscript tools/dod-gate.R` exits 0.

## Coverage

- AC1 → T2
- AC2 → T3, T4
- AC3 → T1, T3, T4
- AC4 → T4
- AC5 → T5
- AC6 → T6

## Tasks

- [x] T1: Measure before coding, on lavaan 0.7.2. For ML, MLR, WLSMV (ordered), and ULSMV,
      make sure that a `max_iter = 2` plant gives exactly `rstarts` (30) non-convergence
      warnings at one level under varimax, oblimin, and geomin. A second rotation pass doubles
      the count. Log the figures in the work log. If any estimator doubles the count, stop and
      amend the rule.
- [x] T2: In `.esem_rotation_args()`, move varimax to the warnings-on branch. Update the comment
      at `R/engine_esem.R:57` and the unit test (AC1).
- [x] T3: In `.esem_fit_one()`, count the non-convergence warnings for every rotation except
      `"none"`. Read `rstarts` from `lavaan::lavInspect(fit, "options")`. If
      `n >= max(rstarts, 1)`, return `status = "nonconverged"` with the count, and show the
      total as `max(rstarts, 1)`. Else return `"ok"` with `n` and `N`. In the assembly loop
      (`R/engine_esem.R:505` and `:525`), word both warnings with the count. The all-fail
      warning says that the kept rotation did not converge. Remove the `# nocov` marks on the
      nonconverged branch that the new tests reach. Update the comment at `:366`.
- [x] T4: Write AC2's partial plants, AC3's four all-fail plants, and AC4's assertions. As in
      the M90 test (`tests/testthat/test-oblique-branches.R:351`), pin `.esem_rotation_args()`
      to the real `lavaan::efa` formals and run `.esem_lapply` as serial `lapply` in every
      plant. The AC3 plants replace the M90 test.
- [x] T5: Edit the `rotation` paragraph of `?ackwards` (`R/ackwards.R:208`). Read the `seed`
      paragraph (`R/ackwards.R:161`) and make it agree. Run `devtools::document()`, and add the
      NEWS.md bullet and the DESIGN §4 line. Run `Rscript tools/check-prose.R` on each edited
      doc file.
- [x] T6: Run `Rscript tools/dod-gate.R`.

## Work log

- 2026-10-06: created by /milestone-plan.
- 2026-10-06: question set: what to plan, with no work named — the ESEM varimax rotation-warning candidate row (row absorbed and removed in the plan commit).
- 2026-10-06: question set: failure handling for every lavaan rotation — count failed starts, keep the level when some fail, truncate at the level before when all fail.
- 2026-10-06: plan gate chose count-and-truncate over warn-and-keep with M90's wording because that wording is false when every start failed and IP7 truncates a non-converged level; falsified by a lavaan version that reports the kept start's convergence directly.
- 2026-10-06: measured on lavaan 0.7.2 — default varimax ESEM on bfi25 gave 0 rotation warnings at k = 2 to 8 (ML) and 2 to 6 (WLSMV); a max_iter = 2 plant gave 30 of 30, max_iter = 30 gave 1 of 30; warn = FALSE gave 0.
- 2026-10-06: criteria audit (full mode, fresh Opus reader) returned 10 findings, all accepted: AC1 names the two f forms and keeps the oblique asserts; AC2 uses a deterministic warning plant instead of a seeded max_iter or a second real call; mock pinning and serial lapply moved to T4 as instrument properties; AC2 warning matched by k; AC3 uses k_max = 4, adds geomin, reads N from options, drops the unreachable rstarts = 0 case to T3; AC4 narrowed to the two named fits and states what the baseline test shows; AC5 DESIGN content named; T5 checks the seed paragraph.
- 2026-10-06: collision sweep — absorbs the [low] candidate "ESEM varimax hides lavaan's rotation non-convergence"; extends M90 (archive); no D-entry rejects it; GitHub inbox has 0 open issues and 0 open PRs.
- 2026-10-06: T1 (lavaan 0.7.2, max_iter = 2 at k = 3): varimax, oblimin, and geomin each gave 30 of 30 warnings under ML, MLR, ULSMV (sim16) and WLSMV, ULSMV (bfi25 ordered); lavInspect converged stayed TRUE; rstarts = 0 gave 1. No doubling, so the counting rule stands.
- 2026-10-06: T2-T4 done. `.esem_rotation_starts()` reads rstarts and treats 0 or an unreadable value as one start, so an unreadable total truncates loudly rather than hiding a failure. The kept-level warning drops lavaan's raw text for the count and says the kept start can be one that did not converge. The model non-convergence branch keeps its `# nocov`, and the rotation branch beside it is now covered.
- 2026-10-06: tests for AC2-AC4 live in the new `tests/testthat/test-esem-rotation-starts.R`, and `.warnings_of()` moved to `helper-data.R` so both files share it. Each loop iteration scopes its plant with `local()`, because a second plant in one test otherwise wraps the first mock. Against master's engine the partial and all-fail tests fail 5/12 and 20/24 expectations. Full suite: 820 tests, 0 failed, 2 skipped (pre-existing On Mac skips).
- 2026-10-06: T5 done. The `rotation` help splits the ESEM rule into its own paragraph and scopes the PCA and EFA sentences to oblique rotation, as before. The `seed` paragraph already says every lavaan rotation starts at random, so it needed no edit. NEWS bullet and DESIGN §4 paragraph added. check-prose clean on R/ackwards.R and NEWS.md.
- 2026-10-06: claim audit: 28 claims read, 5 corrected — NEWS.md, R/engine_esem.R, tests/testthat/test-esem-rotation-starts.R
- 2026-10-06: claim-audit fixes: NEWS now says fits where at least one start converges keep the same levels and values; the truncation warning drops "random" so it holds for rstarts = 0; two code comments and one test name made exact; the AC4 test also asserts no "did not converge" warning. The same reader's one re-read found all five hold.
- 2026-10-06: T6 done. `Rscript tools/dod-gate.R` exit 0: vignette freshness clean, check 0 errors, 0 warnings, 0 notes, coverage 100%, styler and lintr clean, pkgdown index complete. The first gate run failed only on a styler trailing-blank-line fix in test-oblique-branches.R. Status set to review.

## Decisions

- 2026-10-06 (review): ESEM rotation failure rule. M90 kept an ESEM level whatever lavaan's rotation warnings said, and kept varimax silent so that the default path did not change. M098 replaces both at the user's plan-gate choice (work log, question set): every lavaan rotation counts its failed starts, and a level where every start failed truncates under IP7, because the start lavaan keeps then did not converge. Default varimax output changes only for a fit where all 30 starts fail, and no measured fit did (bfi25 k = 2 to 8 ML and 2 to 6 WLSMV, the two AC4 fits, the engines vignette fit three times). D-034's "no current numerical output changes" bound the oblique decision, not later convergence handling, so no D-entry supersedes it. Kept milestone-local, as M90 kept its rule.

## Review

- AC1: fresh run (2026-10-06, HEAD 0b7e7d8) of test-oblique-branches.R: ".esem_rotation_args() turns lavaan's rotation warnings on for every rotation" passes 6 of 6 expect_identical() calls. The four named in AC1 (none and varimax, each with and without `"rotation_args"`) are there, and the oblimin and geomin calls stay. PASS.
- AC2: fresh run of test-esem-rotation-starts.R "a level where some rotation starts fail is kept": 12 expectations, 0 failed. For varimax and oblimin (ML, confirmed `x$meta$estimator` is ML on sim16), k = 3 of k_max = 4, 7 planted warnings: levels 1 to 4 present, exactly one "random starts" warning, and it carries "at k = 3:" and "7 of <N> random starts did not converge" with N read from `x$fits[["3"]]` options and N > 7. Against master's engine this test failed 5 of 12 (implement log). PASS.
- AC3: fresh run of "a level where every rotation start fails ends the hierarchy": 24 expectations, 0 failed, over varimax, oblimin, and geomin on ML (sim16) and varimax on WLSMV (bfi25, polychoric). Each asserts the estimator, `names(x$levels)` = c("1", "2"), one warning containing "rotation did not converge at k = 3", "<N> of <N> starts did not converge", and "the kept rotation did not converge". The test reads N from the level 2 fit's options because level 3 is dropped. The code reads the warning's N from the level 3 fit (`.esem_rotation_starts(fit)` in `.esem_fit_one()`), with the same options. The M90 test is gone from test-oblique-branches.R. Against master this test failed 20 of 24. PASS.
- AC4: fresh run of "default varimax fits whose rotation converges raise no rotation warning" (the two named fits): 3 expectations, 0 failed. It asserts no "random starts", "rotation algorithm", or "did not converge" text. test-baseline-m89.R: 5 tests, 0 failed (esem_sim16 21 of 21). `git diff --quiet master -- tests/testthat/fixtures/baseline-m89.rds` reports no change. PASS.
- AC5: read at HEAD. The `rotation` help (R/ackwards.R, "ESEM counts the failed starts under every lavaan rotation, varimax included") states the keep and truncate rule. The NEWS.md first development bullet states the same rule and says varimax ESEM showed no rotation warning before. The DESIGN §4 "ESEM rotation starts" paragraph states the rule and that N is `rotation.args$rstarts`. PASS.
- AC6: fresh `Rscript tools/dod-gate.R` at HEAD dc0541e exited 0: vignette freshness clean, check 0 errors, 0 warnings, 0 notes, coverage 100%, styler and lintr clean, pkgdown index complete. PASS.
- AC6 (re-run after the review fixes at f1148ca): `Rscript tools/dod-gate.R` exit 0 at HEAD 4a9b3b4: vignette freshness clean, check 0 errors, 0 warnings, 0 notes, coverage 100%, styler and lintr clean, pkgdown index complete. Full suite after the fixes: 820 tests, 0 failed, 2 skipped (On Mac). PASS.
- consistency gate: `cairn_validate.py` exit 0 (16 work-log WARNs, all in M84's file). No principle text changed, so `cairn_impact` is skipped. `devtools::document()` (run by the gate) left `man/` and NAMESPACE with no diff. README untouched. pkgdown index complete. NEWS entry present. check clean.
- spawned: diff-bug, blame-history, prior-review
- diff-bug #1: the start count matches a fixed phrase, so a line break inside lavaan's message would hide a failed start — fix now (whitespace collapsed before matching, planted wrapped warnings in the AC2 test), fixed f1148ca
- diff-bug #2: one warning per failed start was measured only on lavaan 0.7.2, not on the pre-0.7 `rotation.args` path — follow-up (extended the `.esem_rotation_args()` [low] candidate row)
- diff-bug #3: a kept level's failed-start count lives only in the warning, not on the object — follow-up (new [low] candidate row "Record a partial ESEM rotation failure on the object")
- diff-bug #4: the AC3 test reads N from the level 2 fit, not the dropped level 3 fit — reject, false as a defect: every level gets the same lavaan options, and the code reads the warning's N from the level 3 fit (AC3 evidence)
- diff-bug #5: the partial test plants warnings rather than real failed starts — reject, planned change: the plan audit chose the deterministic plant, and AC3's real-lavaan plants show real failures reach the counter under varimax
- diff-bug #6: kept varimax fits now store `warn = TRUE` in their options — reject, planned change: the plan turns the warnings on for varimax, and oblique fits stored it since M90
- diff-bug #7: "1 of 1 starts" and the repeated "did not converge" read awkwardly — reject, style
- diff-bug #8: "30 random starts by default" implies an `ackwards()` setting, and "each level" includes the unrotated k = 1 — fix now (help and NEWS say each rotated level, lavaan's default, not changed by `ackwards()`), fixed f1148ca
- diff-bug #9: NEWS does not say partial-failure varimax fits now warn — reject, false: the bullet says varimax showed no rotation warning before and that a warning now gives the count
- diff-bug #10: precomputed vignettes were not regenerated — reject, false as a defect: the engines vignette's ESEM fit raised 0 rotation warnings in 3 branch runs, and `warn = TRUE` left a seeded fit's `est.std` identical
- diff-bug #11: the unreadable-rstarts fallback is `# nocov` and untested — fix now (nocov removed, unit test for an empty, non-numeric, and erroring options read), fixed f1148ca
- diff-bug #12: the AC boxes were unticked at review time — reject, false: each box is ticked against its evidence line above
- diff-bug #13: the reversal of M90's keep rule has no D-entry — fix now (milestone-local Decisions entry giving the rule, its source, and why D-034 is not superseded), fixed f1148ca
- blame-history #1: M90's NEWS bullet in the same development section still says ESEM always keeps the level — fix now (bullet now says EFA keeps the level and points ESEM to the new bullet), fixed f1148ca
- blame-history #2: M90 kept varimax silent so the default path stayed unchanged, and D-034 says no current numerical output changes — fix now (Decisions entry answers both), fixed f1148ca
- blame-history #3: an oblique level where every start fails now truncates, where M90 kept it — reject, planned change (AC3, the user's plan-gate choice)
- blame-history #4: new `# nocov` on the reachable fallback — fix now (same fix as diff-bug #11), fixed f1148ca
- blame-history #5: the removed M90 test checked no warning at k = 2, and the all-fail test does not — fix now (each all-fail plant asserts no "at k = 2" warning), fixed f1148ca
- blame-history #6: `.warnings_of()` in helper-data.R could shadow a local definition — reject, false: grep finds no other definition
- blame-history #7: the kept-level warning drops lavaan's raw text — reject, requests nothing: the planned count replaces it and the old prefix stays
- prior-review #1: `# nocov` on the fallback regresses M90's AC13 finding — fix now (same fix as diff-bug #11), fixed f1148ca
- prior-review #2: the NEWS claim of unchanged values was not checked — reject, false: `est.std` is identical with `warn = TRUE` and `FALSE` (seeded bfi, 5 factors), and the baseline ESEM test passes
- prior-review #3: "every lavaan rotation" overstates scope at the unrotated k = 1 — reject, false: k = 1 is not rotated (the help now says each rotated level)
