# M098: Counted lavaan rotation failures on every ESEM rotation

- **Status:** planned
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP6, IP7
- **Resolves:** —
- **Surface tier:** user-facing — it changes the warnings and levels that `ackwards(engine = "esem")` returns
- **Branch/PR:** —

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

- [ ] AC1: If `f` contains `"rotation_args"`, `.esem_rotation_args("varimax", f)` returns
      `list(rotation = list("varimax", warn = TRUE))`. Else it returns
      `list(rotation = "varimax", rotation.args = list(warn = TRUE))`.
      `.esem_rotation_args("none", f)` returns `list(rotation = "none")` for `f` with and
      without `"rotation_args"`. A unit test in `tests/testthat/test-oblique-branches.R` asserts
      these four results with `expect_identical()`, and its oblimin and geomin assertions stay.
- [ ] AC2: Some starts fail at a level, but not all. Then `ackwards()` keeps that level, and
      exactly one of its warnings contains both `k = 3` and `"<n> of <N> random starts"`. Here
      `<n>` is the count of lavaan non-convergence warnings at that level, and `<N>` is
      `rotation.args$rstarts` from the lavaan options of that fit. A test plants this case at
      k = 3 of `k_max = 4` for `rotation = "varimax"` and for `rotation = "oblimin"` (ML, sim16).
      The plant forwards the real `lavaan::efa` fit and then raises a fixed count `n` of
      warnings in lavaan's exact text, with `0 < n < N`. The test asserts that level 3 is
      present and that the warning shows the planted `n`.
- [ ] AC3: Every start fails at level k. Then `ackwards()` returns levels `1..k-1` only. One
      warning contains `k = 3`, contains `"<N> of <N>"` with `<N>` read as in AC2, and says
      that the kept rotation did not converge. A test plants `max_iter = 2` at k = 3 of
      `k_max = 4` in the real lavaan call for four fits: varimax, oblimin, and geomin on ML
      (sim16), and varimax on WLSMV (`bfi25`, `cor = "polychoric"`). For each fit it asserts
      that `names(x$levels)` is `c("1", "2")` and the warning text above. These plants replace
      the M90 test "lavaan's rotation warning is shown under an oblique rotation only".
- [ ] AC4: Two named default fits raise no rotation warning. A test asserts that no warning
      that contains `"random starts"` or `"rotation algorithm"` comes from
      `ackwards(sim16, k_max = 4, engine = "esem", seed = 1)` or from
      `ackwards(bfi25, k_max = 4, engine = "esem", cor = "polychoric", seed = 1)`. Separately,
      `tests/testthat/test-baseline-m89.R` passes with its fixture unchanged. That file
      suppresses warnings, so it shows unchanged values and level names only.
- [ ] AC5: The rule is documented in three places. The `rotation` paragraph of `?ackwards`
      says that ESEM reports failed lavaan rotation starts under every rotation, varimax
      included, keeps a level where some starts failed, and ends the hierarchy at the level
      before one where every start failed. A NEWS.md bullet under the development version says
      the same and names varimax as newly covered. DESIGN §4 states the keep and truncate rule
      and that `<N>` is lavaan's `rstarts`.
- [ ] AC6: `Rscript tools/dod-gate.R` exits 0.

## Coverage

- AC1 → T2
- AC2 → T3, T4
- AC3 → T1, T3, T4
- AC4 → T4
- AC5 → T5
- AC6 → T6

## Tasks

- [ ] T1: Measure before coding, on lavaan 0.7.2. For ML, MLR, WLSMV (ordered), and ULSMV,
      make sure that a `max_iter = 2` plant gives exactly `rstarts` (30) non-convergence
      warnings at one level under varimax, oblimin, and geomin. A second rotation pass doubles
      the count. Log the figures in the work log. If any estimator doubles the count, stop and
      amend the rule.
- [ ] T2: In `.esem_rotation_args()`, move varimax to the warnings-on branch. Update the comment
      at `R/engine_esem.R:57` and the unit test (AC1).
- [ ] T3: In `.esem_fit_one()`, count the non-convergence warnings for every rotation except
      `"none"`. Read `rstarts` from `lavaan::lavInspect(fit, "options")`. If
      `n >= max(rstarts, 1)`, return `status = "nonconverged"` with the count, and show the
      total as `max(rstarts, 1)`. Else return `"ok"` with `n` and `N`. In the assembly loop
      (`R/engine_esem.R:505` and `:525`), word both warnings with the count. The all-fail
      warning says that the kept rotation did not converge. Remove the `# nocov` marks on the
      nonconverged branch that the new tests reach. Update the comment at `:366`.
- [ ] T4: Write AC2's partial plants, AC3's four all-fail plants, and AC4's assertions. As in
      the M90 test (`tests/testthat/test-oblique-branches.R:351`), pin `.esem_rotation_args()`
      to the real `lavaan::efa` formals and run `.esem_lapply` as serial `lapply` in every
      plant. The AC3 plants replace the M90 test.
- [ ] T5: Edit the `rotation` paragraph of `?ackwards` (`R/ackwards.R:208`). Read the `seed`
      paragraph (`R/ackwards.R:161`) and make it agree. Run `devtools::document()`, and add the
      NEWS.md bullet and the DESIGN §4 line. Run `Rscript tools/check-prose.R` on each edited
      doc file.
- [ ] T6: Run `Rscript tools/dod-gate.R`.

## Work log

- 2026-10-06: created by /milestone-plan.
- 2026-10-06: question set: what to plan, with no work named — the ESEM varimax rotation-warning candidate row (row absorbed and removed in the plan commit).
- 2026-10-06: question set: failure handling for every lavaan rotation — count failed starts, keep the level when some fail, truncate at the level before when all fail.
- 2026-10-06: plan gate chose count-and-truncate over warn-and-keep with M90's wording because that wording is false when every start failed and IP7 truncates a non-converged level; falsified by a lavaan version that reports the kept start's convergence directly.
- 2026-10-06: measured on lavaan 0.7.2 — default varimax ESEM on bfi25 gave 0 rotation warnings at k = 2 to 8 (ML) and 2 to 6 (WLSMV); a max_iter = 2 plant gave 30 of 30, max_iter = 30 gave 1 of 30; warn = FALSE gave 0.
- 2026-10-06: criteria audit (full mode, fresh Opus reader) returned 10 findings, all accepted: AC1 names the two f forms and keeps the oblique asserts; AC2 uses a deterministic warning plant instead of a seeded max_iter or a second real call; mock pinning and serial lapply moved to T4 as instrument properties; AC2 warning matched by k; AC3 uses k_max = 4, adds geomin, reads N from options, drops the unreachable rstarts = 0 case to T3; AC4 narrowed to the two named fits and states what the baseline test shows; AC5 DESIGN content named; T5 checks the seed paragraph.
- 2026-10-06: collision sweep — absorbs the [low] candidate "ESEM varimax hides lavaan's rotation non-convergence"; extends M90 (archive); no D-entry rejects it; GitHub inbox has 0 open issues and 0 open PRs.

## Decisions

## Review
