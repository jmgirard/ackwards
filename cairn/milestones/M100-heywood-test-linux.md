# M100: Heywood ESEM test green on Linux under lavaan 0.7.3

- **Status:** planned
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP7
- **Resolves:** #111 closes
- **Surface tier:** user-facing — the test ships in the tarball and runs under R CMD check on CRAN, and the fix can change the ESEM engine
- **Branch/PR:** —

## Goal

The ESEM Heywood test passes again on the GitHub Ubuntu runners under lavaan 0.7.3, so the master CI and the weekly upstream watch are green.

## Scope

**In:** Reproduce the failure on Linux. Name the failing level and the condition at that level. Repair the test plant. If the cause is in package code, repair the engine, add a NEWS entry, and add a test that fails before the fix. Keep the test able to fail. Dispatch the watch workflow on the branch.

**Out:** CRAN flavors: all 13 report 0.2.0 OK (read 2026-10-18), so no action. The pre-0.7 `rotation.args` candidate row stays as it is. If the diagnosis shows a lavaan defect, an upstream lavaan report becomes a candidate row.

Context. Since lavaan 0.7.3 (CRAN 2026-10-08), the test "ESEM warns on improper solution (Heywood case) but still builds" (`tests/testthat/test-esem.R:535`) errors on ubuntu-latest R 4.6.1. The error is "The esem engine failed to build at least 2 converged levels (k_eff = 1 ...)". Watch run 38267912233 and master run 39166011797 show it. The test passes on macOS with R 4.6.1, lavaan 0.7.3, and psych 2.6.9. Watch run 37323818672 with lavaan 0.7.2 and psych 2.6.9 passed its tests. The engine ends the hierarchy at a level where the lavaan model fit did not converge. It also ends it where every rotation start warned "did not converge" (`R/engine_esem.R:151-171`, from M098). The field `x$k_max` holds the built depth (`R/ackwards.R:1129`).

## Acceptance criteria

- [ ] AC1: The Heywood test is the `test_that()` block in `tests/testthat/test-esem.R` whose `expect_warning()` matches "Heywood". At the branch head, under lavaan 0.7.3 or later, it passes and does not skip in two environments. The first is macOS local. The second is a Linux environment in which the block as it exists at master commit `46bb2f6` fails with the error "failed to build at least 2 converged levels". That environment is a Docker container or a GitHub Actions `ubuntu-latest` job. The evidence for each run is `as.data.frame(testthat::test_file(...))`, run with `NOT_CRAN=true` and filtered to the block, with `failed` 0 and `skipped` FALSE.
- [ ] AC2: You dispatch `re-oracle.yaml` on the branch with `gh workflow run re-oracle.yaml --ref <branch>`. The run ends with conclusion `success`, and its `headSha` equals the final branch head. Its "Record upstream versions" step prints a `lavaan=` version of 0.7.3 or later.
- [ ] AC3: The repaired plant puts the Heywood case at the deepest requested level, and the built depth stays at 2 or more. On macOS local, the block fails under each of three uncommitted edits to `R/engine_esem.R`. First, an edit removes the Heywood warning, and the block fails at its `expect_warning(..., "Heywood")`. Second, an edit changes the residual test from `<= 0` to `< 0`, and the block fails at the same expectation. Third, an edit keeps the warning but drops the Heywood level from the hierarchy. The failures of the block then include its check that the built depth equals the requested depth. For each edit, the evidence shows that expectation among the failures, not an error or a skip.
- [ ] AC4: One of two states holds. In the first, the branch changes no file under `R/` and leaves `NEWS.md` unchanged, so `git diff --stat master...HEAD -- R/ NEWS.md` prints nothing. In the second, the branch changes at least one file under `R/`. Each behavior-changing hunk under `R/` then comes with a `NEWS.md` bullet under the development version that describes the changed behavior. Each such hunk also comes with a test that fails at the master `R/` and passes at the branch head in the AC1 Linux environment. If the change alters when ESEM ends the hierarchy, a D-entry in `cairn/DECISIONS.md` records it.
- [ ] AC5: `Rscript tools/dod-gate.R` exits 0 on macOS local at the branch head.

## Coverage

- AC1 → T1, T3
- AC2 → T5
- AC3 → T2, T4
- AC4 → T2
- AC5 → T6

## Tasks

- [ ] T1: Reproduce on Linux. Pull `rocker/r-ver:4.6.1`, install lavaan 0.7.3, psych 2.6.9, and the package dependencies, and run the block at `46bb2f6`. If Docker does not reproduce the error, push a temporary branch-scoped GitHub Actions job that runs `testthat::test_file()` on ubuntu-latest (the M66 lesson), and revert it before review. In the reproducing environment, run a diagnostic that names the failing level and its condition: lavaan fit not converged, every rotation start failed, or a fit error. Write the result as a work-log line.
- [ ] T2: Repair from the diagnosis. For a fragile plant, build a plant whose Heywood case is at the deepest requested level and whose lower levels converge on both platforms. For a package-code cause, change `R/engine_esem.R` with a NEWS bullet, a before-and-after test, and a D-entry if the truncation rule changes (AC4). Rewrite the comment above the block so it says why the plant is built as it is.
- [ ] T3: Run the block at the branch head on macOS local and in the T1 Linux environment, and record the filtered `test_file()` results (AC1).
- [ ] T4: Run the three uncommitted probe edits of AC3 on macOS local, one at a time, and record which expectation fails under each. Restore `R/engine_esem.R` after each.
- [ ] T5: Push the branch and dispatch `re-oracle.yaml` on it (AC2). If an unrelated upstream break fails the run, stop at the red-CI stop.
- [ ] T6: Run `Rscript tools/dod-gate.R` on macOS local (AC5).

## Work log

- 2026-10-18: created by /milestone-plan.
- 2026-10-18: tier check: not run — work named at the question set.
- 2026-10-18: question set: what to work on — issue #111 as a milestone.
- 2026-10-18: question set: Linux reproduction — pull rocker/r-ver:4.6.1 (about 1 GB), fall back to temporary CI pushes, and push the branch early to dispatch the watch.
- 2026-10-18: question set: cause in package code — fix it in this milestone with NEWS, a before-and-after test, and a D-entry if the truncation rule changes, with no stop.
- 2026-10-18: question set: assign #111 to the operator's GitHub account — yes.
- 2026-10-18: investigation: lavaan 0.7.2 to 0.7.3 is the only version change between passing watch run 37323818672 and failing run 38267912233. The block passes on macOS with lavaan 0.7.3 and psych 2.6.9 in a scratch library. `formals(lavaan::efa)` is the same in 0.7.2 and 0.7.3.
- 2026-10-18: inbox sweep: 1 open issue (#111, taken as `Resolves: #111 closes`), 0 open PRs.
- 2026-10-18: collision check: no D-entry or planned milestone overlaps. M098 (archived) added the all-starts-failed rule that T1 must rule in or out. The `.esem_rotation_args()` candidate row does not fire because lavaan 0.7.3 keeps the `rotation_args` formal.
- 2026-10-18: lessons harvested: M66 (temporary push trigger for a CI-only job), M83 (measure the rate of an intermittent failure), M89 (seed an ESEM oracle), M098 (`local_mocked_bindings()` in loops).
- 2026-10-18: criteria audit (full mode, fresh Opus reader) returned 13 findings, all fixed toward the narrower promise. The Linux evidence became a run that reproduces the failure, with per-block results. AC2 now binds `headSha`. AC3 now places the Heywood level, names the depth check, and adds the `< 0` probe. AC4 now needs at least one behavior-changing `R/` hunk and a D-entry for a truncation change. The test-comment criterion was dropped as a recording act, and the diagnosis moved to T1's work-log line.
- 2026-10-18: plan gate chose a Docker reproduction over CI-only runs because one CI try takes 10 to 15 minutes; falsified by Docker not reproducing the failure, which T1's CI fallback covers.
