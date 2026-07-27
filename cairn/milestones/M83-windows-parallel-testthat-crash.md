# M83: Windows parallel-testthat crash — measure, diagnose, mitigate

- **Status:** planned
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** —

## Goal

Make the Windows CI job a trustworthy signal again by measuring, diagnosing, and
mitigating the intermittent access violation that kills a testthat parallel worker.

## Scope

**In:** A dispatchable Windows stress workflow that runs the suite repeatedly and
tallies crashes; a measured baseline crash rate on current `master`; bisection
across the plausible triggers (worker count, `EFAtools::CD()`, `psych::polychoric`
→ `mnormt`, runner memory); a mitigation, preferring the CRAN-matching 2-worker
setting where the evidence permits; a measured post-fix rate; and re-verification
of the pending 0.2.0 tarball if the fix reaches shipped content.

**Out:** Filing the bug upstream if it lands in a dependency's compiled code →
ROADMAP candidate, with the reproducer and evidence attached. The two Ubuntu CI
failures of 2026-07-23 → nothing owed: diagnosed at plan time as
`setup-r-dependencies` failing on `there is no package called 'pak'`, GitHub
Action infrastructure rather than package code, and not a test failure. Any
change to the package's own runtime parallelism (`future` plans in `ackwards()` /
`boot_edges()`) → out; this milestone touches test execution only.

## Acceptance criteria

- [ ] AC1: A `workflow_dispatch` workflow runs the test suite on `windows-latest`
      a configurable number of times (default 30) within one job, prints a
      per-iteration pass/crash tally, and exits non-zero if any iteration dies
      with exit code `-1073741819`.
- [ ] AC2: A baseline run of AC1's workflow against unmitigated code records the
      crash count over ≥30 iterations; the tally is quoted in the work log.
- [ ] AC3: For each candidate trigger examined — worker count, `EFAtools::CD()`,
      `polychoric`/`mnormt`, runner memory — the milestone records the iteration
      tally or log excerpt that implicates or eliminates it. A trigger left
      unexamined is named as such rather than omitted.
- [ ] AC4: A mitigation is applied, and AC1's workflow reports zero access
      violations over ≥60 post-fix iterations. This bounds the rate rather than
      proving absence: against the AC2 baseline it is the evidence the milestone
      claims, and the milestone states the bound rather than declaring the crash
      impossible.
- [ ] AC5: The milestone states whether the mitigation changed tarball content.
      If it did, `NEWS.md` carries a user-facing line, and R-hub `atlas` + `nold`
      are re-run green on the mitigated code so the pending 0.2.0 resubmission
      rests on fresh evidence. If it is confined to `.github/`, the milestone
      records that no tarball content changed and no re-verification is owed.
- [ ] AC6: `Rscript tools/dod-gate.R` passes, and the full CI matrix is green on
      the merge commit.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4, T5
- AC5 → T6
- AC6 → T6

## Tasks

- [ ] T1: Author `.github/workflows/windows-stress.yaml` — `workflow_dispatch`
      with an `iterations` input (default 30), looping `R CMD check`'s test phase
      (or `devtools::test()`) on `windows-latest`, capturing each iteration's exit
      code and detecting `-1073741819`. Model the job on `macos-oldrel-check.yaml`
      (M82's single-flavour job).
- [ ] T2: Run it on unmitigated `master` for ≥30 iterations; record the tally.
      Re-run if zero crashes appear, since the observed rate is ~4% and a single
      30-iteration clean sweep is ~29% likely by chance.
- [ ] T3: Bisect with stress variants, one variable at a time:
      `TESTTHAT_CPUS` 4 → 2 → 1; `EFAtools` absent (its `suggest_k()` CD path
      skips); `cor = "polychoric"` paths skipped; runner memory logged per
      iteration. Record each tally against `TESTTHAT_CPUS: 4`
      ([R-CMD-check.yaml:68](.github/workflows/R-CMD-check.yaml:68)).
- [ ] T4: Apply the mitigation the T3 evidence supports, preferring
      `TESTTHAT_CPUS: 2` on the Windows job (what CRAN effectively runs) so
      M48's parallel speedup survives everywhere else.
- [ ] T5: Re-run the stress workflow for ≥60 iterations on the mitigated code.
- [ ] T6: Settle the tarball question per AC5 — `NEWS.md` line and a green R-hub
      `atlas` + `nold` re-run if shipped content changed, otherwise a recorded
      statement that it did not — then run the DoD gate and open the PR.

## Work log

- 2026-07-27: created by /milestone-plan.
- 2026-07-27: plan gate chose a before/after stress-workflow measurement over counting consecutive green CI runs because at the observed ~4% rate 20 clean runs occur by luck 44% of the time; falsified by a post-fix stress run that still crashes.
- 2026-07-27: plan gate chose implementing before the 0.2.0 resubmission over deferring until CRAN accepts, at the user's direction; falsified by a mitigation that changes tarball content and fails re-verification, which would return the release to the pre-M83 tarball.

## Decisions

## Review
