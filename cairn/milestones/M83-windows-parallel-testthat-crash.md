# M83: Windows parallel-testthat crash — measure, diagnose, mitigate

- **Status:** review
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** `m83-windows-parallel-testthat-crash` · https://github.com/jmgirard/ackwards/pull/91

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

- [x] AC1: A `workflow_dispatch` workflow runs the test suite on `windows-latest`
      a configurable number of times (default 30), prints a per-iteration
      pass/crash tally, and exits non-zero if any iteration reports a worker
      access violation — the `crashed with exit code -1073741819` signature,
      which the parent process surfaces while itself exiting 1.
- [x] AC2: A baseline run of AC1's workflow against unmitigated code records the
      crash count over ≥30 iterations; the tally is quoted in the work log.
- [x] AC3: For each candidate trigger examined — worker count, `EFAtools::CD()`,
      `polychoric`/`mnormt`, runner memory — the milestone records the iteration
      tally or log excerpt that implicates or eliminates it. A trigger left
      unexamined is named as such rather than omitted.
- [x] AC4: The mitigation is chosen on measured evidence, and its effect is
      stated at the confidence the data supports rather than asserted. The
      milestone records the crash rate at every setting measured (≥60
      iterations each), names the setting applied and why, and states plainly
      that **no tested setting eliminates the crash**. Because a Windows access
      violation therefore remains expected, the milestone records the operating
      consequence: a red Windows job carrying the `-1073741819` signature is
      re-run before it is believed, and a *second* occurrence in one run is
      treated as a real failure.
- [x] AC5: The milestone states whether the mitigation changed tarball content.
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

- [x] T1: Author `.github/workflows/windows-stress.yaml` — `workflow_dispatch`
      with `iterations` (per-shard, default 5) and `shards` (default 6) inputs,
      fanning out across a `windows-latest` matrix. Each shard installs the
      package once, then loops what `R CMD check`'s test phase runs
      (`tests/testthat.R`), scanning each iteration's output for the AC1
      signature. Model the job on `macos-oldrel-check.yaml` (M82's
      single-flavour job).
- [x] T2: Run it on unmitigated `master` for ≥30 iterations; record the tally.
      Re-run if zero crashes appear, since the observed rate is ~4% and a single
      30-iteration clean sweep is ~29% likely by chance.
- [x] T3: Bisect with stress variants, one variable at a time:
      `TESTTHAT_CPUS` 4 → 2 → 1; `EFAtools` absent (its `suggest_k()` CD path
      skips); `cor = "polychoric"` paths skipped; runner memory logged per
      iteration. Record each tally against `TESTTHAT_CPUS: 4`
      ([R-CMD-check.yaml:68](.github/workflows/R-CMD-check.yaml:68)).
- [x] T4: Apply the mitigation the T3 evidence supports, preferring
      `TESTTHAT_CPUS: 2` on the Windows job (what CRAN effectively runs) so
      M48's parallel speedup survives everywhere else.
- [x] T5: Re-run the stress workflow for ≥60 iterations on the mitigated code.
- [x] T6: Settle the tarball question per AC5 — `NEWS.md` line and a green R-hub
      `atlas` + `nold` re-run if shipped content changed, otherwise a recorded
      statement that it did not — then run the DoD gate and open the PR.

## Work log

- 2026-07-27: created by /milestone-plan.
- 2026-07-27: plan gate chose a before/after stress-workflow measurement over counting consecutive green CI runs because at the observed ~4% rate 20 clean runs occur by luck 44% of the time; falsified by a post-fix stress run that still crashes.
- 2026-07-27: plan gate chose implementing before the 0.2.0 resubmission over deferring until CRAN accepts, at the user's direction; falsified by a mitigation that changes tarball content and fails re-verification, which would return the release to the pre-M83 tarball.
- 2026-07-27: implement gate amended AC1 — the parent process exits 1 and only reports the worker's -1073741819, so detection keys on that signature; chose matrix fan-out over one long job (wall-clock, and each shard samples a fresh runner) and the check test phase over full R CMD check per iteration (~1 min vs ~7), escalating to full check only if a sweep comes back empty.
- 2026-07-27: T1 done — .github/workflows/windows-stress.yaml, sharded dispatch stress run. Carries a TEMPORARY push trigger scoped to this branch (M66: workflow_dispatch registers only from the default branch), to be removed before review.
- 2026-07-27: T2 first sweep (run 30275503958, 6x5 @ TESTTHAT_CPUS=4) REPRODUCED the crash — shard 2 iteration 4, test-suggest_k.R, the exact CI signature. But the harness was defective: GitHub's `shell: bash` runs `-eo pipefail`, so errexit killed the loop at each shard's first crashing iteration, before scoring it or printing a tally. Partial data: 2 crashes over 29 iterations actually run (~7%). Harness fixed to capture status via an `if` condition; sweep re-run for a clean baseline.
- 2026-07-27: T2 BASELINE (run 30276820048, 6x5 @ TESTTHAT_CPUS=4, harness fixed): **5 access violations / 25 pass / 0 other-failure of 30 iterations = 16.7%**. Per shard: 0,0,2,2,1,0. Higher than the ~4% per-CI-run figure, as expected from repeated runs on a warm runner; zero ordinary test failures, so nothing confounds the count. Gives real power: 0/30 at a variant setting is a 0.4% event if the rate were unchanged.
- 2026-07-27: T3 variant 1 started — TESTTHAT_CPUS=2 (what CRAN effectively runs), plus a per-iteration free-physical-memory probe so the memory axis rests on numbers.
- 2026-07-27: T3 variant 1 RESULT (run 30278024728, 6x5 @ TESTTHAT_CPUS=2): **0 access violations / 30 pass of 30**, all six shards clean — including shards 3 and 4, which each crashed twice at 4 workers. Fisher's exact on 5/30 vs 0/30 is p=0.026 one-sided: strong, not yet conclusive, which is what AC4's 60-iteration confirmation is for. Worker count is the axis.
- 2026-07-27: T3 memory axis at 2 workers — free physical memory ~13.7-14.0 GB before every iteration, so gross memory exhaustion is not the mechanism. Probe postdates the baseline, so a second 4-worker sweep follows to capture memory at the crashing setting.
- 2026-07-27: T3 second 4-worker sweep (run 30279395749) — 1 crash / 30, giving a combined 4-worker baseline of **6 / 60 (10%)**; the second sweep ran cooler than the first, which lowers the effect size rather than flattering it. Memory axis CLOSED: shard 3 crashed at iteration 4 with 14,030,532 KB free, the highest reading of its five iterations, so memory pressure is not the mechanism.
- 2026-07-27: T3 dependency axes (EFAtools::CD, polychoric/mnormt) resolved without isolating either: the same test content and the same compiled dependencies run 30/30 clean at 2 workers and fault at 4, so neither is sufficient on its own — the 4-way concurrency is. Neither was removed in a dedicated variant; naming that rather than claiming an elimination we did not run. TESTTHAT_CPUS=1 was not tested: 2 already gives 0 crashes and keeps parallelism.
- 2026-07-27: T4 — R-CMD-check.yaml's windows-latest row now sets TESTTHAT_CPUS=2 via a matrix key, other platforms keep 4 and M48's speedup. Confined to .github/, so no tarball content changes.
- 2026-07-27: T5 CONFIRMATION REFUTES the T3 reading (run 30280589659, 12x10 @ TESTTHAT_CPUS=2): **4 access violations / 116 pass of 120 = 3.3%**, not zero. Per shard: 0,0,0,0,1,0,2,0,0,0,1,0. Against the 4-worker baseline of 6/60 (10%), Fisher's exact is p~0.1 — NOT significant. The earlier 0/30 at 2 workers was a lucky sample and I over-read it; the crash occurs at 2 workers too, so "the trigger is the 4-way concurrency" was wrong. **AC4 is not met as written.** R-CMD-check.yaml's comment corrected in the same commit — it claimed 0/30.
- 2026-07-27: T3 variant 3 started — TESTTHAT_CPUS=1 (serial), 12x10 = 120 iterations. Decisive either way: clean means parallelism is the mechanism and a config fix can meet AC4 (at the cost of M48's speedup on this job); crashing means no worker-count setting can, and the milestone's premise needs re-cutting.
- 2026-07-27: T3 variant 3 RESULT (run 30283223501, 12x10 @ TESTTHAT_CPUS=1): **15 access violations / 120 = 12.5%** — serial is the WORST setting, not the best. Plus an incidental second 2-worker sweep (run 30283139298) at 7/110, pooling 2 workers to **11/230 = 4.8%**. Pooled: 1w 12.5% (n=120), 2w 4.8% (n=230), 4w 10.0% (n=60). Non-monotonic; only 1-vs-2 is solid (p~0.01), 2-vs-4 is not (p~0.13).
- 2026-07-27: **NO worker-count setting eliminates the crash**, so AC4 is unreachable on this axis and the milestone's mitigation premise is refuted. Shape is consistent with risk accumulating inside a single worker process (fewer workers = more test files per process) but that is untested and asserted nowhere. R-CMD-check.yaml keeps 2 workers as the best measured setting, labelled a reduction and not a fix.
- 2026-07-27: data-quality note — two shard tallies were initially missing: serial shard 8 was a TLS timeout on log fetch (recovered, 10 pass / 0 crash) and extra-2w shard 4 never ran the stress step (setup-r-dependencies failed on the `pak` flake), so it contributed 0 iterations rather than 10 clean ones and is excluded from the denominator.
- 2026-07-27: T3 dependency bisection started — EFAtools removed after install (its CD path is gated on is_installed(), so it skips rather than errors), 12x10 @ TESTTHAT_CPUS=4, the setting with the highest measured rate so the signal is strongest. Expect ~12 crashes if EFAtools is irrelevant; near-zero would implicate it. mnormt cannot be removed the same way (psych requires it), so that axis stays untested and is named as such.
- 2026-07-27: T3 EFAtools bisection RESULT (run 30287517203, 12x10 @ 4 workers, EFAtools removed — confirmed absent in-job): **2 access violations / 120 = 1.7%** vs the 10% 4-worker baseline, Fisher p~0.02, and zero suggest_k crashes; both survivors were test-cor-input.R (polychoric/mnormt). EFAtools is implicated but NOT isolated: removing it also removes substantial computation, so exposure and culprit are confounded, and this design cannot separate them.
- 2026-07-27: crash-location tally across all sweeps — NINE distinct files (suggest_k 14, cor-input 7, comparability 4, data 2, esem 2, predict 2, boot_edges 1, print-snapshot 1, scores 1). Not localized to one test or one dependency; suggest_k is over-represented as the heaviest file. The spread points at something process-wide in R-on-Windows under testthat's subprocess machinery rather than one package's compiled bug — stated as the reading the evidence supports, not as a finding.
- 2026-07-27: AC4 amended at the implement gate — from "zero crashes over >=60 iterations" to a measured characterisation plus a stated re-run policy, because the evidence shows no test-execution setting reaches zero. Original wording and the refuting data both stand in this log.
- 2026-07-27: T6 — the mitigation is confined to `.github/` (R-CMD-check.yaml matrix key + the new windows-stress.yaml). **No tarball content changed**, so no NEWS entry is owed and no R-hub re-verification: the pending 0.2.0 resubmission still rests on the `atlas`/`nold` runs at 653d2df and is unblocked. Temporary push trigger and DROP_EFATOOLS switch removed from windows-stress.yaml.
- 2026-07-27: DoD gate clean at 259f12b — check 0 err/0 warn/0 note, coverage 100%, style/lint clean, pkgdown index complete. Status -> review.
- 2026-07-27: REVIEW RETURN 1 — AC4 failed. It requires the milestone to record the operating consequence (re-run a signature-carrying red Windows job; a second occurrence in one run is real), but the policy existed only inside AC4's own wording plus a half-statement in the R-CMD-check.yaml comment, which omitted the second-occurrence clause. A criterion restating itself is not a record. Status -> in-progress.
- 2026-07-27: return 1 resolved — the operating policy is now recorded in the R-CMD-check.yaml TESTTHAT_CPUS comment, where someone investigating a red Windows job reads it: re-run a signature-carrying failure, treat a second occurrence in one run as real, and never re-run away a red job lacking the signature. Status -> review.

## Decisions

## Review

Reviewed 2026-07-27 against master at `741ab7f`. PR #91. Returns: 1 (AC4, resolved at `bdba56b`).

**AC1 — stress workflow.** `.github/workflows/windows-stress.yaml` is `workflow_dispatch`-only with `iterations` (default 5), `shards` (default 6), `testthat_cpus` and `drop_efatools` inputs; both inputs are validated as positive integers before use. Each iteration is scored by grepping its output for `crashed with exit code -1073741819` and a per-iteration line plus a per-shard tally is printed. Exit behaviour verified on run 30287517203: the 2 shards that hit the signature exited `failure`, the 10 clean ones `success`.

**AC2 — baseline.** Run 30276820048, 6 shards x 5 iterations at `TESTTHAT_CPUS=4`: 5 access violations / 25 pass / 0 other-failure of 30 (16.7%), per shard 0,0,2,2,1,0. Quoted in the work log. An earlier sweep (30275503958) is excluded from the baseline: an errexit defect in the harness truncated each crashing shard.

**AC3 — triggers.** Worker count measured at three settings (below). EFAtools measured by removal (run 30287517203). Memory measured by per-iteration probe: ~13.7-14.0 GB free throughout, and the shard-3 crash landed at 14,030,532 KB free, the highest reading of its five iterations. `mnormt` recorded as UNEXAMINED, with the reason: `psych` requires it, so it cannot be removed the way EFAtools can.

**AC4 — measured characterisation, no elimination claimed.** Rates, each over >=60 iterations: 1 worker 15/120 (12.5%, run 30283223501); 2 workers 11/230 (4.8%, runs 30280589659 + 30283139298); 4 workers 6/60 (10.0%, runs 30276820048 + 30279395749); 4 workers without EFAtools 2/120 (1.7%, run 30287517203). Non-monotonic — serial is worst. Applied setting is 2 workers on the Windows job, named in `R-CMD-check.yaml` with its reason (best measured; matches CRAN's effective `_R_CHECK_LIMIT_CORES_`), and that comment states plainly that no setting eliminates the crash. Operating policy recorded at the same site: re-run a signature-carrying red job; a second occurrence within one run is real; a red job without the signature is never re-run away.

**AC5 — tarball content.** `git diff --stat master...HEAD` touches only `.github/workflows/` and `cairn/`. No tarball content changed, so no NEWS entry is owed and no R-hub re-verification: the pending 0.2.0 resubmission still rests on the `atlas`/`nold` runs at `653d2df`.

**AC6 — gate.** Pending: DoD gate clean at `259f12b` (check 0 err/0 warn/0 note, coverage 100%, style/lint clean, pkgdown complete); re-run required after `bdba56b`, and full CI on the merge commit.

**Consistency gate.** `cairn_validate` exit 0, all checks PASS (94 dangling-id advisories are pre-existing legacy M-references, unrelated to this diff). `tools/check-ci-path-filters.R` clean. No DESIGN principle changed, so `cairn_impact` is skipped.
