# M82: macOS oldrel (arm64) CI visibility + skip-filter repair

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP6
- **Branch/PR:** `m82-macos-oldrel-ci-visibility` · [#88](https://github.com/jmgirard/ackwards/pull/88)

## Goal

Make CRAN's `r-oldrel-macos-arm64` check flavour visible to CI before submission, and repair the
stale docs skip-filter that now runs both check matrices on every tracking commit.

## Scope

**In:** a new push-to-master workflow running `R CMD check` on `macos-latest` at `oldrel-1`, with a
step that fails on an x86_64 interpreter; repair of all four `paths-ignore` blocks (two in
`R-CMD-check.yaml`, two in `test-coverage.yaml`) — dead entries removed, cairn tracking paths
added, ledger-read paths deliberately left unmatched.

**Out:**
- macOS x86_64 (`macos-13`) coverage → candidate row; the 0.1.1 failure was arm64-specific and
  CRAN's `r-oldrel-macos-x86_64` passed (commit `5d2fcbd`).
- Standing R-hub macOS pre-submission step in `PROFILE.md`'s release-walk → candidate row; the new
  job covers the same ground without the maintainer having to remember it.
- Running the new job per-PR → declined by the owner at the 2026-07-27 plan gate; push-to-master
  only, accepting that the signal arrives after the merge.
- `cran-comments.md` platform-table refresh → the next release walk (owner-timed).
- `pkgdown.yaml` → untouched; it carries no `paths-ignore` filter by design.
- Any package or test change → none; fork safety was fixed at `5d2fcbd`.

## Acceptance criteria

- [x] AC1 `.github/workflows/macos-oldrel-check.yaml` exists; its `on:` block is exactly
      `push: branches: [master]` plus `workflow_dispatch`, with no `pull_request`, and carries the
      same repaired `paths-ignore` filter AC4 specifies; its single
      `macos-latest` job runs `setup-pandoc@v2`, `setup-r@v2` at `r: 'oldrel-1'`,
      `setup-r-dependencies@v2` (`needs: check`, `extra-packages: any::rcmdcheck`), and
      `check-r-package@v2`.
- [x] AC2 A step placed **before** the check step evaluates
      `grepl("^aarch64-apple-darwin", R.version$platform)` and exits non-zero when false, so an
      x86_64 interpreter on the arm64 runner reddens the job instead of passing green while
      covering the flavour that already passed on CRAN.
- [x] AC3 The job ran green on the milestone branch with that step reporting an
      `aarch64-apple-darwin*` platform, **and** the file that lands on master differs from the
      commit that produced that run only inside the `push: branches:` list. Evidence: the numeric
      run ID (durable after the trigger revert), the `gh run view` conclusion line, the platform
      line quoted from the log, and `git diff <tested-commit>..HEAD -- <the workflow>`.
- [x] AC4 In all five `paths-ignore` blocks (`R-CMD-check.yaml:10`, `:20`;
      `test-coverage.yaml:7`, `:17`; and the new `macos-oldrel-check.yaml`'s
      own): no listed literal path is absent from the repo; each block
      ignores `cairn/ROADMAP.md`, `cairn/LESSONS.md`, `cairn/PROFILE.md`, `cairn/milestones/**`,
      `cairn/legacy/**`, `cairn/reviews/**`; and no pattern in any block matches
      `cairn/DECISIONS.md`, `cairn/DESIGN.md`, or any `cairn/references/**` path, which
      `tools/check-ledger-anchors.R` reads. Evidence: a recorded shell check of all three clauses.
- [x] AC5 The justifying comment above each repaired block states the ledger-read carve-out and why
      blanket `cairn/**` would silently disable that checker.
- [x] AC6 `Rscript tools/dod-gate.R` clean — confirming in particular that the ledger-anchor and
      vignette-freshness guards still pass and no package surface moved.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3, T4
- AC4 → T5, T7
- AC5 → T6
- AC6 → T7

## Tasks

- [x] T1 Author `.github/workflows/macos-oldrel-check.yaml` per AC1. Header comment names the 0.1.1
      `r-oldrel-macos-arm64` ERROR and its fix (`5d2fcbd`), why the flavour was invisible (both
      matrices run macOS at `release` only), and that push-only rather than per-PR was the owner's
      call at the 2026-07-27 plan gate.
- [x] T2 Add the architecture-assertion step ahead of the check step, so a fallback fails fast
      rather than after a full check.
- [x] T3 Temporarily add the milestone branch to the workflow's `push: branches:`; push; confirm the
      run green and the platform line; record run ID, conclusion line, and platform line.
- [x] T4 Revert the temporary trigger (M66 lesson); capture the AC3 diff proving nothing else in the
      file moved since the tested commit.
- [x] T5 Repair the four existing `paths-ignore` blocks: drop `DESIGN.md`, `MILESTONES.md`,
      `ROADMAP.md`, `nested-cv-guide.md`; add the six cairn tracking paths; leave the three
      ledger-read paths unmatched. The new workflow's own block (T1) uses the same list.
- [x] T6 Rewrite the justifying comment above each block per AC5 (currently
      `R-CMD-check.yaml:5-9`).
- [x] T8 (discovered) `tools/check-ci-path-filters.R` — base-R guard asserting AC4's clauses plus
      block existence, wired into `tools/dod-gate.R`, an `R-CMD-check.yaml` fail-fast step, and
      `tests/testthat/test-ci-path-filters.R`. Makes AC4's carve-out mechanical rather than a
      comment nobody re-reads.
- [x] T7 Verify AC4's three clauses mechanically, then run `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-27: T7 + status review. DoD gate exit 0: check 0 err/0 warn/0 note, coverage 100.00%, styler/lintr clean, pkgdown index complete, all three source guards clean. No NEWS.md entry — nothing user-visible ships (.github/, tools/, tests/ only).
- 2026-07-27: T3+T4. Run 30245095260 (commit bb69970) succeeded on macos-latest at R 4.5.3; assertion printed `platform: aarch64-apple-darwin20`; R CMD check reported `Status: OK`. Temporary trigger reverted; `git diff bb69970` over the workflow shows 0 changed lines outside the `branches:` line.
- 2026-07-27: T5+T6+T8. Discovered sub-task T8 (minor amendment): AC5's carve-out was a comment only, so added tools/check-ci-path-filters.R to enforce it. Written in base R after a first yaml-package draft — it runs as a CI fail-fast step before setup-r-dependencies, where only base R exists. Inversion-verified: dead literal, dropped cairn entry, blanket cairn/**, and a deleted block each exit 1 with the right message; clean exits 0.
- 2026-07-27: T1+T2 done in one commit (minor merge — creating the file without the assertion then patching it is make-work). Gated amendment: the new workflow carries the same repaired paths-ignore filter, so AC1 and AC4 now cover five blocks, not four. Assertion verified both ways locally: passes on aarch64-apple-darwin25.4.0, fails on a simulated x86_64-apple-darwin20.
- 2026-07-27: created by /milestone-plan; promotes the 2026-07-27 CI-flavour-visibility candidate row. Owner chose push-to-master over per-PR at the plan gate. Two fresh-context criteria audits ran; the second found `test-coverage.yaml` carries the same stale filter block twice, widening scope from two blocks to four.

## Decisions

## Review

Reviewed 2026-07-27 on branch `m82-macos-oldrel-ci-visibility`, PR #88. master had not moved
since the branch was cut (0 commits ahead/behind), so no re-merge was needed.

### Evidence per criterion

- **AC1** — `on:` block read from the file: `push: branches: [master]` + `workflow_dispatch`, no
  `pull_request`, with the nine-entry repaired filter. Steps in order: `actions/checkout@v4`,
  `setup-pandoc@v2`, `setup-r@v2` (`r-version: 'oldrel-1'`), the assertion, `setup-r-dependencies@v2`
  (`needs: check`, `extra-packages: any::rcmdcheck`), `check-r-package@v2`. `runs-on: macos-latest`.
- **AC2** — assertion at line 70, check step at line 80, so it precedes the check *and* the
  dependency install. Fresh local execution of the same expression: passes on
  `aarch64-apple-darwin25.4.0`, fires on a simulated `x86_64-apple-darwin20`.
- **AC3** — re-proved after the review fixes changed the header comment, rather than reading the
  criterion charitably: run [30247411925](https://github.com/jmgirard/ackwards/actions/runs/30247411925),
  `conclusion=success`, sha `3059906`; log shows R 4.5.3, `platform: aarch64-apple-darwin20`,
  `Status: OK`. `git diff -U0 3059906` over the workflow: **0** changed lines outside `branches:`.
  (Superseded runs: 30245095260 at `bb69970`, `Status: OK`, predating the comment fix;
  30247047807 at `0a2829f`, `Status: 1 NOTE` — see F7.)
- **AC4** — `Rscript tools/check-ci-path-filters.R` exit 0 across 5 blocks. Guard now pinned by
  committed fixtures rather than hand-run inversions: against the pre-fix checker (`f141594`) the
  new tests produce 8 clean failures; against the fixed checker 19 assertions pass.
- **AC5** — both rewritten comments read from the files; each names the ledger-read carve-out.
  R-CMD-check.yaml's retains the original pkgdown-not-filtered note verbatim (confirmed by the
  blame-history lens against commit `9030125`).
- **AC6** — `Rscript tools/dod-gate.R` exit 0: check 0 errors / 0 warnings / 0 notes [142s],
  coverage 100.00% [77s], styler clean, lintr clean, pkgdown reference index complete; all three
  source-checkout guards clean.

### Consistency gate

`cairn_validate` exit 0 — `weight caps`, `roadmap<->disk orphans`, `scaffold present`,
`coverage complete`, `binding criteria` all PASS; the 91 `dangling id tokens` WARNs are
pre-existing legacy M-ids. No `DESIGN.md` principle changed, so `cairn_impact` was skipped.
Profile `consistency-gate` slot: check/coverage/style/lint/pkgdown via the DoD gate above; no NEWS
entry owed — nothing user-visible ships. One `.Rbuildignore` entry was owed and added (F7).

### Recorded inconsistencies in plan-owned text (not edited review-side)

Two stale sentences survive in plan-owned sections. Neither affects what shipped, and review does
not edit plan-owned wording, so both are recorded rather than fixed:

1. **Scope "In:"** still says "all four `paths-ignore` blocks (two in `R-CMD-check.yaml`, two in
   `test-coverage.yaml`)". The gated amendment at implement made it five (the new workflow carries
   its own); AC1 and AC4 were updated, Scope was not.
2. **T1** says the flavour was invisible because "both matrices run macOS at `release` only" — the
   same false claim F6 corrected in the workflow header. `test-coverage.yaml` is ubuntu-only.

### Independent review — three lenses + scorer

Diff-bug [O]: 5 findings. Blame-history [S]: none — it traced every modified line to its
originating commit (`9030125` for the filters, `a187c07` for the cairn/ migration) and found the
change a faithful continuation. Prior-PR [S]: 1 finding; the existence probe
(`gh api .../pulls/comments`) returned `[]`, so the GitHub thread walk was skipped as designed and
the archive was the evidence base.

**Actioned:**

- **F1 (93) — fixed.** The block parser ended a list at the first non-item line, so a comment
  inside it truncated the block and entries below went unexamined; a blanket `cairn/**` written
  under a comment passed clean — the exact mutation clause 3 exists to catch. Reproduced before
  fixing. A block now ends only on a dedent; blanks and comments are skipped; an unreadable line
  inside a block is reported. Regression-tested by the `comment-hidden` fixture.
- **F4 (82) — fixed.** The test file's only checker assertion (`expect_equal(problems,
  character(0))`) would have passed against a stub, leaving F1's territory untested. Six fixture
  workflow trees now pin each failure path.
- **F6 (76) — fixed by exception.** Below the 80 threshold, actioned anyway: the header claimed
  `test-coverage.yaml` "runs macOS at `release`", but it is `runs-on: ubuntu-latest` with no matrix
  and runs macOS at no version. The score reflects doc-only impact, not doubt about the fact; a
  verified-false claim is fixed regardless (M63/M67 lesson).
- **F7 — self-found during the AC3 re-run, fixed.** The new fixtures ship a literal `.github`
  directory, tripping `checking for hidden files and directories`; `Status:` went `OK` → `1 NOTE`.
  Added `^tests/testthat/fixtures/ci-filters$` to `.Rbuildignore`; tarball verified to exclude
  them, and the tests already skip in the built package. `Status: OK` restored.

**Logged, below threshold, not actioned:**

- **F2 (74)** — flow-style `paths-ignore:` lists recorded no block, and two blocks resolving to one
  event name silently overwrote. Incidentally covered by F1's rewrite, which now reports both.
- **F3 (58)** — clause 2 is applied to every block found, not only the five in `REQUIRED_BLOCKS`,
  so a future unrelated workflow with a narrow filter would fail the gate. Inert today (the only
  blocks that exist are the five); stricter than AC4 requires, not in conflict with it.
- **F5 (42)** — clause 1's literal inputs are themselves in the ignore list, so a commit deleting
  one skips CI and the staleness surfaces later on an unrelated commit. Inherent to any
  filter-scoped guard, and mitigated by `dod-gate.R` running the checker unconditionally pre-push.
