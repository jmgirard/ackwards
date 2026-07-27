# M82: macOS oldrel (arm64) CI visibility + skip-filter repair

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP6
- **Branch/PR:** `m82-macos-oldrel-ci-visibility`

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

- [ ] AC1 `.github/workflows/macos-oldrel-check.yaml` exists; its `on:` block is exactly
      `push: branches: [master]` plus `workflow_dispatch`, with no `pull_request`, and carries the
      same repaired `paths-ignore` filter AC4 specifies; its single
      `macos-latest` job runs `setup-pandoc@v2`, `setup-r@v2` at `r: 'oldrel-1'`,
      `setup-r-dependencies@v2` (`needs: check`, `extra-packages: any::rcmdcheck`), and
      `check-r-package@v2`.
- [ ] AC2 A step placed **before** the check step evaluates
      `grepl("^aarch64-apple-darwin", R.version$platform)` and exits non-zero when false, so an
      x86_64 interpreter on the arm64 runner reddens the job instead of passing green while
      covering the flavour that already passed on CRAN.
- [ ] AC3 The job ran green on the milestone branch with that step reporting an
      `aarch64-apple-darwin*` platform, **and** the file that lands on master differs from the
      commit that produced that run only inside the `push: branches:` list. Evidence: the numeric
      run ID (durable after the trigger revert), the `gh run view` conclusion line, the platform
      line quoted from the log, and `git diff <tested-commit>..HEAD -- <the workflow>`.
- [ ] AC4 In all five `paths-ignore` blocks (`R-CMD-check.yaml:10`, `:20`;
      `test-coverage.yaml:7`, `:17`; and the new `macos-oldrel-check.yaml`'s
      own): no listed literal path is absent from the repo; each block
      ignores `cairn/ROADMAP.md`, `cairn/LESSONS.md`, `cairn/PROFILE.md`, `cairn/milestones/**`,
      `cairn/legacy/**`, `cairn/reviews/**`; and no pattern in any block matches
      `cairn/DECISIONS.md`, `cairn/DESIGN.md`, or any `cairn/references/**` path, which
      `tools/check-ledger-anchors.R` reads. Evidence: a recorded shell check of all three clauses.
- [ ] AC5 The justifying comment above each repaired block states the ledger-read carve-out and why
      blanket `cairn/**` would silently disable that checker.
- [ ] AC6 `Rscript tools/dod-gate.R` clean — confirming in particular that the ledger-anchor and
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
- [ ] T3 Temporarily add the milestone branch to the workflow's `push: branches:`; push; confirm the
      run green and the platform line; record run ID, conclusion line, and platform line.
- [ ] T4 Revert the temporary trigger (M66 lesson); capture the AC3 diff proving nothing else in the
      file moved since the tested commit.
- [ ] T5 Repair the four existing `paths-ignore` blocks: drop `DESIGN.md`, `MILESTONES.md`,
      `ROADMAP.md`, `nested-cv-guide.md`; add the six cairn tracking paths; leave the three
      ledger-read paths unmatched. The new workflow's own block (T1) uses the same list.
- [ ] T6 Rewrite the justifying comment above each block per AC5 (currently
      `R-CMD-check.yaml:5-9`).
- [ ] T7 Verify AC4's three clauses mechanically, then run `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-27: T1+T2 done in one commit (minor merge — creating the file without the assertion then patching it is make-work). Gated amendment: the new workflow carries the same repaired paths-ignore filter, so AC1 and AC4 now cover five blocks, not four. Assertion verified both ways locally: passes on aarch64-apple-darwin25.4.0, fails on a simulated x86_64-apple-darwin20.
- 2026-07-27: created by /milestone-plan; promotes the 2026-07-27 CI-flavour-visibility candidate row. Owner chose push-to-master over per-PR at the plan gate. Two fresh-context criteria audits ran; the second found `test-coverage.yaml` carries the same stale filter block twice, widening scope from two blocks to four.

## Decisions

## Review
