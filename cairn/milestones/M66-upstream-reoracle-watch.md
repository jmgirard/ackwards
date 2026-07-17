# M66: Upstream re-oracle watch (weekly CRAN-current CI)

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** IP8, IP9
- **Branch/PR:** m66-upstream-reoracle-watch · [PR #70](https://github.com/jmgirard/ackwards/pull/70)

## Goal

A weekly scheduled CI job reruns the full check against fresh CRAN psych/lavaan and opens or
updates a GitHub issue on failure, mechanizing DESIGN §2's track-CRAN-current posture (lavaan
0.7 broke the suite unnoticed until the 0.1.1 release walk — PR #66).

## Scope

**In:** one new dedicated workflow (`.github/workflows/re-oracle.yaml`): weekly cron +
`workflow_dispatch`, single ubuntu-latest / R-release job, dependencies from plain CRAN (not
RSPM binaries, which can lag a fresh release), `R CMD check`, an issue-on-failure step (open
or update one labeled issue, never duplicates) recording the installed psych/lavaan versions,
and a `dry_fail` dispatch input to exercise the notification path end-to-end.

**Out:** multi-platform matrix (push/PR CI already runs the 5-config matrix; the watch hunts
upstream drift, not platform drift); R-devel (rhub.yaml + the release walk cover it); any
R-code change — this milestone touches `.github/` only. Auto-disable caveat (GitHub suspends
cron workflows after ~60 days of repo inactivity) is documented in the workflow header, not
worked around.

## Acceptance criteria

- [x] `re-oracle.yaml` exists with a weekly `schedule:` trigger and `workflow_dispatch`
      (with a `dry_fail` input); header comment states purpose, the CRAN-not-RSPM choice, and
      the 60-day auto-disable caveat.
- [x] A normal run completes green on the milestone branch (via the temporary branch-push
      trigger — `workflow_dispatch` only registers once the file is on the default branch) and
      its log records the installed psych and lavaan versions.
- [x] A `dry_fail` run opens a labeled issue containing the psych/lavaan versions; a second
      `dry_fail` run updates that same issue (no duplicate); the issue is closed afterward as
      test evidence; the temporary branch-push trigger is removed before review.
- [x] The diff touches only `.github/` (no package surface; `git diff --stat` evidence), so
      no R-level verify is required beyond the standard PR CI being green.

## Coverage

- AC1 → T1
- AC2 → T2, T3
- AC3 → T2, T4
- AC4 → T1, T2

## Tasks

- [x] T1: author `re-oracle.yaml` — cron (Mon 06:00 UTC) + `workflow_dispatch` with boolean
      `dry_fail`; ubuntu-latest, R release, `use-public-rspm: false` with the CRAN cloud
      mirror; `rcmdcheck` step (dry_fail forces a failing step instead); `permissions:`
      `issues: write` scoped to the job.
- [x] T2: issue-on-failure step (`actions/github-script` or `gh` CLI): find open issue by
      label `upstream-watch`; comment if found, create if not; body includes run link + psych/
      lavaan versions.
- [x] T3: version-recording step: print psych/lavaan `packageVersion()` to the log and expose
      as step outputs for T2.
- [x] T4: verify end-to-end on the milestone branch via a temporary `push` trigger scoped to
      the branch (dry_fail defaulted on branch pushes as needed): one green run, two `dry_fail`
      runs (issue create, then update), close the test issue, remove the temporary trigger,
      capture evidence for review; post-merge, one confirming `workflow_dispatch`.

## Work log

- 2026-07-17: created by /milestone-plan (promotes the 2026-07-17 candidate row; DESIGN §2
  track-CRAN-current posture, D-031).
- 2026-07-17: T1–T3 — authored `.github/workflows/re-oracle.yaml` (weekly cron + dispatch
  `dry_fail`, CRAN-not-RSPM, version-recording step, gh-CLI issue-on-failure). Notification
  tool = gh CLI (question gate: preinstalled on ubuntu-latest, no third-party action).
- 2026-07-17: T4 — end-to-end verified on branch via a temporary push trigger, then reverted.
  AC2 green run [29600465220] `success` (17m37s), log recorded psych=2.6.5 / lavaan=0.7.2.
  AC3: force-fail run [29601623972] opened issue #69 (label `upstream-watch`, body with run
  link + versions); re-trigger [29601746921] commented on the *same* #69 (1 open issue, 1
  comment — no duplicate). Closed #69 as test evidence; restored the production workflow (no
  TEMP markers remain). Post-merge confirming `workflow_dispatch` is owner follow-up.
- 2026-07-17: review — all 4 ACs verified fresh; cairn_validate clean; 3-lens fan-out + scorer.
  F1 (concurrency race, 72/100) fixed at owner's election at merge gate: added workflow-level
  `concurrency: { group: re-oracle, cancel-in-progress: false }`.

## Decisions

## Review

**Reviewed 2026-07-17 · PR #70 · branch `m66-upstream-reoracle-watch` vs `master`.**

### Acceptance-criteria evidence (fresh, by command)

- **AC1 (PASS).** `re-oracle.yaml` present; `schedule: cron '0 6 * * 1'` (weekly Mon 06:00
  UTC) + `workflow_dispatch` with boolean `dry_fail` input (`grep` confirmed). Header carries
  `PURPOSE.`, `CRAN, NOT RSPM.`, and `AUTO-DISABLE CAVEAT.` sections.
- **AC2 (PASS).** Re-queried run [29600465220]: `conclusion=success`, `event=push` on the
  branch; log records `psych=2.6.5`, `lavaan=0.7.2` (package passes against current CRAN
  lavaan 0.7.2).
- **AC3 (PASS).** Re-queried GitHub: run [29601623972] (`failure`) opened issue #69 (label
  `upstream-watch`, body references that run + versions); run [29601746921] (`failure`)
  commented on the *same* #69 — issue now CLOSED, **0** open `upstream-watch` issues, no
  duplicate. Body renders clean (no heredoc `EOF` artifact — verified against #69). Production
  workflow has no `TEMP-M66-TEST` / branch-push scaffolding.
- **AC4 (PASS).** `git diff --name-only master..HEAD` = `.github/workflows/re-oracle.yaml`,
  `cairn/ROADMAP.md`, `cairn/milestones/M66-*.md`; zero package-surface files
  (`R/`,`man/`,`tests/`,`vignettes/`,`data`,`DESCRIPTION`,`NAMESPACE`).

### Consistency gate

- `cairn_validate` → exit 0, all checks PASS (78 dangling-id WARNs are advisory pre-existing
  refs to entombed pre-migration milestone IDs; untouched by this diff).
- No DESIGN principle changed (DESIGN.md untouched) → `cairn_impact` skipped.
- R-toolchain `consistency-gate`: **vacuous** — `git diff` proves no package surface, so
  `document()` no-diff / generated-file / README / pkgdown-reference checks are trivially clean;
  a `.Rbuildignore`'d CI workflow is dev infra (not user-visible), so no NEWS entry required;
  `re-oracle.yaml` is under the existing `^\.github$` ignore. Authoritative `R CMD check` =
  PR #70's CI matrix (AC4's plan-sanctioned waiver of local R-level verify).

### Independent fresh-context review (3 lenses + scorer)

- **[O] diff-bug (Opus):** verified dry_fail `if:` complementarity (exactly one of force-fail
  / real-check runs on every event; uses typed `inputs` not stringy `github.event.inputs`),
  `failure()` semantics, heredoc column-0 rendering, no injection, least-privilege
  permissions, `tee`→step-output propagation. One finding (F1, below).
- **[S] blame-history (Sonnet):** change consistent with the intent/history it touches — every
  divergence from sibling workflows (`use-public-rspm: false`, scoped `permissions`, no
  paths-ignore, omitted vignette-freshness step) is plan-called-for, not a silent undo. Its
  heredoc "aside" was a false positive (tested the raw file, not the YAML-stripped script;
  refuted against #69's clean body).
- **[S] prior-PR-comments (Sonnet):** no prior-PR review-comment evidence in this repo
  (lessons live in cairn, not GitHub PRs); no regressions.

**Below-threshold findings — logged, not actioned (1):**
- **F1 (score 72/100, `concurrency`):** the find-then-create issue step has no `concurrency:`
  group; a slow scheduled run overlapping a manual `workflow_dispatch`, both failing, could
  each read `existing=""` and open two issues. Real but low-probability (weekly cron can't
  self-overlap; needs a coincident manual dispatch and two independent failures on a solo
  repo). Not an AC failure — AC3 covers the verified sequential dedup, not concurrent
  in-flight. Scorer 72 < 80 → excluded from actioned list; one-line fix
  (`concurrency: { group: re-oracle }`) available if the owner opts in.
  **Disposition: FIXED at the owner's election at the merge gate (2026-07-17)** —
  workflow-level `concurrency: { group: re-oracle, cancel-in-progress: false }` added
  (serializes runs so the issue read-then-write can't race a duplicate). Workflow-level
  key only; does not touch any T4-verified step path.
