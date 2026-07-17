# M66: Upstream re-oracle watch (weekly CRAN-current CI)

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** IP8, IP9
- **Branch/PR:** m66-upstream-reoracle-watch

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

- [ ] `re-oracle.yaml` exists with a weekly `schedule:` trigger and `workflow_dispatch`
      (with a `dry_fail` input); header comment states purpose, the CRAN-not-RSPM choice, and
      the 60-day auto-disable caveat.
- [ ] A normal run completes green on the milestone branch (via the temporary branch-push
      trigger — `workflow_dispatch` only registers once the file is on the default branch) and
      its log records the installed psych and lavaan versions.
- [ ] A `dry_fail` run opens a labeled issue containing the psych/lavaan versions; a second
      `dry_fail` run updates that same issue (no duplicate); the issue is closed afterward as
      test evidence; the temporary branch-push trigger is removed before review.
- [ ] The diff touches only `.github/` (no package surface; `git diff --stat` evidence), so
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
- [ ] T4: verify end-to-end on the milestone branch via a temporary `push` trigger scoped to
      the branch (dry_fail defaulted on branch pushes as needed): one green run, two `dry_fail`
      runs (issue create, then update), close the test issue, remove the temporary trigger,
      capture evidence for review; post-merge, one confirming `workflow_dispatch`.

## Work log

- 2026-07-17: created by /milestone-plan (promotes the 2026-07-17 candidate row; DESIGN §2
  track-CRAN-current posture, D-031).
- 2026-07-17: T1–T3 — authored `.github/workflows/re-oracle.yaml` (weekly cron + dispatch
  `dry_fail`, CRAN-not-RSPM, version-recording step, gh-CLI issue-on-failure). Notification
  tool = gh CLI (question gate: preinstalled on ubuntu-latest, no third-party action).

## Decisions

## Review
