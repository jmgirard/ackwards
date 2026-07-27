# M82: macOS oldrel (arm64) CI visibility + skip-filter repair

**Status:** done (2026-07-27, PR #88 https://github.com/jmgirard/ackwards/pull/88)

**Goal:** Make CRAN's `r-oldrel-macos-arm64` check flavour visible to CI before submission, and
repair the stale docs skip-filter that ran both check matrices on every tracking commit.

**Outcome:** New `.github/workflows/macos-oldrel-check.yaml` — `R CMD check` on `macos-latest` at
`oldrel-1`, push-to-master + `workflow_dispatch`, never per-PR (owner's call; signal arrives
post-merge). An assertion step fails the job unless `R.version$platform` matches
`^aarch64-apple-darwin`, so a `setup-r` x86_64 fallback cannot pass green while checking
`r-oldrel-macos-x86_64` — the flavour that already passed at 0.1.1. All five `paths-ignore` blocks
repaired (four dead entries dropped, six cairn tracking paths added), leaving `cairn/DECISIONS.md`,
`cairn/DESIGN.md` and `cairn/references/**` unfiltered because `check-ledger-anchors.R` reads them
from source. New base-R `tools/check-ci-path-filters.R` enforces that carve-out as a CI fail-fast
step, in `dod-gate.R`, and over six fixture trees.

**Decisions:** Push-to-master over per-PR (owner, plan gate). Standing R-hub step declined as
subsumed → candidate row, with macOS x86_64. Gated amendments: four filter blocks → five; Scope/T1.

**Review:** 3 lenses + scorer. Actioned F1 (93 — the parser ended a list at the first non-item
line, so a blanket `cairn/**` under a comment passed clean), F4 (82 — tests would have passed
against a stub), F6 (76, sub-threshold but verified false), F7 (self-found — fixtures tripped a
hidden-directory NOTE). Logged: F2 (74), F3 (58), F5 (42). AC3 was re-proved at a fresh run, not
read charitably, after a comment fix broke it as written.
