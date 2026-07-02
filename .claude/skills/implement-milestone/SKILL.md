---
name: implement-milestone
description: Implements the milestone plan already established in this conversation, following CLAUDE.md's dev workflow and definition of done. Use right after planning a milestone in this same session.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS matches the milestone number recorded in CLAUDE.md's "## Current focus"
   section. If it doesn't match, stop and flag the mismatch rather than proceeding.
0b. Confirm you are on the milestone's feature branch (created by /plan-milestone, named
    `m$ARGUMENTS-<slug>`), not on `master`. If you are on `master`, stop and flag — milestone
    work must not be committed directly to `master`.

# Implement Milestone $ARGUMENTS

Implement Milestone $ARGUMENTS using the plan already established earlier in this conversation.

**If no milestone plan is visible in the current context, stop and ask for it — do not invent or assume one.**

Follow `CLAUDE.md`'s dev workflow and definition of done throughout:
- Build in small, reviewable increments: one coherent change → test → commit → repeat. Do not produce one large commit at the end.
- After any roxygen change, run `devtools::document()` and commit the regenerated `NAMESPACE`/docs.
- **When you add or remove an exported object** (a function, dataset, or S3 generic — i.e. the
  regenerated `NAMESPACE` gained/lost an `export()`/`S3method()`, or you added a dataset under
  `data/`), update the `reference:` list in `_pkgdown.yml` to cover it in the same commit. This is
  easy to forget: it does **not** break local `R CMD check`, only the pkgdown GitHub Action ("Topics
  missing from index"). Guard against it explicitly — see the `pkgdown::check_pkgdown()` step in the
  gate below.
- Write or update tests for new behavior.
- **Verification cadence — don't re-run the suite needlessly:**
  - *While iterating / before each commit:* run only the **relevant** test file(s) —
    `devtools::test(filter = "<pattern>")` or `testthat::test_file("tests/testthat/test-<x>.R")` —
    not the whole suite. Branches are squash-merged, so a momentarily-red intermediate commit never
    ships; the full gate below catches everything before the PR.
  - *When you do run the whole suite*, prefix with `TESTTHAT_CPUS=8` — the suite is parallel
    (M48; ~27s with 8 workers vs ~81s serial; testthat defaults to only 2 workers without it).
  - *Run tests in a single pass that surfaces failures **with** their details* — capture the result
    (`res <- devtools::test(...)`; inspect `as.data.frame(res)`) or use a non-silent reporter. Never
    run silently "to see if it's green" and then re-run "to see what broke."
  - *`devtools::check()` already runs the entire test suite* (plus examples and vignette rebuild),
    and `covr::package_coverage()` runs it **again**. So do **not** stack `test()` → `check()` →
    `coverage()` at a gate — that executes the suite ~3×. Run `check()` once; it subsumes `test()`.
  - *Mid-implementation `check()` runs* (when you just want the R CMD check signal) can skip the
    ~74s vignette rebuild with `devtools::check(vignettes = FALSE)`. Do **one** full check
    (vignettes included) before the PR.
  - *Never run two package-touching R processes concurrently* (e.g. a background `check()` while a
    foreground `test()` runs) — they interfere and can produce spurious failures.
  - *Don't run a bare `devtools::load_all()` in its own `Rscript` call* — each `Rscript` is a fresh
    process, so nothing persists; `test()`/`check()` load the package themselves. And avoid
    `cd <repo> && …` compounds (permission-prompt noise) — use absolute paths.
  - *In tests, reuse the `cached()` fit memo* (`tests/testthat/helper-data.R`) instead of refitting
    identical `ackwards()` objects — but never for reproducibility / serial-vs-parallel oracles
    (a cached second call asserts nothing) or fits wrapped in condition expectations.
- **The definition-of-done gate (run once, before the PR):** `Rscript tools/dod-gate.R` — it runs
  `devtools::check()` (0/0/0) → `covr::package_coverage()` (target 100%) → `styler::style_pkg()` →
  `lintr::lint_package()` → `pkgdown::check_pkgdown()` serially in one process with sensible
  `TESTTHAT_CPUS`, printing every failure and exiting non-zero on any. This is the only place the
  full suite must run end-to-end; there is no need for a separate `devtools::test()` because
  `check()` already ran it. `pkgdown::check_pkgdown()` is sub-second and mirrors the pkgdown GHA —
  the script runs it unconditionally (gated on `rlang::is_installed("pkgdown")`; if pkgdown is
  absent, instead eyeball that every `export()` in `NAMESPACE` appears in `_pkgdown.yml`'s
  `reference:` contents). It catches exported topics missing from the reference index, which local
  `R CMD check` does not.
- Respect the ask-first guardrails in `CLAUDE.md`: flag before adding an `Imports` dependency, introducing Rcpp, changing a resolved default, or touching git history/tags.
- If you hit a design ambiguity not covered by `DESIGN.md`, stop and ask rather than deciding unilaterally.
- Update `NEWS.md` for user-visible changes.

**Milestone documentation (do this before the final commit — it is part of the definition of done):**
- **`MILESTONES.md` is the single source of truth for milestone history.** Add an `M$ARGUMENTS (done)`
  entry **in numeric order** (place it immediately after `M{$ARGUMENTS − 1}`; never append out of
  order or skip a number), in the same style as the existing entries (what shipped + the durable
  decisions/amendments + final test count and R CMD check status).
- Add the matching one-line index entry to CLAUDE.md's "## Completed milestones" list (numeric
  order), and update "## Current focus" to record the milestone as completed.
- Record any DESIGN.md contract changes in the relevant numbered section (§1–14) — **not** as a
  second copy of the milestone log (DESIGN.md §15 is only a pointer to `MILESTONES.md`).
- Do **not** duplicate the milestone narrative across files: detail in `MILESTONES.md`, one line in
  CLAUDE.md, user-facing notes in `NEWS.md`, design-contract deltas in DESIGN.md §1–14.

**Finishing the milestone (branch → PR, not direct to `master`):**
- Push the feature branch and open a PR into `master`, then **squash-merge it immediately** once
  the local definition of done is green (`devtools::check()` 0/0/0 + tests + coverage +
  style/lint): `gh pr merge <n> --squash --delete-branch`. Do **not** gate the merge on remote CI,
  and do **not** use `gh pr merge --auto` — `master` is deliberately unprotected, so auto-merge has
  no required check to wait on and silently no-ops. CI still runs on the PR/push as an
  after-the-fact signal (glance at it later; fix-forward if a platform job goes red), but it does
  not block the merge. Don't spin up a synchronous CI watch or a background CI poller "to wait for
  green" — that reintroduces exactly the wait this workflow rejects. After merging, sync local:
  `git switch master && git pull --ff-only`. **If that fast-forward fails**, local `master` held
  commits that were not on the PR's base (typically unpushed direct-to-`master` edits, which get
  bundled into the squash on the remote side). Do **not** force-push or `reset --hard`: stop,
  inspect the divergence (`git log --oneline --left-right origin/master...master`), and reconcile
  deliberately — merging `origin/master` into local `master` (a merge commit, then push) is the
  non-destructive fix and preserves both sides. The clean prevention is /plan-milestone step 8a's
  in-sync check; when that passes, this fast-forward always succeeds. Do not push milestone commits
  straight to `master`.
  - **Exception — release / CRAN-submission milestones** (CRAN-prep, version bump, release): for
    these, **do** wait for the *full* green CI matrix before merging (CRAN runs exactly that
    matrix and rejects platform failures the local macOS `check()` can't see). Only here is a
    synchronous/polled CI wait warranted; flag it explicitly when you take it.
- Trivial, isolated doc-typo fixes may still go directly to `master` at the user's discretion;
  anything touching `R/`, `tests/`, `DESCRIPTION`, or vignettes goes through a PR.

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step you trigger after reviewing the work.