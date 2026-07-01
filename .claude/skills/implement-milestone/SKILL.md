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
- Write or update tests for new behavior; run `devtools::test()` before each commit.
- Run `devtools::check()` before considering any sub-step done; resolve errors/warnings, triage notes.
- Style and lint (`styler::style_pkg()`, `lintr::lint_package()`).
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
  `git switch master && git pull`. Do not push milestone commits straight to `master`.
  - **Exception — release / CRAN-submission milestones** (CRAN-prep, version bump, release): for
    these, **do** wait for the *full* green CI matrix before merging (CRAN runs exactly that
    matrix and rejects platform failures the local macOS `check()` can't see). Only here is a
    synchronous/polled CI wait warranted; flag it explicitly when you take it.
- Trivial, isolated doc-typo fixes may still go directly to `master` at the user's discretion;
  anything touching `R/`, `tests/`, `DESCRIPTION`, or vignettes goes through a PR.

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step you trigger after reviewing the work.