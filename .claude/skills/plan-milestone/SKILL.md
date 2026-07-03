---
name: plan-milestone
description: Scaffolds a milestone plan against CLAUDE.md and DESIGN.md, prompting for milestone-specific details before producing the plan.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

# Plan Milestone $ARGUMENTS

1. Read CLAUDE.md's "## Current focus" section to confirm the currently/previously recorded
   milestone number. If $ARGUMENTS isn't exactly one more than that (or doesn't otherwise make
   sense — e.g., re-planning the same milestone), stop and flag the discrepancy before proceeding.

2. Ask the user for milestone-specific details: what this milestone covers, why it comes next,
   and any constraints, deviations from DESIGN.md, or open questions they already have in mind.
   Wait for their answer — do not guess or proceed on assumptions.

Once they respond:

3. Re-read CLAUDE.md and DESIGN.md fresh. For prior-milestone detail, the single log is
   `MILESTONES.md` (CLAUDE.md carries only a one-line index); deferred/out-of-scope items are in
   DESIGN.md §14 and CLAUDE.md's "Out of scope for now". **If [`ROADMAP.md`](../../../ROADMAP.md)
   has a section covering the milestone being planned** (a scoped candidate, a deferred extension,
   or banked source notes), read it and treat it as the primary brief for steps 2/5 — it holds the
   driving rationale and banked decisions that the CLAUDE.md one-liner will later compress.
4. Confirm this milestone's scope against the deferred items in DESIGN.md §14 and CLAUDE.md's
   "out of scope for now" list (and that it isn't already done per the `MILESTONES.md` log); flag
   if the user's details conflict with either.
5. Propose a concrete plan: files to create/modify, the order of implementation, and testable
   acceptance criteria in the same style as Milestone 1's (see `MILESTONES.md`). If the milestone
   adds a new **exported** object (function, dataset, S3 generic), list `_pkgdown.yml`'s `reference:`
   list among the files to modify and make "new export is in the pkgdown reference index
   (`pkgdown::check_pkgdown()` passes)" an explicit acceptance criterion — a missing entry breaks the
   pkgdown GHA but not local `R CMD check`, so it is easy to miss.
6. Flag any design ambiguity the brief doesn't resolve — do not invent a resolution.
7. Do not write any code. This is planning only; implementation happens via
   /implement-milestone afterward.

8. Once the user explicitly approves the plan:
   a. Create the milestone's feature branch off an up-to-date `master`
      (e.g. `git switch master && git pull && git switch -c m$ARGUMENTS-<short-slug>`).
      **Before branching, confirm local `master` is *in sync* with `origin/master`, not merely
      "up to date".** Run `git fetch origin` then
      `git rev-list --left-right --count origin/master...master` — it must return `0	0`. If the
      right-hand count is > 0, local `master` is *ahead* of origin (unpushed commits — `git pull`
      does **not** flag this, since there is nothing to fetch): **stop and flag.** Push or reconcile
      those commits first. Branching from an ahead-of-origin `master` silently bundles those commits
      into the milestone's squash-merge and forces a post-merge divergence (this bit us once — see
      the recovery note in /implement-milestone's finishing steps).
      All milestone work — planning commit and implementation — lands on this branch, not on
      `master`. `master` is merged into only via a PR, squash-merged as soon as the local
      definition of done is green — not gated on remote CI (see /implement-milestone). Do **not**
      commit milestone work directly to `master`.
   b. On the branch, update CLAUDE.md's "## Current focus" section to replace the previous
      milestone's scope/acceptance criteria with Milestone $ARGUMENTS's, and commit that change
      on its own (not bundled with any implementation commit).