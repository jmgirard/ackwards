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

**CLAUDE.md is the single source of truth for *how* to work.** It is always in context; do not
re-derive its guidance here. Follow, in particular:
- **"Dev workflow"** — the verification cadence (targeted `test(filter=)` while iterating; prefix
  full-suite runs with `TESTTHAT_CPUS=8`; `check()` subsumes `test()`, so don't stack
  `test()`→`check()`→`coverage()`; `check(vignettes = FALSE)` for mid-work checks; the `cached()`
  fit memo and its caveats; no two concurrent package-touching R processes).
- **"Definition of done (every change)"** — including the one-command gate `Rscript tools/dod-gate.R`
  (check 0/0/0 → coverage → style → lint → `pkgdown::check_pkgdown()`), run **once** before the PR.
- **"Git"** and **"Ask-first / guardrails"** — the branch→PR→squash-merge protocol, the
  release-milestone CI exception, and when to flag before acting.

What this skill adds on top of that shared guidance:

- **Build in small, reviewable increments:** one coherent change → relevant test(s) → commit →
  repeat. Do not produce one large commit at the end.
- After any roxygen change, run `devtools::document()` and commit the regenerated `NAMESPACE`/docs
  in the same increment.
- **When you add or remove an exported object** (the regenerated `NAMESPACE` gained/lost an
  `export()`/`S3method()`, or you added a dataset under `data/`), update `_pkgdown.yml`'s
  `reference:` list to cover it **in the same commit**. This does not break local `R CMD check`,
  only the pkgdown GitHub Action — the `pkgdown::check_pkgdown()` step in the DoD gate is the guard.
- Write or update tests for new behavior; update `NEWS.md` for user-visible changes.
- If you hit a design ambiguity not covered by `DESIGN.md`, stop and ask rather than deciding
  unilaterally.

**Milestone documentation (do this before the final commit — it is part of the definition of done).**
The rule is *one detailed home, everything else a pointer* — do not duplicate the milestone
narrative across files:
- **`MILESTONES.md`** (single source of truth): add a detailed `M$ARGUMENTS (done)` entry **in
  numeric order** (immediately after `M{$ARGUMENTS − 1}`; never append out of order or skip a
  number), in the style of existing entries — what shipped + durable decisions/amendments + final
  test count and `R CMD check` status.
- **CLAUDE.md "## Completed milestones"**: add the matching index entry — a **true one-liner**
  (title + a single clause naming the headline artifact; no second sentence — see that section's
  header). Then set "## Current focus" back to the status slot: `In flight: nothing. Last shipped:
  M$ARGUMENTS …`.
- **`NEWS.md`**: user-facing notes only.
- **DESIGN.md §1–14**: record any design-contract change in the relevant numbered section — **not**
  as a second copy of the milestone log (§15 is only a pointer to `MILESTONES.md`).
- **`ROADMAP.md`**: if this milestone had a section there, delete it (it's now shipped).

**Finishing the milestone.** Follow CLAUDE.md's "Git" section: push the feature branch, open a PR
into `master`, and **squash-merge immediately** once the local DoD is green
(`gh pr merge <n> --squash --delete-branch`) — do **not** gate on remote CI or use `--auto`. Then
sync local: `git switch master && git pull --ff-only`. If that fast-forward fails, stop and follow
CLAUDE.md's divergence-recovery note (inspect, then merge `origin/master` in — never force-push or
`reset --hard`). The **release-milestone exception** (wait for the full green CI matrix) is in
CLAUDE.md's Git section — flag it explicitly if this is such a milestone.

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step you trigger after reviewing the work.
