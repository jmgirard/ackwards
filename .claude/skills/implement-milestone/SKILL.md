---
name: implement-milestone
description: Implements the milestone plan already established in this conversation, following CLAUDE.md's dev workflow and definition of done. Use right after planning a milestone in this same session.
disable-model-invocation: true
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS matches the milestone number recorded in CLAUDE.md's "## Current focus"
   section. If it doesn't match, stop and flag the mismatch rather than proceeding.

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

Do not run `/post-milestone-review` yourself at the end — that's a separate, deliberate step you trigger after reviewing the work.