---
name: post-milestone-review
description: Read-only conformance audit of a completed milestone against CLAUDE.md acceptance criteria/invariants and the DESIGN.md contracts. Use after finishing a milestone, before starting the next.
disable-model-invocation: true
allowed-tools: Read, Grep, Glob, Bash
argument-hint: "[milestone-number]"
---

0. Confirm $ARGUMENTS is recorded as done: it has a detailed `M$ARGUMENTS (done)` entry in
   `MILESTONES.md` **and** a matching one-line entry in CLAUDE.md's "## Completed milestones"
   index. If either is missing, stop and flag: either the milestone has not been implemented yet
   (run /implement-milestone first) or the number is wrong. Do NOT check "## Current focus" —
   /implement-milestone clears that section on completion, so a successfully implemented milestone
   will never appear there.

# Post-Milestone Review: Milestone $ARGUMENTS

This is a **read-only audit**. Do not edit, create, or fix any files — report only.

1. Re-read `CLAUDE.md` and `DESIGN.md` fresh; don't rely on conversation memory of what they say.
2. Run the test suite and `devtools::check()` and report actual results, not assumptions.
3. For each acceptance criterion for Milestone $ARGUMENTS in `CLAUDE.md`: state **Met / Partially Met / Not Met**, with the specific file/test as evidence.
4. For each of the seven invariants in `CLAUDE.md`: state whether the implementation upholds it, citing the file/location. Flag silent violations even if no test currently catches them.
5. List untested edge cases implied by `DESIGN.md`/`CLAUDE.md` that the current test suite does not exercise (e.g., non-convergence handling, score-SD standardization, sign-alignment ties, algebra-vs-scores cross-check coverage).
6. Dependency hygiene: anything misplaced between `Imports`/`Suggests` in `DESCRIPTION`; any engine call missing an `rlang::check_installed()` guard.
7. Roxygen check: does documentation explain *why* defaults were chosen (per CLAUDE.md's documentation standard), not just *what*.
8. Milestone-log hygiene: confirm `MILESTONES.md` is the only detailed log and is internally
   consistent — the **numbered** entries run in numeric order with **no gaps** (M1..M$ARGUMENTS all
   present); the CLAUDE.md "## Completed milestones" index has exactly one matching line per
   numbered entry, and each is a genuine **one-liner** (title + single clause — flag any that has
   grown into a paragraph, the specific drift the single-log structure exists to prevent); and the
   narrative is not duplicated into DESIGN.md §15 (which must remain a pointer). The
   "Between-milestone changes" section at the end of `MILESTONES.md` is **exempt** from the numeric
   no-gaps check (it deliberately holds non-numbered PRs); just confirm any code PR merged since the
   prior milestone that isn't this milestone appears there. Flag any missing, out-of-order,
   duplicated, or over-long entry.
9. Output a single triaged list: **Blocking** (must fix before the next milestone) / **Should-fix** / **Nice-to-have**. End with one line: **READY** or **NOT READY** for the next milestone, and why.

Do not rewrite, fix, or refactor anything during this review — that's a separate step.
