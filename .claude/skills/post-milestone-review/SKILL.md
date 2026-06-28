---
name: post-milestone-review
description: Read-only conformance audit of a completed milestone against CLAUDE.md acceptance criteria/invariants and the DESIGN.md contracts. Use after finishing a milestone, before starting the next.
model: opus
disable-model-invocation: true
allowed-tools: Read, Grep, Glob, Bash
argument-hint: [milestone-number]
---

# Post-Milestone Review: Milestone $ARGUMENTS

This is a **read-only audit**. Do not edit, create, or fix any files — report only.

1. Re-read `CLAUDE.md` and `DESIGN.md` fresh; don't rely on conversation memory of what they say.
2. Run the test suite and `devtools::check()` and report actual results, not assumptions.
3. For each acceptance criterion for Milestone $ARGUMENTS in `CLAUDE.md`: state **Met / Partially Met / Not Met**, with the specific file/test as evidence.
4. For each of the seven invariants in `CLAUDE.md`: state whether the implementation upholds it, citing the file/location. Flag silent violations even if no test currently catches them.
5. List untested edge cases implied by `DESIGN.md`/`CLAUDE.md` that the current test suite does not exercise (e.g., non-convergence handling, score-SD standardization, sign-alignment ties, algebra-vs-scores cross-check coverage).
6. Dependency hygiene: anything misplaced between `Imports`/`Suggests` in `DESCRIPTION`; any engine call missing an `rlang::check_installed()` guard.
7. Roxygen check: does documentation explain *why* defaults were chosen (per CLAUDE.md's documentation standard), not just *what*.
8. Output a single triaged list: **Blocking** (must fix before the next milestone) / **Should-fix** / **Nice-to-have**. End with one line: **READY** or **NOT READY** for the next milestone, and why.

Do not rewrite, fix, or refactor anything during this review — that's a separate step.
