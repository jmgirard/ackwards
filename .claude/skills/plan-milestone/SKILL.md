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

3. Re-read CLAUDE.md and DESIGN.md fresh.
4. Confirm this milestone's scope against DESIGN.md §15 and CLAUDE.md's "out of scope for now"
   list; flag if the user's details conflict with either.
5. Propose a concrete plan: files to create/modify, the order of implementation, and testable
   acceptance criteria in the same style as Milestone 1's.
6. Flag any design ambiguity the brief doesn't resolve — do not invent a resolution.
7. Do not write any code. This is planning only; implementation happens via
   /implement-milestone afterward.

8. Once the user explicitly approves the plan, update CLAUDE.md's "## Current focus" section
   to replace the previous milestone's scope/acceptance criteria with Milestone $ARGUMENTS's,
   and commit that change on its own (not bundled with any implementation commit).