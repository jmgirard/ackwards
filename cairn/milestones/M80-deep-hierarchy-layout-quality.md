# M80: Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** m80-deep-hierarchy-layout-quality

## Goal

Reduce edge crossings and label collisions in deep (k>=10) hierarchies by making
`ba_layout`'s **ordinal ordering** pass iterative (alternating up/down barycenter
sweeps) and adding edge-label dodging (candidate G — Forbes's "levels got a bit
bent" at ten levels).

## Scope

**In:**
- Replace the single top-down ordinal pass with a **primary-forest traversal
  ordering**: the primary edges form a forest (one primary parent per node), so a
  depth-first, subtree-contiguous leaf order eliminates primary-tree crossings.
  Chosen keep-best against the single-pass barycenter seed (lexicographic:
  primary crossings, then all crossings), so shallow layouts never regress.
  Deterministic tie-break (`rank(ties.method = "first")`). **Pass 2**
  (primary-child x-assignment) is left unchanged — it supplies the "parent sits
  directly above its primary child" property that AC2 and the existing test lock
  (gate 2026-07-24: ordering only). Amended 2026-07-24: barycenter *sweeps*
  yielded no reduction (single pass is already a local optimum) — see work log.
- Extract a pure, dependency-light `.dodge_edge_labels(lx, ly, threshold)` helper
  and wire it into `autoplot`'s edge-label path (`autoplot.R:592-598`) so
  overlapping midpoint labels (the `.99`/`1.00`/`.64` collisions Forbes saw) are
  offset. Helper is unit-testable without ggplot2 (gate 2026-07-24: pure helper).
- Add a pure `.count_crossings(layout)` helper and lock the layout with
  `testthat::expect_snapshot_value()` coordinate snapshots (shallow + k=10) — not
  vdiffr (gate 2026-07-24).
- Preserve invariants: every level on its own row (`y = -level`); level-1
  anchoring at `x = 0`; parentage/tree faithfulness unchanged.

**Out:**
- vdiffr / rendered-image snapshots — declined; coordinate snapshots cover layout
  regression (gate 2026-07-24).
- Iterating x-positions with barycenter sweeps — declined (gate 2026-07-24); would
  risk the direct-above-child alignment and its test.
- Manual factor-ordering **override** → candidate F (depends on this milestone; it
  overrides the Pass-1 ordinal permutation this milestone rewrites).
- Secondary edges → M79 (shipped; shared surface, independent).

## Acceptance criteria

- [ ] AC1: On the k=10 fixture (`ackwards(forbes2023, k_max = 10, pairs = "all")`),
      `ba_layout` produces **zero** primary-edge crossings (down from the
      single-pass baseline of 3, via `.count_crossings()` on the primary-edge
      subset) and never increases primary or total crossings on shallow fixtures.
- [ ] AC2: Layout stays faithful — each node at `y = -level`, parentage unchanged,
      every parent with primary children still at the mean x of its primary
      children (Pass 2 unchanged), and reordering confined to the Pass-1
      barycenter permutation within levels (no cross-level relabeling).
- [ ] AC3: `.dodge_edge_labels()` offsets overlapping edge labels so no two label
      anchors sit within the threshold distance; verified by a unit test on
      colliding inputs.
- [ ] AC4: Layout is deterministic — `expect_snapshot_value()` coordinate
      snapshots for a shallow fixture (`sim16`, k=5) and the k=10 fixture
      regenerate identically across runs.
- [ ] AC5: `devtools::check()` clean (0 err/0 warn/0 note); full suite green.

## Coverage

- AC1 → T1, T2
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T4

## Tasks

- [x] T1: (test-first) add a pure `.count_crossings(layout)` helper; assert the
      single-pass seed yields 3 primary crossings on the k=10 fixture and the new
      ordering yields 0 (the baseline the fix must beat).
- [x] T2: Replace Pass 1 with `.primary_forest_order()` chosen keep-best against
      the barycenter seed (lexicographic primary-then-all crossings) in
      `ba_layout`; leave Pass 2 x-assignment unchanged; update the roxygen
      docstring. Assert AC1 (0 primary crossings on k=10, no shallow regression)
      and AC2 (faithfulness + parent-above-primary-child).
- [x] T3: Extract `.dodge_edge_labels(lx, ly, threshold)` pure helper; call it from
      the `geom_label` path (`autoplot.R:592-598`); unit-test on colliding coords.
- [ ] T4: Add `expect_snapshot_value()` coordinate snapshots (shallow `sim16` k=5 +
      k=10 fixture); NEWS.md entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-24: re-planned via /milestone-plan after code investigation. Corrected
  the false "existing vdiffr snapshots" premise (none exist — layout is locked by
  structural invariants in `test-layout.R`); resolved four gate questions (see
  Decisions). Scope now: iterate Pass-1 ordering only, pure dodge + crossing-count
  helpers, coordinate snapshots.
- 2026-07-24: **substantive amendment (user-approved)** — the planned "alternating
  up/down barycenter sweeps" yield NO crossing reduction: the single top-down pass
  is already a local optimum (proven by sweeps + adjacent transposition both
  re-converging to it), and in `.assign_x` only level-K order is free. The all-edge
  crossing metric (3762 on k=10) is insensitive to ordering (dominated by the dense
  all-pairs secondary edges). Adopted `.primary_forest_order()` instead: the
  primary edges form a forest, and a subtree-contiguous DFS leaf order takes k=10
  primary crossings 3→0 (before/after rendered and confirmed). Amended Scope In(1)
  + AC1 accordingly; Pass-2 x-placement and D-015 unchanged (still a layered
  barycenter layout, not a tree rendering).

- 2026-07-24: T1+T2 done — added `.count_crossings()`/`.count_crossings_xmap()`,
  `.seed_order()`, `.primary_forest_order()`; rewired `ba_layout` Stage 1 to
  keep-best (seed vs traversal). k=10 primary crossings 3→0, shallow unchanged,
  deterministic. Tests + docstring updated; layout suite green.
- 2026-07-24: T3 done — `.dodge_edge_labels()` (vectorised force-directed
  repulsion, deterministic) wired into `autoplot`'s `show_r` label path; unit
  tests on coincident/band/separated/single inputs + a k=10 `show_r` render.
  Vectorised after an interpreted-loop version rendered the 312-label dense view
  in ~10s (now ~1.9s).

## Decisions

- Layout regression is locked by `expect_snapshot_value()` coordinate snapshots,
  not vdiffr — no new dependency, exact, and the coordinates *are* the layout
  (gate 2026-07-24).
- Edge-label dodging lives in a pure `.dodge_edge_labels()` helper (dependency-
  light, unit-testable without ggplot2), matching D-015's layout/render split
  (gate 2026-07-24).
- AC1 anchors on a `.count_crossings()` helper + a captured single-pass baseline
  constant on the k=10 fixture (gate 2026-07-24).
- Stage-1 ordering touches ordering only; Pass-2 primary-child x-assignment
  is unchanged, preserving the direct-above-child alignment (gate 2026-07-24).
- Ordering mechanism is primary-forest DFS (keep-best vs the barycenter seed),
  not barycenter sweeps, which yielded no reduction (amendment 2026-07-24). Within
  D-015: Pass-2 barycenter placement and all secondary/skip edges are unchanged —
  the forest seeds left-right order only, not a tree rendering.

## Review
