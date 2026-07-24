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
- Make **Pass 1 (ordinal ordering)** iterative: alternating up/down barycenter
  sweeps (Sugiyama-style crossing reduction) until crossings stabilize or an
  iteration cap is hit, replacing the single top-down sweep at `layout.R:47-76`.
  Deterministic tie-break (`rank(ties.method = "first")`, as today). **Pass 2**
  (primary-child x-assignment, `layout.R:78-132`) is left unchanged — it supplies
  the "parent sits directly above its primary child" property that AC2 and the
  existing test lock (gate 2026-07-24: ordering only).
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
      `ba_layout` produces **strictly fewer** edge crossings than the captured
      single-pass baseline, measured by `.count_crossings()`; on shallow fixtures
      crossings are never increased.
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

- [ ] T1: (test-first) add a pure `.count_crossings(layout)` helper; capture the
      current single-pass algorithm's crossing count on the k=10 fixture as a
      documented baseline constant, asserting the fixture actually exhibits
      reducible crossings (if it does not, pick/construct a deep fixture that
      does).
- [ ] T2: Make Pass 1 ordering iterative (alternating up/down barycenter sweeps,
      iteration cap, deterministic tie-break) in `ba_layout` (`layout.R:47-76`);
      leave Pass 2 (`78-132`) unchanged; update the roxygen docstring. Assert AC1
      (new < baseline on k=10) and AC2 (faithfulness + parent-above-primary-child).
- [ ] T3: Extract `.dodge_edge_labels(lx, ly, threshold)` pure helper; call it from
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

## Decisions

- Layout regression is locked by `expect_snapshot_value()` coordinate snapshots,
  not vdiffr — no new dependency, exact, and the coordinates *are* the layout
  (gate 2026-07-24).
- Edge-label dodging lives in a pure `.dodge_edge_labels()` helper (dependency-
  light, unit-testable without ggplot2), matching D-015's layout/render split
  (gate 2026-07-24).
- AC1 anchors on a `.count_crossings()` helper + a captured single-pass baseline
  constant on the k=10 fixture (gate 2026-07-24).
- Iterative sweeps touch Pass-1 ordering only; Pass-2 primary-child x-assignment
  is unchanged, preserving the direct-above-child alignment (gate 2026-07-24).

## Review
