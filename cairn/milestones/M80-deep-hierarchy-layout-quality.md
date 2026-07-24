# M80: Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging

- **Status:** planned
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** —

## Goal

Reduce edge crossings and label collisions in deep (k>=10) hierarchies by replacing `ba_layout`'s single-pass ordering with iterative barycenter sweeps and adding edge-label dodging (candidate G — Forbes's "levels got a bit bent" at ten levels).

## Scope

**In:**
- Replace the single top-down ordinal pass + single bottom-up x-assignment (`layout.R:47-132`) with **alternating up/down barycenter sweeps** (Sugiyama-style crossing reduction) iterated until crossings stabilize or an iteration cap is hit — resolving the order-by-parents vs position-by-children conflict the single pass cannot.
- Add **edge-label dodging** so overlapping midpoint labels (the `.99`/`1.00`/`.64` collisions Forbes saw) are offset.
- Preserve invariants: every level stays on its own horizontal row (`y = -k`); parentage/tree faithfulness unchanged; level-1 anchoring per the existing convention.
- vdiffr snapshots including a k>=10 case; a crossing-count helper asserting reduction vs the single-pass baseline on a Forbes-k=10-shaped fixture (e.g. AMH `k_max = 10`).

**Out:** manual factor-ordering *override* → F (candidate; it overrides this ordinal pass, so it depends on this milestone). Secondary edges → M79 (shared surface, independent).

## Acceptance criteria

- [ ] AC1: On a k>=10 fixture, `ba_layout` produces strictly fewer edge crossings than the single-pass baseline, measured by a crossing-count helper.
- [ ] AC2: Layout stays faithful — each node at `y = -k`, parentage unchanged, reordering limited to the barycenter permutation within levels (no cross-level relabeling).
- [ ] AC3: Overlapping edge labels are dodged — no two labels within a threshold distance share a position.
- [ ] AC4: Existing shallow-hierarchy vdiffr snapshots are unchanged, or updated with justification (layout stays deterministic).
- [ ] AC5: `devtools::test()` clean; `devtools::check()` clean.

## Coverage

- AC1 → T1, T2
- AC2 → T2
- AC3 → T3
- AC4 → T2, T3
- AC5 → T4

## Tasks

- [ ] T1: (test-first) add a crossing-count helper + a baseline assertion on a k>=10 fixture (deep sim or AMH `k_max = 10`).
- [ ] T2: Implement alternating barycenter sweeps in `ba_layout` (`layout.R:47-132`) with an iteration cap; hold the row (`y=-k`)/parentage invariants.
- [ ] T3: Implement edge-label dodging in `autoplot`'s label placement (the `geom_label`/midpoint path).
- [ ] T4: Refresh/justify vdiffr snapshots; NEWS.md entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.

## Decisions

## Review
