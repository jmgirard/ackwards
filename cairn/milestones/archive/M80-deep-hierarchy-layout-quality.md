# M80: Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging

**Status:** done (2026-07-24, PR #85 https://github.com/jmgirard/ackwards/pull/85)

**Goal:** Remove primary-tree edge crossings and edge-label collisions in deep (k>=10) bass-ackwards diagrams (Forbes's "levels got a bit bent" at ten levels).

**Outcome:** `ba_layout()` Stage-1 ordering now scores two candidates keep-best (lexicographic primary-then-all crossings): the historical single top-down barycenter seed (`.seed_order`) and a new `.primary_forest_order()` — a DFS over the primary-parent forest (one primary parent per node) that lays each subtree out contiguously. This drops primary-tree crossings 3->0 on the AMH `k_max=10` example; shallow layouts never regress. Pass-2 primary-child x-placement extracted verbatim to `.assign_x` (unchanged). Added `.count_crossings()`/`.count_crossings_xmap()` (adjacent-band Sugiyama crossing number) and `.dodge_edge_labels()` — a vectorised, deterministic force-directed repulsion wired into `autoplot(show_r=TRUE)`. Coordinate snapshots (`_snaps/layout.md`) lock the layout; all new helpers internal.

**Decisions:** Planned "alternating barycenter sweeps" replaced mid-flight (user-approved) by primary-forest ordering — sweeps yielded no reduction (single pass already a local optimum; all-edge metric insensitive). D-015 intact: forest seeds left-right order only, not a tree rendering. No D-entry (milestone-local).

**Review:** Three lenses — blame-history and prior-review clean; diff-bug (Opus) found no correctness bug, two test-adequacy gaps (scored 85/45). F1 (85) fixed: shallow test now asserts total crossings not increased vs seed. F2 (45) logged. Nothing graduated/retired.
