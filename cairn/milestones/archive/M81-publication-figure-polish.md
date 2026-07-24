# M81: Publication-figure polish — item lists, per-node box sizes, manual factor ordering

**Status:** done (2026-07-24, PR #86 https://github.com/jmgirard/ackwards/pull/86)

**Goal:** Add three opt-in `autoplot()`/`ba_layout()` controls for publication figures: item lists under the deepest-level factors, per-node box sizes, and manual left-to-right factor ordering.

**Outcome:** Three opt-in controls, defaults byte-identical. (a) `autoplot(show_items=, n_items=, item_cut=)` lists the top deepest-level items via `top_items()` (variable labels when present) — `.ba_item_text()` helper, below boxes (vertical) / right (horizontal), margin widened. (b) `node_width`/`node_height` accept a scalar or a named vector keyed by factor id (per-box override, unnamed fall back to 0.8/0.4) — `.resolve_node_size()` helper, per-node `nw`/`nh` mapped in `geom_tile`, edge faces + overlap warning made vector-safe. (c) `order=` on `ba_layout()` (and forwarded from `autoplot()`) pins the deepest level's order via `.resolve_manual_order()`; upper levels re-center over primary children (D-015 intact); skips M80 keep-best (collapses when the leaf order is fixed). NEWS + visualization vignette (3 new figures).

**Decisions:** none (milestone-local). Deepest-level-only `order` semantics chosen at plan gate over all-levels (upper-level independent pinning deferred; Out).

**Review:** three-lens fan-out (diff-bug/blame-history/prior-review) zero actioned findings; scorer no-op. One dropped non-finding: horizontal item-list crowding for long lists (accepted aesthetic tradeoff). DoD gate 0/0/0, coverage 100%. Nothing graduated/retired.
