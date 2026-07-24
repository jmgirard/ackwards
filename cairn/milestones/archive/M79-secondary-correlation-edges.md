# M79: Secondary-correlation edges in the pruned/publication view

**Status:** done (2026-07-24, PR #84 https://github.com/jmgirard/ackwards/pull/84)

**Goal:** In the pruned (`drop_pruned`) view, optionally draw the between-level
correlations the single-strongest-ancestor primary view hides, in a channel
distinct from the sign encoding.

**Outcome:** `autoplot.ackwards()` gains `show_secondary` (default `FALSE`).
Under `drop_pruned = TRUE`, `show_secondary = TRUE` draws every kept cross-level
pair with `|r| >= cut_show` that is not the primary edge — cross-branch second
parents *and* same-lineage skip arcs (a direct skip-level `r` is its own
non-transitive fact, per D-032/D-017). `.drop_pruned_nodes()` now returns
`$secondary`, a set-difference on the all-pairs table it already builds (IP1).
The secondary layer is a separate `geom_segment` at `alpha = 0.4` and
`linewidth = min(0.3, 0.6 * width_val)` (thinner than the primary in every mode),
plain ends, inheriting the sign colour/linetype. Vignette + roxygen + NEWS.

**Decisions:** none (worked under IP1, GP2; no principle added or changed).

**Review:** three-lens — blame-history 0, prior-review 0, diff-bug 2. F1 (85)
fixed: hardcoded secondary `linewidth = 0.3` inverted "thinner" under a user
`edge_linewidth < 0.3`; now scales under the primary width, with a regression
test. F2 (40) logged (secondary unlabeled under `show_r` — intentional).
Gate: check 0/0/0, coverage 100%, lint 0.
