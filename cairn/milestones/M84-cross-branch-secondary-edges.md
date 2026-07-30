# M84: Cross-branch-only secondary edges in the pruned view

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, IP5, IP6
- **Branch/PR:** `m84-cross-branch-secondary-edges`

## Goal

Let the pruned view's secondary edges be restricted to cross-branch links, omitting
same-lineage skip arcs, matching the convention Forbes (2023) Figure 6B uses.

## Scope

**In:** lineage classification of `.drop_pruned_nodes()`'s secondary set by a walk up the
reduced primary-parent chain (IP5 — lineage from edges, never IDs); a `secondary_scope`
argument on `autoplot.ackwards()` defaulting to `"cross_branch"`; a cli warning when the
chosen scope leaves nothing to draw; the Figure 6B fidelity oracle extracted into
`cairn/references/forbes2023.md`; roxygen, NEWS, vignette, and a D-entry for the default.

**Out:** the non-pruned rendering path (`show_skip` draws adjacent + skip edges by design —
unchanged); which pairs count as secondary at all (M79's set-difference — unchanged);
manual ordering on intermediate levels (M81 `order` covers the deepest level only — no
current request); gap-tolerant redundancy chains and the oblique-rotation question →
ROADMAP candidates, both gated on the Forbes design session.

## Acceptance criteria

- [ ] AC1. `.drop_pruned_nodes()` returns its `secondary` table with an added logical
      column `same_lineage`, `TRUE` exactly when the row's `from` node lies on the `to`
      node's primary-parent chain in the reduced kept-only tree built from the returned
      `edges` table. A test on a fit carrying at least one `TRUE` and one `FALSE` row
      asserts the column against an ancestor walk computed in the test.
- [ ] AC2. On the AMH fit (`forbes2023`, `k_max = 10`) pruned to the node set of Forbes
      (2023) Figure 6B, `secondary_scope = "cross_branch"` yields exactly the three
      secondary edges that figure draws — `E1→J7`, `E2→G7`, `G7→J7` in her component
      labels — and no others, each with `same_lineage == FALSE`. (Source: forbes2023.pdf
      Fig. 6B, p. 9; its Note defines dashed lines there as "secondary component
      correlations .3 ≤ |r| < .9".)
- [ ] AC3. On that same fit `secondary_scope = "all"` draws every one of AC2's three edges
      plus ≥1 same-lineage arc Figure 6B omits — so the two scopes differ on Forbes's own
      data — and equals the full M79 set (every non-primary kept cross-level pair with
      `|r| >= cut_show`).
- [ ] AC4. `secondary_scope` is validated at the top of `autoplot.ackwards()` by
      `match.arg()`, so an unrecognized value errors naming the accepted values whatever
      `show_secondary`/`drop_pruned` are set to; `"cross_branch"` is in force when the
      argument is not supplied.
- [ ] AC5. Under `show_secondary = TRUE`, a scope leaving no edge with `|r| >= cut_show`
      emits a cli warning naming the scope and the cutoff, mirroring the primary path's
      existing "No edges with |r| >= cut_show to display"; a test fires it on a fit whose
      cross-branch set is empty (bfi25 `k_max = 4`: 0 of 3 cross-branch rows clear 0.3).
- [ ] AC6. The default change is announced (IP6): `NEWS.md` states under the development
      version that `show_secondary = TRUE` now draws cross-branch links only and names
      `secondary_scope = "all"` as the pre-change behavior; roxygen documents both values
      and why cross-branch is the default, citing D-033.
- [ ] AC7. `vignettes/ackwards-visualization.Rmd.orig` gains a passage rendering both
      scopes on a fit where they differ, `Rscript vignettes/precompute.R` is re-run and the
      regenerated `.Rmd` + `vignettes/assets/` committed, and `Rscript tools/dod-gate.R`
      exits 0.

## Coverage

- AC1 → T2
- AC2 → T1, T5
- AC3 → T5
- AC4 → T3
- AC5 → T4
- AC6 → T6
- AC7 → T7

## Tasks

- [ ] T1. Transcribe Figure 6B's secondary-edge convention into
      `cairn/references/forbes2023.md` from the rendered page image, not `pdftotext`
      (LESSONS 2026-07-19): the three dashed edges, their cross-branch classification
      against the solid primary tree, and the Note's `.3 ≤ |r| < .9` band.
- [ ] T2. Add `same_lineage` to the secondary table in `.drop_pruned_nodes()`
      ([layout.R:191](../../R/layout.R:191)) via a parent-chain walk over `edges_kept`;
      unit tests in `test-layout.R` beside the M79 block (line 428).
- [ ] T3. Add `secondary_scope` to `autoplot.ackwards()` — `match.arg` at the top beside
      `sign_by` ([autoplot.R:311](../../R/autoplot.R:311)), filter in the `show_secondary`
      branch ([autoplot.R:504](../../R/autoplot.R:504)); render tests for both scopes.
- [ ] T4. Add the empty-scope cli warning + regression test.
- [ ] T5. Add the Figure 6B fidelity test to `test-forbes-fidelity.R`, reusing its Forbes→
      ackwards label mapper (line 28) and the `cached()` AMH fit.
- [ ] T6. roxygen + `devtools::document()`; NEWS entry flagging the default change; D-entry.
- [ ] T7. Vignette passage + `Rscript vignettes/precompute.R`, reverting timing-only churn
      in untouched vignettes (LESSONS M61/M75); DoD gate.

## Work log

- 2026-07-30: created by /milestone-plan, from Forbes's 2026-07-30 email reply.
- 2026-07-30: criteria audit ([O], fresh context) returned 6 findings, all disposed at the
  plan gate — AC4 was internally contradictory (an unrecognized value had to both error and
  be silently ignored under `show_secondary = FALSE`); AC2 named `from`/`to` keys absent
  from `layer_data()`; AC2/AC3 would have passed vacuously on the M79 `cached()` fit
  (0 cross-branch rows above `cut_show`); AC3's "reproduces the 0.2.0 layer" had no in-repo
  referent; AC1 pinned neither column name nor polarity; AC5 named the freshness checker
  rather than the action it cannot enforce.
- 2026-07-30: plan gate chose a separate `secondary_scope` argument over folding a
  character mode into `show_secondary` because one dial that is both a flag and a mode
  reads worse in `match.arg` errors and in the roxygen; falsified by any call site needing
  the scope without `show_secondary`.
- 2026-07-30: plan gate chose `"cross_branch"` as the default over `"all"` because Figure 6B
  is the published convention and the argument shipped three days ago in 0.2.0, so
  installed-base cost is near zero; falsified by a user report of a changed figure.
- 2026-07-30: plan chose a parent-chain walk over a `level_to - level_from >= 2` proxy for
  same-lineage; the audit measured the two disagreeing on real fits (bfi25 `k_max = 8`:
  24 same-lineage / 6 cross-branch, not separable by gap); falsified by a fit where every
  gap-2+ secondary row is same-lineage.
- 2026-07-30: plan chose the Figure 6B edge set as the fidelity oracle over stipulate-and-
  self-check; the figure carries exactly three secondary edges and zero same-lineage arcs,
  so the criterion can fail; falsified by the figure proving unreadable at render time.

## Decisions

## Review
