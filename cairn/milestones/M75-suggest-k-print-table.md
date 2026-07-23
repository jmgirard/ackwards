# M75: Beautify suggest_k() print output as an aligned criteria table

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP2
- **Branch/PR:** m75-suggest-k-print-table

## Goal

Render `print.suggest_k()`'s criteria block as a column-aligned monospace table
with a header row and an in-output glyph legend, replacing the per-`k`
concatenated lines.

## Scope

**In:** Redesign the criteria display inside `print.suggest_k()`
(`R/suggest_k.R:541-693`) into an aligned table — one header row (`k` + one
column per *requested* criterion: PA-PC, PA-FA, MAP, VSS-1, VSS-2, CD), numeric
cells (MAP/VSS/CD-RMSE) right-aligned, PA/CD cells shown as retain (`✔`) / not
(`·`), the per-criterion optimal `k` starred, plus a compact legend line
explaining the glyphs. Padding is ANSI-width-aware (colored glyphs). Preserve
every existing behavior: dynamic column set (only requested criteria), the
CD-unavailable footer (matrix-input vs. EFAtools-absent), the recommendations
block, the consensus range / undetermined / no-`Inf` path (M42/m1), and the
cautionary footer. Update the suggest-k vignette prose + re-precompute; note the
change in NEWS.md.

**Out:** `autoplot.suggest_k()` styling (already polished; stays as-is — deferred
to a future demand-gated milestone if wanted). Any change to the returned
object's fields, the criteria semantics (D-013), or the report-only stance
(D-014). No new dependency — rendering stays within `cli` (GP2/GP5).

## Acceptance criteria

- [ ] AC1: `print.suggest_k()` emits the criteria as an aligned table with a
      header row whose columns are exactly the requested criteria (`k` first);
      verified by a new `snap_print()` snapshot in `test-print-snapshot.R` for a
      full-criteria object, reviewed as intentional.
- [ ] AC2: Column alignment is ANSI-width-correct — padding computed on display
      width (e.g. `cli::ansi_nchar`/`cli::ansi_align`), not raw `nchar`; a test
      strips ANSI from the rendered header + data rows and asserts every column's
      cells share the same start/end offsets (numeric columns right-aligned).
- [ ] AC3: The output carries a legend line explaining the glyphs (`✔` = retained
      up to the criterion's ceiling, `*` = that criterion's optimal k); asserted
      present in the snapshot / a `grepl` check.
- [ ] AC4: Behavior preserved — a subset-criteria object shows only those columns
      (snapshot), the M42/m1 undetermined-consensus object still prints
      `undetermined` with no `Inf` (`test-suggest_k.R:183` still passes), and the
      CD-unavailable footer (matrix-input and EFAtools-absent wordings) still
      renders.
- [ ] AC5: `vignettes/ackwards-suggest-k.Rmd.orig` prose (~L212-231) describes the
      new table + legend, re-precomputed via `Rscript vignettes/precompute.R`
      (untouched vignettes/assets reverted per M61), vignette-freshness check
      passes, and NEWS.md records the user-visible output change.
- [ ] AC6: `Rscript tools/dod-gate.R` clean — check 0 err/0 warn/0 note, coverage,
      style, lint, pkgdown.

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T1, T2
- AC4 → T1, T2
- AC5 → T3
- AC6 → T4

## Tasks

- [ ] T1: Tests-first. In `test-print-snapshot.R` add `snap_print()` snapshots for
      (a) a full-criteria `suggest_k` object and (b) a subset-criteria object; add
      an ANSI-stripped alignment assertion (header + rows share column offsets,
      numeric columns right-aligned) and a legend-present check. Confirm the
      existing M42/m1 undetermined test (`test-suggest_k.R:183`) still targets the
      new renderer.
- [ ] T2: Rewrite the criteria block of `print.suggest_k()` (`R/suggest_k.R`) to
      build a header + rows table with ANSI-aware padding (reuse `.ok_glyph`;
      add/borrow a pad helper or `cli::ansi_align`), star optimal cells, emit via
      `cli::cli_verbatim`/`cli_text`, and print the legend. Keep dynamic columns,
      CD footer, recommendations, consensus (incl. undetermined/no-Inf), and
      cautionary notes intact. Make T1 green.
- [ ] T3: Update `ackwards-suggest-k.Rmd.orig` prose to describe the table+legend;
      `Rscript vignettes/precompute.R`; `git checkout --` untouched
      `vignettes/*.Rmd` + `vignettes/assets/`; add a NEWS.md bullet.
- [ ] T4: Run `Rscript tools/dod-gate.R`; resolve any snapshot/coverage/lint fallout.

## Work log

- 2026-07-23: created by /milestone-plan.

## Decisions

## Review
