# M75: Beautify suggest_k() print output as an aligned criteria table

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP2
- **Branch/PR:** m75-suggest-k-print-table · https://github.com/jmgirard/ackwards/pull/79

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

- [x] AC1: `print.suggest_k()` emits the criteria as an aligned table with a
      header row whose columns are exactly the requested criteria (`k` first);
      verified by a new `snap_print()` snapshot in `test-print-snapshot.R` for a
      full-criteria object, reviewed as intentional.
- [x] AC2: Column alignment is ANSI-width-correct — padding computed on display
      width (e.g. `cli::ansi_nchar`/`cli::ansi_align`), not raw `nchar`; a test
      strips ANSI from the rendered header + data rows and asserts every column's
      cells share the same start/end offsets (numeric columns right-aligned).
- [x] AC3: The output carries a legend line explaining the glyphs (`✔` = retained
      up to the criterion's ceiling, `*` = that criterion's optimal k); asserted
      present in the snapshot / a `grepl` check.
- [x] AC4: Behavior preserved — a subset-criteria object shows only those columns
      (snapshot), the M42/m1 undetermined-consensus object still prints
      `undetermined` with no `Inf` (`test-suggest_k.R:183` still passes), and the
      CD-unavailable footer (matrix-input and EFAtools-absent wordings) still
      renders.
- [x] AC5: `vignettes/ackwards-suggest-k.Rmd.orig` prose (~L212-231) describes the
      new table + legend, re-precomputed via `Rscript vignettes/precompute.R`
      (untouched vignettes/assets reverted per M61), vignette-freshness check
      passes, and NEWS.md records the user-visible output change.
- [x] AC6: `Rscript tools/dod-gate.R` clean — check 0 err/0 warn/0 note, coverage,
      style, lint, pkgdown.

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T1, T2
- AC4 → T1, T2
- AC5 → T3
- AC6 → T4

## Tasks

- [x] T1: Tests-first. In `test-print-snapshot.R` add `snap_print()` snapshots for
      (a) a full-criteria `suggest_k` object and (b) a subset-criteria object; add
      an ANSI-stripped alignment assertion (header + rows share column offsets,
      numeric columns right-aligned) and a legend-present check. Confirm the
      existing M42/m1 undetermined test (`test-suggest_k.R:183`) still targets the
      new renderer.
- [x] T2: Rewrite the criteria block of `print.suggest_k()` (`R/suggest_k.R`) to
      build a header + rows table with ANSI-aware padding (reuse `.ok_glyph`;
      add/borrow a pad helper or `cli::ansi_align`), star optimal cells, emit via
      `cli::cli_verbatim`/`cli_text`, and print the legend. Keep dynamic columns,
      CD footer, recommendations, consensus (incl. undetermined/no-Inf), and
      cautionary notes intact. Make T1 green.
- [x] T3: Update `ackwards-suggest-k.Rmd.orig` prose to describe the table+legend;
      `Rscript vignettes/precompute.R`; `git checkout --` untouched
      `vignettes/*.Rmd` + `vignettes/assets/`; add a NEWS.md bullet.
- [x] T4: Run `Rscript tools/dod-gate.R`; resolve any snapshot/coverage/lint fallout.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: T1+T2 done — aligned criteria table + adaptive legend in print.suggest_k(); snapshot + ANSI-alignment tests added. Not-retained glyph kept as grey dash (cli::symbol$bullet aliases to "*" in ASCII terminals, colliding with the optimal star). suggest_k + print-snapshot suites 0 fail/0 warn.
- 2026-07-23: T3 done — vignette prose (ackwards-suggest-k) describes the table+legend; precompute re-run; intro + girard vignettes regenerated (they print(sk), so the new table appears); engines/visualization + autoplot PNGs reverted as gt-ID/render noise (M61). NEWS.md bullet added. Vignette-freshness OK.

## Decisions

## Review

**Fresh evidence (2026-07-23, branch m75-suggest-k-print-table @ 4 ahead / 0 behind master):**

- AC1 — `devtools::test(filter="print-snapshot")`: 22 pass / 0 fail / 0 skip. The
  full-criteria snapshot in `_snaps/print-snapshot.md` shows the header row
  `k  PA-PC  PA-FA  MAP  VSS-1  VSS-2  CD` with one column per requested criterion.
- AC2 — the "renders a column-aligned table with a legend" test passes: it strips
  ANSI, asserts all table lines share one width, and asserts the MAP decimal-point
  character index is constant across every data row (right-alignment proof).
- AC3 — same test asserts a legend line carrying both "retained" and "optimal";
  snapshot shows `✔ retained   * optimal k   - not retained`.
- AC4 — subset snapshot (`criteria = c("map","vss")`) shows only `k/MAP/VSS-1/VSS-2`
  columns and a star-only legend; `test-suggest_k.R:183` (M42/m1 undetermined, no
  `Inf`) passes (3 assertions, 0 fail); full `suggest_k` suite 0 fail/0 warn covers
  the CD-unavailable footers.
- AC5 — `tools/check-vignette-freshness.R`: "Vignette freshness OK". NEWS.md carries
  the "aligned criteria table" bullet (grep=1). suggest-k/intro/girard `.Rmd`
  regenerated; engines/visualization + autoplot PNGs reverted as gt-ID/render noise.
- AC6 — `devtools::check()`: 0 errors / 0 warnings / 0 notes. `covr::package_coverage()`:
  100% overall, 0 uncovered lines in `R/suggest_k.R`. `devtools::document()` no-diff;
  `pkgdown::check_pkgdown()` clean; styler/lintr clean (DoD gate T4).

**Consistency gate:** `cairn_validate` exit 0 ("all checks passed"; 82 pre-existing
advisory warnings on legacy milestone-ID references, none from M75). No principle
change (GP2 worked-under, not modified) → `cairn_impact` skipped. Toolchain
consistency-gate (document no-diff, pkgdown, NEWS entry, freshness, full check) all
clean.
