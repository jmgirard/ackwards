# M75: Beautify suggest_k() print output as an aligned criteria table

**Status:** done (2026-07-23, PR #79 https://github.com/jmgirard/ackwards/pull/79)

**Goal:** Render `print.suggest_k()`'s criteria block as a column-aligned table
with a header row and glyph legend, replacing the per-k concatenated lines.

**Outcome:** `print.suggest_k()` (`R/suggest_k.R`) now prints a header row (`k` +
one column per requested criterion) and one row per k, padded with ANSI-width-aware
`cli::ansi_align`/`ansi_nchar` so coloured glyphs align by display width. Numeric
columns (MAP/VSS) are right-aligned with the optimal-`*` in a reserved trailing
slot (decimals never shift); PA/CD show `✔`/`-` retention, CD also stars its
optimal k; an adaptive legend names only the glyphs shown. Dynamic columns, CD
footers, and the consensus/no-Inf block are unchanged; returned object and values
unchanged (display-only). Added two snapshots + an ANSI-alignment test; updated the
suggest-k vignette + intro/girard + NEWS.

**Decisions:** none promoted. Not-retained glyph is a grey dash, not
`cli::symbol$bullet` (collapses to `*` in ASCII terminals, aliasing the optimal
`*`); autoplot styling out of scope.

**Review:** three lenses + scorer. Diff-bug (Opus) + blame-history (Sonnet): no
defects. Prior-review (Sonnet): one M61 scoped-diff regression (score 90) — cli
`[Nms]` timing churn in regenerated intro/girard, fixed by reverting the 8 timing
lines to master. Two below-80 logged (12, 8). Gate: 0/0/0, coverage 100%, CI green.
