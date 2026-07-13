# M59: De-duplicate console output & plot builders — done

**Outcome:** `print.ackwards`, `print.summary_ackwards`, and `autoplot()` now
render their shared blocks through single-source helpers, so no user-visible
string/count/glyph is maintained twice. Added `.print_ba_header`,
`.print_near_singular(min_eig, detail)`, `.print_honesty_footer(prune_note)`,
`.prune_digest` (renamed from `.summary_prune`, now read by both surfaces) in
`print.R`; `.fmt_pct` + `.ok_glyph` in `utils.R`; `.ba_add_nodes` +
`.ba_finish_theme` shared by `autoplot`'s main & degenerate paths; `autoplot`
fit-caption threshold via `.format_r`. Plus the four planned comment fixes.

**Two deliberate cosmetic changes** (per-surface wording otherwise preserved):
`summary()`'s fit glyph now matches `print()`'s convergence mark
(terminal-adaptive `cli` symbol, was hard-coded `✔`/`✘`); `print()`'s
cumulative-variance percent is fixed 1-dp (`20.0%` not `20%`). Values unchanged.

**Key facts:** no prior snapshot tests → added `test-print-snapshot.R` safety
net (+ `.fmt_pct`/`.ok_glyph` unit tests). Splitting the one-line convergence
`if/else` surfaced a pre-existing untested branch that line-based coverage had
masked; covered it. No DESIGN/DECISIONS change; no new deps. Gate: 0/0/0,
coverage 100%, style/lint/pkgdown clean; 3 review lenses clean. Merged
local-green (non-release). PR: https://github.com/jmgirard/ackwards/pull/59
