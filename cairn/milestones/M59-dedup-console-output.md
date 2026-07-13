<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section. -->
# M59: De-duplicate console output & plot builders

- **Status:** review
- **Priority:** normal
- **Depends on:** M58
- **Principles touched:** —   <!-- works under the cli-only-output and honesty-caveat conventions (D-001, D-008); numbers/strings preserved -->
- **Branch/PR:** m59-dedup-console-output · [PR #59](https://github.com/jmgirard/ackwards/pull/59)

## Goal

Route `print`, `summary`, and `autoplot`'s repeated output blocks through shared helpers so every user-visible string, count, and glyph has one source of truth.

## Scope

**In:**
- Shared helpers for the blocks currently copied between `print.ackwards` and `print.summary_ackwards`: the engine/rotation/basis/n/k header (`print.R:20-28` ≡ `summary.R:53-61`), the load-bearing honesty-caveat footer (`print.R:139-147` ≡ `summary.R:248-256`), and the near-singularity caution (`print.R:77-83` ≈ `summary.R:169-177`).
- One prune-count digest read by **both** surfaces, replacing the two divergent extraction paths (`print.ackwards` reaching into `x$prune$…` inline at `print.R:86-126` vs `summary` routing through `.summary_prune` at `summary.R:196-236`/`304-334`) — the structural-signal count, phi-note, and redundant-node count are currently computed twice.
- `autoplot`: `.ba_add_nodes(p, nodes, …)` node-drawing helper + a shared theme finisher used by both the main path and `.ba_degenerate_plot` (`autoplot.R:553-567` ≡ `631-645`; `606-612` ≡ `654-660`).
- Consistency unifications (snapshot-captured): one convention for the 1-decimal cumulative-variance percent (`print` `round` vs `summary` `sprintf` — `print.R:49` vs `summary.R:78`); one convention for the tick/cross glyphs (`col_green(symbol$tick)` vs raw `✔` — `print.R:50` vs `summary.R:131`); route `autoplot`'s fit-caption threshold (`autoplot.R:794`) through `.format_r()` like the edge labels.
- Comment fixes in touched files: `summary.R:88` stale hardcoded line-ref → symbol reference; `summary.R:70` restating comment removed; one-line rationale added to `layout.R` `.spread_positions` (`:214-235`) and `autoplot` `dodge_x`.

**Out:**
- Core/engine and `ackwards()` dedup → "Remaining dedup pass" candidate.
- No rewording of user-visible content — strings are consolidated, not changed, except the two deliberate glyph/percent unifications, which snapshot tests capture.

## Acceptance criteria

- [x] `print` and `summary` share one header helper, one honesty-caveat footer, one near-singularity caution, and one prune-count digest — no duplicated inline copies remain (grep-verifiable); existing `test-print.R` behavior preserved.
- [x] Cumulative-variance percent and convergence/fit glyphs use a single convention across `print` + `summary`; `autoplot`'s fit-caption threshold routes through `.format_r()`.
- [x] `autoplot`'s main and degenerate paths draw nodes and finish the theme via the shared `.ba_add_nodes()` helper + theme finisher.
- [x] `print`/`summary` snapshot tests are green, with the two deliberate output unifications captured in reviewed, intentional snapshot updates.
- [x] Comment fixes applied: `summary.R:88` stale ref → symbol; `summary.R:70` restating comment gone; rationale comments added to `.spread_positions` and `dodge_x`.
- [x] `Rscript tools/dod-gate.R` clean.

## Coverage

- AC1 → T1, T2
- AC2 → T3
- AC3 → T4
- AC4 → T6
- AC5 → T5
- AC6 → T6

## Tasks

- [x] T1 — Extract `.print_ba_header()`, the honesty-caveat footer, and the near-singularity caution into shared helpers; route both `print` methods through them.
- [x] T2 — Build one prune-count digest (extend/reuse `.summary_prune`) and have both `print.ackwards` and `summary` read it; remove the duplicated structural-signal / phi-note / redundant-count inline code in `print.R`.
- [x] T3 — Unify the cumulative-variance percent format and the tick/cross glyphs across `print` + `summary`; route `autoplot`'s threshold through `.format_r()`.
- [x] T4 — Extract `.ba_add_nodes()` + shared theme finisher; route both the main and `.ba_degenerate_plot` paths through them.
- [x] T5 — Fold in the comment fixes (`summary.R:88`, `summary.R:70`, `.spread_positions`, `dodge_x`).
- [x] T6 — Update/accept `print`/`summary` snapshots to lock the unified output; `Rscript tools/dod-gate.R` clean; commit with tracking.

## Work log

- 2026-07-12: created by /milestone-plan (bucket 2 of the 2026-07-12 codebase audit; sequenced after M58 to avoid utils.R churn overlap).
- 2026-07-12: set in-progress; branch m59-dedup-console-output.
- 2026-07-12: gate — glyph unification resolved to shared `.ok_glyph()` cli symbols, print keeps green/red, summary fit glyphs stay uncolored in grey line (user: "what do you recommend"). Percent → shared `.fmt_pct()` (sprintf 1-dp); prune digest → extend `.summary_prune` (renamed `.prune_digest`). No prior snapshot tests existed → establishing a snapshot net first (test-print-snapshot.R) to lock current output before refactor.
- 2026-07-12: T6(net) — snapshot safety net capturing current print/summary output (pca, efa fit-glyph line, pruning digest, near-singular caution) committed as baseline.
- 2026-07-12: T1 — `.print_ba_header()`, `.print_near_singular(min_eig, detail)`, `.print_honesty_footer(prune_note=NULL)` added to print.R; both surfaces routed through them. Near-singular `detail` and prune note wording differ per surface (≈, not ≡) so are passed in, preserving exact strings. Snapshots + test-print.R green (pure refactor).
- 2026-07-12: T2 — renamed `.summary_prune` → `.prune_digest`, moved to print.R shared helpers; `print.ackwards` now reads redundant/phi-note/structural counts from it instead of reaching into `x$prune$…` inline. Display wording stays per-surface. Snapshots + test-print.R green (pure refactor).
- 2026-07-12: T3 — added `.fmt_pct()` (sprintf 1-dp) + `.ok_glyph()` (cli::symbol tick/cross) to utils.R; print cum + summary cum/per-factor percent route through `.fmt_pct`; print convergence keeps green/red, summary fit-meets glyph now bare cli symbol (was raw ✔/✘). autoplot fit-caption threshold → `.format_r()`. Deliberate glyph change captured in accepted snapshot (✘→✖/x); percent convention locked by new `.fmt_pct` unit test (no snapshot diff on chosen data). New `.fmt_pct`/`.ok_glyph` unit tests. No new deps (dropped a withr reach).
- 2026-07-12: T4 — extracted `.ba_add_nodes(p, nodes, w, h)` + `.ba_finish_theme(p, legend, show_labels, direction)`; main render + `.ba_degenerate_plot` both routed through them. Layout/autoplot tests green (pure refactor).
- 2026-07-12: T5 — comment fixes: dropped `# pre-pulled for readability` (summary.R); stale `engine_pca.R:65` line-ref → symbol ref (`pca_levels()`'s `fit_info`); rationale added to `.spread_positions` (why spread+re-centre) and comparability `dodge_x` (dodge formula). Comment-only.
- 2026-07-12: T6 — DoD gate: check 0/0/0, style/lint/pkgdown clean. Coverage exposed a pre-existing gap: print's non-converged red-cross branch (`else` of the convergence glyph) was never tested — the old single-line `if/else` form masked it under line-based coverage; splitting it for `.ok_glyph` revealed it. Added a non-converged-level print snapshot (patches `converged <- FALSE`); coverage back to 100.00%. NEWS: one cosmetic bullet on the summary glyph. Status → review.

## Decisions

## Review

_Reviewed 2026-07-13 (same-session; evidence re-gathered by command). PR [#59](https://github.com/jmgirard/ackwards/pull/59)._

**Acceptance-criteria evidence (fresh):**
- AC1 — grep: `"Engine"` header, `series of linked solutions` caveat, `Near-singular correlation matrix` wrapper, and `.prune_digest` def each appear exactly once (all in `print.R`; zero in `summary.R`). `test-print.R` (incl. M50 single-rule test) green in the gate suite.
- AC2 — `.fmt_pct`/`.ok_glyph` used by both surfaces (`print.R:128,130,132`; `summary.R:68,76,125`); zero leftover raw `✔/✘/✔/✘` in print/summary; `autoplot.R:789` caption uses `.format_r(cut$threshold)`.
- AC3 — `autoplot.R:553,592` (main) and `:644,653` (degenerate) both route through `.ba_add_nodes()` + `.ba_finish_theme()`.
- AC4 — `test-print-snapshot.R` green in gate; accepted snapshot captures the glyph unification (EFA fit line `✘`→cli cross; non-converged print level shows the cross). Percent unification locked by the `.fmt_pct` unit test (no snapshot diff on the sampled data).
- AC5 — stale `engine_pca.R:65` ref gone; `# pre-pulled for readability` gone; rationale comments present in `layout.R` (`.spread_positions`) and `comparability.R` (`dodge_x`).
- AC6 — `Rscript tools/dod-gate.R`: **GATE PASSED** — check 0/0/0 (102s), coverage **100.00%**, styler clean, lintr clean, pkgdown index complete.

**Consistency gate:** `cairn_validate.py` → all checks pass (incl. coverage-complete, mirror, single in-progress, caps). No principle change (Principles touched: —) → `cairn_impact` skipped. r-package `consistency-gate` (check/style/lint/pkgdown) clean via the DoD gate.

**Independent fresh-context review (3 lenses + scorer):**
- [O] diff-bug (Opus): **No correctness findings.** Verified helper-by-helper output-equivalence (cli glue/`\\`-continuation, `.prune_digest` count mapping, `.ok_glyph` no ANSI nesting, autoplot layer order, `.format_r` threshold identity). One transparency note (dropped per taxonomy, not a defect): the diff carries two sanctioned output changes — glyph **and** print's percent trailing-`.0` — and NEWS mentioned only the glyph.
- [S] blame-history (Sonnet): **No findings.** M50 single-rule preserved; near-singular terse/detailed split intact; D-001 caveat byte-identical; `.summary_prune`→`.prune_digest` rename has zero stale refs.
- [S] prior-PR-comments (Sonnet): **No prior-PR evidence** (repo reviews run through cairn milestones + local-green merges; PRs carry no inline review comments) — lens no-ops.
- Scorer: **no scorable findings** (all three lenses clean) → nothing ≥80 to triage, nothing logged below 80.

**Action taken:** acted on the Opus transparency note (review-side doc fix) — broadened the NEWS bullet to document both console-consistency changes (glyph + print percent), not just the glyph. No code change.
