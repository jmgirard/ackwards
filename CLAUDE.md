# CLAUDE.md — `ackwards`

Operating manual for AI-assisted development of this package. Read `DESIGN.md` (repo root) first
and treat it as the **source of truth** for all design decisions; this file covers *how we work*,
not *what we're building*. When this file and `DESIGN.md` disagree, `DESIGN.md` wins for design and
this file wins for process — and flag the conflict.

## What this is

`ackwards` is an R package implementing Goldberg's (2006) bass-ackwards method and modern
descendants (PCA / EFA / ESEM engines) for hierarchical structural analysis. Extract solutions from
1..k factors, then characterize the hierarchy via between-level factor-score correlations. Full
rationale, contracts, object spec, and resolved defaults are in `DESIGN.md`.

**Note:** Forbes (2023) footnote 3 cites this package (`github.com/jmgirard/ackwards`) as the
reference implementation of the extended bass-ackwards approach. Fidelity to the paper's algorithm
is the baseline contract for anything Forbes-related; additive enrichments are acceptable but the
default output must reproduce Forbes's examples exactly.

## Completed milestones

- **M1 (done):** PCA engine + algebra `compute_edges()` + result object + `print`/`tidy`/`glance`,
  validated against `psych::bassAckward()`.
- **M2 (done):** `ba_layout()` + `autoplot()` (adjacent-level diagram) + `suggest_k()`.
- **M3 (done):** EFA engine (`psych::fa()`, tenBerge weights, algebra path) + materialized-scores
  route + algebra-vs-scores cross-check tests.
- **M4 (done):** ESEM engine (`lavaan::efa()`, WLSMV, tenBerge weights, rotation-aware SEs,
  per-level fit indices) + `cor = "polychoric"` for all engines + `loadings_se` in level contract
  + convergence truncation tested + `estimator` argument.
- **M5 (done):** Forbes extension — `pairs = "all"`, `prune = "redundant"/"artefact"`, Tucker's φ
  chains (DFS enumeration, global retain set), annotated `autoplot()` (skip-level arcs, pruned
  fill), `tidy(what = "nodes")`, `augment.ackwards()` print caveat.
- **M6 (done):** Storage materialization + cfQ cleanup — `scores = TRUE` / `keep_fits = TRUE`
  storage, `augment.ackwards()`, `tidy(what = "scores")`, cfQ hard error, cross-check tests.
- **M7 (done):** Documentation — README.Rmd, intro vignette, pkgdown site, three targeted
  vignettes (engines, ordinal, Forbes extension).
- **M8 (done):** Plot customization — `show_r`/`r_digits`, `mono`, `show_level_labels`/
  `level_label_size`, `node_labels`, `primary_only`, `drop_pruned`/`compress_levels` on
  `autoplot.ackwards()`; private `.drop_pruned_nodes()` helper in `layout.R`.
- **M9 (done):** Visualization round 2 — `show_arrows`, `edge_linewidth`, `legend` on
  `autoplot.ackwards()`; new `ackwards-visualization.Rmd` vignette; Forbes vignette slimmed toward
  the paper; intro vignette trimmed and stale comment fixed.

## Current focus — no active milestone

If a step needs a design decision not covered in `DESIGN.md`, **stop and ask** rather than guessing.

## Invariants — do not violate without flagging

These encode hard-won reasoning from the design phase. Changing them is a design decision, not a
refactor.

1. **One edge path.** All between-level correlations go through `compute_edges()`. Use the algebra
   (`W'RW`, standardized) when scoring is linear; materialize scores only when nonlinear (EAP) or
   when the user asks. **Always** standardize by real score SDs `sqrt(diag(W'RW))` — never assume
   unit variance (Bartlett/oblique scores are not unit-variance).
2. **Keep the cross-check.** Retain the scores route even where algebra is the default, and keep the
   test asserting they agree within tolerance for linear engines.
3. **Light core, heavy opt-in.** The object always carries loadings/variance/fit/weights/edges/
   lineage/`R`/meta. `scores`, raw `fits`, raw `data` are NULL by default and recomputable.
4. **Sign alignment anchors to the primary parent**, not "all positive" (that's impossible).
5. **Lineage lives in edges, never in IDs.** `m{k}f{j}` are stable labels; parentage is in the edge
   structure.
6. **Loud defaults.** Announce consequential auto-choices via cli (e.g., the ordinal-detection
   warning). Advise loudly; never switch basis silently.
7. **Convergence is data, not an error.** A non-converging level warns + is skipped; the object
   still builds to the deepest converged level. Never let one bad level abort the run.

## Resolved defaults (see `DESIGN.md` §9, §14)

Orthogonal CF (`cfT`, κ = 1/p, ≈ varimax) rotation — **only supported rotation; `cfQ`/oblique is
out of scope** (it confounds the cross-level signal) · `cor = "pearson"` with ordinal-detection
warning · `tenBerge` scoring (on the active basis) · WLSMV estimator for ordinal ESEM ·
Forbes extension **off** · `k` required · sign `align = TRUE` · `scores`/`keep_fits` stored =
`FALSE`. Don't change these silently.

## Dependencies (see `DESIGN.md` §12)

Keep `Imports` lean: `stats`, `utils`, `cli`, `rlang`, `generics`. Everything else (`psych`,
`GPArotation`, `lavaan`, `ggplot2`, testing/docs infrastructure) goes in `Suggests`, gated by
`rlang::check_installed()`. **Do not add to `Imports` without flagging it.** **No Rcpp** — profile
first; the heavy compute already lives in compiled deps (§3).

Current `Imports`: `cli`, `generics`, `rlang`, `stats`, `utils`.
Current `Suggests`: `covr`, `ggplot2`, `GPArotation`, `knitr`, `lavaan (>= 0.6-13)`, `psych`,
`rmarkdown`, `testthat (>= 3.0.0)`. `suggest_k()` uses `psych::fa.parallel` and `psych::vss` (no
separate `EGAnet`/`paran` dep). Visualization uses `ggplot2` directly (no `ggraph`/`igraph`/
`tidygraph`). `methods` is **not** imported (no `methods::` usage). `clue` was removed in M5.

## Dev workflow

R >= 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard devtools loop:

```r
devtools::load_all()      # load for interactive testing
devtools::document()      # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()          # run testthat suite
devtools::check()         # full R CMD check
styler::style_pkg()       # format
lintr::lint_package()     # lint
```

Scaffolding helpers: `usethis::use_r()`, `use_test()`, `use_package()`, `use_testthat(3)`,
`use_github_action("check-standard")`. Use testthat 3e, roxygen2 for all exported functions
(document the *why* of each default, runnable `@examples`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- `devtools::check()` clean.
- Styled and linted.
- Public-facing change reflected in NEWS.md and (if user-visible) the relevant `@examples`/vignette.

## Git

- Default branch is **`master`**.
- **Do not touch** the `legacy` branch or the `v0-legacy` tag — they preserve the pre-AI code.
- Small, focused commits with imperative messages (e.g., `Add PCA engine and level contract`).
- Don't force-push `master`. Don't commit data, credentials, or large binaries.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don't invent a design decision.
- Adding an `Imports` dependency, introducing Rcpp, or changing a resolved default → flag for
  approval first.
- Touching `git history`, `legacy`, tags, or anything destructive → confirm first.
- Prefer wrapping established engines (`psych`, `lavaan`, `GPArotation`) over reimplementing
  numerics.

## Out of scope for now

- **EAP scoring** for ordinal ESEM — deferred; linear tenBerge scores cover the common case.
- **Oblique rotation** — `cfQ` hard-errors as unsupported (confounds cross-level signal). No plans to add it.
- **Higher-order SEM / Schmid-Leiman** — out of scope per §2; `ackwards` is correlation-based, not SEM-based.
- **M5 deferred improvements** — structural artefact signals (split-then-merge, <3-item factors,
  orphans), φ as default for non-PCA redundancy, bootstrap CIs on skip-level edges. Logged in
  `DESIGN.md` §14.
