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

One line each; **full detail lives in [`MILESTONES.md`](MILESTONES.md)** (the single source of
truth). Add new milestones there in numeric order as part of the definition of done.

- **M1** — PCA engine + `compute_edges()` algebra + result object + `print`/`tidy`/`glance`
- **M2** — `ba_layout()` + `autoplot()` adjacent-level diagram + `suggest_k()`
- **M3** — EFA engine (tenBerge) + materialized-scores route + algebra-vs-scores cross-check
- **M4** — ESEM engine (lavaan WLSMV) + `cor="polychoric"` + `loadings_se` + `estimator`
- **M5** — Forbes extension (`pairs="all"`, `prune`, Tucker's φ chains, annotated `autoplot()`)
- **M6** — Storage materialization (`keep_scores`/`keep_fits`) + `augment()` + cfQ cleanup
- **M7** — Documentation (README.Rmd, intro + engines/ordinal/forbes vignettes, pkgdown)
- **M8** — Plot customization (`autoplot.ackwards()` args + `.drop_pruned_nodes()`)
- **M9** — Visualization round 2 + `ackwards-visualization.Rmd`
- **M10** — Conformance + robustness (`summary()`, ESEM Heywood warning, spearman+esem warning)
- **M11** — Edge-label polish + `show_r` decoupling (APA `.format_r()`)
- **M12** — Best-practice `suggest_k` (PA-FA, VSS-1/2, CD) + `autoplot.suggest_k()`
- **M13** — Rotation honesty (removed `kappa`/`rotation` args; cfT → varimax)
- **M14** — Dedicated `suggest_k()` vignette
- **M15** — Naming clarity pass (`k`→`k_max`, `method`→`engine`, `scores`→`keep_scores`, …)
- **M16** — Estimator-aware missing-data handling (`missing=` arg)
- **M17** — GitHub 0.1.0 release prep (MIT license, `inst/CITATION`, version bump)
- **M18** — Factor interpretation & label scaffolding (`top_items()`, `label_template()`)
- **M19** — Dedicated interpretation/labeling vignette
- **M20** — CRAN submission readiness + example legibility
- **M21** — Onboarding & usability pass (`psych`→Imports, drop `GPArotation`, `bfi25` dataset)
- **M22** — Correlation-matrix input (PCA/EFA)
- **M23** — Test-coverage hardening (→ 100%)
- **M24** — Vignette communication pass (`gt` comparison tables)
- **M25** — Deferred-items pass (`suggest_k` `criteria=`, artefact signals, φ auto-default)
- **M26** — Faster ESEM on large item sets (cached sample stats + parallel per-level fits)
- **M27** — ESEM fit & SEs as first-class output (glance fit, wide fit table, cutoff flags, loading CIs, fit plot, vignette framing)

## Current focus

**M28 — CD correctness & honesty fix.** Make `suggest_k()`'s Comparison Data (CD) output and
plot faithfully represent `EFAtools::CD`'s sequential one-sided Wilcoxon test (`alpha = 0.3`),
not a minimization.

Two verified defects:
1. **Trailing-zero bug.** `EFAtools::CD` fills `RMSE_eigenvalues` columns only up to `k_cd + 1`;
   the rest stay literal `0`. `suggest_k()` averages/plots all `k_max` columns, so the RMSE curve
   plunges to 0 on the right and `which.min()` lands on a fake zero (e.g. `bfi25`, `k_max = 8`:
   argmin = 8, star at 6). Fix: mask `cd_rmse` columns beyond `min(k_cd + 1, k_max)` to `NA`
   (keep the tested-but-rejected `k_cd + 1` level visible).
2. **Mislabeling.** The plot facet `"CD (RMSE, minimize)"` and the `autoplot.suggest_k` roxygen
   ("first crosses below the comparison-data average") falsely frame CD as minimization. Fix:
   rename the facet to an honest label, correct the roxygen to describe the sequential-test
   stopping rule and that the starred k need not be the visible minimum. (Vignette prose is
   already correct.)

Also: record the two deferred-item assessments in DESIGN.md §14 — **EAP scoring declined**
(attenuates the cross-level signal the method measures; `compute_edges()` seam preserved,
tenBerge-on-polychoric covers the case); **bootstrap CIs on skip-level edges remain deferred**
(own milestone; perf-heavy, Forbes-only).

Acceptance: (AC1) `cd_rmse` has no spurious `0`, entries beyond `k_cd + 1` are `NA`,
`min(cd_rmse, na.rm=TRUE)` is a genuinely computed level. (AC2) CD panel draws only over
computed levels (no plunge to 0), stars `k_cd`, facet no longer says "minimize". (AC3) roxygen +
vignette describe the sequential-test rule; no doc claims CD minimizes RMSE. (AC4) DESIGN §14
records EAP declined-with-rationale and bootstrap CIs still-deferred; CLAUDE.md "out of scope"
reflects the EAP decision. (AC5) `check()` clean, styled, linted, coverage retained.

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

**Varimax** (orthogonal) rotation — hardcoded internal constant since M13; no `rotation` argument;
oblique rotation **out of scope** (it confounds the cross-level signal) · `cor = "pearson"` with ordinal-detection
warning · `tenBerge` scoring (on the active basis) · WLSMV estimator for ordinal ESEM ·
Forbes extension **off** · `k_max` required · sign `align_signs = TRUE` · `keep_scores`/`keep_fits` stored =
`FALSE` · `redundancy_phi`: `NULL` (default) auto-resolves — `"pca"` → no φ filter (exact W′RW
algebra); `"efa"`/`"esem"` → `0.95` (Lorenzo-Seva & ten Berge 2006; factor-score indeterminacy
off-PCA makes `|r|`-only liberal). `NA` is the explicit opt-out. Announce auto-resolve loudly
(Invariant 6). Don't change these silently.

## Dependencies (see `DESIGN.md` §12)

`psych` is in **Imports** (M21) — it is the engine substrate for the default PCA and EFA paths and
for polychoric correlations; placing it in Suggests would require an install prompt for core
functionality. The SEM + plotting + optional-criterion stacks remain in `Suggests`. `GPArotation`
was **removed entirely** (M21) — varimax routes through base `stats::varimax` and GPArotation never
enters `loadedNamespaces()` on any supported path. **Do not add further to `Imports` without
flagging it.** **No Rcpp** — profile first; the heavy compute already lives in compiled deps (§3).

Current `Imports`: `cli`, `generics`, `psych`, `rlang`, `stats`, `utils`.
Current `Suggests`: `covr`, `EFAtools`, `future`, `future.apply`, `ggplot2`, `gt`, `knitr`,
`lavaan (>= 0.6-13)`, `rmarkdown`, `testthat (>= 3.0.0)`. (`future` itself is declared because the
parallel test calls `future::plan()` directly; `future.apply` is the dispatch backend.) `suggest_k()` uses
`psych::fa.parallel(fa="both")` + `psych::vss` (PA-PC, PA-FA, MAP, VSS-1/2) and optionally
`EFAtools::CD()` (gated by `rlang::is_installed()`); no separate `EGAnet`/`paran` dep.
`future.apply` (M26) is the optional parallel backend for the ESEM per-level fits, gated by
`rlang::is_installed()` with a serial `lapply` fallback — users opt in via `future::plan()`; no
`ncores` arg, no `future`/`parallel` in Imports. Visualization uses `ggplot2` directly
(no `ggraph`/`igraph`/`tidygraph`). `methods` is **not** imported (no `methods::` usage). `clue`
was removed in M5.

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
- For a milestone: a detailed entry added to `MILESTONES.md` **in numeric order** + a one-line
  index entry here under "Completed milestones". `MILESTONES.md` is the single source of truth for
  milestone history — never re-log it in DESIGN.md §15 (a pointer) or duplicate it across files.

## Git

- Default branch is **`master`**; it stays green and releasable. Milestone work happens on a
  feature branch (`m{N}-<slug>`) and merges to `master` via a **PR** once CI is green — don't
  commit milestone work (anything touching `R/`, `tests/`, `DESCRIPTION`, vignettes) straight to
  `master`. Trivial isolated doc fixes may go direct at the user's discretion.
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
- **Oblique rotation** — varimax is hardcoded; no `rotation` argument; oblique confounds the cross-level signal. No plans to add it.
- **Higher-order SEM / Schmid-Leiman** — out of scope per §2; `ackwards` is correlation-based, not SEM-based.
- **Bootstrap CIs on skip-level edges** — deferred from M5; still out of scope. Logged in
  `DESIGN.md` §14. (Structural artefact signals and φ-default for non-PCA redundancy were
  reactivated and completed in M25.)
