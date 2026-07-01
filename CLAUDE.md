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
- **M28** — CD correctness & honesty fix (`cd_rmse` trailing-zero bug; "minimize" label/roxygen corrected to sequential-test framing)
- **M29** — Strip milestone numbers from user-facing docs (`NEWS.md` `(M24)` tag removed; regression test guards `NEWS.md`/`README.md`/vignettes)
- **M30** — Citation hygiene (`inst/CITATION` Girard-only; `ackwards()` `@references` gains Forbes; README citation prose corrected)
- **M31** — Correctness & output-honesty sweep (ESEM fit row reports scaled variants under WLSMV/ULSMV/MLR — `p_value`/`CFI`/`TLI`/`RMSEA` + `BIC`; `_meets` cleanup; `cor = "polychoric"` + ML/MLR guard; `fa.parallel`/`seed` doc confirmed correct; intro/suggest_k vignette drift fixed)
- **M32** — API-shape & naming resolutions (`tidy(what="fit")` `index`→`statistic`; `k_max` kept in both `ackwards()`/`suggest_k()`, roxygen-disambiguated; cutoffs pass/fail flag output removed, `.fit_cutoffs()` kept for reference lines; variance reported as 0–1 `proportion`/`cumulative`; plus M31-deferred `$meta$estimator` + `summary()` scaled-fit footnote)
- **M33** — simulated Gaussian dataset (`sim16`: 1000×16 continuous, known 1→2→4 hierarchy, no ordinal-detection warning, guaranteed redundant-chain + artefact signals at `k_max=5`)
- **M34** — Pruning verb: extracted `prune()` as a standalone, pipeable S3 generic off `ackwards()` (five prune args removed from `ackwards()`; canonical `"artifact"` naming with `"artefact"` alias; manual + mixed pruning; edges recomputed fresh inside `prune()`, `x$edges` never mutated)
- **M35** — autoplot & visualization: sign-propagation bugfix (primary-parent edges now always non-negative per DESIGN §7); `autoplot()` `sign_by`/`magnitude_by` configurable, always-legended encodings (`cut_strong` retired; `mono` kept as wrapper); `direction="horizontal"` layout; `colour_*` aliases + `color_edge`; `ggsave` documented (not re-exported); code-coupled viz/intro/README prose
- **M36** — interpretation functions: `augment()` `append` (scores-only) + `id_cols` passthrough; `top_items()` `by=c("factor","item")` + variable-label display (`label (code)`, fit-time capture into `meta$item_labels`, `show_labels`); interpret-vignette prose edits (cut=0.5 lead, `by="item"`/label demos, metatrait/HiTOP naming advice) + intro score de-indexing
- **M37** — engines vignette (doc-only): at-a-glance table fixes (parallel EFA/ESEM, substrate/correlations/estimators rows, χ² symbol, no WLSMV parenthetical); ESEM reframed as continuous (ML/MLR/FIML) *and* ordinal (WLSMV); per-level-fit converse (bad fit weakens incident edges); autoplot `k_max=3` truncation note; sim16-vs-bfi25 framing; runnable `psych::corFiml()` MAR route for PCA/EFA (on continuous sim16, both caveats, forward-ref M38); Missing-data/Performance trims + `library(future)` style; Hu & Bentler citation. Epic renumbered M31–M39 (new code milestone M38 inserted).

## Current focus

M37 is complete (see `MILESTONES.md` for detail). Next up in the M31–M39 documentation/UX epic is
**M38** — a **code** milestone: promote `missing = "fiml"` to a first-class route for PCA/EFA that
auto-routes to `psych::corFiml()` (currently `missing="fiml"` errors for PCA/EFA). Not yet planned;
run `/plan-milestone 38` before starting. Its brief in `ROADMAP.md` carries the banked decisions —
string-valued `n_obs` (`"total"`/`"complete"`/`"effective"`) + the FIML-fit-index literature task
for the defensible default, the guard matrix, the DESIGN §9/§14 sign-off (it reverses a resolved
default), and the corFiml-speed note.

Remaining milestones in the epic:
M38 `missing = "fiml"` for PCA/EFA (corFiml auto-route; **code milestone** — DESIGN sign-off +
tests); M39 narrative & remaining prose (intro, suggest_k, ordinal, forbes, README).

These one-liners are a lossy index. The **full driving rationale, banked decisions, and the raw
pkgdown-review notes** behind M35–M39 live in [`ROADMAP.md`](ROADMAP.md) — read it before running
`/plan-milestone N` for any of them. `MILESTONES.md` remains the source of truth for *completed*
milestones; `ROADMAP.md` is its forward-looking counterpart for *pending* ones.

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

**Efficiency (the suite takes minutes — don't re-run it needlessly).** `check()` already runs the
full test suite *and* examples *and* rebuilds vignettes (~3 min: ~90s tests + ~74s vignettes + ~20s
examples), and `covr::package_coverage()` runs the suite *again* — so `test()` → `check()` →
`coverage()` at one gate executes the suite ~3×. Instead: iterate with **targeted**
`devtools::test(filter = "<x>")` / `testthat::test_file()`; run failing tests **once** in a way that
shows the details (capture `res <-` or use a non-silent reporter — never silent-then-rerun); skip the
vignette rebuild during mid-work checks with `check(vignettes = FALSE)`; and at the final gate run
`check()` **once** (it subsumes `test()`) then `coverage()` once. Never run two package-touching R
processes concurrently.

Scaffolding helpers: `usethis::use_r()`, `use_test()`, `use_package()`, `use_testthat(3)`,
`use_github_action("check-standard")`. Use testthat 3e, roxygen2 for all exported functions
(document the *why* of each default, runnable `@examples`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- **New/removed exported object → pkgdown reference updated.** Whenever `NAMESPACE` gains or loses
  an `export()`/`S3method()`, or a dataset is added under `data/`, update the `reference:` list in
  `_pkgdown.yml` to match, then verify with `pkgdown::check_pkgdown()`. This is the exact check the
  pkgdown GitHub Action runs — an exported topic missing from the reference index fails that
  workflow (and only that workflow; local `R CMD check` won't catch it). Run it at the DoD gate
  unconditionally (it's sub-second and gated on `rlang::is_installed("pkgdown")`); it costs nothing
  and closes a recurring CI-break.
- `devtools::check()` clean (run **once** at the gate — it subsumes `devtools::test()`; see the
  efficiency note under *Dev workflow*). Coverage checked once, not per sub-step.
- Styled and linted.
- Public-facing change reflected in NEWS.md and (if user-visible) the relevant `@examples`/vignette.
- For a milestone: a detailed entry added to `MILESTONES.md` **in numeric order** + a one-line
  index entry here under "Completed milestones". `MILESTONES.md` is the single source of truth for
  milestone history — never re-log it in DESIGN.md §15 (a pointer) or duplicate it across files.

## Git

- Default branch is **`master`**; it stays green and releasable. Milestone work happens on a
  feature branch (`m{N}-<slug>`) and merges to `master` via a **PR**, **squash-merged as soon as
  the local definition of done is green** (`devtools::check()` 0/0/0 + tests + coverage +
  style/lint) — *not* gated on remote CI. `master` is deliberately **not branch-protected** (owner
  decision, 2026-07-01): required checks would gate every merge on the ~8–15 min CI matrix, pure
  latency for a solo pre-CRAN repo where local `check()` already ran. CI still runs on every push
  as an after-the-fact signal but does not block the merge; **don't** use `gh pr merge --auto`
  (silently no-ops without required checks) and **don't** synchronously watch or background-poll
  CI. Don't commit milestone work (anything touching `R/`, `tests/`, `DESCRIPTION`, vignettes)
  straight to `master`. Trivial isolated doc fixes may go direct at the user's discretion — but
  **push them immediately.** An unpushed commit on `master` leaves local ahead of `origin/master`;
  `git pull` won't surface it (nothing to fetch), and the next milestone branch — cut from local
  `master` — will bundle it into that milestone's squash-merge and force a post-merge divergence.
  /plan-milestone step 8a now guards against branching from an ahead-of-origin `master`.
  - **Exception — release / CRAN-submission milestones** (e.g. CRAN-prep, version-bump/release):
    here you **do** wait for the *full* green CI matrix before merging, because CRAN runs exactly
    that matrix (macOS/Windows/Ubuntu × release/devel/oldrel) and will reject platform failures
    the local macOS `check()` can't see. For these, don't merge on local-green alone. This
    exception also reactivates once the package has real users or collaborators (a red `master`
    then blocks others / ships bugs) — treat that transition as the trigger to reconsider branch
    protection.
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

- **EAP scoring** for ordinal ESEM — declined (M28): EAP's shrinkage attenuates the cross-level correlations that bass-ackwards exists to measure; tenBerge covers the common case. Seam preserved in `compute_edges()` but implementing EAP is not planned.
- **Oblique rotation** — varimax is hardcoded; no `rotation` argument; oblique confounds the cross-level signal. No plans to add it.
- **Higher-order SEM / Schmid-Leiman** — out of scope per §2; `ackwards` is correlation-based, not SEM-based.
- **Bootstrap CIs on skip-level edges** — deferred from M5; still out of scope. Logged in
  `DESIGN.md` §14. (Structural artefact signals and φ-default for non-PCA redundancy were
  reactivated and completed in M25.)
