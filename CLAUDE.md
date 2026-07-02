# CLAUDE.md — ackwards

Operating manual for AI-assisted development of this package. Read
`DESIGN.md` (repo root) first and treat it as the **source of truth**
for all design decisions; this file covers *how we work*, not *what
we’re building*. When this file and `DESIGN.md` disagree, `DESIGN.md`
wins for design and this file wins for process — and flag the conflict.

## What this is

`ackwards` is an R package implementing Goldberg’s (2006) bass-ackwards
method and modern descendants (PCA / EFA / ESEM engines) for
hierarchical structural analysis. Extract solutions from 1..k factors,
then characterize the hierarchy via between-level factor-score
correlations. Full rationale, contracts, object spec, and resolved
defaults are in `DESIGN.md`.

**Note:** Forbes (2023) footnote 3 cites this package
(`github.com/jmgirard/ackwards`) as the reference implementation of the
extended bass-ackwards approach. Fidelity to the paper’s algorithm is
the baseline contract for anything Forbes-related; additive enrichments
are acceptable but the default output must reproduce Forbes’s examples
exactly. **This contract is test-backed** (M44):
`tests/testthat/test-forbes-fidelity.R` reproduces the paper’s three
simulation studies against expected values computed with Forbes’s own
reference implementation (OSF `pcwm8`; provenance in the fixture), and
the M44 feasibility study additionally verified the 155-variable AMH
applied example to 3.9e-14 (not committed as a test — the AMH matrix
carries no OSF license; see `ROADMAP.md` unscheduled items).

## Completed milestones

One line each; **full detail lives in
[`MILESTONES.md`](https://jmgirard.github.io/ackwards/MILESTONES.md)**
(the single source of truth). Add new milestones there in numeric order
as part of the definition of done.

- **M1** — PCA engine +
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  algebra + result object + `print`/`tidy`/`glance`
- **M2** —
  [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md) +
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  adjacent-level diagram +
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
- **M3** — EFA engine (tenBerge) + materialized-scores route +
  algebra-vs-scores cross-check
- **M4** — ESEM engine (lavaan WLSMV) + `cor="polychoric"` +
  `loadings_se` + `estimator`
- **M5** — Forbes extension (`pairs="all"`, `prune`, Tucker’s φ chains,
  annotated
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md))
- **M6** — Storage materialization (`keep_scores`/`keep_fits`) +
  [`augment()`](https://generics.r-lib.org/reference/augment.html) + cfQ
  cleanup
- **M7** — Documentation (README.Rmd, intro + engines/ordinal/forbes
  vignettes, pkgdown)
- **M8** — Plot customization
  ([`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  args + `.drop_pruned_nodes()`)
- **M9** — Visualization round 2 + `ackwards-visualization.Rmd`
- **M10** — Conformance + robustness
  ([`summary()`](https://rdrr.io/r/base/summary.html), ESEM Heywood
  warning, spearman+esem warning)
- **M11** — Edge-label polish + `show_r` decoupling (APA `.format_r()`)
- **M12** — Best-practice `suggest_k` (PA-FA, VSS-1/2, CD) +
  [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
- **M13** — Rotation honesty (removed `kappa`/`rotation` args; cfT →
  varimax)
- **M14** — Dedicated
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  vignette
- **M15** — Naming clarity pass (`k`→`k_max`, `method`→`engine`,
  `scores`→`keep_scores`, …)
- **M16** — Estimator-aware missing-data handling (`missing=` arg)
- **M17** — GitHub 0.1.0 release prep (MIT license, `inst/CITATION`,
  version bump)
- **M18** — Factor interpretation & label scaffolding
  ([`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md),
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md))
- **M19** — Dedicated interpretation/labeling vignette
- **M20** — CRAN submission readiness + example legibility
- **M21** — Onboarding & usability pass (`psych`→Imports, drop
  `GPArotation`, `bfi25` dataset)
- **M22** — Correlation-matrix input (PCA/EFA)
- **M23** — Test-coverage hardening (→ 100%)
- **M24** — Vignette communication pass (`gt` comparison tables)
- **M25** — Deferred-items pass (`suggest_k` `criteria=`, artefact
  signals, φ auto-default)
- **M26** — Faster ESEM on large item sets (cached sample stats +
  parallel per-level fits)
- **M27** — ESEM fit & SEs as first-class output (glance fit, wide fit
  table, cutoff flags, loading CIs, fit plot, vignette framing)
- **M28** — CD correctness & honesty fix (`cd_rmse` trailing-zero bug;
  “minimize” label/roxygen corrected to sequential-test framing)
- **M29** — Strip milestone numbers from user-facing docs (`NEWS.md`
  `(M24)` tag removed; regression test guards
  `NEWS.md`/`README.md`/vignettes)
- **M30** — Citation hygiene (`inst/CITATION` Girard-only;
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  `@references` gains Forbes; README citation prose corrected)
- **M31** — Correctness & output-honesty sweep (ESEM fit row reports
  scaled variants under WLSMV/ULSMV/MLR —
  `p_value`/`CFI`/`TLI`/`RMSEA` + `BIC`; `_meets` cleanup;
  `cor = "polychoric"` + ML/MLR guard; `fa.parallel`/`seed` doc
  confirmed correct; intro/suggest_k vignette drift fixed)
- **M32** — API-shape & naming resolutions (`tidy(what="fit")`
  `index`→`statistic`; `k_max` kept in both
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)/[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md),
  roxygen-disambiguated; cutoffs pass/fail flag output removed,
  `.fit_cutoffs()` kept for reference lines; variance reported as 0–1
  `proportion`/`cumulative`; plus M31-deferred `$meta$estimator` +
  [`summary()`](https://rdrr.io/r/base/summary.html) scaled-fit
  footnote)
- **M33** — simulated Gaussian dataset (`sim16`: 1000×16 continuous,
  known 1→2→4 hierarchy, no ordinal-detection warning, guaranteed
  redundant-chain + artefact signals at `k_max=5`)
- **M34** — Pruning verb: extracted
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) as
  a standalone, pipeable S3 generic off
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  (five prune args removed from
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md);
  canonical `"artifact"` naming with `"artefact"` alias; manual + mixed
  pruning; edges recomputed fresh inside
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md),
  `x$edges` never mutated)
- **M35** — autoplot & visualization: sign-propagation bugfix
  (primary-parent edges now always non-negative per DESIGN §7);
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  `sign_by`/`magnitude_by` configurable, always-legended encodings
  (`cut_strong` retired; `mono` kept as wrapper);
  `direction="horizontal"` layout; `colour_*` aliases + `color_edge`;
  `ggsave` documented (not re-exported); code-coupled viz/intro/README
  prose
- **M36** — interpretation functions:
  [`augment()`](https://generics.r-lib.org/reference/augment.html)
  `append` (scores-only) + `id_cols` passthrough;
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  `by=c("factor","item")` + variable-label display (`label (code)`,
  fit-time capture into `meta$item_labels`, `show_labels`);
  interpret-vignette prose edits (cut=0.5 lead, `by="item"`/label demos,
  metatrait/HiTOP naming advice) + intro score de-indexing
- **M37** — engines vignette (doc-only): at-a-glance table fixes
  (parallel EFA/ESEM, substrate/correlations/estimators rows, χ² symbol,
  no WLSMV parenthetical); ESEM reframed as continuous (ML/MLR/FIML)
  *and* ordinal (WLSMV); per-level-fit converse (bad fit weakens
  incident edges); autoplot `k_max=3` truncation note; sim16-vs-bfi25
  framing; runnable
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html) MAR
  route for PCA/EFA (on continuous sim16, both caveats, forward-ref
  M38); Missing-data/Performance trims +
  [`library(future)`](https://future.futureverse.org) style; Hu &
  Bentler citation. Epic renumbered M31–M39 (new code milestone M38
  inserted).
- **M38** — `missing = "fiml"` for PCA/EFA (code):
  `engine = "pca"/"efa"` + `cor = "pearson"` now routes `R` through
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
  (MAR-valid FIML) into the `W'RW` algebra (Invariant-1-clean, one
  corFiml call/run, no new dep), reversing the M16 “FIML errors for
  PCA/EFA” default; `.resolve_missing()` gains a `cor` guard (errors for
  non-Pearson PCA/EFA basis + WLSMV/ULSMV); `n_obs` string `"total"`
  (default, Enders 2010) / `"complete"` selects the approximate EFA
  fit-index N (Zhang & Savalei 2020; point estimates unaffected;
  announced via cli); `"effective"` dropped (no canonical formula);
  DESIGN §9/§14 (items 32–33) sign-off.
- **M39** — narrative & remaining prose (doc-only; final milestone of
  the M31–M39 epic): clarity pass across
  intro/suggest_k/ordinal/forbes/README — `print(sk)`/`print(x)`;
  orthogonal-varimax≡CF(1/p) rationale; SE/CI-NA-under-PCA note;
  `top_items(cut = 0.5)`; dropped `keep_scores` demo; **removed the
  README+intro red-arrow explanation** (no red arrow exists after the
  M35 sign fix); de-indexed README score columns; dropped the Waller
  citation nudge; suggest_k `sim16`-idealized vs `bfi25`-realistic
  contrast (inline-computed); ordinal binary/tetrachoric +
  WLSMV-polychoric + scores-trustworthy clarifications; forbes
  verbatim-heading rewrites + visible `redundancy_r` chunk +
  structural-table gt highlight + Lorenzo-Seva ref; refs alphabetized.
  Three code/viz asks (ordinal `categorical` flag, ordinal
  corr-comparison viz, forbes pruned-level label styling) spun off to
  **M40** (`ROADMAP.md`, DESIGN §14). No `R/`/NAMESPACE/export change.
- **M40** — deferred code/viz asks (final M31–M40 milestone): ordinal
  `categorical` flag **declined** (redundant — `cor = "polychoric"`
  already auto-selects WLSMV; would only add a conflict surface + §9
  change for zero capability; discoverability handled in docs); ordinal
  corr-comparison now a **dodged bar chart** (ten `N1`–`N5` item pairs,
  `fill = basis`, hidden reshape code) replacing the two raw
  `round(x$r)` matrices (also fixed stale N1–N2 figures 0.73→0.79);
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  **italicises a fully-pruned level’s** axis label (new
  `.fully_pruned_levels()`, `fontface` aesthetic through
  `.ba_level_labels()`, both directions) — partially-pruned levels stay
  plain. No new/removed export, no signature or dependency change.
- **M41** — independent Fable review (review-only): statistical core
  verified clean numerically
  (tenBerge/W′RW/signs/Forbes/ESEM-fit/suggest_k); findings — 1 Critical
  (EFA chi/p-value pairing), 6 Major (drop_pruned adjacent-pairs M34
  regression, pre-M38 engines-vignette FIML prose, CD-mechanism
  misstatement, Forbes artifact-zero framing, `cut_strong` remnant,
  untested Forbes-fidelity contract), 11 Minor, 4 enhancements; full
  §9/§14 defaults audit (all sound; one sound-but-misjustified wording);
  report + M42/M43/M44 triage in `ROADMAP.md`. No code change.
- **M42** — review fixes, code: EFA `chi` now the likelihood-ratio
  `STATISTIC` matching `p_value`/RMSEA/TLI (C1); `.drop_pruned_nodes()`
  recomputes all-pairs edges fresh, fixing the M34 `pairs = "adjacent"`
  regression (M1); `print.suggest_k` “undetermined” consensus (m1) +
  PA-cap announcements (m2); `cut_show`/`n_iter` validation (m8);
  ordinal warning names flagged columns (e3, `detect_ordinal()` returns
  names); dead `esem_levels(n_obs)` removed (m7); stale comments fixed
  (m6); EFA-aware fit-plot caption (m10). No export/signature/dependency
  change.
- **M43** — review fixes, docs (doc-only): engines vignette rewritten
  around first-class `missing = "fiml"` for PCA/EFA +
  `fm = "pca"`/sample-size fixes (M2, m4, m5); suggest-k CD mechanism
  corrected + worked-BFI prose inline-computed/drift-proof (M3, m11);
  Forbes artifact section rewritten to report-and-judge with a top-\|φ\|
  table, `cut_strong` remnant removed, chain-retention example corrected
  (M4, M5, m3); `redundancy_phi` PCA rationale corrected to score
  *determinacy* across DESIGN §9/CLAUDE/roxygen/vignette (e1); sim16 doc
  comments modernized (m9). No behavior/export change.
- **M44** — Forbes-fidelity fixture (closes the M41→M44 review arc):
  found the paper’s own OSF project (`pcwm8`: simulations script,
  reference implementation, AMH matrix); head-to-head vs her own
  functions matched edges to 3.9e-14 incl. the full 155-variable AMH
  example; shipped a 3.7 KB license-clean fixture (three
  seed-regenerated simulations + her implementation’s expected outputs)
  and `test-forbes-fidelity.R` (65 assertions: edges/φ/chase
  paths/retention); contract annotated **test-backed**; AMH commit
  deferred (no OSF license — options logged for owner outreach). No
  export/dependency change.
- **M45** — out-of-sample scoring (train/test): fit-time item moments
  stored (`meta$item_means`/`item_sds`; NULL for cor-matrix input);
  [`augment()`](https://generics.r-lib.org/reference/augment.html) gains
  `scaling = c("fit", "sample")` with **`"fit"` default** (training
  moments — one metric across train/test/subsets; `"sample"` = pre-M45
  opt-in and the cor-matrix route; DESIGN §14 item 34); new exported
  **`predict.ackwards(object, newdata, scaling)`** ≡
  `augment(append = FALSE)` (equivalence test-asserted; `_pkgdown.yml`
  updated); intro-vignette train/test subsection. New export, no new
  dependency.
- **M46** — Girard extension (replicability-gated hierarchies): new
  exported
  **[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)**
  — Everett (1983)/Goldberg (1990) split-half factor comparability per
  level per factor (full-sample-anchored labels, greedy-with-removal
  matching, cross-solution correlations through
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  on the pooled R, Tucker’s φ alongside, report-first/nothing
  auto-flagged; PCA/EFA + pearson/spearman, `n_splits = 10`, seeded) +
  `print`/`autoplot` methods; capstone vignette `ackwards-girard`
  (“Replicability-Gated Hierarchies: A Recommended Workflow”) with the
  six-step workflow + common-mistakes section; completes the triad
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  (range) ·
  [`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
  (floor) ·
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  (differentiation); DESIGN §14 item 35; ESEM/polychoric extensions
  deferred to `ROADMAP.md`. In passing: `cli::symbol$phi` glyph fix in
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md).
  New export, no new dependency.

## Current focus

**M46 is complete** (2026-07-01) — the Girard extension:
replicability-gated hierarchies. Every acceptance criterion met:
exported
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
(Everett 1983 / Goldberg 1990 split-half factor comparability,
full-sample-anchored labels, cross-solution correlations through
[`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
on the pooled R — Invariant 1; Tucker’s φ alongside; report-first,
nothing auto-flagged) with print/autoplot methods; self-comparability ≡
1, matching invariance, seed reproducibility, and the Invariant-2
cross-solution oracle all test-asserted; discriminating tests pass
strikingly (sim16’s true 1→2→4 levels ≈ .99 vs. the overextracted m5f5
at .14; bfi25’s floor lands exactly on the Big Five); capstone vignette
`ackwards-girard` + cross-links from intro/suggest-k/forbes + README;
DESIGN §14 item 35; ESEM/polychoric extensions logged in `ROADMAP.md`.
Post-review follow-up (same day) addressed all five findings:
NA-degradation guard in cross-solution matching, package-wide
unknown-`...` rejection (`ackwards`/`suggest_k`/ `comparability`),
`k_requested` in the object + truncation note in
[`print()`](https://rdrr.io/r/base/print.html), the `.fit_half()` error
branch exercised (no `# nocov`), and tests for `n_splits = 1` / spearman
/ truncated-anchor / autoplot-with-NA. Post-follow-up gate: `check()`
0/0/0 (vignettes rebuilt), 1791 pass / 0 fail / 0 skip, coverage 100%,
style/lint clean, `check_pkgdown()` clean. Detail in `MILESTONES.md`
(M46).

**Next up: nothing queued.** `ROADMAP.md` carries only unscheduled ideas
(AMH fidelity extension pending the owner’s Forbes outreach; e2 dual EFA
chi-squares; e4 bootstrap edge CIs;
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
ESEM/polychoric extensions). `MILESTONES.md` remains the source of truth
for *completed* milestones.

## Invariants — do not violate without flagging

These encode hard-won reasoning from the design phase. Changing them is
a design decision, not a refactor.

1.  **One edge path.** All between-level correlations go through
    [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md).
    Use the algebra (`W'RW`, standardized) when scoring is linear;
    materialize scores only when nonlinear (EAP) or when the user asks.
    **Always** standardize by real score SDs `sqrt(diag(W'RW))` — never
    assume unit variance (Bartlett/oblique scores are not
    unit-variance).
2.  **Keep the cross-check.** Retain the scores route even where algebra
    is the default, and keep the test asserting they agree within
    tolerance for linear engines.
3.  **Light core, heavy opt-in.** The object always carries
    loadings/variance/fit/weights/edges/ lineage/`R`/meta. `scores`, raw
    `fits`, raw `data` are NULL by default and recomputable.
4.  **Sign alignment anchors to the primary parent**, not “all positive”
    (that’s impossible).
5.  **Lineage lives in edges, never in IDs.** `m{k}f{j}` are stable
    labels; parentage is in the edge structure.
6.  **Loud defaults.** Announce consequential auto-choices via cli
    (e.g., the ordinal-detection warning). Advise loudly; never switch
    basis silently.
7.  **Convergence is data, not an error.** A non-converging level
    warns + is skipped; the object still builds to the deepest converged
    level. Never let one bad level abort the run.

## Resolved defaults (see `DESIGN.md` §9, §14)

**Varimax** (orthogonal) rotation — hardcoded internal constant since
M13; no `rotation` argument; oblique rotation **out of scope** (it
confounds the cross-level signal) · `cor = "pearson"` with
ordinal-detection warning · `tenBerge` scoring (on the active basis) ·
WLSMV estimator for ordinal ESEM · Forbes extension **off** · `k_max`
required · sign `align_signs = TRUE` · `keep_scores`/`keep_fits` stored
= `FALSE` · `redundancy_phi`: `NULL` (default) auto-resolves — `"pca"` →
no φ filter (component scores are determinate, so `|r|` is the true
between-component correlation); `"efa"`/`"esem"` → `0.95` (Lorenzo-Seva
& ten Berge 2006; factor-score indeterminacy makes `|r|`-only liberal).
`NA` is the explicit opt-out. Announce auto-resolve loudly (Invariant
6). Don’t change these silently.

## Dependencies (see `DESIGN.md` §12)

`psych` is in **Imports** (M21) — it is the engine substrate for the
default PCA and EFA paths and for polychoric correlations; placing it in
Suggests would require an install prompt for core functionality. The
SEM + plotting + optional-criterion stacks remain in `Suggests`.
`GPArotation` was **removed entirely** (M21) — varimax routes through
base [`stats::varimax`](https://rdrr.io/r/stats/varimax.html) and
GPArotation never enters
[`loadedNamespaces()`](https://rdrr.io/r/base/ns-load.html) on any
supported path. **Do not add further to `Imports` without flagging it.**
**No Rcpp** — profile first; the heavy compute already lives in compiled
deps (§3).

Current `Imports`: `cli`, `generics`, `psych`, `rlang`, `stats`,
`utils`. Current `Suggests`: `covr`, `EFAtools`, `future`,
`future.apply`, `ggplot2`, `gt`, `knitr`, `lavaan (>= 0.6-13)`,
`rmarkdown`, `testthat (>= 3.0.0)`. (`future` itself is declared because
the parallel test calls
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
directly; `future.apply` is the dispatch backend.)
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
uses `psych::fa.parallel(fa="both")` +
[`psych::vss`](https://rdrr.io/pkg/psych/man/VSS.html) (PA-PC, PA-FA,
MAP, VSS-1/2) and optionally
[`EFAtools::CD()`](https://rdrr.io/pkg/EFAtools/man/CD.html) (gated by
[`rlang::is_installed()`](https://rlang.r-lib.org/reference/is_installed.html));
no separate `EGAnet`/`paran` dep. `future.apply` (M26) is the optional
parallel backend for the ESEM per-level fits, gated by
[`rlang::is_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
with a serial `lapply` fallback — users opt in via
[`future::plan()`](https://future.futureverse.org/reference/plan.html);
no `ncores` arg, no `future`/`parallel` in Imports. Visualization uses
`ggplot2` directly (no `ggraph`/`igraph`/`tidygraph`). `methods` is
**not** imported (no `methods::` usage). `clue` was removed in M5.

## Dev workflow

R \>= 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard
devtools loop:

``` r

devtools::load_all()      # load for interactive testing
devtools::document()      # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()          # run testthat suite
devtools::check()         # full R CMD check
styler::style_pkg()       # format
lintr::lint_package()     # lint
```

**Efficiency (the suite takes minutes — don’t re-run it needlessly).**
`check()` already runs the full test suite *and* examples *and* rebuilds
vignettes (~3 min: ~90s tests + ~74s vignettes + ~20s examples), and
[`covr::package_coverage()`](http://covr.r-lib.org/reference/package_coverage.md)
runs the suite *again* — so `test()` → `check()` → `coverage()` at one
gate executes the suite ~3×. Instead: iterate with **targeted**
`devtools::test(filter = "<x>")` /
[`testthat::test_file()`](https://testthat.r-lib.org/reference/test_file.html);
run failing tests **once** in a way that shows the details (capture
`res <-` or use a non-silent reporter — never silent-then-rerun); skip
the vignette rebuild during mid-work checks with
`check(vignettes = FALSE)`; and at the final gate run `check()` **once**
(it subsumes `test()`) then `coverage()` once. Never run two
package-touching R processes concurrently.

Scaffolding helpers: `usethis::use_r()`, `use_test()`, `use_package()`,
`use_testthat(3)`, `use_github_action("check-standard")`. Use testthat
3e, roxygen2 for all exported functions (document the *why* of each
default, runnable `@examples`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- **New/removed exported object → pkgdown reference updated.** Whenever
  `NAMESPACE` gains or loses an `export()`/`S3method()`, or a dataset is
  added under `data/`, update the `reference:` list in `_pkgdown.yml` to
  match, then verify with
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html).
  This is the exact check the pkgdown GitHub Action runs — an exported
  topic missing from the reference index fails that workflow (and only
  that workflow; local `R CMD check` won’t catch it). Run it at the DoD
  gate unconditionally (it’s sub-second and gated on
  `rlang::is_installed("pkgdown")`); it costs nothing and closes a
  recurring CI-break.
- `devtools::check()` clean (run **once** at the gate — it subsumes
  `devtools::test()`; see the efficiency note under *Dev workflow*).
  Coverage checked once, not per sub-step.
- Styled and linted.
- Public-facing change reflected in NEWS.md and (if user-visible) the
  relevant `@examples`/vignette.
- For a milestone: a detailed entry added to `MILESTONES.md` **in
  numeric order** + a one-line index entry here under “Completed
  milestones”. `MILESTONES.md` is the single source of truth for
  milestone history — never re-log it in DESIGN.md §15 (a pointer) or
  duplicate it across files.

## Git

- Default branch is **`master`**; it stays green and releasable.
  Milestone work happens on a feature branch (`m{N}-<slug>`) and merges
  to `master` via a **PR**, **squash-merged as soon as the local
  definition of done is green** (`devtools::check()` 0/0/0 + tests +
  coverage + style/lint) — *not* gated on remote CI. `master` is
  deliberately **not branch-protected** (owner decision, 2026-07-01):
  required checks would gate every merge on the ~8–15 min CI matrix,
  pure latency for a solo pre-CRAN repo where local `check()` already
  ran. CI still runs on every push as an after-the-fact signal but does
  not block the merge; **don’t** use `gh pr merge --auto` (silently
  no-ops without required checks) and **don’t** synchronously watch or
  background-poll
  101. Don’t commit milestone work (anything touching `R/`, `tests/`,
       `DESCRIPTION`, vignettes) straight to `master`. Trivial isolated
       doc fixes may go direct at the user’s discretion — but **push
       them immediately.** An unpushed commit on `master` leaves local
       ahead of `origin/master`; `git pull` won’t surface it (nothing to
       fetch), and the next milestone branch — cut from local `master` —
       will bundle it into that milestone’s squash-merge and force a
       post-merge divergence. /plan-milestone step 8a now guards against
       branching from an ahead-of-origin `master`.

  - **Exception — release / CRAN-submission milestones**
    (e.g. CRAN-prep, version-bump/release): here you **do** wait for the
    *full* green CI matrix before merging, because CRAN runs exactly
    that matrix (macOS/Windows/Ubuntu × release/devel/oldrel) and will
    reject platform failures the local macOS `check()` can’t see. For
    these, don’t merge on local-green alone. This exception also
    reactivates once the package has real users or collaborators (a red
    `master` then blocks others / ships bugs) — treat that transition as
    the trigger to reconsider branch protection.
- **Do not touch** the `legacy` branch or the `v0-legacy` tag — they
  preserve the pre-AI code.
- Small, focused commits with imperative messages (e.g.,
  `Add PCA engine and level contract`).
- Don’t force-push `master`. Don’t commit data, credentials, or large
  binaries.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don’t invent a design decision.
- Adding an `Imports` dependency, introducing Rcpp, or changing a
  resolved default → flag for approval first.
- Touching `git history`, `legacy`, tags, or anything destructive →
  confirm first.
- Prefer wrapping established engines (`psych`, `lavaan`, `GPArotation`)
  over reimplementing numerics.

## Out of scope for now

- **EAP scoring** for ordinal ESEM — declined (M28): EAP’s shrinkage
  attenuates the cross-level correlations that bass-ackwards exists to
  measure; tenBerge covers the common case. Seam preserved in
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  but implementing EAP is not planned.
- **Oblique rotation** — varimax is hardcoded; no `rotation` argument;
  oblique confounds the cross-level signal. No plans to add it.
- **Higher-order SEM / Schmid-Leiman** — out of scope per §2; `ackwards`
  is correlation-based, not SEM-based.
- **Bootstrap CIs on skip-level edges** — deferred from M5; still out of
  scope. Logged in `DESIGN.md` §14. (Structural artefact signals and
  φ-default for non-PCA redundancy were reactivated and completed in
  M25.)
