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
default output must reproduce Forbes's examples exactly. **This contract is test-backed** (M44):
`tests/testthat/test-forbes-fidelity.R` reproduces the paper's three simulation studies against
expected values computed with Forbes's own reference implementation (OSF `pcwm8`; provenance in
the fixture), and the M44 feasibility study additionally verified the 155-variable AMH applied
example to 3.9e-14 (not committed as a test — the AMH matrix carries no OSF license; see
`ROADMAP.md` unscheduled items).

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
- **M38** — `missing = "fiml"` for PCA/EFA (code): `engine = "pca"/"efa"` + `cor = "pearson"` now routes `R` through `psych::corFiml()` (MAR-valid FIML) into the `W'RW` algebra (Invariant-1-clean, one corFiml call/run, no new dep), reversing the M16 "FIML errors for PCA/EFA" default; `.resolve_missing()` gains a `cor` guard (errors for non-Pearson PCA/EFA basis + WLSMV/ULSMV); `n_obs` string `"total"` (default, Enders 2010) / `"complete"` selects the approximate EFA fit-index N (Zhang & Savalei 2020; point estimates unaffected; announced via cli); `"effective"` dropped (no canonical formula); DESIGN §9/§14 (items 32–33) sign-off.
- **M39** — narrative & remaining prose (doc-only; final milestone of the M31–M39 epic): clarity pass across intro/suggest_k/ordinal/forbes/README — `print(sk)`/`print(x)`; orthogonal-varimax≡CF(1/p) rationale; SE/CI-NA-under-PCA note; `top_items(cut = 0.5)`; dropped `keep_scores` demo; **removed the README+intro red-arrow explanation** (no red arrow exists after the M35 sign fix); de-indexed README score columns; dropped the Waller citation nudge; suggest_k `sim16`-idealized vs `bfi25`-realistic contrast (inline-computed); ordinal binary/tetrachoric + WLSMV-polychoric + scores-trustworthy clarifications; forbes verbatim-heading rewrites + visible `redundancy_r` chunk + structural-table gt highlight + Lorenzo-Seva ref; refs alphabetized. Three code/viz asks (ordinal `categorical` flag, ordinal corr-comparison viz, forbes pruned-level label styling) spun off to **M40** (`ROADMAP.md`, DESIGN §14). No `R/`/NAMESPACE/export change.
- **M40** — deferred code/viz asks (final M31–M40 milestone): ordinal `categorical` flag **declined** (redundant — `cor = "polychoric"` already auto-selects WLSMV; would only add a conflict surface + §9 change for zero capability; discoverability handled in docs); ordinal corr-comparison now a **dodged bar chart** (ten `N1`–`N5` item pairs, `fill = basis`, hidden reshape code) replacing the two raw `round(x$r)` matrices (also fixed stale N1–N2 figures 0.73→0.79); `autoplot()` **italicises a fully-pruned level's** axis label (new `.fully_pruned_levels()`, `fontface` aesthetic through `.ba_level_labels()`, both directions) — partially-pruned levels stay plain. No new/removed export, no signature or dependency change.
- **M41** — independent Fable review (review-only): statistical core verified clean numerically (tenBerge/W′RW/signs/Forbes/ESEM-fit/suggest_k); findings — 1 Critical (EFA chi/p-value pairing), 6 Major (drop_pruned adjacent-pairs M34 regression, pre-M38 engines-vignette FIML prose, CD-mechanism misstatement, Forbes artifact-zero framing, `cut_strong` remnant, untested Forbes-fidelity contract), 11 Minor, 4 enhancements; full §9/§14 defaults audit (all sound; one sound-but-misjustified wording); report + M42/M43/M44 triage in `ROADMAP.md`. No code change.
- **M42** — review fixes, code: EFA `chi` now the likelihood-ratio `STATISTIC` matching `p_value`/RMSEA/TLI (C1); `.drop_pruned_nodes()` recomputes all-pairs edges fresh, fixing the M34 `pairs = "adjacent"` regression (M1); `print.suggest_k` "undetermined" consensus (m1) + PA-cap announcements (m2); `cut_show`/`n_iter` validation (m8); ordinal warning names flagged columns (e3, `detect_ordinal()` returns names); dead `esem_levels(n_obs)` removed (m7); stale comments fixed (m6); EFA-aware fit-plot caption (m10). No export/signature/dependency change.
- **M43** — review fixes, docs (doc-only): engines vignette rewritten around first-class `missing = "fiml"` for PCA/EFA + `fm = "pca"`/sample-size fixes (M2, m4, m5); suggest-k CD mechanism corrected + worked-BFI prose inline-computed/drift-proof (M3, m11); Forbes artifact section rewritten to report-and-judge with a top-|φ| table, `cut_strong` remnant removed, chain-retention example corrected (M4, M5, m3); `redundancy_phi` PCA rationale corrected to score *determinacy* across DESIGN §9/CLAUDE/roxygen/vignette (e1); sim16 doc comments modernized (m9). No behavior/export change.
- **M44** — Forbes-fidelity fixture (closes the M41→M44 review arc): found the paper's own OSF project (`pcwm8`: simulations script, reference implementation, AMH matrix); head-to-head vs her own functions matched edges to 3.9e-14 incl. the full 155-variable AMH example; shipped a 3.7 KB license-clean fixture (three seed-regenerated simulations + her implementation's expected outputs) and `test-forbes-fidelity.R` (65 assertions: edges/φ/chase paths/retention); contract annotated **test-backed**; AMH commit deferred (no OSF license — options logged for owner outreach). No export/dependency change.
- **M45** — out-of-sample scoring (train/test): fit-time item moments stored (`meta$item_means`/`item_sds`; NULL for cor-matrix input); `augment()` gains `scaling = c("fit", "sample")` with **`"fit"` default** (training moments — one metric across train/test/subsets; `"sample"` = pre-M45 opt-in and the cor-matrix route; DESIGN §14 item 34); new exported **`predict.ackwards(object, newdata, scaling)`** ≡ `augment(append = FALSE)` (equivalence test-asserted; `_pkgdown.yml` updated); intro-vignette train/test subsection. New export, no new dependency.
- **M46** — Girard extension (replicability-gated hierarchies): new exported **`comparability()`** — Everett (1983)/Goldberg (1990) split-half factor comparability per level per factor (full-sample-anchored labels, greedy-with-removal matching, cross-solution correlations through `compute_edges()` on the pooled R, Tucker's φ alongside, report-first/nothing auto-flagged; PCA/EFA + pearson/spearman, `n_splits = 10`, seeded) + `print`/`autoplot` methods; capstone vignette `ackwards-girard` ("Replicability-Gated Hierarchies: A Recommended Workflow") with the six-step workflow + common-mistakes section; completes the triad `suggest_k()` (range) · `comparability()` (floor) · `prune()` (differentiation); DESIGN §14 item 35; ESEM/polychoric extensions deferred to `ROADMAP.md`. In passing: `cli::symbol$phi` glyph fix in `prune()`. New export, no new dependency.
- **M47** — bootstrap edge CIs: new exported **`boot_edges(x, data, n_boot = 1000, conf = 0.95, seed)`** — nonparametric bootstrap SEs + percentile CIs on every between-level edge (resurrects the §14 e4 deferral), a standalone pipeable verb re-supplying `data` (Invariant 3). Per replicate: resample rows → recompute R (fit's cor/missing routine) → refit → **anchor** each level to the full-sample solution (M46 matching + sign orientation) → edges via `compute_edges()` (Invariant 1). Upfront seeded indices → **serial ≡ parallel** (`future.apply`, M26); failed replicates dropped + counted (`n_boot_ok`, Invariant 7). `tidy(what = "edges")` gains `se`/`lo`/`hi`/`n_boot_ok`; `print`/`summary` note coverage; `meta$fm` stored for EFA refits. PCA/EFA + pearson/spearman only (ESEM/polychoric deferred to `ROADMAP.md`; cor-matrix/esem/polychoric objects error). DESIGN §14 item 36 (amends the "reuse `loadings_se`" phrasing; percentile-CI + Fisher-z-oracle rationale). New export, no new dependency.
- **M48** — performance & workflow pass (meta/process; no package-code change): suite-wide `cached()` test-fit memo in `helper-data.R` + parallel testthat (suite 93.6s→81.2s serial→**26.9s** at `TESTTHAT_CPUS=8`; full check 319.8s→241s; 1859 assertions unchanged; boot reproducibility/serial≡parallel oracles and MC sizes deliberately untouched); transcript-mined workflow audit (30k messages, 45 sessions: 249 full-suite runs, ~600 cd-compounds, 308 bare `load_all`s) → **`tools/dod-gate.R`** one-command DoD gate (dogfooded) + CLAUDE.md/implement-milestone cadence text in lockstep; expanded read-only permission allowlist staged for owner review (auto-mode classifier correctly blocked agent self-widening). Report-only: forbes vignette (23.9s) dominates rebuilds.
- **M49** — Initial CRAN release (0.1.0) + a robustness arc from real-data testing: **Phase A** roadmap cleanup (e2 declined §14.40; `label_items()`/third-dataset declines + factor-label 0.2.0 banked §14.41); **Phase B** doc/pedagogy pass (README tighten + labels showcase; intro reframed "basic toolkit"; `top_items()` label format → **`code: label`**; pkgdown regroup; girard Pearson-basis; engines ESEM reframe; interpret real-labels + `labelled::var_label()`; vocab); **robustness** (owner hit a 142-item/n=222 clinical-scale polychoric failure) — new **`correct`** arg (psych continuity correction; `correct = 0` fix + actionable errors §14.42), new exported **`check_items()`** (+ `print`/`[`; internal screen errors on constant, warns on near-constant §14.43), muffled psych's flood + a single **near-singular** warning + durable `meta$near_singular`/`min_eigenvalue` re-surfaced by `print`/`summary` + a **"When to trust the result"** `?ackwards` section (§14.44); **Phase C** release mechanics (NEWS → dated `0.1.0`; cran-comments refresh; README CRAN line + lifecycle→**stable**; DESCRIPTION/CITATION/DOIs). New exports: `check_items` (+ methods), `correct` arg, `meta` conditioning fields. Full-CI-matrix merge; `v0.1.0` retagged. No new dependency.
- **M50** — release polish (code; interleaved before M49's CRAN mechanics): `bfi25` ships 25 public-domain IPIP variable labels (Goldberg 1999; `data-raw/bfi25.R`) → `top_items()` prints `label (code)` free (plain attrs; row-subsetting drops them, so fit direct with `missing = "listwise"`); cli consistency (engine name lowercase everywhere incl. `print.comparability`; `summary()` fixed-decimal percentages + blank-line-separated levels; `top_items()` group spacing; single consolidated pruned footer in `print`/`summary`); roxygen examples split by intent (mechanics → continuous `sim16`, content → `bfi25`+polychoric); `suggest_k()` gains the ordinal-detection warning (Invariant-6 symmetry, screening-context wording pointing at the final `ackwards()` fit). `label_items()` setter + third dataset declined; factor-label pipeline deferred to 0.2.0. No new export, no dependency change.
- **M51** — factor-label pipeline (first 0.2.0-cycle milestone; DESCRIPTION → `0.1.0.9000`): new exported **`set_factor_labels(x, labels)`** (pipeable verb; merges/updates, `NULL` clears, `NA`/`""` removes one, errors on unknown ID) + **`factor_labels(x)`** getter, storing names in `meta$factor_labels` (rides through `prune`/`boot_edges`/`augment`/`predict` for free; DESIGN §14 item 45). Display form `label (id)` in `summary()` (per-level + lineage), `print()` (coverage line), `top_items(by="factor")` headers; `autoplot()` shows the label only (call-time `node_labels` overrides per node); `tidy()` gains `factor_label`/`from_label`/`to_label` **only when set** (ID columns untouched, Invariant 5). Purely additive; interpret + visualization vignettes updated. No new dependency.

## Current focus

**M51 is complete** — Factor-label pipeline (first milestone of the 0.2.0 development cycle;
DESCRIPTION bumped to `0.1.0.9000`). Purely additive per DESIGN §14 item 45: persistent *factor*
labels, kept lexically distinct from M50's *item* labels (`meta$item_labels`). New exports
**`set_factor_labels(x, labels)`** (pipeable verb; merges/updates, `NULL` clears, `NA`/`""` removes
one, errors on unknown ID) + **`factor_labels(x)`** getter, storing names in `meta$factor_labels`
(rides through `prune`/`boot_edges`/`augment`/`predict` for free). Display form `label (id)` in
`summary()` (per-level listing + lineage tree), `print()` (coverage line), and
`top_items(by="factor")` headers; `autoplot()` shows the label only with call-time `node_labels`
overriding per node; `tidy()` adds `factor_label`/`from_label`/`to_label` **only when set** (ID
columns never mutated, Invariant 5). Interpret + visualization vignettes updated; `test-factor-labels.R`
(8 tests). Suite **1978 pass / 0 fail / 0 skip**; coverage **100%**; `R CMD check` **0/0/0**;
styler/lintr/pkgdown clean. Detail in `MILESTONES.md` (M51). **Next up: nothing queued** — see the
remaining 0.2.0 candidates below.

**Owner-only release tail (from M49, still pending):** win-builder/R-hub remote checks and the
interactive `devtools::submit_cran()` email the maintainer. Until CRAN accepts 0.1.0, the
`install.packages("ackwards")` line in the README points at a package not yet on CRAN — update the
release date / any "on CRAN" phrasing when it lands. If CRAN bounces 0.1.0, a patch branches from
the `v0.1.0` tag (the `0.1.0.9000` dev bump on `master` does not block that).

**Remaining 0.2.0 candidates** (`ROADMAP.md`): demand-gated `comparability()`/`boot_edges()`
ESEM/polychoric extensions and the AMH-fidelity item pending the owner's Forbes outreach.
`MILESTONES.md` remains the source of truth for completed milestones.

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
`FALSE` · `redundancy_phi`: `NULL` (default) auto-resolves — `"pca"` → no φ filter (component
scores are determinate, so `|r|` is the true between-component correlation); `"efa"`/`"esem"` →
`0.95` (Lorenzo-Seva & ten Berge 2006; factor-score indeterminacy makes `|r|`-only liberal).
`NA` is the explicit opt-out. Announce auto-resolve loudly (Invariant 6). Don't change these
silently.

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

**Efficiency (don't re-run the suite needlessly).** The suite is parallel since M48
(`Config/testthat/parallel`, slowest files first): ~27s with `TESTTHAT_CPUS=8`, ~81s serial —
**always prefix suite runs with `TESTTHAT_CPUS=8`** (testthat defaults to 2 workers without it).
`check()` still runs the full suite *and* examples *and* rebuilds vignettes (~74s of vignettes —
`ackwards-forbes.Rmd` alone is ~24s), and `covr::package_coverage()` runs the suite *again* — so
`test()` → `check()` → `coverage()` at one gate executes the suite ~3×. Instead: iterate with
**targeted** `devtools::test(filter = "<x>")` / `testthat::test_file()`; run failing tests **once**
in a way that shows the details (capture `res <-` or use a non-silent reporter — never
silent-then-rerun); skip the vignette rebuild during mid-work checks with
`check(vignettes = FALSE)`; and at the final gate run `Rscript tools/dod-gate.R` — it executes the
whole DoD sequence (check → coverage → style → lint → pkgdown) serially in one process with
sensible `TESTTHAT_CPUS`, and prints/exits on any failure. Never run two package-touching R
processes concurrently. Two transcript-mined anti-patterns to avoid (M48): a bare
`devtools::load_all()` in its own `Rscript` call does nothing persistent (each `Rscript` is a fresh
process; `test()`/`check()` load the package themselves), and `cd <repo> && …` compound commands
trigger avoidable permission prompts — use absolute paths. In tests, reuse the `cached()` fit memo
(`tests/testthat/helper-data.R`) instead of refitting identical `ackwards()` objects — but never
for reproducibility/serial-vs-parallel oracles (a cached second call asserts nothing) or fits
wrapped in condition expectations; and treat cached objects as read-only — rebinding a returned
value is copy-on-modify safe, but environment-bearing components (e.g. lavaan fits under
`keep_fits = TRUE`) are shared by reference across cache hits.

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
- The whole gate above is one command: `Rscript tools/dod-gate.R` (M48) — check → coverage →
  style → lint → pkgdown, serial, one process, non-zero exit on any failure.
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
- ~~**Bootstrap CIs on skip-level edges**~~ — **shipped M47** as the `boot_edges()` verb (all
  edges, adjacent and skip-level; PCA/EFA + pearson/spearman). Deferred from M5, reactivated from
  DESIGN §14 e4. (Structural artefact signals and φ-default for non-PCA redundancy were
  reactivated and completed in M25.)
