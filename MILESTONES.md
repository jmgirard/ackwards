# Milestone history — `ackwards`

This is the **single detailed log** of completed milestones. `CLAUDE.md` carries only a
compact one-line index (and the current focus); `DESIGN.md` §15 points here. Keep this file the
source of truth for milestone detail — add each new milestone entry here, in numeric order, as
part of the milestone's definition of done.

Forward-looking / deferred items are not here; they live in `DESIGN.md` §14 (decisions remaining)
and `CLAUDE.md`'s "Out of scope" list. User-facing change notes live in `NEWS.md`.

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
- **M10 (done):** Conformance + robustness — `summary.ackwards()` + `print.summary_ackwards()`
  (previously documented but unimplemented §10 method); ESEM Heywood/improper-solution warning
  (`theta ≤ 0`, parity with EFA engine); `cor="spearman"` + `method="esem"` inconsistency warning;
  DESIGN.md §8 reconciled to list only PA + MAP (EKC/EGA marked out of scope).
- **M11 (done):** Edge-label polish + `show_r` decoupling — APA `.format_r()` helper (strip leading
  zero, pad trailing zeros, suppress `-.00`), `geom_label` with perpendicular offset + white halo +
  `r_label_size` arg; **decouple `show_r` from `drop_pruned`** (default `FALSE` everywhere); Forbes
  vignette updated to two-figure (labeled + unlabeled) treatment; `.lintr` added to `.Rbuildignore`
  (R CMD check fully clean: 0 errors, 0 warnings, 0 notes).
- **M12 (done):** Best-practice `suggest_k` — PA-FA added alongside PA-PC (`psych::fa.parallel(fa =
  "both")`); VSS-1/VSS-2 surfaced from existing `psych::vss()` call; Comparison Data (CD) added via
  `EFAtools::CD()` gated by `rlang::is_installed()` (skips gracefully when absent); new `seed` arg;
  enriched `suggest_k` object (`k_parallel_pc`, `k_parallel_fa`, `k_vss1`, `k_vss2`, `k_cd`,
  `cd_available`, expanded `criteria` table); redesigned `print.suggest_k()` multi-criterion table;
  new `autoplot.suggest_k()` three-panel ggplot2 diagnostic (scree/PA + MAP + VSS); `EFAtools` added
  to Suggests; DESIGN.md §8 and §12 updated.
- **M13 (done):** Rotation honesty — removed dead `kappa` argument; removed `rotation` argument
  entirely (only varimax is valid; exposed as a user arg it implied quartimax/equamax were options);
  renamed "cfT" → "varimax" throughout all three engine internals, the result object, print output,
  README, and docs; `@section Defaults` explains T′=T⁻¹ → W′RW algebra + varimax = what all
  reference papers used; DESIGN.md §4, §9, §14.1, §14.7 updated.
- **M14 (done):** Dedicated `suggest_k()` vignette — new `vignette("ackwards-suggest-k")` ("Choosing
  k: How Many Factors?") with all five criteria (pros/cons, bias direction, engine pairing), argument
  coverage (cor/n_iter/seed/k_max including ordinal→Pearson and PA-non-reproducibility caveats), and
  a worked BFI recommendation; intro vignette Step 1 trimmed to default call + pointer; `_pkgdown.yml`
  lists the new article first under Deep dives; README stale two-criteria description corrected to
  five criteria; DESIGN.md §8 and §15 updated. Post-review fixes: worked-example prose corrected to
  match actual output (PA-PC=5/PA-FA=6/CD=8, CD-outlier explanation); CD table "Conservative" →
  "Accurate in simulation; can over-retain on large, correlated samples"; three new tests covering
  `autoplot.suggest_k()` and `print.suggest_k()` for the `k_parallel_fa=NA` and `cd_available=FALSE`
  branches. (719 tests pass, 1 skip.)
- **M15 (done):** Naming clarity & consistency pass — `k`→`k_max`, `method`→`engine`,
  `scores`→`keep_scores`, `align`→`align_signs` on `ackwards()`; `$method`→`$engine`,
  `$cor_type`→`$cor` on the result object; `method`→`edge_method` on `compute_edges()`. All S3
  methods, tests, 6 vignettes, README, NEWS, CLAUDE.md, and DESIGN.md updated. (724 tests pass,
  1 skip; 0/0/0 R CMD check.)
- **M16 (done):** Estimator-aware missing-data handling — new `missing = c("pairwise","listwise",
  "fiml")` argument on `ackwards()`. Default `"pairwise"` preserves existing behaviour and warns
  when NAs present; for ESEM WLSMV/ULSMV passes `available.cases` to lavaan (full N, honest).
  `"listwise"` reduces to complete cases pre-fit for consistent N. `"fiml"` (ESEM ML/MLR only)
  uses Full Information ML via lavaan and derives edge R from the FIML saturated model. Fixes ESEM
  ML/MLR fit-vs-edges inconsistency for `"listwise"` and `"fiml"` (including a silent bug where the
  h1 extraction fell back to pairwise). Adds `.resolve_missing()` helper; records
  `meta$missing`/`meta$n_complete`. 36 tests in `test-missing.R`; missing-data section added to
  `ackwards-engines.Rmd`; DESIGN.md §9 and §14 updated. (778 tests pass, 1 skip; 0/0/0 R CMD check.)
- **M17 (done):** GitHub 0.1.0 release prep — license switched CC BY 4.0 → MIT (`LICENSE` stub +
  `LICENSE.md` updated, DESCRIPTION updated); version bumped `0.0.0.9000 → 0.1.0`; README MIT badge
  + comment mismatch fixed; `inst/CITATION` added (Goldberg 2006 + package), version field dynamic
  via `meta[["Version"]]`, hand-written Goldberg prose removed from README to avoid duplication;
  README rebuilt with correct two-entry citation output; NEWS.md restructured to curated
  capability-grouped summary (development history dropped — pre-release, all captured in git);
  `_pkgdown.yml` 0.1.0 release URL registered; pkgdown rebuilt cleanly. Post-review: citation
  guard test added to `test-utils.R`. Note: pkgdown 2.2.0 renders all root `.md` files via
  `package_mds()` with no config-based exclusion — CLAUDE.md/DESIGN.md appear as unlinked pages;
  not fixable without moving files. (787 tests pass, 1 skip; 0/0/0 R CMD check; all URLs clean.)
  Tag: owner runs `git tag -a v0.1.0 -m "ackwards 0.1.0" && git push origin v0.1.0`.
- **M18 (done):** Factor interpretation & label scaffolding — `top_items()` (salient per-factor item
  listing, `|loading| >= cut`, grouped cli print) and `label_template()` (node_labels scaffold,
  styles: "id"/"forbes"/"blank", prints editable c(...) literal). Both are pure consumers of the
  existing light core; no new dependencies, no invariant or default changes. Intro vignette and
  visualization vignette updated; pkgdown reference index updated; DESIGN.md §10/§11/§15 updated.
  Post-review hardening: `label_template(style = "forbes")` now guards `k_max > 26` (LETTERS
  exhaustion) with a loud `cli_abort` and documents the constraint; the duplicate `ba_layout()`
  call in the forbes branch was cached; the `sort = FALSE` test asserts order against `tidy()`
  rather than set membership; EFA smoke tests added for both helpers (engine-agnosticism);
  stale `tests/testthat/_problems/` removed.
- **M19 (done):** Dedicated interpretation/labeling vignette — documentation-only. New
  `vignette("ackwards-interpret")` ("Interpreting and Labeling Factors") owns the `top_items()` →
  name → `label_template()` → `autoplot(node_labels=)` workflow plus hierarchy-aware naming (parent
  vs child, blends, reorganizing factors via lineage/edges) and the sign-alignment caveat. Listed in
  pkgdown Deep dives after `ackwards-suggest-k`. Intro Step 5 trimmed to a slim `top_items()` example
  + pointer; visualization vignette keeps the labeling mechanic + cross-ref (naming judgment moved to
  the new vignette). DESIGN.md §15 entry (amends §10/§11). Original milestone was
  documentation-only; post-review hardening added a guard test for the vignette's edges idiom,
  `top_items()`/`label_template()`→`autoplot()` idiom smoke tests, switched the interpret vignette
  to `cor = "polychoric"` (consistency + drops the ordinal warning), and **fixed swapped m5f3/m5f5
  labels** (Conscientiousness/Openness) in the interpret and visualization vignettes.
  (935 tests pass, 1 skip; 0/0/0 R CMD check.)
- **M20 (done):** CRAN submission readiness + example legibility. A *release-readiness* milestone
  (not a DESIGN.md §15 feature milestone) ahead of tagging/submitting 0.1.0. Five waves:
  (1) Statistical/correctness: `.standardize()` (na.rm-aware) replaces `scale()` in
  `.compute_scores()` and `compute_edges()` scores route; `detect_ordinal()` guarded against
  all-`NA` columns; stale "oblique rotations" wording in `compute_edges()` roxygen fixed.
  (2) Example conversions: all `\dontrun{}` removed; fast examples use `requireNamespace()` guards;
  slow/stochastic `suggest_k` blocks use `\donttest{}`.
  (3) Submission metadata: three DOIs added to `DESCRIPTION` (Goldberg 2006, Waller 2007, Forbes
  2023); `NEWS.md` restructured into a single `0.1.0` entry; `cran-comments.md` added.
  (4) Example legibility: `tidy.ackwards()` gains `sort = c("none","strength")` for edges;
  README/vignettes rewrote `order(-abs(...))` → `tidy(sort="strength")`, `identical(round(...))` →
  `all.equal()`, `grep("^\\.m5",...)` → `startsWith(names(...), ".m5")`, double-`rbind` pattern →
  intermediate variable.
  (5) Verify: styler + lintr clean; `R CMD check --as-cran` → 0/0/0; README rebuilt.
  Owner next steps: `devtools::check_win_devel()`, `rhub::rhub_check()` before actual CRAN upload.
  (947 tests pass, 1 skip; 0/0/0 R CMD check.)
- **M21 (done):** Onboarding & usability pass (pre-CRAN). Four parts:
  (A) `psych` Suggests→Imports (engine substrate for default PCA/EFA; never needed an install prompt
  for core use); `GPArotation` removed entirely (varimax routes through `stats::varimax`; verified
  never loaded). All `check_installed("psych")` guards removed; two vestigial `skip_if_not_installed
  ("GPArotation")` test lines deleted. DESIGN.md §3/§12 updated.
  (B) `bfi25` dataset bundled: 1 000 rows sampled from `psych::bfi[, 1:25]` (seed 42, NAs preserved)
  via `data-raw/bfi25.R`; documented in `R/data.R` with `@source` (Revelle/psych/SAPA/IPIP — items
  are public-domain IPIP). All `@examples` + 6 vignettes + README.Rmd migrated to `bfi25`; suggest-k
  worked-example prose regenerated (n=875: PA-PC=5, PA-FA=6, MAP=5, VSS-1=4, VSS-2=5, CD=6;
  consensus 4–6). Oracle snapshot stays on full `psych::bfi[, 1:25]`.
  (C) README "Learn more" table now lists all 7 vignettes (Visualization and Interpreting & labeling
  were missing).
  (D) `autoplot.suggest_k()` adds a 4th "CD (RMSE, minimize)" panel (2×2 grid when CD available;
  unchanged 3-panel single-column otherwise); `suggest_k` object gains `cd_rmse` field. CD vline
  removed from MAP panel. Tests for 4-panel/3-panel branches and `cd_rmse` field.
  Post-review: committed stale README hero figure (bfi25 render); added `test-data.R` (6 assertions
  on `bfi25` shape/cols/class/NAs/range); corrected `data.R` `@source` year 2025→2026.
  (974 tests pass, 1 skip; 0/0/0 R CMD check.)
- **M22 (done):** Correlation-matrix input (PCA/EFA-only). `ackwards()` and `suggest_k()` now
  accept a pre-computed correlation matrix detected automatically (square, symmetric, unit diagonal).
  New `n_obs` arg (required for EFA+R, optional for PCA+R). ESEM gated off (lavaan needs raw data).
  `cor`/`missing` args ignored+warned for R input; `$cor` stored as `NA`; print shows
  `"(user-supplied matrix)"`. `keep_scores=TRUE`/`augment()`/`tidy(what="scores")` all error
  clearly. CD gated off in `suggest_k()` with info note. Edges from R-matrix and raw-data paths are
  identical within floating-point tolerance (same W'RW algebra). `.is_cor_matrix()` +
  `.validate_cor_matrix()` + `.check_maybe_cov_matrix()` helpers in `utils.R`. `meta$input_type`
  field added. Non-ASCII chars replaced across all R files (0/0/0 clean). Post-review: non-PD
  warning tested; covariance-matrix detection added (targeted error via `.check_maybe_cov_matrix()`
  in both functions); `prune="redundant"` + cor_matrix test; `autoplot.suggest_k()` + cor_matrix
  test; invalid `n_obs` tests for `suggest_k()`; `missing(missing)` comment; NEWS folded into 0.1.0.
  44 tests in `test-cor-input.R`, 41 in `test-utils.R`. (1035 tests pass, 1 skip; 0/0/0 check.)
- **M23 (done):** Test-coverage hardening — raised `covr` from **93.37% → 100%** overall
  (every file at 100%). Strategy:
  test first (`# nocov` only for genuinely unreachable defensives). New tests: `suggest_k.R`
  cor-ignored/spearman+CD/CD-dash branches; `label_template.R` k>26 guard; `prune.R` `.tucker_phi`
  all-zeros / `.phi_pairs` adjacent / `print.ackwards` artefact+phi-note; `compute_edges.R`
  algebra explicit path; `engine_esem.R` NULL-SE skip in `tidy(what="loadings_se")`; `summary.R`
  EFA chi/dof row, empty-redundant "(none)", phi-note, empty-lineage. `# nocov` markers on
  unreachable engine defensives: PCA + EFA k=1 sign flip, EFA convergence-fail break + tenBerge
  fallback, ESEM lavaan version guard + convergence fail + std_sol NULL + tenBerge fallback +
  W NULL + warning muffler + null-fit + cov2cor unreachable branch + Pearson fallback +
  fitMeasures error handler + Phi shape fallback; `prune.R` empty-chains guard + null-nodes
  branch; `layout.R` top-down zero-weight + bottom-up orphaned-parent fallbacks; `summary.R`
  PCA eigenvalue miss + null-nodes branch. **No behavior change.** Local verification only
  (no CI coverage workflow). Post-review: rewrote the `print.top_items` empty-shown-level test
  to deterministically hit `interpret.R` `next` (was skipping on bfi25 -> covered the line +
  dropped the spurious 2nd skip); `# nocov`'d the EFA k=1 sign flip to match the PCA analogue;
  added a covr-limitation note explaining the ESEM warning-muffler `# nocov` (live code covr
  cannot instrument, not dead code). (1080 tests pass, 1 skip; 0/0/0 R CMD check.)
- **Test-suite speedup (post-M23, no source/behavior change):** `EFAtools::CD()` dominated test
  time (~8 s/call on the full 2800×25 `psych::bfi` vs ~0.1 s for `fa.parallel` + `vss`).
  `test-suggest_k.R` rewritten to run the CD-bound structural/logic/print/autoplot tests on small
  `bfi25` subsets (CD ~0.8 s, cached), keeping CD/PA/MAP/VSS *real* — only the data is smaller — plus
  one full-`psych::bfi` integration smoke for the `n_obs=2800`/`n_vars=25` assertion. PA-cap test
  moved to `bfi25[, 1:15]`/`k_max=2` (genuine clamp). `devtools::test()` ~80 s → ~46 s; coverage
  held at 100%. Also fixed three leaked secondary warnings so the suite is warning-free: the
  spearman/CD `n_factors_max` notice (`test-suggest_k.R`) and the NA-pairwise + ordinal notices in
  two `test-cor-input.R` `n_obs`-ignored tests (switched to continuous/complete data and within-CD
  `k_max`); added a deterministic default-`k_max` test to preserve the `suggest_k.R:208` branch.
  Comfortably within CRAN's check-time budget. (1083 tests pass, 1 skip; 100% coverage.)
- **M24 (done):** Vignette communication pass (documentation-only). Stacked long-format `kable`
  comparison tables replaced by wide **`gt`** tables (one row per item/edge, one column per
  engine/basis, explicit Δ column) in `ackwards-engines` (loadings: `|EFA|−|PCA|`; edges:
  `EFA−PCA`) and `ackwards-ordinal` (loadings: `poly−pearson`; edges: `poly−pearson`). Factor/sign
  alignment asserted via `stopifnot()` before differencing; engine disagreements on primary parents
  surface as `NA`. `ackwards-forbes`: `prune-nodes` raw `tidy()` print replaced by styled gt table;
  narrated counts and stale skip-edge claim (m3f2→m5f1, now m3f2→m5f2 on bfi25) converted to inline
  `` `r` `` expressions. `skip-edges`/`thresholds` tables migrated to gt for visual consistency.
  All presentation code in `echo = FALSE` — no plumbing visible to readers. Audit of the other four
  vignettes: no changes needed (raw console prints are pedagogically appropriate). `gt` added to
  Suggests. Post-review (six findings): unified every Δ column to a **magnitude** convention
  (`|x|−|y|`) so the directional captions read correctly for negatively-signed cells — fixes a
  sign-trap where the polychoric-strengthened negative edge `m2f2→m3f2` had shown a negative signed
  Δ (ordinal edges + loadings and engines edges switched from signed to `abs`; engines loadings was
  already `abs`); derived the forbes redundant-level prose (`flagged_ids` + per-level `lvl_summary`)
  inline instead of hard-coding "entire k=4 level / two factors at k=2,3"; converted the
  `prune-artefact` raw print to a styled gt table (call still echoed) with inline `n_artefact`; added
  guard tests for the `knitr::kable` fallback branch, the NA primary-parent-disagreement merge, the
  magnitude-delta sign property, and the inline-derivation helpers. (1123 tests pass, 1 skip;
  0/0/0 R CMD check.)
- **M25 (done):** Deferred-items pass — three waves. (1) `suggest_k()` gains `criteria` arg
  (`rlang::arg_match(multiple=TRUE)`); any subset of five criteria may be requested; non-requested
  `k_*` fields → `NA`; shared computation (`fa.parallel(fa="both")` / `vss()` called at most once);
  `print()` / `autoplot()` render only requested criteria; consensus from requested only. (2) `prune="artefact"`
  now populates `x$prune$structural` — per-factor `few_items` / `orphan` / `split_merge` signals
  (Forbes Fig 2); new `min_items = 3L` and `orphan_r = 0.5` args; flag/report only, never auto-prune;
  `print()`/`summary()` report flagged count. (3) `redundancy_phi = NULL` auto-resolves: PCA →
  no φ (exact algebra); EFA/ESEM → `0.95` (Lorenzo-Seva & ten Berge 2006; announced via cli,
  Invariant 6); `NA` is explicit opt-out. Bootstrap CIs on edges remain deferred. `CLAUDE.md`
  Resolved defaults + Out of scope updated; DESIGN §9/§14 updated. Post-review hardening:
  added loud validation for `min_items` (positive integer) / `orphan_r` (`[0,1]`); roxygen *why*
  for both defaults (three-indicator rule; moderate-correlation rationale); deterministic unit
  tests for the previously-untested `split_merge = TRUE` path (mock loadings) and the Wave-3
  φ-decision outcome change (`.find_redundant_chains` flags under no-φ but not under φ>0.95) plus
  an end-to-end EFA subset-property test (auto-φ flagged set ⊆ |r|-only set); documented `criteria=`
  in `ackwards-suggest-k.Rmd`, structural signals + φ auto-default in `ackwards-forbes.Rmd`,
  `criteria=` in DESIGN §8; simplified an unreachable `n_struct` branch and `# nocov`'d the
  EFAtools-absent CD note (covered only when EFAtools missing). Coverage back to 100%.
  (1219 tests pass, 2 skip; 0/0/0 R CMD check.)

- **M26 (done):** ESEM performance for large item sets — two complementary speedups, no
  behaviour change. (1) **Cached sample statistics:** the ESEM engine no longer recomputes
  lavaan's data-derived sample statistics (thresholds, polychoric matrix, asymptotic weight
  matrix NACOV/WLS.V) at every level — they depend only on the data, not `nfactors`, so they are
  harvested once at the anchor level (k=1) via `fit@SampleStats` and reused for every deeper level
  through lavaan's `slotSampleStats=` argument. Verified **bit-identical** (0.00e+00 loading/edge
  diff) vs. the from-raw path. Dominant saving at large p (polychoric+NACOV recompute was O(p²)+ ×
  k_max). (2) **Parallel per-level fits:** `esem_levels()` refactored into a slim per-level worker
  `.esem_fit_one()` (returns the level contract, not the heavy fit, to avoid serialising a
  duplicate NACOV) dispatched through `.esem_lapply()` — `future.apply::future_lapply` when
  installed (gated by `rlang::is_installed()`), serial `lapply` fallback otherwise. No `ncores`
  arg: users set `future::plan()` (sequential default → no behaviour change; `multisession`/
  `multicore` to parallelize). `future.seed = TRUE` (lavaan::efa is mildly RNG-stochastic, ~1e-6);
  reproducible across plans when `seed` supplied (verified 0.00e+00). Invariant 7 preserved: all
  levels fit, then assembly truncates at the first non-converged/failed level and emits all cli
  warnings in deterministic level order (workers never signal conditions). `future.apply` added to
  **Suggests** (flagged + approved). New `@section Performance` on `ackwards()`; "Performance with
  many items" section in `ackwards-engines.Rmd` (incl. EFA+polychoric cheaper-route pointer); NEWS +
  DESIGN §12/§14 updated. Three new tests (both `.esem_lapply` branches via `local_mocked_bindings`;
  serial-vs-`multicore` identity, skipped on Windows / when future absent). Coverage held at 100%.
  (1225 tests pass, 2 skip; 0/0/0 R CMD check.)
- **M27 (done):** ESEM fit & SEs as first-class output — CRAN-readiness pass making per-level fit
  indices and rotation-aware loading SEs usable and correctly framed. Additive surfacing + docs only;
  no invariant or resolved-default change; no new `Imports`; version stays `0.1.0`.
  (1) **`glance()` carries fit**: deepest-converged-level `CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`
  appended as a consistent five-column set across all engines (NA where unavailable: CFI/SRMR for
  EFA, all five for PCA). `.glance_fit()` internal helper.
  (2) **`tidy(what="fit")` gains `format` and `cutoffs`**: `format="wide"` returns one row per
  **non-anchor** level (k >= 2; the saturated 1-factor anchor is dropped to match `summary()` and
  `autoplot(what="fit")`) with index columns (long default byte-identical to previous output);
  `cutoffs=TRUE` appends a `meets` flag against Hu & Bentler (1999) thresholds (CFI/TLI >= .95,
  RMSEA <= .06, SRMR <= .08) via `.fit_cutoffs()` / `.flag_fit()` / `.fit_long_to_wide()` helpers;
  report-only, never gates.
  (3) **Loading CIs folded into `tidy(what="loadings")`**: added `se`, `ci_lower`, `ci_upper`
  columns for all engines (NA for PCA/EFA, populated for ESEM); `conf_level=0.95` argument controls
  width. The now-redundant `tidy(what="loadings_se")` accessor was **removed** (package unreleased,
  no deprecation needed); the internal `lev$loadings_se` level-contract field is unchanged.
  (4) **`autoplot(x, what="fit")`**: two-panel ggplot2 line chart (CFI/TLI top panel; RMSEA/SRMR
  bottom panel; Hu & Bentler reference lines). `what="hierarchy"` (default) unchanged. PCA returns
  informative empty plot.
  (5) **`summary()` fit lines** annotated with `✔`/`✘` pass/fail marks for thresholded
  indices.
  (6) **Honesty framing in `print()` and `summary()` footer**: per-level fit describes the k-factor
  solution at that level; does not validate edges or the hierarchy.
  (7) **Dedicated section "Per-level fit: what it tells you (and what it doesn't)"** in
  `vignette("ackwards-engines")`: level-vs-hierarchy distinction, reporting workflow (wide table +
  cutoffs + fit plot), ESEM cost/benefit, when to care (exploratory vs publication).
  `vignette("ackwards-ordinal")` updated to use `format="wide"` and cross-reference.
  No new `Imports`; no design invariant changed; version stays `0.1.0`.
  Post-review fixes (from /post-milestone-review): removed the redundant
  `tidy(what="loadings_se")` accessor (per the approved plan; the implementation had kept it);
  wide fit table now drops the anchor level for consistency with `summary()`/the fit plot;
  fit-plot panel titles made engine-aware (EFA shows "TLI"/"RMSEA", not "CFI / TLI"); removed dead
  `excl` variable and an unreachable `is.null(cut)` guard in `.ba_fit_plot()`; added tests for the
  `.glance_fit` truncation/empty-fit branches and the empty fit-plot branch, restoring coverage to
  100%. (1285 tests pass, 2 skip; 0/0/0 R CMD check; coverage 100%.)
