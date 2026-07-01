# Milestone history — ackwards

This is the **single detailed log** of completed milestones. `CLAUDE.md`
carries only a compact one-line index (and the current focus);
`DESIGN.md` §15 points here. Keep this file the source of truth for
milestone detail — add each new milestone entry here, in numeric order,
as part of the milestone’s definition of done.

Forward-looking / deferred items are not here; they live in `DESIGN.md`
§14 (decisions remaining) and `CLAUDE.md`’s “Out of scope” list.
User-facing change notes live in `NEWS.md`.

- **M1 (done):** PCA engine + algebra
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md) +
  result object + `print`/`tidy`/`glance`, validated against
  [`psych::bassAckward()`](https://rdrr.io/pkg/psych/man/bassAckward.html).

- **M2 (done):**
  [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md) +
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  (adjacent-level diagram) +
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md).

- **M3 (done):** EFA engine
  ([`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html), tenBerge
  weights, algebra path) + materialized-scores route + algebra-vs-scores
  cross-check tests.

- **M4 (done):** ESEM engine
  ([`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html), WLSMV,
  tenBerge weights, rotation-aware SEs, per-level fit indices) +
  `cor = "polychoric"` for all engines + `loadings_se` in level contract

  - convergence truncation tested + `estimator` argument.

- **M5 (done):** Forbes extension — `pairs = "all"`,
  `prune = "redundant"/"artefact"`, Tucker’s φ chains (DFS enumeration,
  global retain set), annotated
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  (skip-level arcs, pruned fill), `tidy(what = "nodes")`,
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
  print caveat.

- **M6 (done):** Storage materialization + cfQ cleanup — `scores = TRUE`
  / `keep_fits = TRUE` storage,
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md),
  `tidy(what = "scores")`, cfQ hard error, cross-check tests.

- **M7 (done):** Documentation — README.Rmd, intro vignette, pkgdown
  site, three targeted vignettes (engines, ordinal, Forbes extension).

- **M8 (done):** Plot customization — `show_r`/`r_digits`, `mono`,
  `show_level_labels`/ `level_label_size`, `node_labels`,
  `primary_only`, `drop_pruned`/`compress_levels` on
  [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md);
  private `.drop_pruned_nodes()` helper in `layout.R`.

- **M9 (done):** Visualization round 2 — `show_arrows`,
  `edge_linewidth`, `legend` on
  [`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md);
  new `ackwards-visualization.Rmd` vignette; Forbes vignette slimmed
  toward the paper; intro vignette trimmed and stale comment fixed.

- **M10 (done):** Conformance + robustness —
  [`summary.ackwards()`](https://jmgirard.github.io/ackwards/reference/summary.ackwards.md) +
  [`print.summary_ackwards()`](https://jmgirard.github.io/ackwards/reference/print.summary_ackwards.md)
  (previously documented but unimplemented §10 method); ESEM
  Heywood/improper-solution warning (`theta ≤ 0`, parity with EFA
  engine); `cor="spearman"` + `method="esem"` inconsistency warning;
  DESIGN.md §8 reconciled to list only PA + MAP (EKC/EGA marked out of
  scope).

- **M11 (done):** Edge-label polish + `show_r` decoupling — APA
  `.format_r()` helper (strip leading zero, pad trailing zeros, suppress
  `-.00`), `geom_label` with perpendicular offset + white halo +
  `r_label_size` arg; **decouple `show_r` from `drop_pruned`** (default
  `FALSE` everywhere); Forbes vignette updated to two-figure (labeled +
  unlabeled) treatment; `.lintr` added to `.Rbuildignore` (R CMD check
  fully clean: 0 errors, 0 warnings, 0 notes).

- **M12 (done):** Best-practice `suggest_k` — PA-FA added alongside
  PA-PC (`psych::fa.parallel(fa = "both")`); VSS-1/VSS-2 surfaced from
  existing [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html)
  call; Comparison Data (CD) added via
  [`EFAtools::CD()`](https://rdrr.io/pkg/EFAtools/man/CD.html) gated by
  [`rlang::is_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
  (skips gracefully when absent); new `seed` arg; enriched `suggest_k`
  object (`k_parallel_pc`, `k_parallel_fa`, `k_vss1`, `k_vss2`, `k_cd`,
  `cd_available`, expanded `criteria` table); redesigned
  [`print.suggest_k()`](https://jmgirard.github.io/ackwards/reference/print.suggest_k.md)
  multi-criterion table; new
  [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
  three-panel ggplot2 diagnostic (scree/PA + MAP + VSS); `EFAtools`
  added to Suggests; DESIGN.md §8 and §12 updated.

- **M13 (done):** Rotation honesty — removed dead `kappa` argument;
  removed `rotation` argument entirely (only varimax is valid; exposed
  as a user arg it implied quartimax/equamax were options); renamed
  “cfT” → “varimax” throughout all three engine internals, the result
  object, print output, README, and docs; `@section Defaults` explains
  T′=T⁻¹ → W′RW algebra + varimax = what all reference papers used;
  DESIGN.md §4, §9, §14.1, §14.7 updated.

- **M14 (done):** Dedicated
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  vignette — new
  [`vignette("ackwards-suggest-k")`](https://jmgirard.github.io/ackwards/articles/ackwards-suggest-k.md)
  (“Choosing k: How Many Factors?”) with all five criteria (pros/cons,
  bias direction, engine pairing), argument coverage
  (cor/n_iter/seed/k_max including ordinal→Pearson and
  PA-non-reproducibility caveats), and a worked BFI recommendation;
  intro vignette Step 1 trimmed to default call + pointer;
  `_pkgdown.yml` lists the new article first under Deep dives; README
  stale two-criteria description corrected to five criteria; DESIGN.md
  §8 and §15 updated. Post-review fixes: worked-example prose corrected
  to match actual output (PA-PC=5/PA-FA=6/CD=8, CD-outlier explanation);
  CD table “Conservative” → “Accurate in simulation; can over-retain on
  large, correlated samples”; three new tests covering
  [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
  and
  [`print.suggest_k()`](https://jmgirard.github.io/ackwards/reference/print.suggest_k.md)
  for the `k_parallel_fa=NA` and `cd_available=FALSE` branches. (719
  tests pass, 1 skip.)

- **M15 (done):** Naming clarity & consistency pass — `k`→`k_max`,
  `method`→`engine`, `scores`→`keep_scores`, `align`→`align_signs` on
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md);
  `$method`→`$engine`, `$cor_type`→`$cor` on the result object;
  `method`→`edge_method` on
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md).
  All S3 methods, tests, 6 vignettes, README, NEWS, CLAUDE.md, and
  DESIGN.md updated. (724 tests pass, 1 skip; 0/0/0 R CMD check.)

- **M16 (done):** Estimator-aware missing-data handling — new
  `missing = c("pairwise","listwise", "fiml")` argument on
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  Default `"pairwise"` preserves existing behaviour and warns when NAs
  present; for ESEM WLSMV/ULSMV passes `available.cases` to lavaan (full
  N, honest). `"listwise"` reduces to complete cases pre-fit for
  consistent N. `"fiml"` (ESEM ML/MLR only) uses Full Information ML via
  lavaan and derives edge R from the FIML saturated model. Fixes ESEM
  ML/MLR fit-vs-edges inconsistency for `"listwise"` and `"fiml"`
  (including a silent bug where the h1 extraction fell back to
  pairwise). Adds `.resolve_missing()` helper; records
  `meta$missing`/`meta$n_complete`. 36 tests in `test-missing.R`;
  missing-data section added to `ackwards-engines.Rmd`; DESIGN.md §9 and
  §14 updated. (778 tests pass, 1 skip; 0/0/0 R CMD check.)

- **M17 (done):** GitHub 0.1.0 release prep — license switched CC BY 4.0
  → MIT (`LICENSE` stub + `LICENSE.md` updated, DESCRIPTION updated);
  version bumped `0.0.0.9000 → 0.1.0`; README MIT badge

  - comment mismatch fixed; `inst/CITATION` added (Goldberg 2006 +
    package), version field dynamic via `meta[["Version"]]`,
    hand-written Goldberg prose removed from README to avoid
    duplication; README rebuilt with correct two-entry citation output;
    NEWS.md restructured to curated capability-grouped summary
    (development history dropped — pre-release, all captured in git);
    `_pkgdown.yml` 0.1.0 release URL registered; pkgdown rebuilt
    cleanly. Post-review: citation guard test added to `test-utils.R`.
    Note: pkgdown 2.2.0 renders all root `.md` files via `package_mds()`
    with no config-based exclusion — CLAUDE.md/DESIGN.md appear as
    unlinked pages; not fixable without moving files. (787 tests pass, 1
    skip; 0/0/0 R CMD check; all URLs clean.) Tag: owner runs
    `git tag -a v0.1.0 -m "ackwards 0.1.0" && git push origin v0.1.0`.

- **M18 (done):** Factor interpretation & label scaffolding —
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  (salient per-factor item listing, `|loading| >= cut`, grouped cli
  print) and
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  (node_labels scaffold, styles: “id”/“forbes”/“blank”, prints editable
  c(…) literal). Both are pure consumers of the existing light core; no
  new dependencies, no invariant or default changes. Intro vignette and
  visualization vignette updated; pkgdown reference index updated;
  DESIGN.md §10/§11/§15 updated. Post-review hardening:
  `label_template(style = "forbes")` now guards `k_max > 26` (LETTERS
  exhaustion) with a loud `cli_abort` and documents the constraint; the
  duplicate
  [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)
  call in the forbes branch was cached; the `sort = FALSE` test asserts
  order against
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) rather than
  set membership; EFA smoke tests added for both helpers
  (engine-agnosticism); stale `tests/testthat/_problems/` removed.

- **M19 (done):** Dedicated interpretation/labeling vignette —
  documentation-only. New
  [`vignette("ackwards-interpret")`](https://jmgirard.github.io/ackwards/articles/ackwards-interpret.md)
  (“Interpreting and Labeling Factors”) owns the
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  → name →
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  → `autoplot(node_labels=)` workflow plus hierarchy-aware naming
  (parent vs child, blends, reorganizing factors via lineage/edges) and
  the sign-alignment caveat. Listed in pkgdown Deep dives after
  `ackwards-suggest-k`. Intro Step 5 trimmed to a slim
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  example

  - pointer; visualization vignette keeps the labeling mechanic +
    cross-ref (naming judgment moved to the new vignette). DESIGN.md §15
    entry (amends §10/§11). Original milestone was documentation-only;
    post-review hardening added a guard test for the vignette’s edges
    idiom,
    [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)/[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)→[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
    idiom smoke tests, switched the interpret vignette to
    `cor = "polychoric"` (consistency + drops the ordinal warning), and
    **fixed swapped m5f3/m5f5 labels** (Conscientiousness/Openness) in
    the interpret and visualization vignettes. (935 tests pass, 1 skip;
    0/0/0 R CMD check.)

- **M20 (done):** CRAN submission readiness + example legibility. A
  *release-readiness* milestone (not a DESIGN.md §15 feature milestone)
  ahead of tagging/submitting 0.1.0. Five waves:

  1.  Statistical/correctness: `.standardize()` (na.rm-aware) replaces
      [`scale()`](https://rdrr.io/r/base/scale.html) in
      `.compute_scores()` and
      [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
      scores route; `detect_ordinal()` guarded against all-`NA` columns;
      stale “oblique rotations” wording in
      [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
      roxygen fixed.
  2.  Example conversions: all `\dontrun{}` removed; fast examples use
      [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html)
      guards; slow/stochastic `suggest_k` blocks use `\donttest{}`.
  3.  Submission metadata: three DOIs added to `DESCRIPTION` (Goldberg
      2006, Waller 2007, Forbes 2023); `NEWS.md` restructured into a
      single `0.1.0` entry; `cran-comments.md` added.
  4.  Example legibility:
      [`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
      gains `sort = c("none","strength")` for edges; README/vignettes
      rewrote `order(-abs(...))` → `tidy(sort="strength")`,
      `identical(round(...))` →
      [`all.equal()`](https://rdrr.io/r/base/all.equal.html),
      `grep("^\\.m5",...)` → `startsWith(names(...), ".m5")`,
      double-`rbind` pattern → intermediate variable.
  5.  Verify: styler + lintr clean; `R CMD check --as-cran` → 0/0/0;
      README rebuilt. Owner next steps: `devtools::check_win_devel()`,
      `rhub::rhub_check()` before actual CRAN upload. (947 tests pass, 1
      skip; 0/0/0 R CMD check.)

- **M21 (done):** Onboarding & usability pass (pre-CRAN). Four parts:

  1.  `psych` Suggests→Imports (engine substrate for default PCA/EFA;
      never needed an install prompt for core use); `GPArotation`
      removed entirely (varimax routes through
      [`stats::varimax`](https://rdrr.io/r/stats/varimax.html); verified
      never loaded). All `check_installed("psych")` guards removed; two
      vestigial `skip_if_not_installed ("GPArotation")` test lines
      deleted. DESIGN.md §3/§12 updated.
  2.  `bfi25` dataset bundled: 1 000 rows sampled from
      `psych::bfi[, 1:25]` (seed 42, NAs preserved) via
      `data-raw/bfi25.R`; documented in `R/data.R` with `@source`
      (Revelle/psych/SAPA/IPIP — items are public-domain IPIP). All
      `@examples` + 6 vignettes + README.Rmd migrated to `bfi25`;
      suggest-k worked-example prose regenerated (n=875: PA-PC=5,
      PA-FA=6, MAP=5, VSS-1=4, VSS-2=5, CD=6; consensus 4–6). Oracle
      snapshot stays on full `psych::bfi[, 1:25]`.
  3.  README “Learn more” table now lists all 7 vignettes (Visualization
      and Interpreting & labeling were missing).
  4.  [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
      adds a 4th “CD (RMSE, minimize)” panel (2×2 grid when CD
      available; unchanged 3-panel single-column otherwise); `suggest_k`
      object gains `cd_rmse` field. CD vline removed from MAP panel.
      Tests for 4-panel/3-panel branches and `cd_rmse` field.
      Post-review: committed stale README hero figure (bfi25 render);
      added `test-data.R` (6 assertions on `bfi25`
      shape/cols/class/NAs/range); corrected `data.R` `@source` year
      2025→2026. (974 tests pass, 1 skip; 0/0/0 R CMD check.)

- **M22 (done):** Correlation-matrix input (PCA/EFA-only).
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  and
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  now accept a pre-computed correlation matrix detected automatically
  (square, symmetric, unit diagonal). New `n_obs` arg (required for
  EFA+R, optional for PCA+R). ESEM gated off (lavaan needs raw data).
  `cor`/`missing` args ignored+warned for R input; `$cor` stored as
  `NA`; print shows `"(user-supplied matrix)"`.
  `keep_scores=TRUE`/[`augment()`](https://generics.r-lib.org/reference/augment.html)/`tidy(what="scores")`
  all error clearly. CD gated off in
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  with info note. Edges from R-matrix and raw-data paths are identical
  within floating-point tolerance (same W’RW algebra).
  `.is_cor_matrix()` + `.validate_cor_matrix()` +
  `.check_maybe_cov_matrix()` helpers in `utils.R`. `meta$input_type`
  field added. Non-ASCII chars replaced across all R files (0/0/0
  clean). Post-review: non-PD warning tested; covariance-matrix
  detection added (targeted error via `.check_maybe_cov_matrix()` in
  both functions); `prune="redundant"` + cor_matrix test;
  [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md) +
  cor_matrix test; invalid `n_obs` tests for
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md);
  `missing(missing)` comment; NEWS folded into 0.1.0. 44 tests in
  `test-cor-input.R`, 41 in `test-utils.R`. (1035 tests pass, 1 skip;
  0/0/0 check.)

- **M23 (done):** Test-coverage hardening — raised `covr` from **93.37%
  → 100%** overall (every file at 100%). Strategy: test first (`# nocov`
  only for genuinely unreachable defensives). New tests: `suggest_k.R`
  cor-ignored/spearman+CD/CD-dash branches; `label_template.R` k\>26
  guard; `prune.R` `.tucker_phi` all-zeros / `.phi_pairs` adjacent /
  `print.ackwards` artefact+phi-note; `compute_edges.R` algebra explicit
  path; `engine_esem.R` NULL-SE skip in `tidy(what="loadings_se")`;
  `summary.R` EFA chi/dof row, empty-redundant “(none)”, phi-note,
  empty-lineage. `# nocov` markers on unreachable engine defensives:
  PCA + EFA k=1 sign flip, EFA convergence-fail break + tenBerge
  fallback, ESEM lavaan version guard + convergence fail + std_sol
  NULL + tenBerge fallback + W NULL + warning muffler + null-fit +
  cov2cor unreachable branch + Pearson fallback + fitMeasures error
  handler + Phi shape fallback; `prune.R` empty-chains guard +
  null-nodes branch; `layout.R` top-down zero-weight + bottom-up
  orphaned-parent fallbacks; `summary.R` PCA eigenvalue miss +
  null-nodes branch. **No behavior change.** Local verification only (no
  CI coverage workflow). Post-review: rewrote the `print.top_items`
  empty-shown-level test to deterministically hit `interpret.R`
  [`next`](https://rdrr.io/r/base/Control.html) (was skipping on bfi25
  -\> covered the line + dropped the spurious 2nd skip); `# nocov`’d the
  EFA k=1 sign flip to match the PCA analogue; added a covr-limitation
  note explaining the ESEM warning-muffler `# nocov` (live code covr
  cannot instrument, not dead code). (1080 tests pass, 1 skip; 0/0/0 R
  CMD check.)

- **Test-suite speedup (post-M23, no source/behavior change):**
  [`EFAtools::CD()`](https://rdrr.io/pkg/EFAtools/man/CD.html) dominated
  test time (~8 s/call on the full 2800×25
  [`psych::bfi`](https://rdrr.io/pkg/psych/man/bfi.html) vs ~0.1 s for
  `fa.parallel` + `vss`). `test-suggest_k.R` rewritten to run the
  CD-bound structural/logic/print/autoplot tests on small `bfi25`
  subsets (CD ~0.8 s, cached), keeping CD/PA/MAP/VSS *real* — only the
  data is smaller — plus one
  full-[`psych::bfi`](https://rdrr.io/pkg/psych/man/bfi.html)
  integration smoke for the `n_obs=2800`/`n_vars=25` assertion. PA-cap
  test moved to `bfi25[, 1:15]`/`k_max=2` (genuine clamp).
  `devtools::test()` ~80 s → ~46 s; coverage held at 100%. Also fixed
  three leaked secondary warnings so the suite is warning-free: the
  spearman/CD `n_factors_max` notice (`test-suggest_k.R`) and the
  NA-pairwise + ordinal notices in two `test-cor-input.R`
  `n_obs`-ignored tests (switched to continuous/complete data and
  within-CD `k_max`); added a deterministic default-`k_max` test to
  preserve the `suggest_k.R:208` branch. Comfortably within CRAN’s
  check-time budget. (1083 tests pass, 1 skip; 100% coverage.)

- **M24 (done):** Vignette communication pass (documentation-only).
  Stacked long-format `kable` comparison tables replaced by wide
  **`gt`** tables (one row per item/edge, one column per engine/basis,
  explicit Δ column) in `ackwards-engines` (loadings: `|EFA|−|PCA|`;
  edges: `EFA−PCA`) and `ackwards-ordinal` (loadings: `poly−pearson`;
  edges: `poly−pearson`). Factor/sign alignment asserted via
  [`stopifnot()`](https://rdrr.io/r/base/stopifnot.html) before
  differencing; engine disagreements on primary parents surface as `NA`.
  `ackwards-forbes`: `prune-nodes` raw
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) print
  replaced by styled gt table; narrated counts and stale skip-edge claim
  (m3f2→m5f1, now m3f2→m5f2 on bfi25) converted to inline `` `r` ``
  expressions. `skip-edges`/`thresholds` tables migrated to gt for
  visual consistency. All presentation code in `echo = FALSE` — no
  plumbing visible to readers. Audit of the other four vignettes: no
  changes needed (raw console prints are pedagogically appropriate).
  `gt` added to Suggests. Post-review (six findings): unified every Δ
  column to a **magnitude** convention (`|x|−|y|`) so the directional
  captions read correctly for negatively-signed cells — fixes a
  sign-trap where the polychoric-strengthened negative edge `m2f2→m3f2`
  had shown a negative signed Δ (ordinal edges + loadings and engines
  edges switched from signed to `abs`; engines loadings was already
  `abs`); derived the forbes redundant-level prose (`flagged_ids` +
  per-level `lvl_summary`) inline instead of hard-coding “entire k=4
  level / two factors at k=2,3”; converted the `prune-artefact` raw
  print to a styled gt table (call still echoed) with inline
  `n_artefact`; added guard tests for the
  [`knitr::kable`](https://rdrr.io/pkg/knitr/man/kable.html) fallback
  branch, the NA primary-parent-disagreement merge, the magnitude-delta
  sign property, and the inline-derivation helpers. (1123 tests pass, 1
  skip; 0/0/0 R CMD check.)

- **M25 (done):** Deferred-items pass — three waves. (1)
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  gains `criteria` arg (`rlang::arg_match(multiple=TRUE)`); any subset
  of five criteria may be requested; non-requested `k_*` fields → `NA`;
  shared computation (`fa.parallel(fa="both")` / `vss()` called at most
  once); [`print()`](https://rdrr.io/r/base/print.html) /
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  render only requested criteria; consensus from requested only. (2)
  `prune="artefact"` now populates `x$prune$structural` — per-factor
  `few_items` / `orphan` / `split_merge` signals (Forbes Fig 2); new
  `min_items = 3L` and `orphan_r = 0.5` args; flag/report only, never
  auto-prune;
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  report flagged count. (3) `redundancy_phi = NULL` auto-resolves: PCA →
  no φ (exact algebra); EFA/ESEM → `0.95` (Lorenzo-Seva & ten Berge
  2006; announced via cli, Invariant 6); `NA` is explicit opt-out.
  Bootstrap CIs on edges remain deferred. `CLAUDE.md` Resolved
  defaults + Out of scope updated; DESIGN §9/§14 updated. Post-review
  hardening: added loud validation for `min_items` (positive integer) /
  `orphan_r` (`[0,1]`); roxygen *why* for both defaults (three-indicator
  rule; moderate-correlation rationale); deterministic unit tests for
  the previously-untested `split_merge = TRUE` path (mock loadings) and
  the Wave-3 φ-decision outcome change (`.find_redundant_chains` flags
  under no-φ but not under φ\>0.95) plus an end-to-end EFA
  subset-property test (auto-φ flagged set ⊆ \|r\|-only set); documented
  `criteria=` in `ackwards-suggest-k.Rmd`, structural signals + φ
  auto-default in `ackwards-forbes.Rmd`, `criteria=` in DESIGN §8;
  simplified an unreachable `n_struct` branch and `# nocov`’d the
  EFAtools-absent CD note (covered only when EFAtools missing). Coverage
  back to 100%. (1219 tests pass, 2 skip; 0/0/0 R CMD check.)

- **M26 (done):** ESEM performance for large item sets — two
  complementary speedups, no behaviour change. (1) **Cached sample
  statistics:** the ESEM engine no longer recomputes lavaan’s
  data-derived sample statistics (thresholds, polychoric matrix,
  asymptotic weight matrix NACOV/WLS.V) at every level — they depend
  only on the data, not `nfactors`, so they are harvested once at the
  anchor level (k=1) via `fit@SampleStats` and reused for every deeper
  level through lavaan’s `slotSampleStats=` argument. Verified
  **bit-identical** (0.00e+00 loading/edge diff) vs. the from-raw path.
  Dominant saving at large p (polychoric+NACOV recompute was O(p²)+ ×
  k_max). (2) **Parallel per-level fits:** `esem_levels()` refactored
  into a slim per-level worker `.esem_fit_one()` (returns the level
  contract, not the heavy fit, to avoid serialising a duplicate NACOV)
  dispatched through `.esem_lapply()` —
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  when installed (gated by
  [`rlang::is_installed()`](https://rlang.r-lib.org/reference/is_installed.html)),
  serial `lapply` fallback otherwise. No `ncores` arg: users set
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  (sequential default → no behaviour change; `multisession`/ `multicore`
  to parallelize). `future.seed = TRUE` (lavaan::efa is mildly
  RNG-stochastic, ~1e-6); reproducible across plans when `seed` supplied
  (verified 0.00e+00). Invariant 7 preserved: all levels fit, then
  assembly truncates at the first non-converged/failed level and emits
  all cli warnings in deterministic level order (workers never signal
  conditions). `future.apply` added to **Suggests** (flagged +
  approved). New `@section Performance` on
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md);
  “Performance with many items” section in `ackwards-engines.Rmd`
  (incl. EFA+polychoric cheaper-route pointer); NEWS + DESIGN §12/§14
  updated. Three new tests (both `.esem_lapply` branches via
  `local_mocked_bindings`; serial-vs-`multicore` identity, skipped on
  Windows / when future absent). Coverage held at 100%. (1225 tests
  pass, 2 skip; 0/0/0 R CMD check.)

- **M27 (done):** ESEM fit & SEs as first-class output — CRAN-readiness
  pass making per-level fit indices and rotation-aware loading SEs
  usable and correctly framed. Additive surfacing + docs only; no
  invariant or resolved-default change; no new `Imports`; version stays
  `0.1.0`.

  1.  **[`glance()`](https://generics.r-lib.org/reference/glance.html)
      carries fit**: deepest-converged-level `CFI`, `TLI`, `RMSEA`,
      `SRMR`, `BIC` appended as a consistent five-column set across all
      engines (NA where unavailable: CFI/SRMR for EFA, all five for
      PCA). `.glance_fit()` internal helper.
  2.  **`tidy(what="fit")` gains `format` and `cutoffs`**:
      `format="wide"` returns one row per **non-anchor** level (k \>= 2;
      the saturated 1-factor anchor is dropped to match
      [`summary()`](https://rdrr.io/r/base/summary.html) and
      `autoplot(what="fit")`) with index columns (long default
      byte-identical to previous output); `cutoffs=TRUE` appends a
      `meets` flag against Hu & Bentler (1999) thresholds (CFI/TLI \>=
      .95, RMSEA \<= .06, SRMR \<= .08) via `.fit_cutoffs()` /
      `.flag_fit()` / `.fit_long_to_wide()` helpers; report-only, never
      gates.
  3.  **Loading CIs folded into `tidy(what="loadings")`**: added `se`,
      `ci_lower`, `ci_upper` columns for all engines (NA for PCA/EFA,
      populated for ESEM); `conf_level=0.95` argument controls width.
      The now-redundant `tidy(what="loadings_se")` accessor was
      **removed** (package unreleased, no deprecation needed); the
      internal `lev$loadings_se` level-contract field is unchanged.
  4.  **`autoplot(x, what="fit")`**: two-panel ggplot2 line chart
      (CFI/TLI top panel; RMSEA/SRMR bottom panel; Hu & Bentler
      reference lines). `what="hierarchy"` (default) unchanged. PCA
      returns informative empty plot.
  5.  **[`summary()`](https://rdrr.io/r/base/summary.html) fit lines**
      annotated with `✔`/`✘` pass/fail marks for thresholded indices.
  6.  **Honesty framing in
      [`print()`](https://rdrr.io/r/base/print.html) and
      [`summary()`](https://rdrr.io/r/base/summary.html) footer**:
      per-level fit describes the k-factor solution at that level; does
      not validate edges or the hierarchy.
  7.  **Dedicated section “Per-level fit: what it tells you (and what it
      doesn’t)”** in
      [`vignette("ackwards-engines")`](https://jmgirard.github.io/ackwards/articles/ackwards-engines.md):
      level-vs-hierarchy distinction, reporting workflow (wide table +
      cutoffs + fit plot), ESEM cost/benefit, when to care (exploratory
      vs publication).
      [`vignette("ackwards-ordinal")`](https://jmgirard.github.io/ackwards/articles/ackwards-ordinal.md)
      updated to use `format="wide"` and cross-reference. No new
      `Imports`; no design invariant changed; version stays `0.1.0`.
      Post-review fixes (from /post-milestone-review): removed the
      redundant `tidy(what="loadings_se")` accessor (per the approved
      plan; the implementation had kept it); wide fit table now drops
      the anchor level for consistency with
      [`summary()`](https://rdrr.io/r/base/summary.html)/the fit plot;
      fit-plot panel titles made engine-aware (EFA shows “TLI”/“RMSEA”,
      not “CFI / TLI”); removed dead `excl` variable and an unreachable
      `is.null(cut)` guard in `.ba_fit_plot()`; added tests for the
      `.glance_fit` truncation/empty-fit branches and the empty fit-plot
      branch, restoring coverage to 100%. (1285 tests pass, 2 skip;
      0/0/0 R CMD check; coverage 100%.)

- **M28 (done):** CD correctness & honesty fix — two verified defects in
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
  Comparison Data (CD) output resolved. No invariant or resolved-default
  change; no new `Imports`; version stays `0.1.0`.

  1.  **Trailing-zero bug fixed.**
      [`EFAtools::CD`](https://rdrr.io/pkg/EFAtools/man/CD.html) fills
      its `RMSE_eigenvalues` matrix only up to column `k_cd + 1`;
      remaining columns stayed at the matrix’s zero initialisation.
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
      was averaging all `k_max` columns and plotting them, so the RMSE
      curve plunged to 0 at higher k and
      [`which.min()`](https://rdrr.io/r/base/which.min.html) landed on a
      spurious zero (e.g. `bfi25`, `k_max = 8`: argmin = 8, star at 6).
      Fix: after extracting `cd_rmse`, columns `(k_cd + 2) … k_max` are
      set to `NA_real_`; column `k_cd + 1` (the tested-but-rejected
      level) is kept visible. The
      [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
      CD data frame filters NA rows before building plot data, so
      `ggplot2` sees no NA values and draws no misleading zero-RMSE
      tail.
  2.  **“minimize” label corrected.**
      [`EFAtools::CD`](https://rdrr.io/pkg/EFAtools/man/CD.html) uses a
      sequential one-sided Wilcoxon test (default α = 0.30): retain a
      level while adding it *significantly* reduces RMSE; stop at the
      first non-significant improvement. The starred k is the last
      retained level and need not be the visible minimum. The plot facet
      was renamed from `"CD (RMSE, minimize)"` to
      `"CD (RMSE; sequential test)"`. The
      [`autoplot.suggest_k()`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
      roxygen was corrected to describe the stopping rule and to note
      that the star ≠ guaranteed minimum. The
      `vignettes/ackwards-suggest-k.Rmd` “What it does” and plot-panel
      description sections were updated to match.
  3.  **Deferred-item assessments recorded.** EAP scoring **declined**
      (DESIGN §14 item 13 updated): EAP’s shrinkage attenuates
      cross-level correlations, the primary signal bass-ackwards
      measures;
      [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
      seam preserved but implementation not planned. Bootstrap CIs on
      skip-level edges **remain deferred** with rationale added (DESIGN
      §14 Forbes-extension bullet). CLAUDE.md “Out of scope” updated to
      reflect the EAP decision. Tests: new tests covering no spurious
      zeros in `cd_rmse`, correct NA masking, no NA in plot data, the
      renamed panel label, the star at `k_cd`, and a divergent-fixture
      test asserting the star follows `k_cd` even when the RMSE minimum
      sits at a deeper level (the sequential-test case the relabel
      exists to clarify); plus the updated panel-level string in a prior
      test. Post-review (from /post-milestone-review): added the
      divergent-fixture star test (the one nice-to-have flagged; the
      vector-form EFAtools masking branch remains pre-existing `nocov`,
      untestable on a machine with current EFAtools). (1294 tests pass,
      2 skip; 0/0/0 R CMD check; coverage 100%.)

- **M29 (done):** Strip milestone numbers from user-facing
  documentation. Audit of `vignettes/*.Rmd`, roxygen/`man/*.Rd`, and
  `README.Rmd`/`README.md` found them already free of internal milestone
  tags (e.g. `(M24)`); the only leak was in `NEWS.md`, which ships with
  the package and renders as the pkgdown changelog. Fixed: `NEWS.md`’s
  “Vignette comparison tables reworked for legibility (M24)” heading
  reworded to drop the tag. Internal `R/*.R` code comments (`# (M16)`,
  `# (M26)`, `# (M5)` — developer traceability back to `MILESTONES.md`)
  were deliberately left alone; they are not user-facing and were out of
  scope per the approved plan. New regression test
  `tests/testthat/test-docs-no-milestone-refs.R` asserts no `(M\d+)` tag
  appears in `NEWS.md`, `README.md`, `README.Rmd` (source as well as
  rendered/shipped file), or `vignettes/*.Rmd`; it skips gracefully
  (rather than failing) when those files aren’t reachable from the
  test’s working directory, which happens under a full installed-package
  check cycle rather than `devtools::test()`/source-tree testing. This
  is a deliberate source-tree hygiene check — the docs are not part of
  the installed package, so it is enforced by `devtools::test()` and the
  test-coverage CI job (both run in the source tree) rather than
  shipping doc copies into `inst/` solely to run it under R CMD check.
  Post-review (from /post-milestone-review): broadened the guard to scan
  `README.Rmd` alongside `README.md`, and made DESIGN.md §15’s pointer
  range-free (it had read “M1–M26”, a number that re-stales every
  milestone; §15 remains a pointer to `MILESTONES.md`, not a log). No
  invariant or resolved-default change; no new `Imports`; no DESIGN.md
  contract change; version stays `0.1.0`. (1297 tests pass, 2 skip;
  0/0/0 R CMD check; coverage 100%.)

- **M30 (done):** Citation hygiene — put the right citation in the right
  role throughout the package. No invariant or resolved-default change;
  no new `Imports`; no DESIGN.md contract change; version stays `0.1.0`.

  1.  **`inst/CITATION`: Goldberg entry removed.** It carried two
      [`bibentry()`](https://rdrr.io/r/utils/bibentry.html)s — a
      Goldberg
  2.  *Article* and a Girard *Manual* — so `citation("ackwards")`
      returned both and conflated “how do I cite this software” with
      “what’s the method’s source paper.” Now it holds only the Girard
      *Manual* entry; Goldberg’s method paper remains cited in
      `DESCRIPTION` and roxygen `@references`, not as a software
      co-citation.
  3.  **[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
      `@references` gains Forbes (2023).** It previously listed only
      Goldberg
  4.  and Waller (2007), but
      [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
      implements all three algorithmic sources: the original method
      (Goldberg), the exact `W'RW` edge algebra (Waller), and the
      extended method (Forbes; `pairs = "all"`, `prune =`). Kim &
      Eaton (2015) and Forbush et al. (2024), named in `@details` prose
      as application examples, were deliberately **not** promoted to
      `@references` — they used the method, they aren’t a source of it.
  5.  **Audit sweep found the rest already correct.**
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
      Forbes citation (an overextraction-caution reference, not an
      implementation claim), and the per-vignette `## References`
      sections (`ackwards-intro`/`ackwards-engines` → Goldberg/Waller,
      `ackwards-forbes` → Goldberg/Forbes/Tucker, `ackwards-suggest-k` →
      Forbes, `ackwards-ordinal` → Pearson) were all already scoped
      correctly; no changes needed there.
  6.  **README Citation section rebuilt.** The “please cite both the
      method and the package” prose no longer matches a single-entry
      [`citation()`](https://rdrr.io/r/utils/citation.html); reworded to
      “please cite the package” plus explicit guidance to also cite
      Goldberg (2006) for the method, Waller (2007) if relying on the
      exact algebra, and Forbes (2023) if using the extended method.
      Rebuild incidentally picked up an unrelated but genuine staleness
      fix already shipped in
      [`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
      (M27’s per-level fit-indices honesty note) that had never been
      re-knitted into `README.md`. Tests: updated the existing
      `citation("ackwards")` test (`test-utils.R`) for the single-entry
      state (was asserting two entries incl. a Goldberg year check);
      added `tests/testthat/test-docs-citations.R`, a source-tree
      regression test (same style as `test-docs-no-milestone-refs.R`)
      asserting `man/ackwards.Rd`’s `\references` block carries the
      Goldberg, Waller, and Forbes DOIs. Post-review (from
      /post-milestone-review): addressed the three nice-to-haves. (a)
      The `test-docs-citations.R` `\references` extraction now bounds to
      the block’s closing brace (line-based: first lone `}` after
      `\references{`) instead of running to end-of-file, so a DOI that
      ever moved into `@details` prose could no longer falsely satisfy
      it. (b) Added a `## References` section citing Forbes (2023) to
      `vignettes/ackwards-interpret.Rmd` and
      `vignettes/ackwards-visualization.Rmd`, which name “Forbes (2023)”
      in prose (letter convention / publication style) but previously
      had no local reference entry. (c) Added a source-text guard to
      `test-docs-citations.R` asserting `inst/CITATION` names Girard and
      not Goldberg, complementing the runtime
      [`citation()`](https://rdrr.io/r/utils/citation.html) check
      (catches a reintroduced Goldberg bibentry even where the package
      cannot be loaded). (1303 tests pass, 2 skip; 0/0/0 R CMD check;
      coverage 100%.)

- **M31 (done):** Correctness & output-honesty sweep — first milestone
  of the M31–M38 pkgdown-review epic, sequenced correctness-first (fix
  output bugs before rewriting the vignettes that display them). No
  breaking renames; no invariant or resolved-default change; no new
  `Imports`; no DESIGN.md contract change; version stays `0.1.0`.

  1.  **ESEM `p_value` NA, root-caused and fixed.** Live-repro on
      `bfi25` with `cor = "polychoric"` (WLSMV) showed `p_value` NA at
      every level. Root cause: lavaan’s naive chi-square test has no
      valid reference distribution under WLSMV/ULSMV — lavaan’s own
      `summary(fit, fit.measures = TRUE)` literally labels that column
      “P-value (Unknown)” and reports NA; only the
      mean-and-variance-adjusted “scaled” test (`pvalue.scaled`) has a
      genuine null distribution for these limited-information
      estimators. `engine_esem.R`’s `.esem_fit_one()` now falls back to
      `chisq.scaled`/`df.scaled`/`pvalue.scaled` whenever the naive
      `pvalue` is NA and the scaled variant exists — a no-op for ML/MLR,
      whose naive p-value is already valid. A genuinely saturated level
      (`dof = 0`, e.g. a small item set at high k) still reports NA on
      both tests — there is no chi-square test to perform on a model
      that fits perfectly by construction; this is the “saturated low-k
      level” case the milestone brief anticipated, and it coexists with
      (rather than replaces) the WLSMV-naive-test cause found on the
      non-saturated `bfi25` case.
  2.  **BIC promoted to a first-class ESEM fit index.** Previously
      absent from `fit_info` entirely (not even as NA) —
      `tidy(what = "fit")` silently had no BIC row for ESEM and
      [`glance()`](https://generics.r-lib.org/reference/glance.html)
      only happened to show NA because the key was missing. Now always
      requested: real value under `"ML"`/`"MLR"` (proper
      log-likelihood), genuine `NA` under `"WLSMV"`/`"ULSMV"` (no proper
      log-likelihood for a limited-information estimator — inapplicable,
      not a bug). Loadings SE/CI were already correctly NA-safe for
      PCA/EFA in `tidy(what = "loadings")`; no change needed there.
  3.  **`_meets` cleanup.**
      `tidy(what = "fit", cutoffs = TRUE, format = "wide")`’s pivot
      (`.fit_long_to_wide()`) generated a `{index}_meets` column for
      every index present, including `chi`/`dof`/`p_value`/`BIC`, which
      have no defined threshold and were therefore always NA. Now
      restricted to indices in `.fit_cutoffs()`
      (`CFI`/`TLI`/`RMSEA`/`SRMR`). `format = "long"` is unaffected
      (`meets` stays NA for those rows, as already documented).
  4.  **Guard: `cor = "polychoric"` + `estimator = "ML"`/`"MLR"`.**
      Previously unguarded: polychoric correlations mark every item
      `ordered` for lavaan, and lavaan itself errors on ML/MLR with
      ordered indicators (“estimator ML for ordered data is not
      supported yet”) — but the failure surfaced many calls deep as a
      per-level ESEM warning, then a generic
      [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
      abort that misdiagnosed the cause as multicollinearity.
      [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
      now errors immediately with a targeted message naming the
      incompatibility and the fix. `"WLSMV"`/`"ULSMV"` with a continuous
      `cor` is deliberately left unguarded — lavaan runs it as a valid
      (if atypical) continuous WLS/ADF estimator, confirmed by a passing
      test.
  5.  **`fa.parallel`/`set.seed` reproducibility — confirmed correct,
      already documented, no code bug.**
      [`psych::fa.parallel()`](https://rdrr.io/pkg/psych/man/fa.parallel.html)
      genuinely does not respond to
      [`set.seed()`](https://rdrr.io/r/base/Random.html);
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
      roxygen, `DESIGN.md` §8, and
      [`vignette("ackwards-suggest-k")`](https://jmgirard.github.io/ackwards/articles/ackwards-suggest-k.md)
      already state this honestly (the `seed` argument governs the CD
      step only). Swept all `seed` mentions across vignettes/man pages
      for an implied-but-false PA-reproducibility claim; found none.
  6.  **Doc-bugs mirroring output, fixed in `ackwards-intro.Rmd` and
      `ackwards-suggest-k.Rmd`.** A hardcoded cumulative-variance jump
      (“22.9% → 34.7%”) had drifted from the code’s actual live output
      (23.2% → 35.5%) — replaced with inline `` `r ...` `` expressions
      so it cannot drift again. The lineage walkthrough had `m4f1`’s
      primary children backwards (claimed `m5f2`/`m5f4`; the edge table
      shows `m5f1`/`m5f4`, confirmed via
      `tidy(x, what = "edges", primary_only = TRUE)` and a live
      [`augment()`](https://generics.r-lib.org/reference/augment.html)
      correlation check). The k-by-k diagram narrative misattributed
      which traits differentiate at which level (claimed
      Agreeableness/Extraversion split at k = 4; the live
      [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
      output shows Conscientiousness/Openness split at k = 4, and
      Agreeableness/ Extraversion don’t split until k = 5) — both
      corrected to match a live run. The “no warning this time” claim
      after switching to `cor = "polychoric"` was checked and is still
      accurate. Added a clarifying note to `ackwards-suggest-k.Rmd`
      explaining why the printed “Recommendations” block shows six lines
      for five criteria (`"vss"` is one `criteria` entry sharing one
      [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html) call, but
      reports two numbers, VSS-1 and VSS-2). Tests: extended
      `test-esem.R` with fit-index-naming (`BIC` added), wide-format
      `_meets` column-set, WLSMV-scaled-p-value, and estimator-guard
      coverage. No new test files. (1323 tests pass, 2 skip; 0/0/0 R CMD
      check; coverage 100%.) Post-review (from /post-milestone-review):
      the review’s Should-fix — that the ESEM fit row mixed a scaled
      chi-square/p-value with *naive* CFI/TLI/RMSEA under WLSMV/ULSMV
      (the naive incremental indices are badly optimistic for ordinal
      data; Xia & Yang 2019: naive CFI ~0.89 vs scaled ~0.71 on BFI) —
      was addressed. Owner decision: report the **scaled** variants
      (over robust) for full internal consistency with the reported
      scaled test and to match the Mplus output the reference clinical
      workflows (Kim & Eaton 2015; Forbush 2024) use. `engine_esem.R`
      now applies one unified `.fm_scaled()` rule to the whole row
      (`chi`/`dof`/`p_value` + `CFI`/`TLI`/`RMSEA`): prefer the
      `<index>.scaled` name whenever the estimator emits one, else the
      naive value. This also fixed a latent MLR bug — MLR previously
      reported the *naive* ML chi-square/p-value (its naive p-value is
      non-NA, so the old p-value-only fallback never fired), defeating
      the purpose of robust ML; the whole MLR row is now scaled. `SRMR`
      (no scaled variant) and `BIC` (no scaling concept; NA under
      WLSMV/ULSMV) are unchanged. Roxygen (`tidy`/`glance`) and
      `NEWS.md` document the scaled reporting. The review’s nice-to-have
      test gaps were closed: added ULSMV (scaled p-value + NA BIC), MLR
      (asserts the row equals lavaan’s `*.scaled`, real BIC), and
      continuous-`cor`+WLSMV (populated fit row) tests to `test-esem.R`.
      Deliberately **deferred to M32** (the API-shape/`meta` pass):
      storing the effective `estimator` in `$meta` and adding a
      [`summary()`](https://rdrr.io/r/base/summary.html) footnote that
      names the scaled reporting — a schema change better bundled with
      M32’s meta/column decisions than bolted on here. Landed via
      follow-up branch `m31-followup-fit-variants` → PR (not reopening
      the merged M31 PR). (1340 tests pass, 2 skip; 0/0/0 R CMD check;
      coverage 100%.)
