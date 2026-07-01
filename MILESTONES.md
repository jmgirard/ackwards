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

- **M32 (done):** API-shape & naming resolutions — second milestone of
  the M31–M38 pkgdown-review epic. Four owner-reviewed naming/API
  decisions, all decided *and* implemented (no decide-only/defer split),
  plus one M31-deferred item folded in mid-milestone. All changes are
  breaking with no deprecation path (pre-CRAN, no users).

  1.  **`tidy(x, what = "fit")` long-format key column `index` →
      `statistic`.** `index` read like a row position; it actually held
      fit-index names for EFA/ESEM and eigenvalue positions for PCA.
      Ripples updated: `autoplot.R` (`.ba_fit_plot()`, `what = "fit"`),
      `summary.R` (per-level fit line, PCA eigenvalue lookup). Wide
      format (`format = "wide"`) column names are unaffected (they come
      from the values of this key, not its name).
  2.  **`k_max` naming collision, resolved by keeping the shared name.**
      [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)’s
      `k_max` (extraction depth) and
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
      `k_max` (max factors/components evaluated) are genuinely the same
      dial applied at different stages of the `suggest_k() → ackwards()`
      workflow; renaming either would create vocabulary drift without
      removing real ambiguity. Resolved via roxygen only — each
      function’s `@param k_max` now states its own meaning and
      cross-references the other’s. No API change.
  3.  **`tidy(what = "fit")` cutoff-flag output removed.** The
      `cutoffs = TRUE` argument and the `meets`/`{statistic}_meets`
      columns it produced (`.flag_fit()`) are dropped — a pass/fail
      boolean quietly endorsed Hu & Bentler (1999) thresholds the
      package elsewhere treats as conventional and contested (continuing
      the M28/M31 output-honesty trajectory). `.fit_cutoffs()` is
      retained internally: `autoplot(what = "fit")` still draws
      threshold reference lines and
      [`summary()`](https://rdrr.io/r/base/summary.html) still annotates
      inline with a check/cross mark — both are visual guides, not a
      returned column a caller could mistake for a computed judgment.
      Because `cutoffs` is no longer a formal parameter, passing it is
      now silently absorbed by
      [`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)’s
      `...` (same as any unrecognized argument) rather than erroring.
  4.  **Variance reported as a 0-1 proportion.**
      `tidy(x, what = "variance")` columns
      `variance_pct`/`cumulative_pct` (0-100) →
      `proportion`/`cumulative` (0-1). Aligns
      [`tidy()`](https://generics.r-lib.org/reference/tidy.html) with
      the engine’s internal `variance` slot (already 0-1 —
      `print.ackwards` reads it directly) and with broom/psych
      convention. Percent formatting moves to the display layer:
      [`print()`](https://rdrr.io/r/base/print.html) and
      [`summary()`](https://rdrr.io/r/base/summary.html) compute `* 100`
      at render time rather than storing a percent in tidy data.
  5.  **M31-deferred: effective ESEM estimator recorded in `$meta`.**
      M31’s log explicitly deferred this (“better bundled with M32’s
      meta/column decisions than bolted on here”) — caught during M32
      planning by re-reading the full M31 `MILESTONES.md` entry (not
      just the CLAUDE.md one-liner) and confirmed with the owner before
      folding it in. `x$meta$estimator` now stores the effective
      estimator after auto-selection
      (`"ML"`/`"MLR"`/`"WLSMV"`/`"ULSMV"`; `NA` for PCA/EFA, including
      the `cor_matrix`-input branch, which errors for `engine = "esem"`
      before reaching this point).
      [`summary()`](https://rdrr.io/r/base/summary.html) gains a
      one-line grey footnote naming lavaan’s scaled-variant fit
      reporting whenever the effective estimator is
      `"WLSMV"`/`"ULSMV"`/`"MLR"` (silent for `"ML"`, which has no
      scaled variant). Tests capture
      [`summary()`](https://rdrr.io/r/base/summary.html)’s printed
      output via `capture.output(..., type = "message")` — cli writes to
      the message/stderr stream in non-interactive sessions, not stdout,
      so plain
      [`capture.output()`](https://rdrr.io/r/utils/capture.output.html)
      silently captures nothing. Updated: `R/tidy.R`, `R/autoplot.R`,
      `R/summary.R`, `R/ackwards.R`, `R/suggest_k.R`; roxygen
      regenerated (`man/tidy.ackwards.Rd`, `man/ackwards.Rd`,
      `man/suggest_k.Rd`); vignette prose fixed in `ackwards-intro.Rmd`
      (variance columns) and `ackwards-engines.Rmd` (cutoffs example
      removed, reframed around
      [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)/[`summary()`](https://rdrr.io/r/base/summary.html)
      reference lines); `DESIGN.md` §14 gained items 22–26 and §10’s
      tidy-column reference was refreshed; `NEWS.md` documents all five
      changes. Tests: `test-print.R`, `test-efa.R`, `test-esem.R`
      updated for the renamed/removed columns; `test-esem.R` gained
      `$meta$estimator` coverage (ML default, explicit MLR,
      polychoric-auto WLSMV

  - explicit ULSMV, NA-for-PCA/EFA) and
    [`summary()`](https://rdrr.io/r/base/summary.html)
    scaled-fit-footnote coverage (present for MLR/WLSMV, absent for ML);
    `test-layout.R` asserts `autoplot(what = "fit")` draws the Hu &
    Bentler `geom_hline` reference lines at the charted thresholds
    (post-review nice-to-haves). No new test files. (1346 tests pass, 2
    skip; 0/0/0 R CMD check; coverage 100%.) Post-review
    (`/post-milestone-review`): clean — all acceptance criteria met (the
    one “`cutoffs=` errors” criterion resolved to a silent `...`-absorb
    no-op, an accepted pre-CRAN decision — no deprecation shim needed
    with no users; consistent with the M34 clean-move precedent). The
    three nice-to-have test-hardening items the review raised were added
    (the extra estimator/footnote/reference-line assertions above); no
    Blocking or Should-fix findings.

- **M33 (done):** simulated Gaussian dataset (foundation) — first data
  milestone of the M31–M38 pkgdown-review epic. Ships `sim16`: a 1
  000×16, fully continuous dataset (`data-raw/sim16.R`, `set.seed(42)`,
  base-R Cholesky sampling from an oblique 4-factor population model —
  no `MASS` dependency) alongside `bfi25`. **Population model.** 4 true
  group factors (`f1`–`f4`, 4 items each, loading `0.75`, no
  cross-loadings); factor correlations `0.45` within a metatrait
  (`f1`-`f2`, `f3`-`f4`) and `0.15` between metatraits (`{f1,f2}` =
  metatrait 1, `{f3,f4}` = metatrait 2); uniformly `0.4375` uniquenesses
  by the design’s symmetry. **Two drivers, both realized without a
  contrived mechanism.** (1) *Clean Pearson showcase*: the data are
  continuous (~1000 distinct values/column), so
  `ackwards(sim16, engine = "pca"/"efa")` never triggers the
  ordinal-detection warning that `bfi25`‘s Likert items do. (2)
  *Guaranteed Tucker’s φ / redundancy finding for the Forbes example*:
  because the population has **exactly 4** factors, `k_max = 5` has no
  real 5th dimension to find. The originally planned mechanism (a
  deliberately planted near-duplicate “twin” sub-factor) proved
  unnecessary once prototyped — requesting a 5th factor from a genuinely
  4-factor population organically produces both an orphan factor with
  zero primary-loading items (flagged `few_items` + `orphan` under
  `prune = "artefact"`, default `min_items = 3`) and, because the true
  (non-splitting) factors persist essentially unchanged from `k = 3`
  onward, redundant parent-child chains (`|r| >= .9` and, by default,
  Tucker’s φ `>= .95`) under `prune = "redundant"`. This is a more
  honest illustration of Forbes’ redundancy concept than the originally
  planned twin-factor design, so the plan was revised during
  implementation (recorded here, not silently). **Verified ground
  truth.** `ackwards(sim16, engine = "efa")`: `k=1` general factor
  across all 16 items; `k=2` splits exactly along the metatrait line
  (`i1`-`i8` vs. `i9`-`i16`); `k=4` recovers the 4 true group factors
  exactly. `suggest_k(sim16)` reaches unanimous 6-criteria consensus
  (PA-PC, PA-FA, MAP, VSS-1, VSS-2, CD) of `k = 4`. **No exported
  simulator and no vignette changes** — a saved dataset + reproducible
  generator only, per the approved plan; the vignette restructuring that
  consumes `sim16` (clean-Pearson-first, guaranteed-φ Forbes example) is
  M37/M38. **Deviation from the approved plan (flagged, not silent):**
  no `DESIGN.md` §14 entry. `bfi25` (M21) set the precedent that dataset
  additions with no API/behavioral contract change are documented in
  roxygen (`R/data.R`) and here, not in `DESIGN.md`, which tracks design
  *decisions* rather than data assets. Files: `data-raw/sim16.R`
  (generator), `data/sim16.rda`, `R/data.R` (roxygen doc block with the
  full population model and ground-truth table), `man/sim16.Rd`
  (generated). `.Rbuildignore` gained `ROADMAP.md` (introduced during
  `/plan-milestone 33`, alongside the existing
  `CLAUDE.md`/`DESIGN.md`/`MILESTONES.md` exclusions). Tests:
  `test-data.R` gained six `sim16` tests — shape/type/no-NA; no
  ordinal-detection warning under `pca`/`efa`; exact hierarchy recovery
  at `k=1/2/4`;
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  6-criteria consensus; guaranteed redundant-chain + artefact signals at
  `k_max=5`; and bit-identical regeneration from `data-raw/sim16.R`
  (parses the script, drops the `usethis::use_data()` call via a
  [`deparse()`](https://rdrr.io/r/base/deparse.html) string match rather
  than a literal `usethis::` token, to avoid an “unstated dependency in
  tests” `R CMD check` warning). No new test files. (1372 tests pass, 2
  skip; 0/0/0 R CMD check; coverage 100%.) Post-review
  (`/post-milestone-review`) follow-up, landed via branch
  `m33-followup-realism` → PR (not reopening the merged M33 PR):
  resolved the review’s test-robustness / teaching-realism
  nice-to-haves. **Owner decision on realism:** keep `sim16` clean
  rather than muddying it — its comparative advantage over `bfi25` is a
  *known, cleanly recoverable* structure, and `bfi25` already carries
  the criterion-disagreement lesson (its six
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  criteria span `k = 4`–`6` vs. `sim16`’s unanimous `4`). `sim16` is now
  explicitly framed as the *idealized* case in
  [`?sim16`](https://jmgirard.github.io/ackwards/reference/sim16.md),
  `NEWS.md`, and the M37/M38 `ROADMAP.md` notes, to be contrasted with
  `bfi25` in the vignettes (not presented as typical). Test changes to
  `test-data.R`: (1) the `suggest_k` assertion was loosened from pinning
  all six criteria to `4` exactly (brittle — PA-PC/PA-FA/CD resample and
  can wobble ±1 across platforms/RNG) to a deterministic MAP anchor plus
  a majority-consensus check; (2) added an assertion that the k=5
  artefact factor is a primary loader for *zero* items and is exactly
  the factor flagged `few_items` + `orphan` (validating the
  [`?sim16`](https://jmgirard.github.io/ackwards/reference/sim16.md)
  claim). The `DESIGN.md` §14 non-entry is confirmed intentional (the
  M21 `bfi25` precedent: data assets are documented in roxygen +
  `MILESTONES.md`, not `DESIGN.md`). The Invariant-2 algebra-vs-scores
  cross-check on `sim16` was deferred to M37 (recorded in `ROADMAP.md`),
  where the engines vignette naturally materializes scores. (1370 tests
  pass, 2 skip; 0/0/0 R CMD check; coverage 100%.) Second follow-up
  (branch `m33-followup-testperf-pkgdown` → PR): fixed a pkgdown gap M33
  introduced — `sim16` was added to `man/` but not to `_pkgdown.yml`’s
  reference index, so
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  failed on an undocumented topic; added it to the “Data” section (now
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  → “No problems found”). Also cached the shared polychoric `bfi25` fits
  in `test-vignette-m24.R` (`.vfit` memo, mirroring `test-suggest_k.R`’s
  `.get_sk`): 14
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  calls → 7 distinct fits, cutting that file from ~15.1s to ~8.0s (a
  general test-suite speedup surfaced while profiling for the workflow’s
  new efficiency guidance).

- **M34 (done):** pruning verb — extracted
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) as
  a standalone, pipeable S3 generic off
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
  the second milestone of the M31–M38 documentation/UX epic. **API
  extraction.** The five Forbes-extension pruning args (`prune` →
  renamed `rules`, `redundancy_r`, `redundancy_phi`, `min_items`,
  `orphan_r`) leave
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  entirely and live on
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md):
  `ackwards(...) |> prune(...)`. Clean move, no deprecation shim
  (pre-CRAN, no users) —
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  now explicitly rejects the five removed args if passed (rather than
  silently absorbing them via `...`, which would be a masked-argument
  footgun, not the intended clean break; Invariant 6). **Generic, not a
  plain function.** `prune <- function(x, ...) UseMethod("prune")` +
  [`prune.ackwards()`](https://jmgirard.github.io/ackwards/reference/prune.md),
  so it coexists with the `prune` generics already defined by
  recursive-partitioning packages
  (e.g. [`rpart::prune`](https://rdrr.io/pkg/rpart/man/prune.rpart.html))
  regardless of package load order. No new class or dispatch surface is
  needed beyond that:
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  returns the same `ackwards` object with `$prune` populated (replacing
  any prior pruning), so `print`/`summary`/`tidy`/`glance`/
  `augment`/`autoplot` all work unchanged. **Edges recomputed fresh,
  `x$edges` untouched.**
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  calls `compute_edges(levels = x$levels, R = x$r, pairs = "all")`
  internally on every call (cheap `W'RW` algebra, not re-extraction —
  DESIGN.md §3) and passes the result to the existing
  `.find_redundant_chains()`/`.compute_structural_signals()` helpers via
  a lightweight
  `list(k_max, levels, lineage, edges = list(matrices = ...))` view,
  rather than reading `x$edges` directly. This means: (1) `x$edges` is
  never mutated by pruning (Invariant 1: one edge path); (2)
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  produces identical results regardless of whether the object was fit
  with `pairs = "adjacent"` (the default) or `"all"` — including the
  endpoint-`r` enrichment for chains spanning more than one level, which
  needs skip-level edges the fit-time object may not carry.
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)‘s
  old `pairs` auto-upgrade-to-`"all"`-when-pruning behavior is removed
  along with it; `pairs` is now a pure display/storage setting on
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
  decoupled from pruning. **Manual + mixed pruning.** `manual =`
  (character vector of factor labels) flags nodes directly — standalone
  (`prune(x, manual = c("m4f3"))`, no auto rule needed) or unioned onto
  an auto rule (`prune(x, "redundant", manual = c(...))`). Unknown
  labels error, listing the valid ones. On overlap, the auto rule’s
  `prune_reason` wins (more informative) over `"manual"`. **Naming:
  canonical `"artifact"` (US spelling), `"artefact"` accepted as an
  alias** (owner preference; nod to Commonwealth spelling and to Forbes’
  own usage), normalized internally. Existing code passing `"artefact"`
  keeps working. `"tucker"` was considered and rejected as an
  alias/rename during planning: the mode surfaces more than Tucker’s φ
  (also the `few_items`/`orphan`/`split_merge` structural signals), so
  naming it after the statistic would mislabel the umbrella.
  **[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  gain a “Manual” pruning line** reporting the count and ids of
  explicitly-flagged nodes, independent of which auto rule (if any) was
  also requested. Files: `R/prune.R` (generic + method + refactored
  internal dispatcher — `.find_redundant_chains`, `.phi_pairs`,
  `.compute_structural_signals` unchanged in signature/behavior),
  `R/ackwards.R` (args + validation + auto-upgrade + auto-resolve blocks
  removed; new rejection guard added), `R/print.R`/`R/summary.R`
  (artifact spelling + Manual line), `R/autoplot.R`/`R/data.R` (roxygen
  examples updated to the piped form), `DESIGN.md` (§2, §5-6, §9, §10,
  §11, §14 items 27-31; item 17 struck through and superseded),
  `vignettes/ackwards-forbes.Rmd` (code chunks converted to the piped
  API; the thresholds-sweep chunk now re-prunes one shared fit instead
  of refitting four times, demonstrating the no-re-extraction design —
  deeper prose/formatting notes from `ROADMAP.md` remain M38’s to
  resolve), `ROADMAP.md` (M34 section deleted per its own maintenance
  rule; M38’s dependency note updated). Tests: `test-prune.R` gained
  coverage for standalone/mixed manual pruning, unknown-label errors,
  re-pruning without re-extraction (`x$levels`/`x$r`/`x$edges` identical
  across calls on the same fitted object), the `"artefact"` alias
  canonicalizing to `"artifact"`,
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  rejecting the removed args, and a regression guard confirming
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  recomputes all-pairs edges (including the endpoint-`r` enrichment)
  even when the object was fit with the default `pairs = "adjacent"`.
  Every other test file with a fit-time `ackwards(..., prune = ...)`
  call (`test-cor-input.R`, `test-data.R`, `test-vignette-m24.R`,
  `test-print.R`, `test-layout.R`) was migrated to the piped form. (1400
  tests pass, 2 skip; 0/0/0 R CMD check; coverage 100%.) Post-review
  (`/post-milestone-review`) follow-up, landed via branch
  `m34-followup-review` → PR: the review returned **READY** with no
  Blocking or Should-fix code findings; this follow-up cleared the
  Nice-to-haves plus one prose consistency item. (1) Normalized the two
  remaining internal `artefact` comments in `R/prune.R` and the two
  user-facing `artefact`/`artefactual` prose mentions in
  `README.Rmd`/`README.md` to the canonical `artifact` spelling (the
  alias input path is untouched — `"artefact"` still works). (2)
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  now de-duplicates `manual` labels
  ([`unique()`](https://rdrr.io/r/base/unique.html)) so `$prune$manual`
  is dup-free; flagging was already idempotent. (3) Added tests:
  clearing prior pruning (`prune(x)` on an already-pruned object →
  `$prune` NULL, class preserved), `manual` de-duplication + class
  preservation, and
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) on
  an ESEM object with a polychoric basis (exercises the polychoric-`R`
  path through
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  inside
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  and re-asserts `x$edges` is untouched). The README
  worked-[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)-example
  the review flagged as an *opportunity* remains M38’s
  (narrative/Quick-Start restructuring is M38’s declared scope; the
  corrected prose already describes the feature accurately). (1412 tests
  pass, 2 skip; 0/0/0 R CMD check; coverage 100%.)

- **M35 (done):** autoplot & visualization — the third milestone of the
  M31–M38 documentation/UX epic. Combines one correctness fix with an
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  encoding overhaul plus code-coupled prose. **Sign-propagation bugfix
  (the lead item).** `.align_signs()` chose each factor’s flip from the
  *raw* adjacent edge, ignoring the parent’s own already-applied flip
  (`signs[[k-1]]`). When a parent was itself flipped (e.g. `m2f2` on
  `bfi25`, whose sign-aligned loading colSum is negative so its edge to
  `m1f1` reads positive), every child whose primary parent was that node
  inherited the −1 on the display side — so the *primary* edge
  `m2f2 → m3f2` rendered −0.99. That violated DESIGN §7 (“orient each
  factor so its correlation with its **primary parent** is positive,
  *propagating top-down*”): the code propagated against the unflipped
  parent. Fix (one line in `R/utils.R`): compute the child flip against
  the parent’s *aligned* sign —
  `sk[j] <- if (signs[[k-1L]][parents[j]] * E[parents[j], j] >= 0) 1L else -1L`.
  Now every primary-parent edge is non-negative across pca/efa/esem;
  only genuinely-negative *secondary* edges are red (Invariant 4; DESIGN
  §7 “non-primary edges may legitimately be negative”). This is a
  conformance fix, not a resolved-default change — it restores the
  documented intent. Loading and edge signs shift wherever a parent was
  flipped, so README/vignette rendered outputs regenerated; no snapshot
  hardcoded the affected signs, so the ripple was mechanical.
  **Configurable, always-legended encodings.** New
  `sign_by = c("color","linetype","both","none")` (default `"color"`)
  and `magnitude_by = c("linewidth","none")` (default `"linewidth"`) let
  the user assign sign and magnitude to specific aesthetics; every
  mapped aesthetic now carries its own legend (the old color-mode
  `linetype` = strong/weak encoding was legend-less). `sign_by = "both"`
  maps sign to color *and* linetype, drawing negatives as a distinct
  `twodash` so the redundant pairing stays legible in greyscale, and
  shares the `"Direction"` scale name across both so ggplot2 merges them
  into one legend key. Implemented with
  [`utils::modifyList()`](https://rdrr.io/r/utils/modifyList.html) on a
  base `aes()` so only active channels are mapped; inactive channels
  pass to each geom as constants. **`cut_strong` retired**
  (soft-deprecated: still accepted, warns via cli, no effect) because it
  double-encoded magnitude via linetype — the source of the
  legend-less-linetype confusion. `mono = TRUE` kept as a thin
  convenience wrapper (= `sign_by = "linetype"` + black `color_edge`).
  **Layout orientation.** New `direction = c("vertical","horizontal")`;
  `"horizontal"` transposes the layered layout left-to-right (level 1 at
  left) for wide slides/posters. Edge face-attachment, node coordinates,
  and the level-axis labels (which move to the bottom margin) are all
  orientation-aware and compose with skip arcs, `mono`, `sign_by`, and
  the `drop_pruned` path. Fixes the intro/README text that mis-described
  the vertical layout as left-to-right. **British-spelling aliases.**
  `colour_pos`/`colour_neg`/`colour_edge`/`colour_pruned` accepted and
  normalized to the canonical American `color_*` (favoring `color`,
  mirroring `artifact`/`artefact`); new `color_edge` sets the single
  edge color when sign is not color-encoded. **`ggsave` documented, not
  re-exported** (re-export would drag `ggplot2` from Suggests into
  Imports — dependency guardrail): a roxygen “Saving plots” section plus
  a vignette section call
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  directly. **Code-coupled prose.** Visualization vignette rewritten
  around `sign_by`/`magnitude_by` (each legended), `cut_strong` section
  dropped, layout-orientation and Saving-plots sections added, the
  redundant `suggest_k` diagnostic section trimmed to a pointer, the
  “annotate with \|r\|” heading corrected to “r” (edges show signed r),
  `r_digits = 1L` → `r_digits = 1` (drop the integer-`L` convention
  readers may not know), and prose switched to American “color”. README
  step 3 and the intro vignette corrected (levels stack top-to-bottom,
  not left-to-right) with a `direction` note; README.md regenerated
  (`m2f2 → m3f2` now +0.99). Broader narrative rewrites remain M38’s.
  Files: `R/utils.R` (sign fix), `R/autoplot.R` (encoding args,
  orientation, aliases, ggsave roxygen),
  `tests/testthat/test-compute_edges.R` (`.align_signs` unit test +
  cross-engine primary-edge-non-negative regression),
  `tests/testthat/test-layout.R` (encoding/alias/mono/ direction tests;
  `cut_strong` deprecation test), `NEWS.md`, `README.Rmd`/`README.md`/
  `man/figures/README-plot-1.png`,
  `vignettes/ackwards-visualization.Rmd`,
  `vignettes/ackwards-intro.Rmd`, `DESIGN.md` (§7 note; §11 autoplot
  contract), `ROADMAP.md` (M35 section deleted per its own maintenance
  rule). (1447 tests pass, 2 skip; 0/0/0 R CMD check; coverage 100%.)
  Post-review (`/post-milestone-review`) follow-up, landed via branch
  `m35-followup-review` → PR: the review returned **READY** with no
  Blocking or Should-fix findings; this follow-up cleared the three
  Nice-to-haves. (1) Locked two visual behaviors that were previously
  only smoke-verified: a test asserting `sign_by = "both"` titles *both*
  the colour and linetype scales `"Direction"` (the mechanism by which
  ggplot2 merges them into one legend key), and a
  `direction = "horizontal"`

  - `show_r` composition test. (2) Broadened the encoding tests: a
    `sign_by = "none"` + `magnitude_by = "none"` case (asserts no
    colour/linetype/linewidth scales and that it still builds), and
    coverage of the remaining three British aliases — `colour_neg` (read
    off the manual colour scale’s palette, so it is checked even with no
    negative edge drawn), `colour_edge`, and `colour_pruned` (via
    deterministic `manual =` pruning). (3) Added a roxygen sentence
    explaining *why* colour is the default `sign_by` channel
    (pre-attentive sign reading; leaves linetype free). (1457 tests
    pass, 2 skip; 0/0/0 R CMD check; coverage 100%.)

- **M36 (done):** interpretation functions — the fourth milestone of the
  M31–M38 documentation/UX epic. Additive only: no new exported objects
  (pkgdown reference index unchanged), no new dependency (labels read
  via base [`attr()`](https://rdrr.io/r/base/attr.html) —
  **labelled**/**haven** never enter Imports/Suggests).
  **[`augment()`](https://generics.r-lib.org/reference/augment.html)
  scores-only output.** New `append` argument (default `TRUE` = current
  behaviour, scores appended to `data`); `append = FALSE` returns only
  the `.m{k}f{j}` score columns. New `id_cols` names passthrough columns
  (e.g. a subject id) to carry through under `append = FALSE` for a
  rejoin that survives filtering; `NULL` (default) yields bare scores.
  The assembly was refactored around a single `base` carrier (full
  `data` / `.obs` index / `id_cols` subset / empty frame) so row order
  and count are always preserved —
  `cbind(data, augment(x, data, append = FALSE))` reproduces the
  appended output exactly. Guards: `id_cols` requires `append = FALSE`
  **and** `data` (errors otherwise, since with `append = TRUE` every
  column is already kept and with `data = NULL` there are no source
  columns), an absent `id_cols` column errors, and `append` must be a
  scalar logical.
  **[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  grouping + variable labels.** New `by = c("factor", "item")`:
  `"factor"` (default) is the existing group-items-under-each-factor
  view; `"item"` inverts it to list, per item, the factors it loads on
  (cross-loading view). The row-builder was restructured to first
  assemble the full filtered long table (retaining
  `.factor_ord`/`.item_ord` index columns) and then apply per-group
  `sort` + `n` via [`order()`](https://rdrr.io/r/base/order.html) +
  `stats::ave(..., FUN = seq_along)`, so both orientations and the
  `sort = FALSE` “original order” mode work regardless of `by`. New
  label support:
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  captures each data.frame column’s `"label"` attribute (the attribute
  **labelled**/**haven** write) into `x$meta$item_labels` before
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) strips it (new
  internal `.capture_item_labels()`, `NULL` for matrix/cor-matrix input
  or when nothing is labelled);
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  then prints items as `label (code)` with a per-item bare-code fallback
  (new internal `.format_item_label()`). `show_labels = FALSE` forces
  bare codes; `$data` gains a `label` column when labels are available.
  The `top_items` object grew two fields (`by`, `item_labels`).
  **Vignettes.** Interpret vignette (`ackwards-interpret.Rmd`): lead
  example moved to `cut = 0.5` with the duplicate-cut chunk dropped and
  adjustability moved to prose; the `$data`-internals paragraph and the
  “why polychoric” setup comment removed; cross-loadings section re-cast
  around `by = "item"`; a new “Showing item wording” section
  demonstrates the label interface; a new “Borrowing names for the upper
  levels” subsection gives naming advice via Big Five metatraits
  (Stability/Plasticity) and HiTOP spectra (with DeYoung 2006 and Kotov
  et al. 2017 references); the “Where to go next” section demoted to a
  blockquote note. Intro vignette (`ackwards-intro.Rmd`): the positional
  `names(scored)[26:40]` / `scored[, c(...)]` score indexing replaced
  with `augment(append = FALSE)` (mechanical touch only; broader intro
  narrative remains M38’s). **Decisions banked.** `id_cols` subsumes
  bare-scores (it *is* the escape hatch, not a separate mode) rather
  than adding a passthrough-vs-bare toggle; labels captured at fit time
  (so display “just works” and travels with the object) rather than
  passed to
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  each call; print format is `label (code)` (code retained for
  cross-referencing `tidy(what = "loadings")`); `by = "item"` orders by
  level then strongest loading within level. Files: `R/ackwards.R`
  (label capture into meta), `R/utils.R` (`.capture_item_labels()`),
  `R/augment.R` (`append`/`id_cols`), `R/interpret.R`
  (`by`/`show_labels`/`.format_item_label()`),
  `tests/testthat/test-utils.R`, `tests/testthat/test-scores.R`,
  `tests/testthat/test-interpret.R`, `NEWS.md`, `DESIGN.md` (§10
  `top_items` signature; object-spec `meta$item_labels`),
  `vignettes/ackwards-interpret.Rmd`, `vignettes/ackwards-intro.Rmd`,
  `ROADMAP.md` (M36 section deleted per its own maintenance rule). (1496
  tests pass, 2 skip; 0/0/0 R CMD check; coverage 100%.) Post-review
  (`/post-milestone-review`) follow-up, landed via branch
  `m36-followup-review` → PR: the review returned **READY** with no
  Blocking or Should-fix findings; this follow-up cleared the three
  Nice-to-haves. (1) Locked the print *content* (previously only
  `expect_no_error`-smoke-tested) with
  [`cli::cli_fmt()`](https://cli.r-lib.org/reference/cli_fmt.html)
  capture: asserts
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  renders `label (code)` in the item **header** under `by = "item"` and
  in the item **body** under `by = "factor"`, and that
  `show_labels = FALSE` prints the bare code and never the label
  text. (2) Broadened the `augment(id_cols = )` coverage: a matrix input
  carrying an extra numeric `id` column (previously only data.frame was
  exercised), an `id_cols` that also names a model item (carried through
  verbatim *and* still scored), and an NA-scored row (missing an item)
  confirmed retained in place with its `id` and `NA` scores — the
  positional-rejoin guarantee. (3) Added a
  [`?augment`](https://generics.r-lib.org/reference/augment.html)
  roxygen sentence spelling out that `append = FALSE` returns NA-scored
  rows in place so positional alignment holds even with `NA` scores.
  Code change was documentation-only (the roxygen note); the rest are
  tests. (1509 tests pass, 2 skip; 0/0/0 R CMD check; coverage 100%.)

- **M37 (done):** engines vignette — the fifth milestone of the M31–M39
  documentation/UX epic (the epic was renumbered M31–M38 → M31–M39
  during M37 planning; see below). **Doc-only:** the entire change is
  `vignettes/ackwards-engines.Rmd` — no new/removed exports (pkgdown
  reference index unchanged), no `R/` changes, no new dependency (the
  new example uses
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html), and
  `psych` is already in Imports). **At-a-glance table.** Parallelised
  the EFA/ESEM “what it models” rows (both now “Common (latent)
  variance”); added *engine-substrate*
  ([`psych::principal()`](https://rdrr.io/pkg/psych/man/principal.html)
  / [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) / `lavaan`),
  *correlations*, and *estimators* rows; rendered χ² as a proper symbol
  (was “chi”); dropped the confusing “(WLSMV for ordinal)” parenthetical
  from the Loading-SEs row (now just “Yes”). **ESEM reframed as a
  general engine.** The section now states ESEM fits continuous items
  (ML/MLR, the default, + FIML for missing data) *and* ordinal items
  (WLSMV), correcting the prior impression that ESEM ⇒ ordinal-only; the
  two-capability list grew to three (SEs; full ML incl. FIML; WLSMV).
  **Per-level fit honesty.** Added the converse of the existing “good
  fit doesn’t bless the edges” point: a poorly-fitting level weakens the
  edges *incident to it* — its factors are less well defined — even
  though the edge correlation itself is computed faithfully; edges
  between two well-fitting levels are on firmer ground. **autoplot depth
  note.** Added a paragraph that the running examples deliberately stop
  at `k_max = 3` for build speed/legibility, that
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  points to ~5 for the BFI, and that the truncated fit plot shows the
  *mechanics* of reading a trajectory, not the recommended depth.
  **sim16 vs bfi25 framing.** In the cross-engine convergence
  discussion, added a note that clean convergence is the best case
  (idealized `sim16`), not the norm (bfi25’s
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  criteria span 4–6), with a pointer to the suggest_k vignette for the
  full contrast. **corFiml MAR route (the substantive enrichment).** New
  runnable example showing how to bring FIML-based missing-data handling
  into the PCA/EFA engines via the M22 correlation-matrix seam:
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
  estimates the correlation matrix, which is passed to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  Demoed on the **continuous** `sim16` (with injected NAs) —
  deliberately *not* the ordinal `bfi25`, since corFiml assumes
  multivariate normality. Two caveats carried explicitly: (1) `n_obs` is
  the user’s call — total N is mildly anti-conservative for the EFA fit
  indices, complete-case N conservative, and the loading/edge **point
  estimates are unaffected by `n_obs`**; (2) MVN → continuous-only,
  ordinal items should use `engine = "esem"` + `cor = "polychoric"`.
  Added a missing-data “which option” table row and a forward-reference
  to M38 (which will promote this to a first-class `missing = "fiml"`
  route for PCA/EFA); the Correlation-matrix-input section links to it.
  **Trims.** The Missing-data and Performance sections were compressed
  to a lead paragraph + summary, demoting the fine per-engine semantics
  to
  [`?ackwards`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  The parallel-fit chunk switched to
  [`library(future); plan(multisession)`](https://future.futureverse.org)
  style with a note that `plan()` comes from **future** and is *not*
  re-exported by `future.apply`. **Citations.** Added the missing Hu &
  Bentler (1999) reference (cited in-text for the fit cutoffs; the
  `cutoffs` arg was kept in M32) and alphabetised the References list.
  **Epic renumber (planning-time decision).** During M37 planning the
  owner elected to slot a new *code* milestone — **M38**
  `missing = "fiml"` for PCA/EFA (auto-route via
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)) —
  between the engines vignette and the prose milestone; the former “M38
  — Narrative & remaining prose” became **M39**. `cor = "fiml"` was
  rejected as a category error (corFiml returns a Pearson matrix; FIML
  belongs on the `missing=` axis). Both are logged in `ROADMAP.md`.
  Files: `vignettes/ackwards-engines.Rmd`, `NEWS.md`, `CLAUDE.md`
  (Current focus + Completed index), `ROADMAP.md` (M38 insert + M39
  renumber + corFiml speed note), `MILESTONES.md`. (1509 tests pass, 2
  skip; 0/0/0 R CMD check; coverage 100%.) Post-review
  (`/post-milestone-review`) follow-up, landed via branch
  `m37-followup-review` → PR: the review returned **NOT READY** on one
  Should-fix and its clean-check re-run surfaced a pre-existing blocking
  bug; this follow-up cleared both. (1) **EFA “Estimators” table cell
  corrected.**
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  restricts `fm` via `arg_match(fm, c("minres", "ml", "pa"))`, but the
  new M37 table cell read “OLS / minres / ML” — naming the invalid
  `fm = "ols"` (which errors) and omitting principal-axis (`pa`). Now
  reads `minres` (OLS) / `ml` / `pa` via psych’s `fm=`.

  2.  **Pre-existing
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
      crash fixed** (`R/suggest_k.R:352`). The CD-unavailable branch
      called `cli::cli_inform("i" = "…")` — a bare named bullet
      argument, which leaves `cli_inform()`’s `message` parameter
      missing and errors with *“argument "message" is missing, with no
      default”*. So
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
      crashed on **any machine without `EFAtools`** (the branch is
      reached only when EFAtools is absent), which in turn broke the
      `intro`, `suggest-k`, and `visualization` vignette builds (all
      call
      [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md))
      and would fail CRAN’s no-Suggests check. It was masked at the M37
      gate because the dev machine has EFAtools installed (the covering
      test at `test-suggest_k.R` L518
      `skip_if(is_installed("EFAtools"))` skips there). Fix: wrap the
      bullet in `c(...)`. A multiline-aware scan confirmed this was the
      only bare-bullet `cli_*()` call in `R/`; coverage is unaffected
      (the line sits inside an existing `# nocov` block). Not an
      M37-introduced bug, but M37’s clean re-check is what exposed it.
      Verified without EFAtools: `R CMD check` **Status OK** (0/0/0,
      `_R_CHECK_FORCE_SUGGESTS_=false`); all vignettes re-build OK; full
      suite **1484 pass / 0 fail / 0 error / 8 skip**. (The “1509 pass,
      2 skip” figure above is the EFAtools-*installed* count; the delta
      is the EFAtools-gated tests — the suite total is
      environment-dependent on which Suggests are present.)

- **M38 (done):** `missing = "fiml"` for PCA/EFA — the sixth milestone
  of the M31–M39 epic and its first **code** milestone. Promotes FIML
  from an ESEM-only option to a first-class route for the PCA and EFA
  engines, reversing the M16 contract that `missing = "fiml"` errors for
  those engines. **The route.** Under `engine = "pca"/"efa"` with
  `cor = "pearson"`, `missing = "fiml"` now estimates the correlation
  matrix via
  [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
  (full-information ML; multivariate-normal, MAR-valid) and feeds it to
  the existing `W'RW` between-level algebra. Invariant-1-clean: it just
  supplies a better `R` to the one edge path — no new edge path, no new
  dependency (`psych` already Imports), and one `corFiml()` call per run
  by construction (PCA/EFA compute `R` once and reuse it across all
  levels, unlike the per-level ESEM fits M26 had to hoist). Grew out of
  the M37 doc-planning observation that a user could already smuggle
  FIML in through the M22 correlation-matrix seam
  (`ackwards(psych::corFiml(x), …)`); M38 makes it discoverable and owns
  the `n_obs` tradeoff. New internal helper `.corfiml_R()` wraps the
  call with the same non-PD
  [`psych::cor.smooth()`](https://rdrr.io/pkg/psych/man/cor.smooth.html)
  fallback the polychoric path uses. **Guard matrix.**
  `.resolve_missing()` gained a `cor` argument: FIML is valid for
  `{pca, efa}`-pearson and (unchanged) ESEM-ML/MLR; it **errors** for a
  non-Pearson PCA/EFA basis (Spearman, polychoric — corFiml is MVN-only)
  and for WLSMV/ULSMV. A `cor = "fiml"` alternative was rejected during
  M37 planning as a category error (corFiml returns a Pearson matrix;
  FIML is an estimation-under-missingness method that belongs on the
  `missing=` axis, not the measurement-level `cor=` axis). **`n_obs`
  string option.** On the raw-data FIML PCA/EFA path, `n_obs` accepts
  `"total"` (default — every row contributing to the FIML likelihood,
  matching the convention a genuine FIML analysis reports, Enders 2010,
  and giving cross-engine parallelism with ESEM-ML FIML) or `"complete"`
  (complete-case N, a conservative lower bound). Point estimates
  (loadings/edges) are unaffected by the choice; only the EFA fit
  indices depend on N, and those are *approximate* under this two-step
  (FIML matrix → normal-theory EFA) procedure regardless of N (Zhang &
  Savalei 2020), so the route announces itself and the caveat via cli
  (Invariant 6). `"effective"` was considered and **dropped** at plan
  time: no canonical formula exists for a corFiml→EFA route, so it would
  be a package-invented convention masquerading as a standard. A string
  `n_obs` is rejected off this path; correlation-matrix input still
  requires a numeric N. Score materialisation is unchanged — corFiml
  estimates `R` but does not impute item responses, so
  `keep_scores = TRUE` still yields `NA` score rows for incomplete cases
  (the existing `.compute_scores()` note already covers this
  generically). **DESIGN sign-off.** §9 (`missing` + `n_obs` rows) and
  §14 (new resolved items 32–33) updated to record the reversed default
  and the `n_obs` semantics; roxygen for `missing`/`n_obs` gained the
  *why* plus Enders (2010) and Zhang & Savalei (2020) references. No
  new/removed exports (helpers are internal; pkgdown reference index
  unchanged). Files: `R/utils.R` (`.resolve_missing()` + `cor` guard,
  `.corfiml_R()`), `R/ackwards.R` (route, `n_obs` string handling, cli
  announce, roxygen), `tests/testthat/test-missing.R`, `DESIGN.md`,
  `NEWS.md`, `CLAUDE.md` (Current focus + Completed index),
  `MILESTONES.md`. Verified: `R CMD check` **0/0/0**; full suite **1509
  pass / 0 fail / 8 skip** (EFAtools absent on the dev machine — 8
  EFAtools/doc-context skips; higher with EFAtools installed); `R/`
  coverage **100%** (`R/ackwards.R` and `R/utils.R` both 100%; the only
  package-wide uncovered lines are the EFAtools-gated `suggest_k.R`
  block, which is covered when EFAtools is present); lint clean;
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  clean. (Two pre-existing defensive early-returns in `.is_cor_matrix()`
  / `.check_maybe_cov_matrix()` were given direct unit tests here to
  keep `R/utils.R` at 100% after the M38 additions.) Post-review
  (`/post-milestone-review`) follow-up, landed via branch
  `m38-followup-review` → PR: the review returned **READY** with no
  Blocking/Should-fix items; this cleared the two nice-to-haves.

  1.  **Cross-check exclusion note** (`DESIGN.md` §14): the
      algebra-vs-scores oracle now documents that it does not certify
      the `missing = "fiml"` PCA/EFA path either (not just
      `cor = "polychoric"`) — the algebra uses the `corFiml()` matrix
      while the scores route standardizes raw NA-bearing data, so the
      bases diverge under missingness by design; the oracle tests run
      only on complete-data linear engines, so there is no false-failure
      risk. (2) **`pca` + `fiml` + `n_obs = "complete"` test** added
      (`test-missing.R`) confirming the effective-N `switch` is
      engine-agnostic and that the R matrix/edges are identical across
      `n_obs` choices for PCA (whose fit does not use N). The optional
      advisory-on-complete-data idea was deliberately **not**
      implemented — it would fire noise on legitimate always-FIML
      workflows for no correctness benefit.
