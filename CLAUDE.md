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

## Current focus

**M21 (in progress):** Onboarding & usability pass (pre-CRAN; 0.1.0 still unreleased). A
usability/polish milestone, **not** a DESIGN.md §15 feature milestone. Touches nothing on the
out-of-scope list. Correlation-matrix input was raised but **deferred to M22** (large; PCA/EFA-only,
needs an `n_obs` arg + R-detection + engine gating). Four parts, order A → B → C → D:

**(A) `psych` → Imports; drop `GPArotation` entirely.** The default PCA engine currently sits behind
`rlang::check_installed("psych")` (`engine_pca.R:11`, `engine_efa.R:11`, the polychoric path, and
`suggest_k`) — core functionality behind an install prompt. Move `psych` Suggests→Imports and remove
every `check_installed("psych")` guard; **keep** the `lavaan`/`ggplot2`/`EFAtools` guards. `GPArotation`
is unused — varimax (orthogonal) routes through base `stats::varimax`; verified it never enters
`loadedNamespaces()` on the `psych::pca`/`psych::fa` varimax path — so **remove it from Suggests** and
delete the two vestigial `skip_if_not_installed("GPArotation")` lines (`test-interpret.R:85`,
`test-label_template.R:123`). New `Imports`: `cli, generics, psych, rlang, stats, utils`. Update
DESIGN.md §3/§12 and the CLAUDE.md Dependencies section. **This reverses the documented lean-`Imports`
stance for `psych` only — owner-approved (psych is the engine substrate for the default path; the SEM
and plotting stacks stay opt-in).**
  *Acceptance:* `ackwards(bfi25, k_max = 5)` runs with no install prompt when
  lavaan/ggplot2/EFAtools/GPArotation are absent but psych is present; `grep check_installed R/` shows
  no psych guard while lavaan/ggplot2/EFAtools guards remain; `R CMD check` 0/0/0.

**(B) Bundle `bfi25` example dataset.** `data-raw/bfi25.R` samples ~1000 rows from
`psych::bfi[, 1:25]` with a fixed seed, **preserving the natural NAs** (so `test-missing.R` and the
ordinal-detection warning still have data to bite on) → `data/bfi25.rda` (xz). `R/data.R` documents
`@format`, `@source` (Revelle / psych / SAPA / IPIP; items are public-domain IPIP), and `@details`
explaining the derivation (seed, n, first 25 columns). `DESCRIPTION` gains `LazyData: true`;
`.Rbuildignore` ignores `^data-raw$`. Migrate **all** `@examples` and **all** vignettes from
`psych::bfi[, 1:25]` → `bfi25`, dropping psych data-guards; **regenerate any numeric prose against
bfi25's actual output** (esp. the suggest-k worked example — re-verify like the M14 post-review pass).
**Exception:** the PCA-vs-`psych::bassAckward()` fidelity snapshot test stays on full
`psych::bfi[, 1:25]` (canonical oracle).
  *Acceptance:* `dim(bfi25)` is ~1000 × 25; `?bfi25` documents format/source/derivation;
  examples + vignettes reference `bfi25`; the oracle snapshot is unchanged and passing; `R CMD check`
  clean (no data-compression NOTE).

**(C) README "Learn more" completeness.** README.Rmd lists only 5 of 7 vignettes. Add the missing
**Visualization** and **Interpreting & labeling factors** rows; `devtools::build_readme()`.
  *Acceptance:* all 7 vignettes listed; links resolve to the pkgdown article paths.

**(D) CD panel in `autoplot.suggest_k()`.** `EFAtools::CD()` returns `RMSE_eigenvalues`
(`N_samples × n_factors_max`); its column means are CD's native RMSE-vs-k curve. Store
`cd_rmse = colMeans(cd_out$RMSE_eigenvalues)` (length `n_factors_max`; `NULL` when CD unavailable) in
the `suggest_k` object. Add a 4th panel **"CD (RMSE, minimize)"** plotting `cd_rmse` vs k with `k_cd`
marked, **only when `cd_available`**, and move the CD vline out of the MAP panel into it. Facet **2×2**
(`facet_wrap(ncol = 2)`) when CD is present; keep the current single-column layout when CD is absent
(no dangling empty cell). Tests for the 4-panel / 3-panel branches and the `cd_rmse` field. DESIGN.md
§15 gets the M21 entry + a §8/§12 CD-panel note; NEWS.md bullets land under the **existing unreleased
0.1.0 entry**; refresh the CD figure in the suggest-k or visualization vignette.
  *Acceptance:* with EFAtools present, `autoplot(suggest_k(bfi25, k_max = 8))` renders a 2×2 grid
  including the CD RMSE panel with `k_cd` marked; without EFAtools, 3 panels render with no error; the
  object carries `cd_rmse`.

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
`FALSE`. Don't change these silently.

## Dependencies (see `DESIGN.md` §12)

`psych` is in **Imports** (M21) — it is the engine substrate for the default PCA and EFA paths and
for polychoric correlations; placing it in Suggests would require an install prompt for core
functionality. The SEM + plotting + optional-criterion stacks remain in `Suggests`. `GPArotation`
was **removed entirely** (M21) — varimax routes through base `stats::varimax` and GPArotation never
enters `loadedNamespaces()` on any supported path. **Do not add further to `Imports` without
flagging it.** **No Rcpp** — profile first; the heavy compute already lives in compiled deps (§3).

Current `Imports`: `cli`, `generics`, `psych`, `rlang`, `stats`, `utils`.
Current `Suggests`: `covr`, `EFAtools`, `ggplot2`, `knitr`, `lavaan (>= 0.6-13)`,
`rmarkdown`, `testthat (>= 3.0.0)`. `suggest_k()` uses `psych::fa.parallel(fa="both")` +
`psych::vss` (PA-PC, PA-FA, MAP, VSS-1/2) and optionally `EFAtools::CD()` (gated by
`rlang::is_installed()`); no separate `EGAnet`/`paran` dep. Visualization uses `ggplot2` directly
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
- **Oblique rotation** — varimax is hardcoded; no `rotation` argument; oblique confounds the cross-level signal. No plans to add it.
- **Higher-order SEM / Schmid-Leiman** — out of scope per §2; `ackwards` is correlation-based, not SEM-based.
- **M5 deferred improvements** — structural artefact signals (split-then-merge, <3-item factors,
  orphans), φ as default for non-PCA redundancy, bootstrap CIs on skip-level edges. Logged in
  `DESIGN.md` §14.
