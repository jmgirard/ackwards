# ackwards (development)

## `bfi25` ships item labels; cleaner console output

- The `bfi25` dataset now carries each item's public-domain IPIP stem (Goldberg,
  1999) as a `label` attribute, so `ackwards()` captures it at fit time and
  `top_items()` prints the wording as `code: label` (e.g.
  `E4: Make friends easily`) with no setup. Fit the dataset directly to keep
  the labels: base row-subsetting such as `na.omit(bfi25)` drops plain
  attributes, so use `missing = "listwise"` (or another `missing` option) rather
  than pre-filtering rows.
- `suggest_k()` now warns when the input looks ordinal, matching `ackwards()`
  and `comparability()`. Because `suggest_k()` screens on the Pearson/Spearman
  basis by design, the advice points at the final `ackwards()` fit
  (`cor = "polychoric"`), not at `suggest_k()` itself.
- Console output is more consistent: the engine name renders in lowercase
  everywhere (including `print(comparability(...))`); `summary()` per-level
  percentages use a fixed number of decimals with trailing zeros and blank
  lines separate the level blocks; `top_items()` separates its factor/item
  groups; and a pruned object's `print()`/`summary()` footer is a single
  consolidated note instead of stacked rules.

## Screen items before analysis — `check_items()`

- New `check_items()` verb: a per-item pre-analysis screen that flags the
  columns that break or *silently* degrade a factor analysis — constant (zero
  variance) items, near-constant items (one response category dominates), sparse
  response categories (a polychoric-basis concern), and heavy missingness — with
  one row per item and a printed summary plus guidance. It only reports; it
  never changes your data.
- `ackwards()` runs the same screen internally: it now **errors** on a constant
  item (naming it) instead of letting `psych::polychoric()` silently delete it
  and crash later with `subscript out of bounds`, and **warns** on a
  near-constant item (which could previously produce a plausible-looking but
  meaningless factor with no warning at all).

## Polychoric continuity correction (`correct`)

- `ackwards()` gains a `correct` argument (default `0.5`, psych's own default),
  forwarded to `psych::polychoric()` on the PCA/EFA polychoric path. Set
  `correct = 0` when `psych::polychoric()` fails on real-world ordinal data —
  its own error suggests exactly this, and it typically happens when an item has
  a near-empty (singleton) response category, or items with unequal category
  counts produce a sparse cross-cell. Previously there was no way to pass it, so
  such data could not be fit without pre-computing the matrix by hand.
- The polychoric failure message now names this remedy (`correct = 0`, collapse
  rare categories, or drop the item) instead of passing through psych's opaque
  error, and a polychoric matrix that comes back with undefined (`NA`/`NaN`)
  entries is caught with a clear, item-naming error rather than a base-R
  "missing value where TRUE/FALSE needed" crash.
- When the polychoric matrix is **near-singular** (rank-deficient), `ackwards()`
  now raises a single clear warning instead of letting `psych::fa()` flood the
  console with per-level "determinant … zero" notes, and records it durably:
  `meta$near_singular` / `meta$min_eigenvalue` are set, and `print()`/`summary()`
  re-surface a caution that per-level fit indices (especially CFI) and factor
  scores may be unreliable — so the signal survives on a saved or shared object,
  not just at fit time.
- `?ackwards` gains a **"When to trust the result"** section that tiers every
  diagnostic the function raises as fatal / caution / informational, so it is
  clear which warnings mean "fix before trusting" versus "proceed with care".

## Bootstrap confidence intervals on edges

- New `boot_edges()` verb: nonparametric bootstrap standard errors and
  percentile confidence intervals for every between-level edge. `ackwards()`
  reports each edge as a point estimate; `boot_edges(x, data)` quantifies its
  sampling uncertainty, which matters most where a hard threshold consumes the
  estimate — `prune()`'s `|r| >= redundancy_r` redundancy rule and the Forbes
  (2023) practice of reading off the strongest all-pairs edge. Each replicate
  resamples rows with replacement, recomputes the correlation matrix with the
  fit's basis and missing-data routine, refits the hierarchy, and — crucially —
  **anchors** every replicate level to the full-sample solution (the same
  matching + sign-orientation machinery `comparability()` uses) before its
  edges are computed, so factor label-switching and sign-flipping across
  replicates do not corrupt the pooled distributions. All resample indices are
  drawn upfront from `seed`, so results are reproducible and identical whether
  replicates run serially or in parallel (via `future.apply` under a
  `future::plan()`). Supports the PCA and EFA engines on the Pearson or
  Spearman basis; correlation-matrix, ESEM, and polychoric objects are not
  supported (each errors with a pointer).
- After `boot_edges()`, `tidy(x, what = "edges")` gains `se`, `lo`, `hi`, and
  `n_boot_ok` columns, and `print(x)` / `summary(x)` note the interval
  coverage. The intervals are per-edge error bars: they make sampling
  uncertainty visible but do **not** correct the selection bias of scanning
  many edges for the strongest one (documented in the object and the girard
  vignette).

## Split-half factor comparability (replicability gate)

- New `comparability()` verb: Everett's (1983) factor comparability
  coefficients, extended to every level of a bass-ackwards hierarchy. For
  each of `n_splits` random split-halves, solutions at every level are fit
  independently in each half; each half-solution's factors are matched to the
  full-sample solution's factors (so results are reported under the same
  `m{k}f{j}` labels as `ackwards()`), and the coefficient for a factor is the
  correlation between its two matched half-solution scores — both halves'
  scoring weights applied to the full sample, exactly the gate Goldberg's own
  lab applied to its hierarchies (Goldberg, 1990). Tucker's congruence of the
  matched loading columns is reported alongside. Everything is reported and
  nothing is auto-flagged; conventional .90/.95 benchmarks appear only as
  visual guides. Supports the PCA and EFA engines on the Pearson or Spearman
  basis; reproducible via `seed`.
- New `print.comparability()` (per-level median/min with the weakest factor)
  and `autoplot.comparability()` (per-split and median coefficients for both
  statistics, with benchmark reference lines).
- New vignette `vignette("ackwards-girard")` — *Replicability-Gated
  Hierarchies: A Recommended Workflow* — laying out the recommended
  end-to-end analysis: screen the depth range (`suggest_k()`), gate depth on
  replicability (`comparability()`), fit (`ackwards()`), separate persistence
  from differentiation (`prune()`), interpret with guardrails, and validate
  downstream out of sample (`predict()`). Includes a common-mistakes section.
- `comparability` objects record `k_requested` alongside the effective
  `k_max`, and `print()` notes when the full-sample fit truncated below the
  request. Cross-solution matrices carrying `NA` cells (pathological pairwise
  missingness) now degrade to `NA` coefficients for the level instead of
  erroring mid-run.

## Stricter argument checking

- `ackwards()`, `suggest_k()`, and `comparability()` now error on unknown
  arguments passed through `...`. Previously a misspelled argument — e.g.
  `kmax = 6` or `nsplits = 20` — was silently absorbed and the function ran
  with the default instead. (`ackwards()`'s pointer for the arguments that
  moved to `prune()` is unchanged.)

## Out-of-sample scoring (cross-validation)

- New `predict.ackwards()` method: apply a fitted hierarchy to new
  observations — e.g. the held-out test split of a cross-validation design —
  without retraining. `predict(x, newdata)` returns the factor-score columns
  for every level; it is the idiomatic front door to
  `augment(x, data = newdata, append = FALSE)`.
- `augment()` gains a `scaling` argument controlling which item means/SDs
  standardize the supplied data before the stored weights are applied.
  The new default, `scaling = "fit"`, uses the **fit-time (training) moments**
  — now stored in every raw-data fit — so scores for new observations,
  subsets, and the training data itself all share the training metric.
  `scaling = "sample"` retains the previous behavior (standardize by the
  supplied data's own moments) and is the option for objects fit from a
  correlation matrix, which carry no raw-data moments. Supplying `scaling`
  without `data` is an error (stored scores are returned exactly as computed
  at fit time).
- **Behavior change:** under the new default, `augment(x, data =)` on a
  *subset* of the training data (or any new data) produces slightly different
  values than before — previously the supplied data's own means/SDs were used,
  which made a score depend on which other rows were scored alongside it.
  Scoring the full training set is unchanged. Pass `scaling = "sample"` for
  the old behavior.

## Validation

- The package's fidelity to Forbes's (2023) extended bass-ackward method is
  now test-backed: a fixture test reproduces the paper's three simulation
  studies against expected values computed with Forbes's own reference
  implementation (public OSF materials), covering the between-level
  correlations, loading congruences, redundancy chase paths, and `prune()`'s
  retention decisions.

## Statistical correctness

- **EFA fit rows now report a consistent chi-square/p-value pair.** The `chi`
  statistic for `engine = "efa"` is now the likelihood-ratio chi-square
  (`psych::fa()`'s `STATISTIC`) — the statistic that `p_value`, `RMSEA`, and
  `TLI` are all derived from. Previously `chi` was psych's residual-based
  *empirical* chi-square while `p_value` belonged to the likelihood statistic,
  so quoting the pair (e.g. from `tidy(what = "fit")` or `summary()`)
  misreported. Reported `chi` values will change (they are typically larger);
  `dof`, `p_value`, `RMSEA`, `TLI`, and `BIC` are unchanged.

## Visualization

- Fixed a regression (introduced when pruning became a standalone `prune()`
  verb) where `autoplot(drop_pruned = TRUE)` on an object fitted with the
  default `pairs = "adjacent"` lost every edge into a node below a fully-pruned
  level. The pruned view now recomputes all-pairs edges internally, so it
  renders identically regardless of the fit-time `pairs` setting.
- The `autoplot(what = "fit")` caption now names only the Hu & Bentler
  thresholds actually plotted (EFA panels show TLI/RMSEA only, so CFI/SRMR
  cutoffs are no longer listed there).

## Usability

- The ordinal-detection warning now **names the flagged columns** (truncated
  past eight) instead of reporting only that "one or more columns" look
  ordinal, so the `cor = "polychoric"` advice is actionable for mixed data.
- `print()` on a `suggest_k` object now reports "Consensus: undetermined" when
  every requested criterion returned `NA` (previously it warned and printed an
  `Inf` range).
- `suggest_k()` now announces when a parallel-analysis suggestion exceeds the
  evaluated ceiling (`k_max`) instead of silently reporting the capped value as
  the recommendation.
- `ackwards()` and `autoplot()` now validate `cut_show` (a single number in
  [0, 1]) and `suggest_k()` validates `n_iter` (a single positive integer),
  erroring informatively instead of failing deep inside downstream code.

- `autoplot()` now italicises the axis label of a **fully-pruned level** (one
  whose factors are *all* flagged by `prune()`), denoting its status alongside
  the existing grey node fill. Partially-pruned levels keep a plain label. This
  applies in the normal diagram; under `drop_pruned = TRUE` a fully-pruned
  level's nodes are removed entirely.

## Documentation corrections

- The engines vignette's missing-data section now documents the built-in
  `missing = "fiml"` route for PCA/EFA (it previously described FIML as
  ESEM-only and pointed to a manual `psych::corFiml()` workaround), and its
  engine-choice table no longer recommends the invalid `fm = "pca"` setting.
- The "Choosing k" vignette now describes the Comparison Data method's
  mechanism correctly — comparison datasets reproduce the observed correlation
  structure under a k-factor model as well as the items' marginal
  distributions — and its worked BFI recommendation derives every number from
  the fitted object, so the prose always matches the printed table even though
  parallel-analysis results vary across builds.
- The Forbes vignette's artifact section was rewritten around its actual
  report-and-judge semantics: `prune(x, "artifact")` never flags factors
  automatically, so the section now demonstrates reading the Tucker's-phi
  table (`x$prune$phi`) instead of presenting a by-construction zero count as
  an empirical finding. The redundancy-chain example now states the full
  retention rule (a chain reaching the deepest level keeps only its bottom
  node), and a stale reference to the retired `cut_strong` argument was
  removed.
- The rationale for the `redundancy_phi` default was corrected in the
  documentation: PCA needs no congruence guard because component scores are
  *determinate*, not because the between-level algebra is exact (the algebra
  is equally exact for tenBerge-scored EFA).

## Ordinal vignette

- The "Ordinal Data" vignette now compares Pearson and polychoric item
  correlations with a dodged bar chart (the ten `N1`–`N5` item pairs, coloured
  by basis) instead of two raw correlation-matrix printouts, which were hard to
  compare visually.

## FIML for PCA and EFA

- `missing = "fiml"` is now a first-class option for `engine = "pca"` and
  `"efa"` (previously it errored for those engines). On the Pearson basis it
  estimates the correlation matrix via `psych::corFiml()` — full-information
  maximum likelihood, MAR-valid — and feeds it to the usual between-level
  algebra. This promotes the previously manual `ackwards(psych::corFiml(x), …)`
  pattern to a built-in route. It requires `cor = "pearson"` (corFiml estimates
  a multivariate-normal matrix) and still errors for a non-Pearson PCA/EFA basis
  and for WLSMV/ULSMV.
- New string form of `n_obs` for this path: `"total"` (default — every row
  contributing to the FIML likelihood, the standard FIML convention) or
  `"complete"` (complete-case N, a conservative lower bound). This governs only
  the EFA fit indices, which are approximate under this two-step route; the
  loadings and edges are unaffected by the choice. The route and this caveat are
  announced via a message.

## Bug fixes

- `suggest_k()` no longer errors on machines without the optional **EFAtools**
  package. The "CD requires EFAtools" info message was malformed and raised
  *"argument \"message\" is missing"* instead of informing, crashing any
  `suggest_k()` call (and the vignettes that use it) when EFAtools was absent.

## Engines vignette overhaul

- The "Choosing an Engine" vignette now frames ESEM as a general engine for
  **continuous** (ML/MLR, including FIML for missing data) *and* ordinal (WLSMV)
  data, rather than implying it is ordinal-only. The at-a-glance comparison table
  gains engine-substrate, correlation-type, and estimator rows.
- New guidance on getting **MAR-valid missing-data handling into the PCA and EFA
  engines**: estimate a FIML correlation matrix with `psych::corFiml()` and pass
  it to `ackwards()`, with caveats on the `n_obs` choice and the multivariate-
  normality (continuous-only) assumption.
- Clarifies that a poorly-fitting level weakens the edges *incident to it* (its
  factors are less well defined) even though the edge correlation itself is
  computed faithfully; and that the running examples stop at `k_max = 3` for build
  speed, not because the BFI hierarchy ends there.

## Vignette and README prose pass

- Documentation-only clarity pass across the introduction, choosing-k, ordinal,
  Forbes, and README pages: clearer print calls (`print(sk)`/`print(x)` so
  objects aren't buried in output), an explanation of why bass-ackwards uses
  orthogonal varimax rotation (and that it equals the "CF-VARIMAX" some papers
  report), a note that loading standard errors and confidence intervals are `NA`
  under the PCA/EFA engines and populated only by ESEM, and reference lists put
  in alphabetical order.
- The choosing-k vignette gains a side-by-side contrast of an idealized dataset
  (`sim16`, where the criteria agree) against a realistic one (`bfi25`, where
  they span a band), so the clean case is not mistaken for the norm.
- The ordinal vignette now states plainly that ordinal factor scores are
  trustworthy for downstream use (the standardization caveat is about scaling,
  not validity), that a two-category polychoric correlation *is* the tetrachoric
  correlation (no separate step), and that WLSMV genuinely operates on the
  polychoric basis.

## `augment()` scores-only output

- New `append` argument: `augment(x, data, append = FALSE)` returns only the
  `.m{k}f{j}` score columns, convenient for feeding scores straight into `cor()`,
  `lm()`, or a clustering call. The default `append = TRUE` is unchanged
  (scores appended to `data`).
- New `id_cols` argument names columns of `data` (e.g. a subject id) to carry
  through alongside the scores when `append = FALSE`, so scores can be rejoined
  after filtering. Row order and count are always preserved.

## `top_items()`: item grouping and variable labels

- New `by = "item"` inverts the listing to show, for each item, the factors it
  loads on — making cross-loadings easy to read. The default `by = "factor"` is
  unchanged.
- Variable labels are now supported: when the data carries a `"label"` column
  attribute (as **labelled** and **haven** set) at fit time, `ackwards()`
  captures it and `top_items()` displays items as `code: label`, with a
  per-item fallback to the bare code. Set `show_labels = FALSE` to force codes.

## Sign alignment: primary-parent edges are now always non-negative

`ackwards()` sign-aligns each factor to its primary parent. Previously a factor
whose *own* parent had been flipped could still display a **negative** primary
edge (e.g. `m2f2 -> m3f2` on `bfi25`), because the flip was computed against the
parent's unflipped orientation. Sign now propagates correctly top-down, so every
primary-parent edge is non-negative and only genuinely-negative *secondary* edges
appear red. Loading and edge signs may change wherever a parent was flipped.

## `autoplot()`: configurable, always-legended edge encodings

- New `sign_by` chooses how edge sign is shown: `"color"` (default),
  `"linetype"`, `"both"` (redundant colour + a distinct double-dash for
  negatives), or `"none"`. New `magnitude_by` (`"linewidth"` default, or
  `"none"`) controls the `|r|` encoding. Every mapped aesthetic now gets its own
  legend — nothing is silently encoded.
- New `direction = "horizontal"` draws the hierarchy left-to-right (level 1 at
  left) for wide slides and posters; the default `"vertical"` is unchanged.
- Edge-colour arguments accept British-spelling aliases (`colour_pos`,
  `colour_neg`, `colour_edge`, `colour_pruned`); a new `color_edge` sets the
  single edge colour when sign is not colour-encoded.
- `cut_strong` is **deprecated** and ignored (it double-encoded magnitude via
  linetype); it warns if supplied. `mono = TRUE` is retained as a convenience
  wrapper for `sign_by = "linetype"` with black edges.

## Pruning is now a standalone `prune()` verb (⚠ breaking, no deprecation)

Forbes-extension pruning is extracted out of `ackwards()` into a standalone,
pipeable S3 generic, `prune()`:

```r
x <- ackwards(bfi25, k_max = 5)
xp <- prune(x, "redundant")
```

- The five pruning-related arguments — `prune` (renamed `rules`),
  `redundancy_r`, `redundancy_phi`, `min_items`, `orphan_r` — are **removed
  from `ackwards()` entirely** (clean move; the package is pre-CRAN with no
  users, so there is no deprecation shim). Calling `ackwards()` with any of
  these now errors with a pointer to `prune()`, rather than silently
  discarding them.
- Because extraction (expensive) and pruning (cheap, interpretive) are now
  separate steps, you can re-prune with new thresholds without
  re-extracting: `prune(x, "redundant", redundancy_r = 0.95)` reuses the
  already-fit object.
- New `manual =` argument flags user-named factors directly — standalone
  (`prune(x, manual = c("m4f3"))`) or unioned onto an auto rule
  (`prune(x, "redundant", manual = c("m4f3"))`).
- The canonical rule name is now **`"artifact"`** (US spelling); `"artefact"`
  is still accepted as an alias and normalizes to `"artifact"`.
- `ackwards()`'s `pairs` argument no longer auto-upgrades to `"all"` when
  pruning is requested — `prune()` recomputes its own all-pairs edges on
  demand, so it works correctly regardless of the fit-time `pairs` setting,
  and `x$edges` is never modified by pruning.

## New dataset: `sim16`

A simulated, fully continuous 1 000×16 teaching dataset (`data(sim16)`) with
a known 1 → 2 → 4 bass-ackwards hierarchy, alongside `bfi25`. Two purposes:

- Showcases the default `cor = "pearson"` extraction path cleanly — unlike
  `bfi25`'s Likert items, `sim16` never triggers the ordinal-detection
  warning.
- Deliberately exercises overextraction: the population has exactly 4
  factors, so requesting `k_max = 5` produces a guaranteed redundant chain
  (`prune(x, "redundant")`) and an orphan/few-items artifact signal
  (`prune(x, "artifact")`) for the Forbes/redundancy examples to teach
  against.

`sim16` is deliberately the *idealized* case — its clean signal makes all
`suggest_k()` criteria agree on `k = 4`, in contrast to `bfi25` (real data),
where they span `k = 4`–`6`. The two are complementary teaching foils: `sim16`
for recovering a known structure, `bfi25` for reasoning under criterion
disagreement.

See `?sim16` for the full generative model and ground-truth hierarchy.

## `tidy()` API-shape cleanup: renamed/removed columns

Three `tidy.ackwards()` output changes, all breaking (no deprecation; the
package is pre-CRAN with no users):

- `tidy(x, what = "fit")` — the long-format key column `index` is renamed
  **`statistic`** (it held fit-index names for EFA/ESEM and eigenvalue
  positions for PCA; `index` read like a row position).
- `tidy(x, what = "fit")` — the `cutoffs = TRUE` argument, and the `meets`/
  `{statistic}_meets` columns it produced, are **removed**. A pass/fail
  boolean quietly endorsed Hu & Bentler (1999) thresholds the package
  elsewhere treats as conventional and contested. Conventional-threshold
  reference lines/annotations remain available via `autoplot(what = "fit")`
  and `summary()`.
- `tidy(x, what = "variance")` — `variance_pct`/`cumulative_pct` (0-100) are
  renamed **`proportion`**/**`cumulative`** and now report a 0-1 scale,
  matching the engine's internal variance representation and broom/psych
  convention. `print()` and `summary()` still display percentages.

`k_max` is retained, unchanged, in both `ackwards()` (extraction depth) and
`suggest_k()` (max factors/components evaluated); the two functions' roxygen
now cross-reference the distinct meanings sharing that name.

## Effective ESEM estimator recorded in `$meta`; `summary()` footnote

`x$meta$estimator` now stores the effective ESEM estimator after
auto-selection (`"ML"`, `"MLR"`, `"WLSMV"`, or `"ULSMV"`; `NA` for PCA/EFA).
`summary()` adds a one-line footnote naming lavaan's mean-and-variance-adjusted
("scaled") fit-index reporting whenever the effective estimator is
`"WLSMV"`/`"ULSMV"`/`"MLR"` (not shown for `"ML"`, which has no scaled
variant) — making the scaled-reporting decision visible in the printed
summary, not just in `tidy()`/roxygen documentation.

## Vignette corrections: `ackwards-intro` and `ackwards-suggest-k`

`vignette("ackwards-intro")` hardcoded a cumulative-variance jump ("22.9% →
34.7%") that had drifted from the code's actual output (23.2% → 35.5%); it is
now computed inline from the fitted object so it cannot drift again. The
lineage walkthrough had `m4f1`'s primary children backwards (claimed
`m5f2`/`m5f4`; the edge table shows `m5f1`/`m5f4`), and the diagram narrative
misattributed which traits differentiate at which level (Conscientiousness/
Openness split at k = 4, not Agreeableness/Extraversion, which split at
k = 5) — both corrected to match a live run. `vignette("ackwards-suggest-k")`
now clarifies why the printed "Recommendations" block shows six lines for
five criteria (`"vss"` reports both VSS-1 and VSS-2).

## Guard against `cor = "polychoric"` with an incompatible ESEM estimator

`ackwards(engine = "esem", cor = "polychoric", estimator = "ML")` (or
`"MLR"`) now errors immediately with a clear explanation instead of failing
many calls deep inside the per-level ESEM fit and surfacing a misleading
"failed to build at least 2 converged levels... check your data for
multicollinearity" abort. Polychoric correlations mark every item `ordered`
for lavaan, and lavaan itself does not support ML/MLR estimation on ordered
indicators. `"WLSMV"`/`"ULSMV"` with a continuous `cor` remains allowed (a
valid, if atypical, continuous WLS/ADF estimator).

## ESEM fit indices: scaled variants, honest p-value, and BIC

For ESEM under a **scaled-test estimator** (`"WLSMV"`/`"ULSMV"` for ordinal
data; `"MLR"` for continuous), the entire fit row — `chi`/`dof`/`p_value`
**and** `CFI`/`TLI`/`RMSEA` — now reports lavaan's mean-and-variance-adjusted
("scaled") variant, so every quantity in the row shares one scaling framework.
Three consequences:

* **WLSMV/ULSMV p-value.** The naive chi-square test has no valid reference
  distribution for these limited-information estimators (lavaan's own
  `summary()` labels its p-value "Unknown"/`NA`); the scaled test does, and is
  now what `p_value` reports. A genuinely saturated level (`dof = 0`) still
  reports `NA` — there is no test to perform on a model that fits perfectly by
  construction.
* **WLSMV/ULSMV CFI/TLI/RMSEA.** These previously reported the *naive* indices,
  which are known to be badly optimistic for ordinal data (Xia & Yang, 2019) —
  e.g. on the BFI the naive CFI (~0.89) sat right next to a scaled p-value,
  overstating fit. They now report the scaled variants (~0.71 in that example),
  consistent with the reported test.
* **MLR.** Previously the fit row showed the *naive* ML statistics for an MLR
  fit, defeating the purpose of robust ML; it now reports the scaled
  (Yuan-Bentler) test and indices.

`"ML"` has no scaled variant and continues to report its naive values (the
correct ones for ML). `SRMR` has no scaled variant and is reported as-is.

`BIC` is now a first-class ESEM fit index (previously silently absent):
populated under `"ML"`/`"MLR"` (a proper log-likelihood exists), and `NA`
under `"WLSMV"`/`"ULSMV"` (no proper log-likelihood for these estimators —
genuinely inapplicable, not a bug). `tidy(x, what = "fit")` and `glance(x)`
both surface it consistently across engines.

## `tidy(what = "fit", cutoffs = TRUE)`: no more always-NA `_meets` columns

`format = "wide"` no longer emits a `{index}_meets` column for indices with no
defined threshold (`chi`, `dof`, `p_value`, `BIC`, PCA eigenvalues) — the pivot
previously generated one for every index, and those columns were always `NA`.
`format = "long"` is unaffected (`meets` is still `NA` for those rows, as
documented).

## Citation hygiene

`citation("ackwards")` now returns a single software entry (Girard) instead
of also listing Goldberg (2006) as a package author — the method paper is
cited in `DESCRIPTION` and roxygen `@references`, not as a software
co-citation. `ackwards()`'s `@references` now also lists Forbes (2023)
alongside Goldberg (2006) and Waller (2007), since `ackwards()` implements
all three (original method, exact edge algebra, and the extended method).
The README's Citation section is updated to guide users to the right
paper(s) for what they used.

## Documentation: no internal milestone numbers

User-facing docs (`NEWS.md`, `README.md`, vignettes) no longer reference internal
development-milestone tags — they are meaningless outside this repo's own process.
A regression test now guards against reintroducing them.

## CD criterion: bug fix and honest framing

`suggest_k()` now correctly handles `EFAtools::CD`'s output. `EFAtools::CD`
fills its `RMSE_eigenvalues` matrix only up to the level it actually tested;
trailing columns were literal zeros, causing the plotted RMSE curve to dive
spuriously to zero at higher k and `which.min()` to land on a fake value.
The `cd_rmse` vector now has those unfilled positions set to `NA`, so the curve
terminates at the last genuinely computed level.

The CD plot panel is relabeled from `"CD (RMSE, minimize)"` to
`"CD (RMSE; sequential test)"` to reflect how `EFAtools::CD` actually works: a
sequential one-sided Wilcoxon test (default α = 0.30) that stops at the first
non-significant RMSE improvement. The starred k is the last retained factor and
need not be the visible minimum of the curve. The `autoplot.suggest_k()` and
`suggest_k()` documentation has been updated accordingly.

## Per-level fit indices and loading SEs as first-class output

**`glance()` now carries fit indices** from the deepest converged level. A
consistent five-column set (`CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`) is present
for all engines; indices unavailable for a given engine are `NA` (e.g. `CFI`
and `SRMR` for EFA; all five for PCA).

**`tidy(what = "fit")` gains two new arguments:**
- `format = "wide"`: returns one row per non-anchor level (k >= 2) with one
  column per index — the natural shape for reporting tables. The saturated
  1-factor anchor is dropped, matching `summary()` and `autoplot(what = "fit")`.
- `cutoffs = TRUE`: appends a `meets` column flagging each index against Hu &
  Bentler (1999) conventional thresholds (CFI/TLI >= .95, RMSEA <= .06,
  SRMR <= .08). Indices without a defined threshold (chi, BIC, eigenvalues)
  get `NA`. Thresholds are report-only and never gate any behaviour.

**`tidy(what = "loadings")` now includes `se`, `ci_lower`, and `ci_upper`
columns** for all engines. For ESEM these are populated from the rotation-aware
loading SEs; for PCA and EFA they are `NA`. A `conf_level` argument (default
`0.95`) controls the interval width. This replaces the need to compute CIs by
hand, and supersedes the previous `tidy(what = "loadings_se")` accessor, which
has been **removed** (the package is unreleased, so no deprecation cycle is
needed).

**`autoplot(x, what = "fit")`** produces a two-panel ggplot2 chart of per-level
fit indices (CFI/TLI in the top panel; RMSEA/SRMR in the bottom panel) with
Hu & Bentler (1999) reference lines. PCA objects return an informative empty
plot. `what = "hierarchy"` (default) is unchanged.

**`summary()` fit lines** now append a check (✔) or cross (✘) to each
index compared against its Hu & Bentler (1999) threshold.

**`print()` and `summary()` footer** extended: per-level fit indices describe
how well a k-factor model fits the items at that level — they do not validate
the edges or the hierarchy itself.

**New vignette section** in `vignette("ackwards-engines")`: "Per-level fit:
what it tells you (and what it doesn't)" — covers the level-vs-hierarchy
distinction, reporting with the new wide table and fit plot, the ESEM
cost/benefit tradeoff, and when to care about fit in exploratory vs publication
workflows.

## Faster ESEM on large item sets

The ESEM engine no longer recomputes lavaan's sample statistics (thresholds, the
polychoric correlation matrix, and the asymptotic weight matrix) at every level.
These depend only on the data, so they are now computed once at the first level and
reused for every deeper level via lavaan's `slotSampleStats=`. Solutions are
identical; the redundant recomputation — which dominated run time with hundreds of
items — is removed.

The independent per-level model fits can now also run in parallel. When the
optional `future.apply` package is installed, `ackwards()` dispatches the ESEM
levels through the `future` framework. The default plan is sequential (no behaviour
change); call `future::plan(future::multisession, workers = N)` before `ackwards()`
to parallelize. Results are reproducible across plans when `seed` is supplied.
`future.apply` added to `Suggests`. See the "Performance with many items" section of
`vignette("ackwards-engines")`.

## Tucker's φ default for non-PCA redundancy (⚠ resolved-default change)

`redundancy_phi = NULL` (the default) now **auto-resolves** based on the extraction engine:

* `engine = "pca"` — no φ filter (the closed-form W′RW algebra is exact; `|r|` alone is
  sufficient, as in the original Waller 2007 method).
* `engine = "efa"` or `"esem"` — automatically applies `redundancy_phi = 0.95` (Lorenzo-Seva
  & ten Berge, 2006). Factor-score indeterminacy off-PCA makes `|r|`-only redundancy too
  liberal; the conjunctive φ criterion is the conservative default.

A cli message is emitted whenever auto-resolution applies (Invariant 6: loud defaults). The
resolved value is stored in `x$prune$redundancy_phi`.

**Existing non-PCA calls using `prune = "redundant"` will become more conservative** — some
previously flagged redundant chains may no longer meet the φ criterion. To restore the old
behaviour, pass `redundancy_phi = NA` (explicit opt-out; no φ filter regardless of engine).
Explicit numeric values (e.g., `redundancy_phi = 0.8`) override on any engine.

## Structural artefact signals (Forbes extension)

`ackwards(..., prune = "artefact")` now computes structural signals alongside Tucker's φ
table. Two new arguments on `ackwards()`:

* `min_items = 3L` — minimum primary-item count per factor; factors with fewer items are
  flagged `few_items = TRUE`.
* `orphan_r = 0.5` — minimum adjacent-level `|r|` required for a factor to be considered
  "connected"; factors below this threshold across all adjacent levels are flagged
  `orphan = TRUE`.

A third signal, `split_merge`, is `TRUE` when a factor at level k draws its primary items
from multiple distinct primary factors at level k−1 (items that were in separate groups at
the shallower level merged into a single factor at the deeper level).

All signals are **flag-only** (artefact identification requires researcher judgment;
Forbes, 2023). Results are stored in `x$prune$structural` (one row per factor × level);
`print()` and `summary()` report the count of flagged factors. Both `min_items` (positive
integer) and `orphan_r` (number in `[0, 1]`) are validated with a clear error. See
`vignette("ackwards-forbes")` for a worked example.

## `suggest_k()` criterion selection

`suggest_k()` gains a `criteria` argument that controls which retention criteria are
computed. Any subset of `c("pa_pc", "pa_fa", "map", "vss", "cd")` may be requested;
the default runs all five (identical to prior behaviour). Criteria not requested are
skipped entirely (no computation), delivering a real speed win:

* `"pa_pc"` / `"pa_fa"` share one `psych::fa.parallel()` call — both skipped if neither
  is requested.
* `"map"` / `"vss"` share one `psych::vss()` call — both skipped if neither is
  requested.
* `"cd"` runs `EFAtools::CD()` only when explicitly requested **and** `EFAtools` is
  installed; if requested but unavailable an informational message is emitted.

`k_*` fields for non-requested criteria are `NA_integer_`; the `criteria` data frame
uses `NA` for non-run columns (stable schema). The returned object gains a
`criteria_requested` field. `print()` and `autoplot()` render only the requested
criteria; the consensus range is computed from requested criteria only.

# ackwards 0.1.0

First public release. Licensed under MIT.

## Extraction engines

Three engines are available via the `engine` argument to `ackwards()`:

* **PCA** (`engine = "pca"`) — the original Goldberg (2006) bass-ackwards method.
  Fastest; never fails to converge; exact W′RW algebra for between-level edges.
* **EFA** (`engine = "efa"`) — exploratory factor analysis via `psych::fa()` with
  tenBerge scoring weights and the same closed-form algebra path.
* **ESEM** (`engine = "esem"`) — EFA estimated in a SEM framework via `lavaan::efa()`
  with WLSMV ordinal estimation, per-level fit indices (CFI, TLI, RMSEA, SRMR), and
  rotation-aware standard errors. Required for the clinical/HiTOP workflow.

All engines use **varimax** (orthogonal) rotation — the only rotation compatible with the
W′RW algebra. Oblique rotation is unsupported by design.

## Choosing k — `suggest_k()`

`suggest_k()` runs five complementary retention criteria and returns a consensus range:

* **PA-PC** and **PA-FA** — parallel analysis on the PC and FA eigenvalue bases
  (`psych::fa.parallel(fa = "both")`).
* **MAP** (Velicer 1976) and **VSS-1/VSS-2** (Revelle & Rocklin 1979) via `psych::vss()`.
* **Comparison Data** (Ruscio & Roche 2012) via `EFAtools::CD()` — skipped gracefully when
  `EFAtools` is not installed.

`autoplot.suggest_k()` produces a ggplot2 diagnostic: a four-panel 2×2 grid (scree/PA,
MAP, VSS, CD RMSE) when `EFAtools` is installed, or a three-panel single-column layout
(scree/PA, MAP, VSS) otherwise. The returned object also carries `cd_rmse` (column means
of CD's RMSE eigenvalue matrix) for use in the CD panel.

## Between-level edges — exact W′RW algebra

Between-level factor-score correlations are computed via closed-form `W′RW` algebra for all
linear engines (PCA, EFA, ESEM with tenBerge weights) — no scores need to be materialized.
Standardization uses real score SDs from `sqrt(diag(W′RW))` rather than assuming unit
variance. A materialized-scores fallback is available and a cross-check test verifies the
two routes agree within tolerance.

## Ordinal / Likert data

Pass `cor = "polychoric"` to use polychoric correlations as the basis for PCA or EFA.
For ESEM, lavaan's own polychoric matrix is used for edge computation (consistent with the
fitted model). The WLSMV estimator is the default for `cor = "polychoric"` with ESEM.
`ackwards()` detects likely-ordinal columns (≤ 7 distinct integer values) and emits a
suppressible warning when `cor = "pearson"` is set.

## Missing data — `missing` argument

Three strategies via `missing` on `ackwards()`:

* `"pairwise"` (default) — pairwise-complete correlations; warns when NAs are present.
* `"listwise"` — complete cases only, consistent N across fit and edges.
* `"fiml"` — Full Information ML via lavaan (ESEM ML/MLR only); edges derived from the
  FIML saturated model.

## Forbes (2023) extension

* `pairs = "all"` computes between-level correlations for every level pair, revealing
  skip-level associations.
* `prune = "redundant"` flags chains where adjacent `|r| ≥ 0.9` (Forbes retention rule).
* `prune = "artefact"` computes Tucker's φ for researcher inspection; no auto-flagging.
* `autoplot()` renders pruned views (`drop_pruned = TRUE`) and skip-level arcs.

## Visualization — `autoplot.ackwards()`

A ggplot2 layered diagram with extensive customization: edge color by sign, width/alpha by
`|r|`, solid/dashed by strength; monochrome mode (`mono`); APA-formatted correlation labels
(`show_r`, `r_digits`, `r_label_size`); custom node labels (`node_labels`); primary-only
mode (`primary_only`); level axis labels (`show_level_labels`); arrowhead toggle
(`show_arrows`); fixed edge width (`edge_linewidth`); legend toggle (`legend`); Forbes-style
pruned view (`drop_pruned`, `compress_levels`).

## Factor interpretation

* **`top_items(x, level, cut, n, sort)`** — lists the salient items for each factor,
  filtered to `|loading| >= cut` (default `0.3`) and sorted by descending absolute loading.
  Arguments `level`, `n`, and `sort` allow subsetting levels, capping items per factor, and
  controlling order. Returns a `top_items` S3 object with a grouped cli print method.

* **`label_template(x, style)`** — generates the named-character-vector scaffold for
  `autoplot(x, node_labels = ...)`. Styles: `"id"` (default), `"forbes"` (level-letter +
  within-level index: `"A1"`, `"B1"`, `"B2"`, …), `"blank"`. Factor IDs are in canonical
  layout order; an editable `c(...)` literal is also printed for copy-paste.

* See `vignette("ackwards-interpret")` for the end-to-end naming workflow: reading a factor
  with `top_items()`, hierarchy-aware naming (parent vs. child, blends), sign-alignment
  caveat, and the `top_items()` → `label_template()` → `autoplot(node_labels = ...)` round-trip.

## Dependencies

* **`psych` moved to `Imports`** — the default PCA and EFA engines require psych; placing it in
  Suggests meant a mandatory install prompt for core functionality. The SEM (`lavaan`) and
  plotting (`ggplot2`) stacks remain in Suggests.
* **`GPArotation` removed** — varimax rotation routes through base `stats::varimax`; GPArotation
  was never loaded on any supported path.
* **`bfi25` example dataset** bundled — a 1 000-row, 25-item subset of the SAPA/IPIP Big Five
  data (sampled from `psych::bfi`). Used throughout examples and vignettes so they run without
  reaching into psych's namespace. See `?bfi25` for provenance and `@source`.

## Correlation-matrix input

`ackwards()` and `suggest_k()` now accept a pre-computed **correlation matrix** in place of
raw item data, detected automatically from the matrix shape (square, symmetric, unit diagonal).

* **Engine gating** — only `"pca"` and `"efa"` are supported; `"esem"` errors clearly (lavaan
  requires raw data for WLSMV estimation and per-level fit indices).
* **`n_obs` argument** — new optional argument on both functions. Required for `engine = "efa"`
  with a matrix (psych needs N for chi-square/RMSEA/TLI); optional for `engine = "pca"` (stored
  as `NA` when omitted). Ignored (with a warning) when raw data are supplied.
* **Edge correctness** — edges from `ackwards(R, ...)` match `ackwards(data, ...)` exactly for
  the same correlation matrix: both routes use the closed-form `W'RW` algebra.
* **`cor` and `missing` arguments** — ignored (with a warning if set) for matrix input. `$cor`
  is stored as `NA_character_` and printed as `"(user-supplied matrix)"`.
* **Score paths blocked** — `keep_scores = TRUE`, `augment()`, and `tidy(what = "scores")`
  error clearly: individual-level scoring requires row-level item responses.
* **CD skipped in `suggest_k()`** — Comparison Data resamples raw item distributions and is
  skipped with an info note when a matrix is supplied.
* **Validation** — square, numeric, finite, symmetric, unit diagonal, `|r| <= 1`, no NA;
  synthesises `V1..Vp` dimnames when absent; warns (does not auto-smooth) when not
  positive-definite; errors with a clear message when a covariance matrix is detected.

## Tidy interface and scoring

* `print.ackwards()` — compact cli summary card.
* `summary.ackwards()` — per-level variance, fit indices, lineage list, pruning notes.
* `tidy(x, what = "edges" / "loadings" / "loadings_se" / "variance" / "nodes" / "scores")` —
  long-format tidy tibbles. For `what = "edges"`, `primary_only = TRUE` returns just each
  factor's primary-parent edge (the lineage tree) and `sort = "strength"` orders by `|r|`.
  `what = "loadings_se"` returns the rotation-aware loading standard errors (ESEM only).
* `glance(x)` — one-row model-level summary.
* `augment(x, data = ...)` — appends factor scores to a data frame; recomputes from stored
  weights when scores were not materialized at fit time.

## Documentation

Vignette comparison tables reworked for legibility:

* `ackwards-engines` and `ackwards-ordinal` — stacked long-format `kable` tables replaced
  by **wide gt tables**: one row per item/edge, one column per engine/basis, plus an explicit
  Δ column that is the teaching point (EFA/polychoric attenuation). Factor and sign
  correspondence is asserted before differencing; edge tables surface primary-parent
  disagreements as NA rather than hiding them.
* `ackwards-forbes` — `prune-nodes` raw `tidy()` print replaced by a styled gt table;
  narrated counts (`n_adj`, `n_all`, `n_redundant`, top skip-edge) are now computed via
  inline R expressions so they can never go stale.
* `gt` added to Suggests (vignette table formatting only; never loaded by core).

## Test coverage

Test suite reaches **100% line coverage**. Genuinely unreachable defensive branches
(engine fallbacks, lavaan error/warning handlers) are excluded via `# nocov` rather than
brittle environment-dependent tests.

