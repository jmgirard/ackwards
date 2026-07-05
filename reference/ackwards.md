# Bass-ackwards hierarchical structural analysis

Extracts factor/component solutions at levels 1 through `k`, then
characterises the hierarchy by computing between-level factor-score
correlations. The "hierarchy" is descriptive: edges are score
correlations, not a fitted higher-order SEM.

## Usage

``` r
ackwards(
  data,
  k_max,
  engine = "pca",
  cor = "pearson",
  fm = "minres",
  estimator = NULL,
  missing = "pairwise",
  n_obs = NULL,
  align_signs = TRUE,
  keep_scores = FALSE,
  keep_fits = FALSE,
  seed = NULL,
  pairs = "adjacent",
  cut_show = 0.3,
  correct = 0.5,
  ...
)
```

## Arguments

- data:

  A data frame or numeric matrix of observed variables (items in
  columns, observations in rows). Alternatively, a pre-computed
  **correlation matrix** may be supplied (a square, symmetric, numeric
  matrix with unit diagonal). When a correlation matrix is supplied,
  `engine` must be `"pca"` or `"efa"` (ESEM requires raw data), and the
  `missing` and `cor` arguments are ignored. See the *Correlation-matrix
  input* section below.

- k_max:

  Maximum number of factors/components to extract. Required; use
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  if uncertain. Sets the *depth* of the hierarchy: levels 1 through
  `k_max` are all extracted and retained. (Note:
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)'s
  `k_max` means something related but distinct – the max number of
  factors/components it *evaluates* when recommending a depth, not a
  depth itself; see that function's docs.)

- engine:

  Extraction engine: `"pca"` (default), `"efa"`, or `"esem"`. `"esem"`
  uses [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) with
  rotation-aware SEs and per-level fit indices; recommended for the
  clinical/HiTOP workflow (Kim & Eaton, 2015; Forbush et al., 2024).
  Requires lavaan \>= 0.6-13.

- cor:

  Correlation basis: `"pearson"` (default), `"spearman"`, or
  `"polychoric"`. For PCA/EFA, `"polychoric"` computes a polychoric
  correlation matrix via
  [`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  (requires psych). For ESEM, it triggers WLSMV estimation via lavaan. A
  cli warning is emitted when ordinal-looking columns are detected and
  `cor` is not `"polychoric"`.

- fm:

  Factor extraction method passed to
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html); only used when
  `engine = "efa"`. One of `"minres"` (default, robust OLS), `"ml"`
  (maximum likelihood, gives chi-square fit but converges less reliably
  at deep levels), or `"pa"` (principal axis). Ignored for
  `engine = "pca"`.

- estimator:

  Estimation method for the ESEM engine. `NULL` (default) auto-selects:
  `"WLSMV"` when `cor = "polychoric"`, `"ML"` otherwise. Pass explicitly
  to override: `"ULSMV"` (unweighted WLS), `"MLR"` (robust ML).
  `cor = "polychoric"` with `estimator = "ML"`/`"MLR"` errors (lavaan
  itself does not support ML/MLR on ordered indicators); `"WLSMV"`/
  `"ULSMV"` with a continuous `cor` is allowed (a valid, if atypical,
  continuous WLS/ADF estimator). Ignored for PCA and EFA engines. The
  effective value (after auto-selection) is recorded in
  `x$meta$estimator` (`NA` for PCA/EFA).

- missing:

  How to handle missing item responses. One of:

  - `"pairwise"` (default) – use all available observations pairwise.
    For PCA/EFA this feeds `stats::cor(use = "pairwise.complete.obs")`.
    For ESEM with WLSMV/ULSMV (ordinal), lavaan uses `available.cases`,
    which computes polychoric thresholds and correlations from all rows
    that contribute to each pair – MCAR-valid and uses the full N. For
    ESEM with ML/MLR (continuous), lavaan uses listwise deletion
    internally while edges are computed from a pairwise correlation
    matrix; this minor inconsistency is documented in `$meta`. A warning
    is emitted when incomplete rows are detected.

  - `"listwise"` – only complete rows are used. Reduces data to
    [`stats::complete.cases()`](https://rdrr.io/r/stats/complete.cases.html)
    before fitting, so the correlation matrix, the engine fit, and the
    edges are all consistent. `n_obs` in the result reflects the reduced
    N.

  - `"fiml"` – Full Information Maximum Likelihood. For
    `engine = "esem"` (with `estimator = "ML"`/`"MLR"`), passes
    `missing = "fiml"` to
    [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) and
    derives edge correlations from lavaan's FIML-estimated saturated
    model. For `engine = "pca"`/`"efa"` (M38), the correlation matrix is
    estimated via
    [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html) and
    fed to the usual `W'RW` algebra; this requires `cor = "pearson"`
    (corFiml estimates a multivariate-normal matrix) and the route is
    announced via a cli message. Errors for WLSMV/ULSMV and for a
    non-Pearson PCA/EFA basis. Note: FIML improves estimation under
    missingness but does not impute item responses; score
    materialisation (`keep_scores = TRUE`) still produces `NA` rows for
    incomplete observations. See `n_obs` for the fit-index sample size
    on the PCA/EFA path.

- n_obs:

  Number of observations, or (on the raw-data FIML path) a string
  selecting which N feeds the fit indices.

  - **Correlation-matrix input:** a positive integer. Required for
    `engine = "efa"`
    ([`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) needs N for
    chi-square / RMSEA / TLI); optional for `"pca"` (stored as
    `NA_integer_` if omitted, disabling N-dependent fit statistics).

  - **Raw data:** N is normally taken from `nrow(data)` and a numeric
    `n_obs` is ignored (with a warning). The exception is
    `missing = "fiml"` with `engine = "pca"`/`"efa"` (M38): because
    [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
    estimates the correlation matrix from incomplete rows, `n_obs` may
    be `"total"` (default – every row contributing to the FIML
    likelihood, matching the FIML convention; Enders, 2010) or
    `"complete"` (complete-case N, a conservative lower bound). Point
    estimates (loadings, edges) do not depend on this choice; only the
    EFA fit indices do, and those are *approximate* under this two-step
    (FIML matrix into normal-theory EFA) route regardless of N (Zhang &
    Savalei, 2020). A string `n_obs` is accepted only on this path.

- align_signs:

  Logical; sign-align factors to primary-parent lineage? Default `TRUE`.

- keep_scores:

  Logical; store factor scores in the result? Default `FALSE`
  (recomputable via
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)).
  When `TRUE`, per-observation scores are stored in `x$scores` as a
  named list of `n x k_j` matrices, one per level, standardized by real
  score SDs (see
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)).

- keep_fits:

  Logical; store raw engine fit objects? Default `FALSE`. When `TRUE`,
  the per-level fit objects (psych or lavaan) are stored in `x$fits` as
  a named list indexed by level.

- seed:

  Integer seed for stochastic engines (not used by PCA but captured for
  reproducibility metadata). Default `NULL`.

- pairs:

  Which level pairs to compute edges for: `"adjacent"` (default, classic
  Goldberg – only consecutive levels) or `"all"` (Forbes extension –
  every pair of levels, including skip-level correlations).
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  recomputes its own all-pairs edges on demand regardless of this
  setting, so pruning does not require `pairs = "all"` here.

- cut_show:

  Edges with `|r| >= cut_show` are flagged `above_cut` in
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) output.
  Default `0.3`.

- correct:

  Continuity correction passed to
  [`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  on the PCA/EFA polychoric path (`engine = "pca"`/`"efa"` with
  `cor = "polychoric"`). Default `0.5` (psych's own default), which adds
  that value to zero cells before estimating thresholds. **Set
  `correct = 0`** if
  [`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  fails on your data (its error suggests exactly this) – this typically
  happens when an item has a near-empty response category or when items
  with unequal category counts produce a sparse cross-cell. Ignored on
  other paths: ESEM computes its own polychoric correlations inside
  lavaan, and the Pearson/Spearman bases do not use it.

- ...:

  Reserved for future arguments.

## Value

An object of class `"ackwards"`. See
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
and
[`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
for output methods.

## Defaults and why

- **`engine = "pca"`** – the original Goldberg (2006) method; fastest;
  never fails to converge; the Waller (2007) algebra is exact for
  components.

- **`rotation = "varimax"`** – the `T'=T^-1` property of orthogonal
  rotation enables the closed-form `W'RW` edge algebra and keeps
  within-level factors uncorrelated so cross-level edges reflect only
  the hierarchical signal. Matches Goldberg (2006), Kim & Eaton (2015),
  and Forbush et al. (2024). Varimax is the only supported rotation;
  oblique rotation would confound the between-level signal that is the
  method's core output.

- **`cor = "pearson"`** – no silent basis switching. If your items look
  ordinal (\<= 7 distinct integer values), a cli warning will suggest
  `cor = "polychoric"`, which is available for all three engines.

- **`align_signs = TRUE`** – unaligned signs make the output unreadable.
  Anchor: m1f1 is oriented toward the positive manifold; each subsequent
  factor is flipped so its edge to its primary parent is positive.

- **`keep_scores = FALSE` / `keep_fits = FALSE`** – memory and privacy.
  Scores are O(n x Sigmak) and often sensitive; raw engine fits can be
  large. Both are recomputable from the stored `r` matrix.

## Performance (ESEM, large item sets)

The ESEM engine fits a separate `lavaan` model at every level
1..`k_max`. For ordinal data (`cor = "polychoric"`, WLSMV) the costly
sample statistics lavaan derives from the raw data – thresholds, the
polychoric correlation matrix, and the asymptotic weight matrix – depend
only on the data, not on the number of factors, so they are **computed
once** at the first level and **reused** for every deeper level
(identical solutions, much less work). This matters most when you have
many items (hundreds), where recomputing those statistics at each level
dominated the run time.

The per-level model fits are mutually independent and are dispatched
through the future framework when future.apply is installed. By default
the plan is sequential (no behaviour change). To run the levels in
parallel, set a plan once before calling `ackwards()`:


      future::plan(future::multisession, workers = 4)  # or multicore on Unix
      x <- ackwards(items, k_max = 8, engine = "esem", cor = "polychoric")

Parallelism pays off when the per-level fits are heavy (large `p`,
several levels); for small problems the worker startup cost can outweigh
it. Results are reproducible across plans when `seed` is supplied. PCA
and EFA already compute their correlation matrix once and are
unaffected.

## Correlation-matrix input

When `data` is a pre-computed correlation matrix (square, symmetric,
unit diagonal), `ackwards()` runs entirely from that matrix using the
`W'RW` algebra – no raw item responses are needed. This is useful when
you have a published correlation table or a polychoric matrix computed
externally.

Constraints and behaviour when a correlation matrix is supplied:

- **Engine:** only `"pca"` and `"efa"` are supported. `"esem"` requires
  raw data (for lavaan's own polychoric computation, WLSMV estimation,
  and per-level fit indices) and will error clearly.

- **`n_obs`:** required for `"efa"` (psych needs N for chi-square /
  RMSEA / TLI); optional for `"pca"` (stored as `NA` if omitted).

- **`cor` argument:** ignored – the basis is already determined by the
  matrix you supply. A warning is emitted if you set `cor` explicitly.

- **`missing` argument:** ignored – missingness was handled when
  computing the matrix. A warning is emitted if you set `missing`
  explicitly.

- **Factor scores:** `keep_scores = TRUE` will error;
  [`augment()`](https://generics.r-lib.org/reference/augment.html) and
  `tidy(what = "scores")` will also error because individual-level
  scores require row-level item responses.

- **`$cor` field:** stored as `NA_character_`; printed as
  `"(user-supplied matrix)"`.

## When to trust the result

`ackwards()` raises diagnostics as it fits. They fall into three tiers
by what they mean for whether you should trust and report the solution:

**Fatal – fix before trusting.** The result is undefined or rests on a
broken correlation matrix:

- A **constant item** (no variance) errors – drop it (see
  [`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md)).

- **[`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  fails** – usually a near-empty response category; set `correct = 0` or
  collapse rare categories.

- A **level fails to converge** – the hierarchy is truncated to the
  deepest level that did; do not interpret beyond it.

- A **near-singular correlation matrix** (smallest eigenvalue `< 1e-4`,
  recorded in `meta$near_singular` / `meta$min_eigenvalue` and
  re-surfaced by
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html))
  means per-level fit indices and factor scores are unreliable and the
  loadings/edges rest on a rank-deficient matrix. (For EFA the
  residual-based fallback inflates `TLI`/`RMSEA`; for ESEM `CFI` comes
  back `NA`.) Trim redundant items, use `missing = "listwise"`, or – on
  the polychoric basis – collapse sparse categories or set
  `correct = 0`.

**Caution – interpret carefully.** The solution exists but may be
unstable:

- A **Heywood case** (communality `> 1` / negative uniqueness) at a
  level – check that level's loadings make substantive sense; consider
  fewer factors.

- A **near-constant item** (one response category dominates) can drive a
  meaningless factor – inspect it with
  [`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md).

- **Ordinal data on a Pearson basis** attenuates correlations – fine for
  a quick look or
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  screening, but report the final model on `cor = "polychoric"`.

**Informational – usually fine.** Proceed, just be aware: the
pairwise-missing note, a merely *sparse* (rare-but-present) response
category, and the ordinal-detection warning when you did intend Pearson.

## References

Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research
in Personality*, 40(4), 347–358.
[doi:10.1016/j.jrp.2006.01.001](https://doi.org/10.1016/j.jrp.2006.01.001)

Waller, N. G. (2007). A general method for computing hierarchical
component structures by Goldberg's bass-ackwards method. *Journal of
Research in Personality*, 41(4), 745–752.
[doi:10.1016/j.jrp.2006.08.005](https://doi.org/10.1016/j.jrp.2006.08.005)

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

Enders, C. K. (2010). *Applied Missing Data Analysis*. Guilford Press.

Zhang, X., & Savalei, V. (2020). Examining the effect of missing data on
RMSEA and CFI under normal theory full-information maximum likelihood.
*Structural Equation Modeling*, 27(2), 219–239.
[doi:10.1080/10705511.2019.1642111](https://doi.org/10.1080/10705511.2019.1642111)

## See also

[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
(Forbes-extension redundancy/artifact flagging, piped off the result of
this function)

## Examples

``` r
# sim16 is continuous with a known 1 -> 2 -> 4 hierarchy, so the default
# pearson basis is appropriate. For ordinal items (e.g. bfi25), fit on the
# polychoric basis instead -- see `cor = "polychoric"` and the ordinal
# vignette.
x <- ackwards(sim16, k_max = 4)
print(x)
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 4
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 28.2% variance
#> ✔ k = 2: 2 factors, 46.5% variance
#> ✔ k = 3: 3 factors, 57.5% variance
#> ✔ k = 4: 4 factors, 67.7% variance
#> 
#> ── Edges ──
#> 
#> 9 of 20 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
tidy(x)
#>    from   to level_from level_to            r is_primary above_cut
#> 1  m1f1 m2f1          1        2  0.707296113       TRUE      TRUE
#> 2  m1f1 m2f2          1        2  0.706917399       TRUE      TRUE
#> 3  m2f1 m3f1          2        3 -0.011342395      FALSE     FALSE
#> 4  m2f1 m3f2          2        3  0.697092981       TRUE      TRUE
#> 5  m2f1 m3f3          2        3  0.716891014       TRUE      TRUE
#> 6  m2f2 m3f1          2        3  0.995593231       TRUE      TRUE
#> 7  m2f2 m3f2          2        3  0.074615767      FALSE     FALSE
#> 8  m2f2 m3f3          2        3 -0.056803218      FALSE     FALSE
#> 9  m3f1 m4f1          3        4 -0.016458058      FALSE     FALSE
#> 10 m3f1 m4f2          3        4 -0.004812835      FALSE     FALSE
#> 11 m3f1 m4f3          3        4  0.663913488       TRUE      TRUE
#> 12 m3f1 m4f4          3        4  0.747612767       TRUE      TRUE
#> 13 m3f2 m4f1          3        4  0.940647736       TRUE      TRUE
#> 14 m3f2 m4f2          3        4  0.031735215      FALSE     FALSE
#> 15 m3f2 m4f3          3        4  0.262765620      FALSE     FALSE
#> 16 m3f2 m4f4          3        4 -0.212435736      FALSE     FALSE
#> 17 m3f3 m4f1          3        4  0.051257892      FALSE     FALSE
#> 18 m3f3 m4f2          3        4  0.970778427       TRUE      TRUE
#> 19 m3f3 m4f3          3        4 -0.171584098      FALSE     FALSE
#> 20 m3f3 m4f4          3        4  0.159752219      FALSE     FALSE
glance(x)
#>   engine rotation     cor k_max n_obs deepest_converged n_edges CFI TLI RMSEA
#> 1    pca  varimax pearson     4  1000                 4      20  NA  NA    NA
#>   SRMR BIC
#> 1   NA  NA

# Correlation-matrix input (PCA engine; n_obs optional)
R <- cor(sim16)
x_R <- ackwards(R, k_max = 4)
#> ℹ `n_obs` not supplied; stored as `NA`.
#> ℹ Fit statistics requiring N (chi-square, RMSEA, TLI) are unavailable. Pass
#>   `n_obs = <N>` to enable them.
print(x_R)
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: (user-supplied matrix)
#> n: NA
#> k (max): 4
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 28.2% variance
#> ✔ k = 2: 2 factors, 46.5% variance
#> ✔ k = 3: 3 factors, 57.5% variance
#> ✔ k = 4: 4 factors, 67.7% variance
#> 
#> ── Edges ──
#> 
#> 9 of 20 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
```
