# Bass-ackwards hierarchical structural analysis

Extracts factor or component solutions at levels 1 through `k`. A factor
(or component) is a summary variable that stands in for a group of items
that move together. It then characterises the hierarchy by computing the
correlations between the factor scores of different levels. A factor
score is each person's estimated standing on a factor. The "hierarchy"
is descriptive: edges are score correlations, not a fitted higher-order
SEM.

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
  `engine` must be `"pca"` or `"efa"`, because ESEM (exploratory
  structural equation modeling) needs raw data. The `missing` and `cor`
  arguments are then ignored. See the *Correlation-matrix input* section
  below.

- k_max:

  Maximum number of factors or components to extract. It is required, so
  use
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  if uncertain. It sets the *depth* of the hierarchy: levels 1 through
  `k_max` are all extracted and retained. (The `k_max` of
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  means something related but different. It is the largest number of
  factors that function *evaluates* when recommending a depth, not a
  depth itself. See that function's docs.)

- engine:

  Extraction engine: `"pca"` (default), `"efa"`, or `"esem"`. `"esem"`
  uses [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) with
  rotation-aware SEs and per-level fit indices. ESEM is exploratory
  structural equation modeling, and it is recommended for the clinical
  or HiTOP workflow (Kim & Eaton, 2015, Forbush et al., 2024). It
  requires lavaan \>= 0.6-13.

- cor:

  Correlation basis: `"pearson"` (default), `"spearman"`, or
  `"polychoric"`. For PCA and EFA, `"polychoric"` computes a polychoric
  correlation matrix via
  [`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  (requires psych). Polychoric correlations estimate the correlation
  between the continuous traits assumed to underlie ordered responses.
  For ESEM, it triggers WLSMV estimation via lavaan. A cli warning is
  emitted when ordinal-looking columns are detected and `cor` is not
  `"polychoric"`. Ordinal items have a few ordered categories, such as a
  1 to 5 rating.

- fm:

  Factor extraction method passed to
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html), and used only
  when `engine = "efa"`. One of `"minres"` (the default, an ordinary
  least squares fit that converges reliably), `"ml"` (maximum
  likelihood, gives chi-square fit but converges less reliably at deep
  levels), or `"pa"` (principal axis). Ignored for `engine = "pca"`.

- estimator:

  Estimation method for the ESEM engine. `NULL` (default) auto-selects:
  `"WLSMV"` when `cor = "polychoric"`, `"ML"` otherwise. Pass explicitly
  to override: `"ULSMV"` (unweighted WLS), `"MLR"` (maximum likelihood
  with standard errors that tolerate non-normal data).
  `cor = "polychoric"` with `estimator = "ML"`/`"MLR"` errors, because
  lavaan itself does not support ML or MLR on ordered indicators.
  `"WLSMV"`/`"ULSMV"` with a continuous `cor` is allowed, and gives a
  valid, if atypical, continuous WLS or ADF estimator. Ignored for PCA
  and EFA engines. The effective value (after auto-selection) is
  recorded in `x$meta$estimator` (`NA` for PCA and EFA).

- missing:

  How to handle missing item responses. One of:

  - `"pairwise"` (the default) uses all available observations pairwise.
    For PCA and EFA this feeds
    `stats::cor(use = "pairwise.complete.obs")`. For ESEM with WLSMV or
    ULSMV (ordinal), lavaan uses `available.cases`, which computes
    polychoric thresholds and correlations from all rows that contribute
    to each pair. That is valid under MCAR and uses the full N. For ESEM
    with ML or MLR (continuous), lavaan uses listwise deletion
    internally while edges are computed from a pairwise correlation
    matrix. This small inconsistency is documented in `$meta`. A warning
    is emitted when incomplete rows are detected.

  - `"listwise"` uses only complete rows. It reduces the data to
    [`stats::complete.cases()`](https://rdrr.io/r/stats/complete.cases.html)
    before fitting, so the correlation matrix, the engine fit, and the
    edges are all consistent. `n_obs` in the result reflects the reduced
    N.

  - `"fiml"` is full information maximum likelihood. For
    `engine = "esem"` (with `estimator = "ML"`/`"MLR"`), it passes
    `missing = "fiml"` to
    [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html) and
    derives edge correlations from lavaan's FIML-estimated saturated
    model. For `engine = "pca"`/`"efa"` (M38), the correlation matrix is
    estimated via
    [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html) and
    fed to the usual `W'RW` algebra. That route requires
    `cor = "pearson"`, because corFiml estimates a multivariate-normal
    matrix, and the route is announced via a cli message. It errors for
    WLSMV and ULSMV, and for a non-Pearson PCA or EFA basis. FIML
    improves estimation under missingness but does not impute item
    responses. So score materialisation (`keep_scores = TRUE`) still
    produces `NA` rows for incomplete observations. See `n_obs` for the
    fit-index sample size on the PCA and EFA path.

- n_obs:

  Number of observations. On the raw-data FIML path it may instead be a
  string selecting which N feeds the fit indices. FIML is full
  information maximum likelihood, which uses every observed value
  without dropping incomplete rows.

  - **Correlation-matrix input:** a positive integer. It is required for
    `engine = "efa"`, because
    [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html) needs N for
    chi-square, RMSEA, and TLI. It is optional for `"pca"` (stored as
    `NA_integer_` if omitted). For PCA it is recorded in the result
    metadata and feeds the N-based sampling-adequacy checks only. PCA
    level fit is eigenvalue-based and computes no N-dependent fit
    statistics.

  - **Raw data:** N is normally taken from `nrow(data)` and a numeric
    `n_obs` is ignored (with a warning). The exception is
    `missing = "fiml"` with `engine = "pca"`/`"efa"` (M38).
    [`psych::corFiml()`](https://rdrr.io/pkg/psych/man/corFiml.html)
    estimates the correlation matrix from incomplete rows, so `n_obs`
    may be `"total"` or `"complete"`. The default is `"total"`, meaning
    every row that contributes to the FIML likelihood, which matches the
    FIML convention (Enders, 2010). `"complete"` is the complete-case N,
    a conservative lower bound. Point estimates do not depend on this
    choice. Those are the loadings (the correlation between each item
    and a factor) and the edges. Only the fit indices of EFA
    (exploratory factor analysis) do. Those indices are *approximate*
    whatever N you pick, because of this two-step route: a FIML matrix
    fed into normal-theory EFA (Zhang & Savalei, 2020). A string `n_obs`
    is accepted only on this path.

- align_signs:

  Logical. Sign-align factors to primary-parent lineage? Default `TRUE`.

- keep_scores:

  Logical. Store factor scores in the result? Default `FALSE`
  (recomputable via
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)).
  When `TRUE`, per-observation scores are stored in `x$scores` as a
  named list of `n x k_j` matrices, one per level, standardized by real
  score SDs (see
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)).

- keep_fits:

  Logical. Store raw engine fit objects? Default `FALSE`. When `TRUE`,
  the per-level fit objects (psych or lavaan) are stored in `x$fits` as
  a named list indexed by level.

- seed:

  Integer seed for stochastic engines (not used by PCA but captured for
  reproducibility metadata). Default `NULL`.

- pairs:

  Which level pairs to compute edges for. `"adjacent"` is the default
  and is the classic Goldberg choice of consecutive levels only. `"all"`
  is the Forbes extension, which takes every pair of levels, including
  skip-level correlations. The
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  function recomputes its own all-pairs edges on demand whatever this
  setting is, so pruning does not require `pairs = "all"` here.

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
  fails on your data, and its error suggests exactly this. That failure
  typically happens when an item has a near-empty response category, or
  when items with unequal category counts produce a sparse cross-cell.
  It is ignored on other paths: ESEM computes its own polychoric
  correlations inside lavaan, and the Pearson/Spearman bases do not use
  it.

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

- **`engine = "pca"`** is the original Goldberg (2006) method. PCA
  (principal component analysis) is the fastest engine and never fails
  to converge, and the Waller (2007) algebra is exact for components.

- **`rotation = "varimax"`** keeps the within-level factors mutually
  uncorrelated (orthogonal). A rotation re-orients the factors without
  changing how well they fit, and varimax pushes each item toward one
  factor. So each between-level edge reflects only the cross-level
  relationship. An oblique rotation's correlated within-level factors
  would leak into the between-level edges and confound the between-level
  signal that is the method's core output. The closed-form `W'RW` edge
  algebra is itself exact for any fixed linear scoring, orthogonal or
  not, so the choice is interpretive, not a numerical necessity. It
  matches Goldberg (2006), Kim & Eaton (2015), and Forbush et al.
  (2024). Varimax is the only supported rotation.

- **`cor = "pearson"`** means no silent basis switching. If your items
  look ordinal (\<= 7 distinct integer values), a cli warning will
  suggest `cor = "polychoric"`, which is available for all three
  engines. Ordinal items have a few ordered categories, such as a 1 to 5
  rating. Polychoric correlations estimate the correlation between the
  continuous traits assumed to underlie those ordered responses.

- **`align_signs = TRUE`** because unaligned signs make the output
  unreadable. The anchor is m1f1, which is oriented toward the positive
  manifold. Each subsequent factor is flipped so its edge to its primary
  parent is positive.

- **`keep_scores = FALSE` / `keep_fits = FALSE`** for memory and
  privacy. Scores are O(n x Sigmak) and often sensitive, and raw engine
  fits can be large. Both are recomputable from the stored `r` matrix.

## Performance (ESEM, large item sets)

The ESEM engine fits a separate `lavaan` model at every level
1..`k_max`. For ordinal data (`cor = "polychoric"`, WLSMV) lavaan
derives costly sample statistics from the raw data. Those are the
thresholds, the polychoric correlation matrix, and the asymptotic weight
matrix. They depend only on the data, not on the number of factors. So
they are **computed once** at the first level and **reused** for every
deeper level, which gives identical solutions for much less work. This
matters most when you have many items (hundreds), where recomputing
those statistics at each level dominated the run time.

The per-level model fits are mutually independent and are dispatched
through the future framework when future.apply is installed. By default
the plan is sequential (no behaviour change). To run the levels in
parallel, set a plan once before calling `ackwards()`:


      future::plan(future::multisession, workers = 4)  # or multicore on Unix
      x <- ackwards(items, k_max = 8, engine = "esem", cor = "polychoric")

Parallelism pays off when the per-level fits are heavy (large `p`,
several levels). For small problems the worker startup cost can outweigh
it. Results are reproducible across plans when `seed` is supplied. PCA
and EFA already compute their correlation matrix once and are
unaffected.

## Correlation-matrix input

When `data` is a pre-computed correlation matrix (square, symmetric,
unit diagonal), `ackwards()` runs entirely from that matrix using the
`W'RW` algebra, so no raw item responses are needed. This is useful when
you have a published correlation table or a polychoric matrix computed
externally.

Constraints and behaviour when a correlation matrix is supplied:

- **Engine:** only `"pca"` and `"efa"` are supported. `"esem"` requires
  raw data (for lavaan's own polychoric computation, WLSMV estimation,
  and per-level fit indices) and will error clearly.

- **`n_obs`:** required for `"efa"`, because psych needs N for
  chi-square, RMSEA, and TLI. It is optional for `"pca"`, where it is
  stored as `NA` if omitted and is used only for the N-based
  sampling-adequacy checks and the result metadata. PCA computes no
  N-dependent fit statistics.

- **`cor` argument:** ignored, because the basis is already determined
  by the matrix you supply. A warning is emitted if you set `cor`
  explicitly.

- **`missing` argument:** ignored, because missingness was handled when
  computing the matrix. A warning is emitted if you set `missing`
  explicitly.

- **Factor scores:** `keep_scores = TRUE` will error.
  [`augment()`](https://generics.r-lib.org/reference/augment.html) and
  `tidy(what = "scores")` will also error because individual-level
  scores require row-level item responses.

- **`$cor` field:** stored as `NA_character_`, and printed as
  `"(user-supplied matrix)"`.

## When to trust the result

`ackwards()` raises diagnostics as it fits. They fall into three tiers
by what they mean for whether you should trust and report the solution:

**Fatal, so fix before trusting.** The result is undefined or rests on a
broken correlation matrix:

- A **constant item** (no variance) errors. Drop it (see
  [`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md)).

- **[`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
  fails**, usually because of a near-empty response category. Set
  `correct = 0` or collapse rare categories.

- A **level fails to converge**. The hierarchy is truncated to the
  deepest level that did converge, and you should not interpret beyond
  it.

- A **near-singular correlation matrix** (smallest eigenvalue `< 1e-4`,
  recorded in `meta$near_singular` / `meta$min_eigenvalue` and
  re-surfaced by
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html))
  means per-level fit indices and factor scores are unreliable. The
  loadings and edges then rest on a rank-deficient matrix. (For EFA the
  residual-based fallback inflates `TLI` and `RMSEA`, and for ESEM `CFI`
  comes back `NA`.) Trim redundant items, or use `missing = "listwise"`.
  On the polychoric basis you can instead collapse sparse categories or
  set `correct = 0`.

**Caution, so interpret carefully.** The solution exists but may be
unstable:

- A **Heywood case** at a level, meaning a communality `> 1` or a
  negative uniqueness. Check that level's loadings make substantive
  sense, and consider fewer factors.

- A **near-constant item** (one response category dominates) can drive a
  meaningless factor. Inspect it with
  [`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md).

- **Ordinal data on a Pearson basis** attenuates correlations. That is
  fine for a quick look or
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  screening, but report the final model on `cor = "polychoric"`.

**Informational, and usually fine.** Proceed, just be aware of the
pairwise-missing note, a merely *sparse* (rare-but-present) response
category, and the ordinal-detection warning when you did intend Pearson.

## References

Goldberg, L. R. (2006). Doing it all Bass-Ackwards: The development of
hierarchical factor structures from the top down. *Journal of Research
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

Kaiser, H. F. (1958). The varimax criterion for analytic rotation in
factor analysis. *Psychometrika*, 23(3), 187–200.
[doi:10.1007/BF02289233](https://doi.org/10.1007/BF02289233)

Grice, J. W. (2001). Computing and evaluating factor scores.
*Psychological Methods*, 6(4), 430–450.
[doi:10.1037/1082-989X.6.4.430](https://doi.org/10.1037/1082-989X.6.4.430)

Beauducel, A., Hilger, N., & Kuhl, T. (2024). The trade-off between
factor score determinacy and the preservation of inter-factor
correlations. *Educational and Psychological Measurement*, 84(2),
289–313.
[doi:10.1177/00131644231171137](https://doi.org/10.1177/00131644231171137)

Williams, A. L., Conway, C. C., Olino, T. M., Revelle, W., Zinbarg, R.
E., & HiTOP Utility Workgroup. (2025). Testing criterion validity in
hierarchical models of psychopathology: Comparison of latent-variable
and factor-score approaches. *Clinical Psychological Science*, 13(1),
128–145.
[doi:10.1177/21677026231225414](https://doi.org/10.1177/21677026231225414)

## See also

[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md),
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md),
which is the Forbes-extension flagging of redundancy (one factor adding
nothing over another) and of artifacts, piped off the result of this
function

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
#> ℹ PCA level fit is eigenvalue-based and does not use N; supplying `n_obs = <N>`
#>   records it in the result metadata and enables the N-based sampling-adequacy
#>   checks.
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
