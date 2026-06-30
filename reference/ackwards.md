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
  prune = "none",
  redundancy_r = 0.9,
  redundancy_phi = NULL,
  min_items = 3L,
  orphan_r = 0.5,
  cut_show = 0.3,
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
  if uncertain. Sets the depth of the hierarchy: levels 1 through
  `k_max` are all extracted.

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
  to override: `"ULSMV"` (unweighted WLS), `"MLR"` (robust ML). Ignored
  for PCA and EFA engines.

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

  - `"fiml"` – Full Information Maximum Likelihood. Passes
    `missing = "fiml"` to
    [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html); edge
    correlations are derived from lavaan's FIML-estimated saturated
    model, ensuring consistency. **Only valid for `engine = "esem"` with
    `estimator = "ML"` or `"MLR"`** – errors for PCA, EFA, and
    WLSMV/ULSMV. Note: FIML improves factor estimation under missingness
    but does not impute item responses; score materialisation
    (`keep_scores = TRUE`) still produces `NA` rows for incomplete
    observations.

- n_obs:

  Number of observations. Used only when `data` is a pre-computed
  correlation matrix. Ignored when raw data are supplied (N is
  determined from `nrow(data)`). For `engine = "efa"`, `n_obs` is
  required – [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html)
  needs N to compute fit indices (chi-square, RMSEA, TLI). For
  `engine = "pca"`, `n_obs` is optional; edges use only the `W'RW`
  algebra and never touch N (if `NULL`, stored as `NA_integer_` and fit
  statistics requiring N are unavailable).

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
  every pair of levels). `"all"` reveals associations that span multiple
  levels and is required for redundancy pruning. Setting `prune` to
  anything other than `"none"` automatically upgrades this to `"all"`
  with a message.

- prune:

  Character vector controlling Forbes-extension pruning. Default
  `"none"` (no pruning). Options:

  - `"redundant"` – identify chains of factors connected by
    primary-parent links with `|r| >= redundancy_r` (and optionally
    `phi > redundancy_phi`). Applies Forbes's (2023) retention rule:
    keep the bottom node when the chain reaches level `k` (most
    specific); keep the top node otherwise. Pruning is *flag-only*:
    flagged nodes stay in the object with `pruned = TRUE` and
    `prune_reason = "redundant"` in `x$prune$nodes`.

  - `"artefact"` – compute Tucker's congruence coefficient (phi) for all
    cross-level factor pairs and store in `x$prune$phi` for researcher
    inspection. No factors are auto-flagged; artefact identification
    requires judgment (Forbes, 2023; Wicherts et al., 2016).

- redundancy_r:

  Scalar in `(0, 1]`. Adjacent primary-parent `|r|` threshold for
  redundancy chains. Default `0.9` (Forbes, 2023).

- redundancy_phi:

  Scalar in `(0, 1]`, `NULL` (default, auto), or `NA` (explicit
  opt-out). When `NULL`:

  - `engine = "pca"` — no phi filter (the W'RW algebra is exact; phi
    adds nothing that \|r\| does not already capture).

  - `engine = "efa"` or `"esem"` — automatically set to `0.95`
    (Lorenzo-Seva & ten Berge, 2006). Factor-score indeterminacy off-PCA
    means \|r\|-alone is liberal; the conjunctive phi criterion is the
    conservative default. A cli message announces the resolved value
    (Invariant 6). Pass `NA` to disable phi filtering regardless of
    engine (matches the old `NULL` behaviour). Pass a numeric value to
    override on any engine.

- min_items:

  Minimum number of items for which a factor must be the primary loader
  (highest `|loading|`). Factors with fewer than `min_items` primary
  items are flagged `few_items = TRUE` in `x$prune$structural`. Only
  used when `prune = "artefact"`. Default `3L` – a factor defined by one
  or two items is under-identified and frequently an extraction artefact
  rather than a replicable construct (the classic "three-indicator
  rule"; Forbes, 2023, Fig. 2).

- orphan_r:

  Threshold for the `orphan` structural signal. A factor whose maximum
  **adjacent-level** `|r|` (to the immediately shallower and deeper
  levels) falls below `orphan_r` is flagged `orphan = TRUE` in
  `x$prune$structural` – it does not connect to the neighbouring
  solutions and so does not replicate across the hierarchy. Only used
  when `prune = "artefact"`. Default `0.5` – a moderate correlation; a
  factor that shares less than a quarter of its variance with every
  neighbour is a structural outlier worth inspecting.

- cut_show:

  Edges with `|r| >= cut_show` are flagged `above_cut` in
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) output.
  Default `0.3`.

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

## References

Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research
in Personality*, 40(4), 347–358.
[doi:10.1016/j.jrp.2006.01.001](https://doi.org/10.1016/j.jrp.2006.01.001)

Waller, N. G. (2007). A general method for computing hierarchical
component structures by Goldberg's bass-ackwards method. *Journal of
Research in Personality*, 41(4), 745–752.
[doi:10.1016/j.jrp.2006.08.005](https://doi.org/10.1016/j.jrp.2006.08.005)

## See also

[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)

## Examples

``` r
x <- ackwards(bfi25, k_max = 5)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> Warning: ! One or more columns look like ordinal/Likert items (… "<= 7" distinct integer
#>   values).
#> ℹ Results use a "pearson" basis. Consider `cor = "polychoric"` for ordinal
#>   data.
#> This warning is displayed once per session.
print(x)
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 20.6% variance
#> ✔ k = 2: 2 factors, 31.8% variance
#> ✔ k = 3: 3 factors, 40.1% variance
#> ✔ k = 4: 4 factors, 47.3% variance
#> ✔ k = 5: 5 factors, 53.3% variance
#> 
#> ── Edges ──
#> 
#> 14 of 40 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
tidy(x)
#>    from   to level_from level_to            r is_primary above_cut
#> 1  m1f1 m2f1          1        2  0.878066713       TRUE      TRUE
#> 2  m1f1 m2f2          1        2  0.478538241       TRUE      TRUE
#> 3  m2f1 m3f1          2        3  0.905019446       TRUE      TRUE
#> 4  m2f1 m3f2          2        3  0.084915431      FALSE     FALSE
#> 5  m2f1 m3f3          2        3  0.416808316       TRUE      TRUE
#> 6  m2f2 m3f1          2        3 -0.032217716      FALSE     FALSE
#> 7  m2f2 m3f2          2        3 -0.963373570       TRUE      TRUE
#> 8  m2f2 m3f3          2        3  0.266220556      FALSE     FALSE
#> 9  m3f1 m4f1          3        4  0.996745750       TRUE      TRUE
#> 10 m3f1 m4f2          3        4  0.006661471      FALSE     FALSE
#> 11 m3f1 m4f3          3        4  0.042095931      FALSE     FALSE
#> 12 m3f1 m4f4          3        4  0.068421245      FALSE     FALSE
#> 13 m3f2 m4f1          3        4 -0.016233327      FALSE     FALSE
#> 14 m3f2 m4f2          3        4  0.988176397       TRUE      TRUE
#> 15 m3f2 m4f3          3        4 -0.018059508      FALSE     FALSE
#> 16 m3f2 m4f4          3        4  0.151386067      FALSE     FALSE
#> 17 m3f3 m4f1          3        4 -0.070214243      FALSE     FALSE
#> 18 m3f3 m4f2          3        4 -0.061809300      FALSE     FALSE
#> 19 m3f3 m4f3          3        4  0.861695922       TRUE      TRUE
#> 20 m3f3 m4f4          3        4  0.498728091       TRUE      TRUE
#> 21 m4f1 m5f1          4        5  0.806394128       TRUE      TRUE
#> 22 m4f1 m5f2          4        5  0.001363543      FALSE     FALSE
#> 23 m4f1 m5f3          4        5  0.006259133      FALSE     FALSE
#> 24 m4f1 m5f4          4        5  0.589479673       TRUE      TRUE
#> 25 m4f1 m5f5          4        5 -0.046916846      FALSE     FALSE
#> 26 m4f2 m5f1          4        5 -0.030991419      FALSE     FALSE
#> 27 m4f2 m5f2          4        5  0.998569303       TRUE      TRUE
#> 28 m4f2 m5f3          4        5 -0.015523208      FALSE     FALSE
#> 29 m4f2 m5f4          4        5  0.040546852      FALSE     FALSE
#> 30 m4f2 m5f5          4        5  0.003723118      FALSE     FALSE
#> 31 m4f3 m5f1          4        5 -0.105831031      FALSE     FALSE
#> 32 m4f3 m5f2          4        5  0.006424971      FALSE     FALSE
#> 33 m4f3 m5f3          4        5  0.984798330       TRUE      TRUE
#> 34 m4f3 m5f4          4        5  0.135975178      FALSE     FALSE
#> 35 m4f3 m5f5          4        5  0.021012193      FALSE     FALSE
#> 36 m4f4 m5f1          4        5  0.191253503      FALSE     FALSE
#> 37 m4f4 m5f2          4        5  0.010261396      FALSE     FALSE
#> 38 m4f4 m5f3          4        5  0.025504732      FALSE     FALSE
#> 39 m4f4 m5f4          4        5 -0.185238698      FALSE     FALSE
#> 40 m4f4 m5f5          4        5  0.963510734       TRUE      TRUE
glance(x)
#>   engine rotation     cor k_max n_obs deepest_converged n_edges
#> 1    pca  varimax pearson     5  1000                 5      40

# Correlation-matrix input (PCA engine; n_obs optional)
R <- cor(bfi25, use = "pairwise.complete.obs")
x_R <- ackwards(R, k_max = 5)
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
#> k (max): 5
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 20.6% variance
#> ✔ k = 2: 2 factors, 31.8% variance
#> ✔ k = 3: 3 factors, 40.1% variance
#> ✔ k = 4: 4 factors, 47.3% variance
#> ✔ k = 5: 5 factors, 53.3% variance
#> 
#> ── Edges ──
#> 
#> 14 of 40 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```
