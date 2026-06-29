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
  align_signs = TRUE,
  keep_scores = FALSE,
  keep_fits = FALSE,
  seed = NULL,
  pairs = "adjacent",
  prune = "none",
  redundancy_r = 0.9,
  redundancy_phi = NULL,
  cut_show = 0.3,
  ...
)
```

## Arguments

- data:

  A data frame or numeric matrix of observed variables (items in
  columns, observations in rows). How missing values are handled depends
  on the `missing` argument; see below.

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

  - `"pairwise"` (default) — use all available observations pairwise.
    For PCA/EFA this feeds `stats::cor(use = "pairwise.complete.obs")`.
    For ESEM with WLSMV/ULSMV (ordinal), lavaan uses `available.cases`,
    which computes polychoric thresholds and correlations from all rows
    that contribute to each pair — MCAR-valid and uses the full N. For
    ESEM with ML/MLR (continuous), lavaan uses listwise deletion
    internally while edges are computed from a pairwise correlation
    matrix; this minor inconsistency is documented in `$meta`. A warning
    is emitted when incomplete rows are detected.

  - `"listwise"` — only complete rows are used. Reduces data to
    [`stats::complete.cases()`](https://rdrr.io/r/stats/complete.cases.html)
    before fitting, so the correlation matrix, the engine fit, and the
    edges are all consistent. `n_obs` in the result reflects the reduced
    N.

  - `"fiml"` — Full Information Maximum Likelihood. Passes
    `missing = "fiml"` to
    [`lavaan::efa()`](https://rdrr.io/pkg/lavaan/man/efa.html); edge
    correlations are derived from lavaan's FIML-estimated saturated
    model, ensuring consistency. **Only valid for `engine = "esem"` with
    `estimator = "ML"` or `"MLR"`** — errors for PCA, EFA, and
    WLSMV/ULSMV. Note: FIML improves factor estimation under missingness
    but does not impute item responses; score materialisation
    (`keep_scores = TRUE`) still produces `NA` rows for incomplete
    observations.

- align_signs:

  Logical; sign-align factors to primary-parent lineage? Default `TRUE`.

- keep_scores:

  Logical; store factor scores in the result? Default `FALSE`
  (recomputable via
  [`augment.ackwards()`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)).
  When `TRUE`, per-observation scores are stored in `x$scores` as a
  named list of `n × k_j` matrices, one per level, standardized by real
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
  Goldberg — only consecutive levels) or `"all"` (Forbes extension —
  every pair of levels). `"all"` reveals associations that span multiple
  levels and is required for redundancy pruning. Setting `prune` to
  anything other than `"none"` automatically upgrades this to `"all"`
  with a message.

- prune:

  Character vector controlling Forbes-extension pruning. Default
  `"none"` (no pruning). Options:

  - `"redundant"` — identify chains of factors connected by
    primary-parent links with `|r| >= redundancy_r` (and optionally
    `phi > redundancy_phi`). Applies Forbes's (2023) retention rule:
    keep the bottom node when the chain reaches level `k` (most
    specific); keep the top node otherwise. Pruning is *flag-only*:
    flagged nodes stay in the object with `pruned = TRUE` and
    `prune_reason = "redundant"` in `x$prune$nodes`.

  - `"artefact"` — compute Tucker's congruence coefficient (φ) for all
    cross-level factor pairs and store in `x$prune$phi` for researcher
    inspection. No factors are auto-flagged; artefact identification
    requires judgment (Forbes, 2023; Wicherts et al., 2016).

- redundancy_r:

  Scalar in `(0, 1]`. Adjacent primary-parent `|r|` threshold for
  redundancy chains. Default `0.9` (Forbes, 2023).

- redundancy_phi:

  Scalar in `(0, 1]` or `NULL` (default). If non-`NULL`, Tucker's φ must
  *also* exceed this threshold for a link to be included in a redundancy
  chain (conjunctive with `redundancy_r`). `NULL` means only
  `redundancy_r` is used. Recommended: `0.95` (Lorenzo-Seva & ten Berge,
  2006).

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

- **`engine = "pca"`** — the original Goldberg (2006) method; fastest;
  never fails to converge; the Waller (2007) algebra is exact for
  components.

- **`rotation = "varimax"`** — the `T'=T^-1` property of orthogonal
  rotation enables the closed-form `W'RW` edge algebra and keeps
  within-level factors uncorrelated so cross-level edges reflect only
  the hierarchical signal. Matches Goldberg (2006), Kim & Eaton (2015),
  and Forbush et al. (2024). Varimax is the only supported rotation;
  oblique rotation would confound the between-level signal that is the
  method's core output.

- **`cor = "pearson"`** — no silent basis switching. If your items look
  ordinal (≤ 7 distinct integer values), a cli warning will suggest
  `cor = "polychoric"`, which is available for all three engines.

- **`align_signs = TRUE`** — unaligned signs make the output unreadable.
  Anchor: m1f1 is oriented toward the positive manifold; each subsequent
  factor is flipped so its edge to its primary parent is positive.

- **`keep_scores = FALSE` / `keep_fits = FALSE`** — memory and privacy.
  Scores are O(n × Σk) and often sensitive; raw engine fits can be
  large. Both are recomputable from the stored `r` matrix.

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
if (FALSE) { # \dontrun{
x <- ackwards(psych::bfi[, 1:25], k_max = 5)
print(x)
tidy(x)
glance(x)
} # }
```
