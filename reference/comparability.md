# Split-half factor comparability

Measures how well each factor at each level of a bass-ackwards hierarchy
**replicates** across random split-halves of the sample – Everett's
(1983) factor comparability coefficients, the depth gate Goldberg's own
lab applied to its hierarchies (Goldberg, 1990) and the direct
instrument for the overextraction caution in
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md):
non-replicable structure concentrates in the deeper levels of an
overextracted hierarchy (Forbes, 2023).

## Usage

``` r
comparability(
  data,
  k_max,
  engine = "pca",
  cor = "pearson",
  fm = "minres",
  n_splits = 10L,
  seed = NULL,
  ...
)
```

## Arguments

- data:

  A data frame or numeric matrix of observed variables (items in
  columns, observations in rows). **Raw data only** – splitting needs
  rows, so a correlation matrix is not accepted. Missing values are
  handled pairwise throughout (as in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)'s
  default).

- k_max:

  Maximum number of factors/components to evaluate – normally the same
  value (or one or two above it) you intend to pass to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).
  Required.

- engine:

  Extraction engine: `"pca"` (default) or `"efa"`. `"esem"` is not yet
  supported here (fitting `2 * n_splits` lavaan hierarchies needs its
  own performance treatment); for ESEM workflows, run `comparability()`
  with `engine = "efa"` as a structural screen.

- cor:

  Correlation basis: `"pearson"` (default) or `"spearman"`. As with
  [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md),
  `"polychoric"` is not supported (estimating polychoric matrices in
  every half-sample is slow and unstable); users analysing ordinal data
  should screen replicability on the Pearson basis and fit the final
  model with `cor = "polychoric"` in
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

- fm:

  Factor extraction method passed to
  [`psych::fa()`](https://rdrr.io/pkg/psych/man/fa.html); only used when
  `engine = "efa"`. One of `"minres"` (default), `"ml"`, or `"pa"`.

- n_splits:

  Number of random split-half replicates. Default `10L` (repeated random
  splits, following Goldberg's practice – a single split can mislead by
  luck of the draw). Each replicate fits `2 * k_max` solutions, so the
  default costs 20 hierarchy fits; PCA and EFA are fast enough that this
  is typically a few seconds.

- seed:

  Integer seed for reproducible splits. `NULL` (default) uses the
  current RNG state.

- ...:

  Reserved for future arguments.

## Value

An object of class `"comparability"`. Print it for a per-level summary;
call
[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
on it for a diagnostic plot. The list contains:

- coefficients:

  Data frame with one row per split x level x factor: `split`, `level`,
  `factor` (full-sample `m{k}f{j}` label), `r` (score comparability),
  `phi` (Tucker's congruence of the matched loading columns). `NA` when
  the level did not converge in one of the halves.

- summary:

  Data frame with one row per level x factor: `level`, `factor`,
  `r_median`, `r_min`, `phi_median`, `phi_min` (across splits), and
  `n_splits_ok` (splits in which both halves converged).

- k_max:

  Deepest level evaluated. Can be lower than the `k_max` you asked for
  when the full-sample fit truncated (non-convergence at deep levels);
  the original request is kept in `k_requested`.

- k_requested, n_splits, n_half, engine, cor, fm, n_obs, n_vars, seed:

  Metadata.

## Details

For each of `n_splits` random half-splits, solutions at every level
`1..k_max` are fit independently in each half. Each half-solution's
factors are matched to the **full-sample** solution's factors (so
coefficients are reported under the same `m{k}f{j}` labels you get from
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)),
and the comparability coefficient for a factor is the correlation
between its two matched half-solution scores, computed on the pooled
correlation matrix via the same `W'RW` algebra used for between-level
edges – applying both halves' scoring weights to the full sample,
exactly Everett's procedure. Tucker's congruence coefficient (phi)
between the matched half-solution loading columns is reported alongside:
comparability asks whether the two halves' *scores* agree; phi asks
whether their *loading patterns* agree.

## Interpreting the output

Coefficients near 1 mean the factor re-emerges in independent
half-samples; a factor whose comparability is low is
sample-idiosyncratic and should not anchor substantive interpretation.
Goldberg's lab treated roughly .95 as comfortable and .90 as a floor;
these are conventions, not tests, so `comparability()` reports every
coefficient and flags nothing. The deepest level at which all factors
replicate is a natural **hierarchy floor** for
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)'s
`k_max`; see
[`vignette("ackwards-girard")`](https://jmgirard.github.io/ackwards/articles/ackwards-girard.md)
for the full workflow.

A level that fails to converge in a half-sample yields `NA` coefficients
for that split (convergence is data, not an error); the number of usable
splits per factor is reported in `summary$n_splits_ok` and a message
summarises any shortfall.

## References

Everett, J. E. (1983). Factor comparability as a means of determining
the number of factors and their rotation. *Multivariate Behavioral
Research*, 18(2), 197–218.
[doi:10.1207/s15327906mbr1802_5](https://doi.org/10.1207/s15327906mbr1802_5)

Goldberg, L. R. (1990). An alternative "description of personality": The
Big-Five factor structure. *Journal of Personality and Social
Psychology*, 59(6), 1216–1229.
[doi:10.1037/0022-3514.59.6.1216](https://doi.org/10.1037/0022-3514.59.6.1216)

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

## See also

[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
for the plausible depth *range* (eigenstructure),
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) for
factors that perpetuate without differentiating (redundancy), and
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
for the extraction itself.

## Examples

``` r
# \donttest{
cmp <- comparability(bfi25, k_max = 5, n_splits = 5, seed = 1)
#> Warning: ! 125 rows have missing values; correlations are computed pairwise.
#> ℹ Use `missing = "listwise"` for consistent complete-case analysis.
#> ℹ Fitting 5 split-half replicates (pca, k = 1-5)...
#> ✔ Fitting 5 split-half replicates (pca, k = 1-5)... [679ms]
#> 
cmp
#> 
#> ── Split-Half Factor Comparability (ackwards) ──────────────────────────────────
#> Engine: pca
#> Basis: pearson
#> n: 1,000 (500 per half)
#> Splits: 5
#> Levels: 1-5
#> 
#> ── Comparability by level (median across splits) ──
#> 
#> k = 1: median r 1.00, min r 1.00 (m1f1)
#> k = 2: median r .97, min r .96 (m2f2)
#> k = 3: median r .98, min r .97 (m3f3)
#> k = 4: median r .98, min r .94 (m4f4)
#> k = 5: median r .97, min r .93 (m5f4)
#> ────────────────────────────────────────────────────────────────────────────────
#> Per-factor detail (incl. Tucker's φ) in `$summary`; per-split values in
#> `$coefficients`.
#> Conventional benchmarks: ≥ .95 comfortable, ≥ .90 floor (Everett, 1983;
#> Goldberg, 1990) -- conventions, not tests. Interpret levels whose factors all
#> replicate.
cmp$summary
#>    level factor  r_median     r_min phi_median   phi_min n_splits_ok
#> 1      1   m1f1 0.9971703 0.9927152  0.9902854 0.9835375           5
#> 2      2   m2f1 0.9854606 0.9797007  0.9873563 0.9659103           5
#> 3      2   m2f2 0.9628712 0.9489709  0.9608554 0.9381347           5
#> 4      3   m3f1 0.9846760 0.9558333  0.9841539 0.9534736           5
#> 5      3   m3f2 0.9875921 0.9714452  0.9776666 0.9603330           5
#> 6      3   m3f3 0.9662627 0.9103907  0.9517414 0.9081449           5
#> 7      4   m4f1 0.9850333 0.7238259  0.9807810 0.8407559           5
#> 8      4   m4f2 0.9925570 0.9824683  0.9840938 0.9660301           5
#> 9      4   m4f3 0.9764837 0.9430105  0.9599242 0.9056327           5
#> 10     4   m4f4 0.9440768 0.6169249  0.9266956 0.6122198           5
#> 11     5   m5f1 0.9646207 0.9038571  0.9610989 0.9349538           5
#> 12     5   m5f2 0.9914512 0.9876653  0.9879514 0.9741450           5
#> 13     5   m5f3 0.9906372 0.9857630  0.9775060 0.9588205           5
#> 14     5   m5f4 0.9316262 0.8307607  0.9488651 0.8176100           5
#> 15     5   m5f5 0.9664983 0.9605977  0.9461381 0.9291701           5
# }
```
