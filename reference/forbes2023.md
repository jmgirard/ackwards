# Assessing Mental Health symptom correlation matrix (Forbes 2023 applied example)

The 155 x 155 Spearman correlation matrix among 155 mental-health
symptom variables that forms the applied example in Forbes (2023). It is
a real, deep hierarchy: `ackwards(forbes2023, k_max = 10)` unfolds a
general factor of psychopathology at the top down to 10 fine-grained
components, the worked example that motivates the Forbes extension
(`pairs = "all"`,
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)).

## Usage

``` r
forbes2023
```

## Format

A 155 x 155 numeric matrix of Spearman correlations: symmetric, unit
diagonal, correlations in roughly `0.01`–`0.94`. Row and column names
are the 155 symptom-variable labels (e.g. `Impulsivity`, `Blurting`).

## Source

Forbes's OSF project for the 2023 paper (file `corSpearman_AMH.csv`),
<https://osf.io/pcwm8/>, redistributed here under its Creative Commons
Attribution 4.0 International (CC-BY 4.0) license (see the package's
`LICENSE.note`). The underlying data are from the Assessing Mental
Health study (Forbes et al., 2021).

## Details

Where `sim16` and `bfi25` are teaching foils, `forbes2023` is a
fidelity/reproduction dataset: a large, messy, published case bundled so
the package can reproduce the exact applied example analyzed in Forbes's
paper.

The correlations come from the Assessing Mental Health (AMH) study – the
Australian general-population sample (N = 3,175) of Forbes et al. (2021)
– spanning symptoms of 18 DSM disorders. Being a correlation matrix, it
carries no per-variable sample size; supply `n_obs = 3175` to
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
if you want EFA/ESEM fit statistics scaled to the original sample.

This matrix reproduces Forbes's published results exactly: the package
regression test `test-forbes-fidelity.R` runs
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
on this exported `forbes2023` and matches her reference implementation's
between-level correlations to `1.3e-14` across all 45 level-pairs at
`k_max = 10`.

To regenerate this dataset, run `source("data-raw/forbes2023.R")` from
the package root (downloads the source CSV from OSF).

## References

Forbes, M. K. (2023). Improving hierarchical models of individual
differences: An extension of Goldberg's bass-ackwards method.
*Psychological Methods*.
[doi:10.1037/met0000546](https://doi.org/10.1037/met0000546)

Forbes, M. K., Sunderland, M., Rapee, R. M., Batterham, P. J., Calear,
A. L., Carragher, N., Ruggero, C., Zimmerman, M., Baillie, A. J., Lynch,
S. J., Mewton, L., Slade, T., & Krueger, R. F. (2021). A detailed
hierarchical model of psychopathology: From individual symptoms up to
the general factor of psychopathology. *Clinical Psychological Science*,
9(2), 139–168.
[doi:10.1177/2167702620954799](https://doi.org/10.1177/2167702620954799)

## Examples

``` r
dim(forbes2023)
#> [1] 155 155
forbes2023[1:3, 1:3]
#>                        Impulsivity Reckless_behaviour__i_  Blurting
#> Impulsivity              1.0000000              0.5415115 0.6172762
#> Reckless_behaviour__i_   0.5415115              1.0000000 0.4783731
#> Blurting                 0.6172762              0.4783731 1.0000000
# \donttest{
# The Forbes (2023) applied example: a 10-level hierarchy from 155 symptoms.
x <- ackwards(forbes2023, k_max = 10, pairs = "all")
#> ℹ `n_obs` not supplied; stored as `NA`.
#> ℹ Fit statistics requiring N (chi-square, RMSEA, TLI) are unavailable. Pass
#>   `n_obs = <N>` to enable them.
x
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: (user-supplied matrix)
#> n: NA
#> k (max): 10
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 31.4% variance
#> ✔ k = 2: 2 factors, 36.8% variance
#> ✔ k = 3: 3 factors, 41.1% variance
#> ✔ k = 4: 4 factors, 44.9% variance
#> ✔ k = 5: 5 factors, 48.2% variance
#> ✔ k = 6: 6 factors, 50.8% variance
#> ✔ k = 7: 7 factors, 52.9% variance
#> ✔ k = 8: 8 factors, 54.6% variance
#> ✔ k = 9: 9 factors, 56.2% variance
#> ✔ k = 10: 10 factors, 57.6% variance
#> 
#> ── Edges ──
#> 
#> 312 of 1320 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
# }
```
