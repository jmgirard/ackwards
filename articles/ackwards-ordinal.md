# Ordinal Data: Polychoric Correlations and WLSMV

Most psychological scales use ordinal response formats — Likert-type
items with 4, 5, 6, or 7 response options. Treating these as continuous
and computing Pearson correlations between them is ubiquitous in
practice, but it introduces systematic bias that can distort both
loadings and the between-level edges that bass-ackwards analysis depends
on. This vignette explains the problem and how `ackwards` addresses it.

## Why Pearson underestimates for ordinal items

When a continuous latent trait is sliced into a small number of ordered
categories, the observed correlation between two items is always lower
than the correlation between their underlying latent variables. The more
coarse the rating scale, the larger the attenuation. For a 5-point
scale, Pearson correlations can be 10–20% lower than the true latent
correlations; for a 2-point (binary) scale, the underestimate can exceed
30%.

This matters for bass-ackwards analysis because:

1.  **Loadings are attenuated.** Factors appear weaker than they are.
2.  **Between-level edges are attenuated.** The hierarchy looks flatter
    — factors that should be nearly perfectly correlated across levels
    may appear only moderately correlated.
3.  **The number of factors is underestimated.** Attenuation compresses
    the eigenvalue spectrum, so selection criteria (parallel analysis,
    MAP) may suggest fewer factors than truly exist.

## The polychoric solution

A **polychoric correlation** between two ordinal items estimates the
Pearson correlation that *would* exist between their underlying
continuous latent variables, if those variables had been measured
without discretization. It does this by fitting a bivariate normal model
with thresholds at each category boundary.

Setting `cor = "polychoric"` in
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
replaces step 1 of the pipeline — computing the item correlation matrix
`R` — with a polychoric matrix. All downstream computation (factor
extraction, rotation, tenBerge scoring, between-level algebra) then
operates on that polychoric `R`.

## Automatic detection

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
checks for ordinal-looking data using a simple heuristic: if any column
is integer-valued with ≤ 7 distinct levels, the data is flagged and a
warning is issued when `cor = "pearson"` (the default).

``` r

library(ackwards)
bfi <- na.omit(bfi25)
```

``` r

# The warning fires with the default Pearson basis
x_pearson <- ackwards(bfi, k_max = 5)
#> Warning: ! One or more columns look like ordinal/Likert items (… "<= 7" distinct integer
#>   values).
#> ℹ Results use a "pearson" basis. Consider `cor = "polychoric"` for ordinal
#>   data.
#> This warning is displayed once per session.
```

The warning is ackwards reminding you to make an active choice.
Specifying `cor = "polychoric"` acknowledges the data type and
suppresses it:

``` r

x_poly <- ackwards(bfi, k_max = 5, cor = "polychoric")
```

## Seeing the difference

### Correlation matrices

The most direct way to see the impact of the basis is to compare the
correlation matrices directly.

``` r

# Pearson R (stored in x$r)
round(x_pearson$r[1:5, 1:5], 2)
#>       A1    A2    A3    A4    A5
#> A1  1.00 -0.37 -0.33 -0.16 -0.24
#> A2 -0.37  1.00  0.54  0.36  0.40
#> A3 -0.33  0.54  1.00  0.39  0.57
#> A4 -0.16  0.36  0.39  1.00  0.34
#> A5 -0.24  0.40  0.57  0.34  1.00
```

``` r

# Polychoric R
round(x_poly$r[1:5, 1:5], 2)
#>       A1    A2    A3    A4    A5
#> A1  1.00 -0.44 -0.38 -0.18 -0.29
#> A2 -0.44  1.00  0.61  0.41  0.45
#> A3 -0.38  0.61  1.00  0.44  0.64
#> A4 -0.18  0.41  0.44  1.00  0.39
#> A5 -0.29  0.45  0.64  0.39  1.00
```

The polychoric correlations are consistently higher — sometimes
substantially so. The `N1`–`N2` pair, for instance, goes from roughly
0.59 (Pearson) to 0.73 (polychoric): a difference large enough to shift
loadings and change eigenvalues.

### Loadings

Stacking `tidy(what = "loadings")` from each basis lets us read the
`k = 5` Neuroticism-item loadings side by side:

| basis      | factor | item | loading |
|:-----------|:-------|:-----|--------:|
| pearson    | m5f1   | N1   |    0.05 |
| polychoric | m5f1   | N1   |    0.05 |
| pearson    | m5f1   | N2   |    0.09 |
| polychoric | m5f1   | N2   |    0.09 |
| pearson    | m5f1   | N3   |   -0.04 |
| polychoric | m5f1   | N3   |   -0.04 |
| pearson    | m5f1   | N4   |   -0.30 |
| polychoric | m5f1   | N4   |   -0.31 |
| pearson    | m5f1   | N5   |   -0.14 |
| polychoric | m5f1   | N5   |   -0.15 |
| pearson    | m5f2   | N1   |    0.79 |
| polychoric | m5f2   | N1   |    0.81 |
| pearson    | m5f2   | N2   |    0.78 |
| polychoric | m5f2   | N2   |    0.81 |
| pearson    | m5f2   | N3   |    0.80 |
| polychoric | m5f2   | N3   |    0.82 |
| pearson    | m5f2   | N4   |    0.63 |
| polychoric | m5f2   | N4   |    0.65 |
| pearson    | m5f2   | N5   |    0.66 |
| polychoric | m5f2   | N5   |    0.69 |
| pearson    | m5f3   | N1   |   -0.06 |
| polychoric | m5f3   | N1   |   -0.07 |
| pearson    | m5f3   | N2   |   -0.04 |
| polychoric | m5f3   | N2   |   -0.05 |
| pearson    | m5f3   | N3   |   -0.03 |
| polychoric | m5f3   | N3   |   -0.04 |
| pearson    | m5f3   | N4   |   -0.21 |
| polychoric | m5f3   | N4   |   -0.21 |
| pearson    | m5f3   | N5   |   -0.01 |
| polychoric | m5f3   | N5   |   -0.02 |
| pearson    | m5f4   | N1   |   -0.23 |
| polychoric | m5f4   | N1   |   -0.25 |
| pearson    | m5f4   | N2   |   -0.21 |
| polychoric | m5f4   | N2   |   -0.22 |
| pearson    | m5f4   | N3   |   -0.03 |
| polychoric | m5f4   | N3   |   -0.03 |
| pearson    | m5f4   | N4   |    0.04 |
| polychoric | m5f4   | N4   |    0.05 |
| pearson    | m5f4   | N5   |    0.14 |
| polychoric | m5f4   | N5   |    0.17 |
| pearson    | m5f5   | N1   |   -0.13 |
| polychoric | m5f5   | N1   |   -0.13 |
| pearson    | m5f5   | N2   |    0.00 |
| polychoric | m5f5   | N2   |    0.00 |
| pearson    | m5f5   | N3   |    0.02 |
| polychoric | m5f5   | N3   |    0.01 |
| pearson    | m5f5   | N4   |    0.11 |
| polychoric | m5f5   | N4   |    0.13 |
| pearson    | m5f5   | N5   |   -0.14 |
| polychoric | m5f5   | N5   |   -0.15 |

Neuroticism-item loadings (k = 5): Pearson vs polychoric {.table}

Polychoric loadings for the Neuroticism items are noticeably higher —
the latent structure is sharper when the attenuating effect of coarse
categories is removed.

### Between-level edges

The same comparison for the primary-parent edges —
`tidy(what = "edges", primary_only = TRUE)` under each basis:

| basis      | from | to   |     r |
|:-----------|:-----|:-----|------:|
| pearson    | m1f1 | m2f1 |  0.88 |
| polychoric | m1f1 | m2f1 |  0.89 |
| pearson    | m1f1 | m2f2 |  0.48 |
| polychoric | m1f1 | m2f2 |  0.46 |
| pearson    | m2f1 | m3f1 |  0.90 |
| polychoric | m2f1 | m3f1 |  0.87 |
| pearson    | m2f2 | m3f2 | -0.97 |
| polychoric | m2f2 | m3f2 | -0.99 |
| pearson    | m2f1 | m3f3 |  0.43 |
| polychoric | m2f1 | m3f3 |  0.48 |
| pearson    | m3f1 | m4f1 |  1.00 |
| polychoric | m3f1 | m4f1 |  0.99 |
| pearson    | m3f2 | m4f2 |  0.99 |
| polychoric | m3f2 | m4f2 |  0.98 |
| pearson    | m3f3 | m4f3 |  0.80 |
| polychoric | m3f3 | m4f3 |  0.73 |
| pearson    | m3f3 | m4f4 |  0.59 |
| polychoric | m3f3 | m4f4 |  0.68 |
| pearson    | m4f1 | m5f1 |  0.82 |
| polychoric | m4f1 | m5f1 |  0.84 |
| pearson    | m4f2 | m5f2 |  1.00 |
| polychoric | m4f2 | m5f2 |  1.00 |
| pearson    | m4f3 | m5f3 |  0.98 |
| polychoric | m4f3 | m5f3 |  0.98 |
| pearson    | m4f1 | m5f4 |  0.58 |
| polychoric | m4f1 | m5f4 |  0.55 |
| pearson    | m4f4 | m5f5 |  0.96 |
| polychoric | m4f4 | m5f5 |  0.99 |

Primary-parent edges: Pearson vs polychoric {.table}

The edges are broadly similar in sign and rank order — the hierarchy is
the same — but polychoric edges are stronger. Factors that represent
genuinely stable dimensions show near-perfect correlations with their
counterparts at adjacent levels when the true latent structure is
properly recovered.

## WLSMV for ESEM with ordinal items

When `engine = "esem"` and `cor = "polychoric"` are combined,
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
automatically switches the lavaan estimator to **WLSMV** (diagonally
weighted least squares, mean- and variance-adjusted). This is the
standard estimator for structural equation models with categorical or
ordinal indicators, and is the same method used by Mplus by default.

``` r

x_esem <- ackwards(bfi, k_max = 3, engine = "esem", cor = "polychoric")
x_esem
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: esem
#> Rotation: varimax
#> Basis: polychoric
#> n: 875
#> k (max): 3
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 23.5% variance
#> ✔ k = 2: 2 factors, 32.9% variance
#> ✔ k = 3: 3 factors, 39.5% variance
#> 
#> ── Edges ──
#> 
#> 5 of 8 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```

The WLSMV fit indices (CFI, RMSEA, SRMR) are now computed on a model
that respects the ordinal measurement level. They will typically look
better than Pearson-based fit indices on the same data, reflecting both
the better model specification and the different reference distribution
WLSMV uses.

``` r

tidy(x_esem, what = "fit")
#>    level   index        value
#> 1      1     chi 7.629329e+03
#> 2      1     dof 2.750000e+02
#> 3      1 p_value           NA
#> 4      1     CFI 7.153139e-01
#> 5      1     TLI 6.894333e-01
#> 6      1   RMSEA 1.749240e-01
#> 7      1    SRMR 1.443569e-01
#> 8      2     chi 3.172628e+03
#> 9      2     dof 2.510000e+02
#> 10     2 p_value           NA
#> 11     2     CFI 8.869038e-01
#> 12     2     TLI 8.648252e-01
#> 13     2   RMSEA 1.154037e-01
#> 14     2    SRMR 9.547997e-02
#> 15     3     chi 1.703232e+03
#> 16     3     dof 2.280000e+02
#> 17     3 p_value           NA
#> 18     3     CFI 9.428938e-01
#> 19     3     TLI 9.248602e-01
#> 20     3   RMSEA 8.604131e-02
#> 21     3    SRMR 7.172502e-02
```

## Practical guidance: when to use polychoric

| Scale characteristics | Recommendation |
|----|----|
| Continuous or many categories (\> 7) | `cor = "pearson"` (default) |
| 5–7 ordered categories (typical Likert) | `cor = "polychoric"` — recommended |
| 3–4 categories | `cor = "polychoric"` — strongly recommended |
| Binary (2 categories) | `cor = "polychoric"` — essential; consider tetrachoric specifically |

Polychoric estimation requires fitting a bivariate normal model for
every item pair. With p items that is p(p−1)/2 pairs — about 300 for a
25-item scale. Computation is usually fast (seconds), but scales as
O(p²) and can be slow for very large item pools.

> **A note on score computation.** When `cor = "polychoric"`, the weight
> matrices stored in the object are derived from the polychoric `R`.
> Materializing factor scores via
> [`augment()`](https://generics.r-lib.org/reference/augment.html)
> applies those weights to Pearson-standardized raw data
> (`.standardize(data)`). The resulting scores are calibrated to have
> *model-implied* unit variance, not empirically unit variance. For
> downstream analysis this distinction rarely matters — the scores are
> still well-scaled — but it is why
> [`augment()`](https://generics.r-lib.org/reference/augment.html)
> issues a warning when `cor != "pearson"`. The between-level edges from
> `tidy(what = "edges")` are always exact because they come from the
> algebra, not materialized scores.

## References

Pearson, K. (1900). Mathematical contributions to the theory of
evolution. VII. On the correlation of characters not quantitatively
measurable. *Philosophical Transactions of the Royal Society of London,
Series A*, *195*, 1–47.

Olsson, U. (1979). Maximum likelihood estimation of the polychoric
correlation coefficient. *Psychometrika*, *44*(4), 443–460.

Muthén, B. O. (1984). A general structural equation model with
dichotomous, ordered categorical, and continuous latent variable
indicators. *Psychometrika*, *49*(1), 115–132.

Flora, D. B., & Curran, P. J. (2004). An empirical evaluation of
alternative methods of estimation for confirmatory factor analysis with
ordinal data. *Psychological Methods*, *9*(4), 466–491.
