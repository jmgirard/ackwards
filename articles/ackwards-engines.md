# Choosing an Engine: PCA, EFA, and ESEM

**ackwards** supports three factor extraction engines. They share the
same downstream machinery — the same rotation, the same tenBerge scoring
weights, the same between-level correlation algebra — but differ in
their statistical model and what they report. This vignette explains
when each one is appropriate and what the differences look like in
practice.

## The three engines at a glance

|  | `"pca"` | `"efa"` | `"esem"` |
|----|----|----|----|
| **What it models** | Total item variance | Common (latent) variance | Common variance via lavaan |
| **Communalities** | All 1.0 (by definition) | Estimated from data | Estimated from data |
| **Fit indices** | Eigenvalues only | chi, RMSEA, TLI, BIC | CFI, TLI, RMSEA, SRMR + chi |
| **Loading SEs** | No | No | Yes (WLSMV for ordinal) |
| **Speed** | Fast | Moderate | Slowest |
| **Best for** | Exploration, large k | Latent-factor inference | Model evaluation, ordinal data |

All three produce the same labels (`m{k}f{j}`), the same
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
[`glance()`](https://generics.r-lib.org/reference/glance.html) /
[`augment()`](https://generics.r-lib.org/reference/augment.html)
interface, and comparable between-level edges for well-structured data.
The hierarchy they reveal is usually the same; the statistical
guarantees differ.

## Setup

``` r

library(ackwards)
bfi <- na.omit(psych::bfi[, 1:25])
```

We use the BFI-25 with polychoric correlations throughout so that
differences in output reflect the engine, not the correlation basis.

## PCA: components from total variance

PCA extracts **principal components** — linear combinations of the
observed variables that capture maximum variance, including measurement
error. Every item is modeled with communality 1.0: the components
account for 100% of each item’s variance. This is not a true latent
variable model; it is a data reduction method.

In the bass-ackwards context, PCA is the natural default. It is fast,
always converges, and produces eigenvalues that can guide the choice of
k. Waller (2007) showed that the between-level algebra (`W'RW`) holds
exactly for components, making the edges algebraically exact rather than
approximated from materialized scores.

``` r

x_pca <- ackwards(bfi, k = 3, cor = "polychoric")
x_pca
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: cfT
#> Basis: polychoric
#> n: 2,436
#> k (max): 3
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 22.9% variance
#> ✔ k = 2: 2 factors, 34.7% variance
#> ✔ k = 3: 3 factors, 43.9% variance
#> 
#> ── Edges ──
#> 
#> 5 of 8 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```

The “fit” for PCA is just the eigenvalue of each component — the amount
of variance it captures. There are no chi-square tests, no RMSEA, no
model rejection.

``` r

tidy(x_pca, what = "fit")
#>   level           index    value
#> 1     1 eigenvalue.m1f1 5.725295
#> 2     2 eigenvalue.m2f1 5.725295
#> 3     2 eigenvalue.m2f2 2.960038
#> 4     3 eigenvalue.m3f1 5.725295
#> 5     3 eigenvalue.m3f2 2.960038
#> 6     3 eigenvalue.m3f3 2.293698
```

## EFA: factors from common variance

EFA extracts **latent factors** that model only the variance shared
among items. Each item retains a unique variance (communality \< 1.0)
that the factors do not explain. This is the classical common-factor
model, and it is more appropriate than PCA when you believe the items
are fallible indicators of latent constructs rather than the constructs
themselves.

``` r

x_efa <- ackwards(bfi, k = 3, method = "efa", cor = "polychoric")
x_efa
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: efa
#> Rotation: cfT
#> Basis: polychoric
#> n: 2,436
#> k (max): 3
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 19.9% variance
#> ✔ k = 2: 2 factors, 30% variance
#> ✔ k = 3: 3 factors, 36.9% variance
#> 
#> ── Edges ──
#> 
#> 5 of 8 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```

EFA produces genuine goodness-of-fit indices. These tell you whether the
k factors are sufficient to reproduce the observed correlation matrix
within sampling error.

``` r

tidy(x_efa, what = "fit")
#>    level   index        value
#> 1      1     chi 1.338581e+04
#> 2      1     dof 2.750000e+02
#> 3      1 p_value 0.000000e+00
#> 4      1   RMSEA 1.437845e-01
#> 5      1     TLI 3.418359e-01
#> 6      1     BIC 1.198013e+04
#> 7      2     chi 7.111653e+03
#> 8      2     dof 2.510000e+02
#> 9      2 p_value 0.000000e+00
#> 10     2   RMSEA 1.213312e-01
#> 11     2     TLI 5.312110e-01
#> 12     2     BIC 7.294886e+03
#> 13     3     chi 4.096129e+03
#> 14     3     dof 2.280000e+02
#> 15     3 p_value 0.000000e+00
#> 16     3   RMSEA 1.074290e-01
#> 17     3     TLI 6.323810e-01
#> 18     3     BIC 4.860084e+03
```

The RMSEA values here are large (\> 0.10), indicating that 1–3 factors
do not fully account for the BFI item correlations — unsurprising,
because the true structure is 5 factors. Fit improves steadily from k =
1 to k = 3, which is exactly the kind of evidence bass-ackwards analysis
is designed to make visible.

### How close are EFA and PCA loadings?

For clean, continuous data with moderate-to-strong factor structure, EFA
and PCA loadings are highly correlated but not identical. EFA loadings
are systematically somewhat smaller because they model only the common
variance; PCA inflates loadings by fitting noise alongside signal.

``` r

l_pca <- tidy(x_pca, what = "loadings")
l_efa <- tidy(x_efa, what = "loadings")

# k = 3, first 5 items per factor
l_pca$engine <- "pca"
l_efa$engine <- "efa"
both <- rbind(l_pca, l_efa)
both3 <- both[both$level == 3, ]

# Representative items: highest-loading per factor in the PCA solution
anchors <- c("N1", "N2", "E1", "E2", "C1", "C2")
comp <- both3[both3$item %in% anchors, c("engine", "factor", "item", "loading")]
comp[order(comp$factor, comp$item), ]
#>     engine factor item      loading
#> 81     pca   m3f1   C1  0.089198696
#> 231    efa   m3f1   C1  0.100385156
#> 82     pca   m3f1   C2  0.121469270
#> 232    efa   m3f1   C2  0.116475023
#> 86     pca   m3f1   E1 -0.604309937
#> 236    efa   m3f1   E1 -0.539690270
#> 87     pca   m3f1   E2 -0.649293041
#> 237    efa   m3f1   E2 -0.617496322
#> 91     pca   m3f1   N1 -0.067953439
#> 241    efa   m3f1   N1 -0.075321454
#> 92     pca   m3f1   N2 -0.083392143
#> 242    efa   m3f1   N2 -0.086681968
#> 106    pca   m3f2   C1 -0.049859560
#> 256    efa   m3f2   C1 -0.017670030
#> 107    pca   m3f2   C2  0.024921135
#> 257    efa   m3f2   C2  0.049742317
#> 111    pca   m3f2   E1  0.042380567
#> 261    efa   m3f2   E1  0.052189787
#> 112    pca   m3f2   E2  0.298022534
#> 262    efa   m3f2   E2  0.273511715
#> 116    pca   m3f2   N1  0.793836214
#> 266    efa   m3f2   N1  0.763790131
#> 117    pca   m3f2   N2  0.793788172
#> 267    efa   m3f2   N2  0.763387518
#> 131    pca   m3f3   C1  0.657166570
#> 281    efa   m3f3   C1  0.614744266
#> 132    pca   m3f3   C2  0.633600363
#> 282    efa   m3f3   C2  0.608718517
#> 136    pca   m3f3   E1  0.023066417
#> 286    efa   m3f3   E1  0.007142178
#> 137    pca   m3f3   E2 -0.095262712
#> 287    efa   m3f3   E2 -0.108191392
#> 141    pca   m3f3   N1 -0.075697895
#> 291    efa   m3f3   N1 -0.095752244
#> 142    pca   m3f3   N2 -0.002610903
#> 292    efa   m3f3   N2 -0.029911910
```

EFA loadings for the same items are consistently a few points lower —
the PCA loadings include some noise variance that EFA partitions into
uniqueness. The factor structure (which items define which factor) is
unchanged.

## ESEM: EFA with full model diagnostics

ESEM (exploratory structural equation modeling, Asparouhov & Muthén,
2009) fits the same common-factor model as EFA but uses **lavaan** as
the engine. This unlocks two capabilities that EFA cannot provide:

1.  **Standard errors for every loading**, enabling confidence intervals
    and significance tests.
2.  **The WLSMV estimator** for ordinal data, which is the appropriate
    maximum-likelihood-adjacent estimator for categorical indicators.
    When `cor = "polychoric"` is set with `method = "esem"`, WLSMV is
    used automatically.

``` r

x_esem <- ackwards(bfi, k = 3, method = "esem", cor = "polychoric")
x_esem
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: esem
#> Rotation: cfT
#> Basis: polychoric
#> n: 2,436
#> k (max): 3
#> 
#> ── Levels ──
#> 
#> ✔ k = 1: 1 factor, 23.1% variance
#> ✔ k = 2: 2 factors, 32% variance
#> ✔ k = 3: 3 factors, 38.8% variance
#> 
#> ── Edges ──
#> 
#> 5 of 8 edges have |r| ≥ 0.3
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations.
```

ESEM fit indices include CFI and SRMR in addition to RMSEA and TLI,
giving a richer picture of model adequacy.

``` r

tidy(x_esem, what = "fit")
#>    level   index        value
#> 1      1     chi 2.018685e+04
#> 2      1     dof 2.750000e+02
#> 3      1 p_value           NA
#> 4      1     CFI 7.108666e-01
#> 5      1     TLI 6.845817e-01
#> 6      1   RMSEA 1.724408e-01
#> 7      1    SRMR 1.403309e-01
#> 8      2     chi 9.149066e+03
#> 9      2     dof 2.510000e+02
#> 10     2 p_value           NA
#> 11     2     CFI 8.707941e-01
#> 12     2     TLI 8.455707e-01
#> 13     2   RMSEA 1.206595e-01
#> 14     2    SRMR 9.666487e-02
#> 15     3     chi 4.900412e+03
#> 16     3     dof 2.280000e+02
#> 17     3 p_value           NA
#> 18     3     CFI 9.321535e-01
#> 19     3     TLI 9.107282e-01
#> 20     3   RMSEA 9.173893e-02
#> 21     3    SRMR 7.325403e-02
```

### Loading standard errors

The unique output from ESEM is `loadings_se` — a matrix of standard
errors for every loading at every level. Access it via the level object
directly:

``` r

# SE for k = 3 loadings
se3 <- x_esem$levels[["3"]]$loadings_se
round(se3[c("N1", "N2", "E1", "E2", "C1", "C2"), ], 3)
#>     m3f1  m3f2  m3f3
#> N1 0.015 0.008 0.016
#> N2 0.015 0.008 0.017
#> E1 0.016 0.019 0.019
#> E2 0.012 0.017 0.019
#> C1 0.017 0.018 0.015
#> C2 0.018 0.019 0.014
```

These SEs allow you to construct confidence intervals: loading ± 1.96 ×
SE gives an approximate 95% CI. For the BFI with \> 2,000 participants
the SEs are small; with smaller samples they become important for
judging which loadings are meaningfully non-zero.

## How much do the edges differ?

The primary output of bass-ackwards analysis is the between-level edges.
For well-structured, continuous data, all three engines should agree
closely on the hierarchy.

``` r

e_pca <- tidy(x_pca, what = "edges")
e_efa <- tidy(x_efa, what = "edges")
e_pca$engine <- "pca"
e_efa$engine <- "efa"

# Primary-parent edges at each level transition
primary_pca <- e_pca[e_pca$is_primary, c("engine", "from", "to", "r")]
primary_efa <- e_efa[e_efa$is_primary, c("engine", "from", "to", "r")]
rbind(primary_pca, primary_efa)[order(rbind(primary_pca, primary_efa)$to), ]
#>    engine from   to          r
#> 1     pca m1f1 m2f1  0.8737068
#> 11    efa m1f1 m2f1  0.8900109
#> 2     pca m1f1 m2f2  0.4864530
#> 21    efa m1f1 m2f2  0.4535937
#> 3     pca m2f1 m3f1  0.8157785
#> 31    efa m2f1 m3f1  0.8514122
#> 7     pca m2f2 m3f2 -0.9975486
#> 71    efa m2f2 m3f2 -0.9980294
#> 5     pca m2f1 m3f3  0.5767812
#> 51    efa m2f1 m3f3  0.5208236
```

The r values are very close between engines: the hierarchy that PCA
reveals is essentially the same hierarchy that EFA reveals. This
convergence across methods is reassuring — it suggests the structure is
real and not an artifact of the extraction method.

When the engines disagree on edges, that is itself informative: it
usually indicates factors whose definition depends on whether you
account for measurement error (EFA/ESEM) or not (PCA).

## Choosing an engine

| Situation | Recommendation |
|----|----|
| Exploratory, large k, unknown structure | Start with `"pca"` |
| Latent-variable theory, want to test model fit | `"efa"` |
| Ordinal items + model fit + loading SEs | `"esem"` with `cor = "polychoric"` |
| Replicating Goldberg (2006) or [`psych::bassAckward()`](https://rdrr.io/pkg/psych/man/bassAckward.html) | `"pca"`, `fm = "pca"` |
| Publication with formal model evaluation | `"esem"` |

A practical workflow: start with PCA to get a feel for the hierarchy and
choose k. Switch to EFA or ESEM to confirm and report. If PCA and
EFA/ESEM edges agree, you have robust evidence for the hierarchy; if
they disagree, investigate why.

## References

Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research
in Personality*, *40*(4), 347–358.

Waller, N. G. (2007). A general method for computing hierarchical
component structures by Bass-Ackward factor analysis. *Journal of
Research in Personality*, *41*(4), 745–752.

Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation
modeling. *Structural Equation Modeling*, *16*(3), 397–438.
