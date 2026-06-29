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

x_pca <- ackwards(bfi, k_max = 3, cor = "polychoric")
x_pca
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: pca
#> Rotation: varimax
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

x_efa <- ackwards(bfi, k_max = 3, engine = "efa", cor = "polychoric")
x_efa
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: efa
#> Rotation: varimax
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

The table below pairs the `k = 3` loadings from each engine — built by
stacking `tidy(x_pca, what = "loadings")` and the EFA equivalent — for
six representative items (two each from the Neuroticism, Extraversion,
and Conscientiousness families).

| engine | factor | item | loading |
|:-------|:-------|:-----|--------:|
| pca    | m3f1   | C1   |    0.09 |
| efa    | m3f1   | C1   |    0.10 |
| pca    | m3f1   | C2   |    0.12 |
| efa    | m3f1   | C2   |    0.12 |
| pca    | m3f1   | E1   |   -0.60 |
| efa    | m3f1   | E1   |   -0.54 |
| pca    | m3f1   | E2   |   -0.65 |
| efa    | m3f1   | E2   |   -0.62 |
| pca    | m3f1   | N1   |   -0.07 |
| efa    | m3f1   | N1   |   -0.08 |
| pca    | m3f1   | N2   |   -0.08 |
| efa    | m3f1   | N2   |   -0.09 |
| pca    | m3f2   | C1   |   -0.05 |
| efa    | m3f2   | C1   |   -0.02 |
| pca    | m3f2   | C2   |    0.02 |
| efa    | m3f2   | C2   |    0.05 |
| pca    | m3f2   | E1   |    0.04 |
| efa    | m3f2   | E1   |    0.05 |
| pca    | m3f2   | E2   |    0.30 |
| efa    | m3f2   | E2   |    0.27 |
| pca    | m3f2   | N1   |    0.79 |
| efa    | m3f2   | N1   |    0.76 |
| pca    | m3f2   | N2   |    0.79 |
| efa    | m3f2   | N2   |    0.76 |
| pca    | m3f3   | C1   |    0.66 |
| efa    | m3f3   | C1   |    0.61 |
| pca    | m3f3   | C2   |    0.63 |
| efa    | m3f3   | C2   |    0.61 |
| pca    | m3f3   | E1   |    0.02 |
| efa    | m3f3   | E1   |    0.01 |
| pca    | m3f3   | E2   |   -0.10 |
| efa    | m3f3   | E2   |   -0.11 |
| pca    | m3f3   | N1   |   -0.08 |
| efa    | m3f3   | N1   |   -0.10 |
| pca    | m3f3   | N2   |    0.00 |
| efa    | m3f3   | N2   |   -0.03 |

PCA vs EFA loadings for anchor items (k = 3) {.table}

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
    When `cor = "polychoric"` is set with `engine = "esem"`, WLSMV is
    used automatically.

``` r

x_esem <- ackwards(bfi, k_max = 3, engine = "esem", cor = "polychoric")
x_esem
#> 
#> ── Bass-Ackwards Analysis (ackwards) ───────────────────────────────────────────
#> Engine: esem
#> Rotation: varimax
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

The unique output from ESEM is the rotation-aware **standard error** of
every loading. `tidy(what = "loadings_se")` returns them in the same
long format as `tidy(what = "loadings")`, with an `se` column in place
of `loading`:

``` r

se <- tidy(x_esem, what = "loadings_se")
head(se)
#>   level factor item         se
#> 1     1   m1f1   A1 0.01750773
#> 2     1   m1f1   A2 0.01561537
#> 3     1   m1f1   A3 0.01356797
#> 4     1   m1f1   A4 0.01675980
#> 5     1   m1f1   A5 0.01228831
#> 6     1   m1f1   C1 0.01675244
```

These SEs allow you to construct confidence intervals: loading ± 1.96 ×
SE gives an approximate 95% CI. For the BFI with \> 2,000 participants
the SEs are small; with smaller samples they become important for
judging which loadings are meaningfully non-zero.
(`tidy(what = "loadings_se")` errors for PCA and EFA objects, which
carry no loading standard errors.)

## How much do the edges differ?

The primary output of bass-ackwards analysis is the between-level edges.
For well-structured, continuous data, all three engines should agree
closely on the hierarchy.

The table stacks each engine’s primary-parent edges —
`tidy(what = "edges", primary_only = TRUE)` — at every level transition:

| engine | from | to   |     r |
|:-------|:-----|:-----|------:|
| pca    | m1f1 | m2f1 |  0.87 |
| efa    | m1f1 | m2f1 |  0.89 |
| pca    | m1f1 | m2f2 |  0.49 |
| efa    | m1f1 | m2f2 |  0.45 |
| pca    | m2f1 | m3f1 |  0.82 |
| efa    | m2f1 | m3f1 |  0.85 |
| pca    | m2f2 | m3f2 | -1.00 |
| efa    | m2f2 | m3f2 | -1.00 |
| pca    | m2f1 | m3f3 |  0.58 |
| efa    | m2f1 | m3f3 |  0.52 |

Primary-parent edges: PCA vs EFA {.table}

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

## Missing data

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
accepts a `missing` argument with three options.

**`"pairwise"` (default).** Use all available observations pairwise. The
exact behaviour depends on the engine and estimator:

- **PCA/EFA:** `stats::cor(use = "pairwise.complete.obs")` — uses all
  rows that contribute to each variable pair. MCAR-valid; N = total
  rows.
- **ESEM WLSMV/ULSMV (ordinal):** lavaan `missing = "available.cases"` —
  computes polychoric thresholds and correlations from every row
  contributing to each pair. MCAR-valid; uses the full N, not just
  complete cases. This is the honest interpretation of “pairwise” for
  ordinal data.
- **ESEM ML/MLR (continuous):** lavaan uses listwise deletion internally
  while edge correlations are computed from a separately-computed
  pairwise correlation matrix. This minor inconsistency (fit statistics
  at complete-case N, edges at full N) is documented in `$meta`. Prefer
  `"listwise"` or `"fiml"` when missingness is substantial.

A warning is emitted whenever incomplete rows are detected with this
option.

``` r

# Default: pairwise (warns if NAs present)
x <- ackwards(data_with_nas, k_max = 4)
```

**`"listwise"`.** Data are reduced to complete cases before *all*
downstream steps — correlation matrix, engine fitting, and edges — so
the three quantities are fully consistent. `x$n_obs` and
`x$meta$n_complete` both reflect the reduced sample size. Valid for all
three engines.

``` r

# Consistent complete-case analysis
x <- ackwards(data_with_nas, k_max = 4, missing = "listwise")
x$n_obs         # complete-case N
x$meta$missing  # "listwise"
```

**`"fiml"`.** Full Information Maximum Likelihood — available only for
`engine = "esem"` with `estimator = "ML"` or `"MLR"`. FIML uses
information from all rows (including those with partial data) when
estimating loadings and fit. Edge correlations are derived from lavaan’s
FIML-estimated saturated model (the h1 unrestricted model), ensuring
that fits and edges use the same information. Note: FIML improves
*estimation* but does not impute item responses; score materialisation
(`keep_scores = TRUE`) still yields `NA` for incomplete rows.

``` r

# FIML for ESEM with continuous data
x <- ackwards(data_with_nas, k_max = 4, engine = "esem",
              estimator = "ML", missing = "fiml")
```

`"fiml"` errors clearly for unsupported combinations:

``` r

# Errors: PCA and EFA are correlation-based, not raw-data likelihood
ackwards(data, k_max = 4, engine = "pca", missing = "fiml")

# Errors: WLSMV is limited-information WLS, no FIML extension
ackwards(data, k_max = 4, engine = "esem", cor = "polychoric",
         missing = "fiml")
```

### Which option to use?

| Situation | Recommendation |
|----|----|
| Continuous data, little missingness | `"pairwise"` (default) |
| Ordinal data + WLSMV, any missingness | `"pairwise"` (uses `available.cases` — MCAR-valid, full N) |
| Want consistent fit statistics and edges (continuous ML/MLR) | `"listwise"` |
| ESEM ML/MLR, meaningful missingness, want all rows used in estimation | `"fiml"` |
| MAR-valid with ordinal (not yet built-in) | MI via `lavaan.mi` or `mirt` |

## References

Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research
in Personality*, *40*(4), 347–358.

Waller, N. G. (2007). A general method for computing hierarchical
component structures by Bass-Ackward factor analysis. *Journal of
Research in Personality*, *41*(4), 745–752.

Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation
modeling. *Structural Equation Modeling*, *16*(3), 397–438.
