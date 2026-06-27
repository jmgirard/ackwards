
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# ackwards

> Bass-ackwards hierarchical structural analysis in R

**ackwards** implements Goldberg’s (2006) bass-ackwards method and its
modern extensions for mapping the hierarchical structure of multivariate
data.

The core insight is simple: instead of asking *how many* factors best
describe your data, ask *how the solutions at different levels of
resolution are related to one another*. A 1-factor solution captures the
broadest shared variance; a 5-factor solution captures narrower, more
specific dimensions. Bass-ackwards analysis traces how broad factors
split into narrow ones — and which narrow factors are redundant
re-combinations of broader ones.

The package supports three extraction engines (PCA, EFA, ESEM),
polychoric correlations for ordinal data, and Forbes’s (2023) extension
that reveals skip-level connections and flags redundant or artefactual
factors.

## Installation

``` r
# install.packages("pak")
pak::pak("jmgirard/ackwards")
```

## Quick start

We use the 25-item Big Five Inventory (`psych::bfi`) as a running
example. Because BFI items are recorded on a 6-point ordinal scale, we
set `cor = "polychoric"`.

### Step 1 — Suggest a range of k

`suggest_k()` runs parallel analysis and Velicer’s MAP criterion to help
you choose an upper bound for the hierarchy depth.

``` r
library(ackwards)

bfi <- na.omit(psych::bfi[, 1:25])
sk  <- suppressWarnings(suggest_k(bfi))
#> Parallel analysis suggests that the number of factors =  NA  and the number of components =  5
sk
```

Both criteria agree: k = 5 fits the known Big Five structure of this
instrument.

### Step 2 — Fit the hierarchy

`ackwards()` fits factor models at every level from 1 to k and computes
the between-level factor-score correlations that define the hierarchy.

``` r
x <- ackwards(bfi, k = 5, cor = "polychoric")
x
```

### Step 3 — Visualize

`autoplot()` draws the hierarchical diagram. Each column is a level (k =
1 at left, k = 5 at right); arrows show which narrow factors inherit
from which broad factor. Solid arrows indicate strong connections (\|r\|
≥ 0.6 by default); dashed arrows show weaker ones.

``` r
autoplot(x)
```

<img src="man/figures/README-plot-1.png" alt="" width="100%" />

The five-factor level cleanly splits into the Big Five. The single broad
factor at k = 1 (roughly *general positive character*) differentiates
first into positive vs. negative affect (k = 2), then into successively
narrower traits.

### Step 4 — Inspect and score

`tidy()` extracts any part of the object into a tidy data frame;
`augment()` appends factor scores to your data for downstream analysis.

``` r
# Five strongest adjacent-level edges
edges <- tidy(x, what = "edges")
edges <- edges[order(-abs(edges$r)), ]
head(edges[edges$is_primary, c("from", "to", "r")], 8)
#>    from   to          r
#> 9  m3f1 m4f1  0.9990773
#> 7  m2f2 m3f2 -0.9975486
#> 33 m4f3 m5f3  0.9972853
#> 26 m4f2 m5f1  0.9964532
#> 40 m4f4 m5f5  0.9930407
#> 14 m3f2 m4f2  0.9789921
#> 1  m1f1 m2f1  0.8737068
#> 3  m2f1 m3f1  0.8157785
```

``` r
# Append scores for all 5 levels to the original data frame
scored <- augment(x, data = bfi)
#> Warning: ! Factor scores are standardized using model-implied SDs from a "polychoric"
#>   correlation matrix.
#> ℹ The raw projection uses `scale(data)` (Pearson z-scores), but `score_var`
#>   comes from the "polychoric" R.
#> ℹ Empirical score SDs will differ from 1.0. For non-Pearson analyses,
#>   between-level edges from `tidy()` are the authoritative associations.
dim(scored)   # original 25 items + 1+2+3+4+5 = 15 score columns
#> [1] 2436   40
names(scored)[26:40]
#>  [1] ".m1f1" ".m2f1" ".m2f2" ".m3f1" ".m3f2" ".m3f3" ".m4f1" ".m4f2" ".m4f3"
#> [10] ".m4f4" ".m5f1" ".m5f2" ".m5f3" ".m5f4" ".m5f5"
```

## Learn more

| Vignette | Topic |
|----|----|
| [Introduction](vignettes/ackwards-intro.Rmd) | Full PCA walkthrough: `suggest_k` → `ackwards` → inspect → plot → score |
| Engines & rotation *(coming soon)* | When to choose EFA or ESEM over PCA |
| Ordinal data *(coming soon)* | Polychoric correlations and WLSMV estimation |
| Forbes extension *(coming soon)* | Skip-level edges, redundancy pruning, `pairs = "all"` |

## Citation

If you use **ackwards** in your research, please cite both the method
and the package:

Goldberg, L. R. (2006). Doing it all bass-ackwards: The development of
hierarchical factor structures from the top down. *Journal of Research
in Personality*, *40*(4), 347–358.
<https://doi.org/10.1016/j.jrp.2006.01.001>

``` r
citation("ackwards")
```
