# Simulated continuous bass-ackwards teaching example (16 items, known hierarchy)

A 1 000-row, fully continuous dataset simulated from a population model
with a known 1 -\> 2 -\> 4 bass-ackwards hierarchy, for showcasing the
default `cor = "pearson"` extraction path without the ordinal-detection
warning that `bfi25`'s Likert items trigger.

## Usage

``` r
sim16
```

## Format

A data frame with 1 000 rows and 16 numeric columns (`i1`–`i16`,
continuous, no missing values).

## Source

Simulated; see `data-raw/sim16.R` for the full generative model.

## Details

**Population model.** `Sigma = Lambda %*% Phi %*% t(Lambda) + Psi`, an
oblique common-factor model with 4 true group factors, sampled via a
Cholesky factorization (base R; no `MASS` dependency):

|        |             |               |
|--------|-------------|---------------|
| Factor | Items       | Metatrait     |
| `f1`   | `i1`–`i4`   | 1 (with `f2`) |
| `f2`   | `i5`–`i8`   | 1 (with `f1`) |
| `f3`   | `i9`–`i12`  | 2 (with `f4`) |
| `f4`   | `i13`–`i16` | 2 (with `f3`) |

All items load `0.75` on their true factor (no cross-loadings). Factor
correlations: `0.45` within a metatrait (`f1`-`f2`, `f3`-`f4`), `0.15`
between metatraits. Uniquenesses are `1 - communality` (uniform `0.4375`
by the symmetry of the design above).

**Ground-truth hierarchy** (verified against
`ackwards(engine = "efa")`): `k=1` recovers a single general factor
across all 16 items; `k=2` splits along the metatrait line (`i1`-`i8`
vs. `i9`-`i16`); `k=4` recovers the 4 true group factors exactly; all
six
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
recommendations (five criteria – VSS reports at complexities 1 and 2)
reach a consensus of `k = 4`.

**Idealized by design.** The planted signal is strong and clean, so all
six
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
recommendations converge on `k = 4` – deliberately the *easy* case, for
building intuition about what recovering a known hierarchy looks like.
Real data rarely agree this cleanly: on `bfi25` the same criteria span
`k = 4`–`6`. The two datasets are complementary teaching foils – `sim16`
for "watch the method recover a structure we planted," `bfi25` for
"reason about a hierarchy when the criteria disagree." Present `sim16`'s
consensus as the ideal, not the norm.

**Deliberate overextraction artifact at k=5.** The population has
exactly 4 factors, so requesting a 5th finds no real dimension: EFA
produces an orphan factor with zero primary-loading items. With
`prune(x, "artifact")` (default `min_items = 3`, `orphan_r = 0.5`), that
factor is flagged both `few_items` and `orphan`. Because the true (non-
splitting) factors persist essentially unchanged from `k=3` onward,
their parent-child score correlations approach 1 and are flagged by
`prune(x, "redundant")` (`|r| >= .9` and, under the EFA auto-default,
Tucker's phi `> .95`). This is a textbook overextraction artifact,
included so the Forbes/redundancy examples have a guaranteed finding to
teach against (unlike `bfi25`, which does not reliably trigger one).

To regenerate this dataset, run `source("data-raw/sim16.R")` from the
package root (`set.seed(42)`).

## Examples

``` r
dim(sim16)
#> [1] 1000   16
head(sim16)
#>           i1         i2          i3         i4         i5          i6
#> 1  1.3709584  2.6935162  1.65649737  1.0043474  0.4985536  0.62179404
#> 2 -0.5646982  0.1157001 -0.37601948 -0.8080455 -0.9279296  0.23675384
#> 3  0.3631284  1.0068595 -0.83720071 -0.1617170 -0.2746798  0.12758795
#> 4  0.6328626  0.6676658 -1.07970658 -0.7959827  0.3055024  0.01789392
#> 5  0.4042683 -0.5960341 -1.06548809  0.4971529 -1.9685936 -1.31694755
#> 6 -0.1061245 -0.5536923  0.04465981 -0.8170172 -1.0854028 -0.72248740
#>           i7         i8         i9        i10        i11         i12        i13
#> 1  0.7158967  1.6662461  0.1762268 -0.5591902  0.8260793 -0.05320752  0.2350657
#> 2 -1.2592888 -0.1193249 -0.3315514 -0.6307922 -0.5144422 -0.06145734  0.5210552
#> 3 -0.7035972  0.5190214 -1.5114876 -3.4055614 -1.9790895 -2.95739348 -0.6546594
#> 4  0.5725013 -0.2194903  0.4623181 -0.3135998 -0.2535567  0.70600027 -0.3650701
#> 5 -1.2225681 -0.4665789 -1.3092214 -0.5390024 -0.1509965 -0.08517385 -0.3197751
#> 6 -0.8549671 -0.7564596  0.3796133  0.1051229 -1.1088126 -0.14489947 -0.5835958
#>           i14        i15        i16
#> 1 -1.25462666 -2.0439388 -1.4869628
#> 2  0.08671852  0.1475108 -0.4135835
#> 3  0.72176685 -1.5012947 -1.4849784
#> 4  0.67296942  1.4489548 -0.2946622
#> 5 -0.30006027 -0.8133962  0.8467458
#> 6  1.48728557  0.1890597  2.4643392
```
