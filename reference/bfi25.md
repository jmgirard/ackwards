# Big Five Inventory – 25-item IPIP example dataset

A 1 000-row subset of the 25 IPIP Big Five personality-marker items from
Revelle's `psych` package
([`psych::bfi`](https://rdrr.io/pkg/psych/man/bfi.html)). Included so
that package examples and vignettes run without reaching into `psych`'s
namespace directly.

## Usage

``` r
bfi25
```

## Format

A data frame with 1 000 rows and 25 integer columns (Likert responses
scored 1–6, with some `NA`s reflecting genuine missing values in the
original survey). The 25 items span the Big Five personality domains:

|           |                   |
|-----------|-------------------|
| Columns   | Domain            |
| `A1`–`A5` | Agreeableness     |
| `C1`–`C5` | Conscientiousness |
| `E1`–`E5` | Extraversion      |
| `N1`–`N5` | Neuroticism       |
| `O1`–`O5` | Openness          |

## Source

Revelle, W. (2026). *psych: Procedures for Psychological, Psychometric,
and Personality Research*. Northwestern University, Evanston, Illinois.
R package. <https://CRAN.R-project.org/package=psych>

Data collected via the SAPA (Synthetic Aperture Personality Assessment)
project (<https://www.sapa-project.org/>). Items are public-domain IPIP
(International Personality Item Pool) markers developed by Goldberg
(1999):

Goldberg, L. R. (1999). A broad-bandwidth, public domain, personality
inventory measuring the lower-level facets of several five-factor
models. In I. Mervielde, I. Deary, F. De Fruyt, & F. Ostendorf (Eds.),
*Personality Psychology in Europe*, Vol. 7, pp. 7–28. Tilburg University
Press.

## Details

Derived from `psych::bfi[, 1:25]` by sampling 1 000 rows with
`set.seed(42)`. Missing values are preserved as in the original data.
The full dataset
([`psych::bfi`](https://rdrr.io/pkg/psych/man/bfi.html)) contains
responses from 2 800 participants collected via the SAPA project.

To regenerate this dataset, run `source("data-raw/bfi25.R")` from the
package root.

## Examples

``` r
dim(bfi25)
#> [1] 1000   25
head(bfi25)
#>   A1 A2 A3 A4 A5 C1 C2 C3 C4 C5 E1 E2 E3 E4 E5 N1 N2 N3 N4 N5 O1 O2 O3 O4 O5
#> 1  6  4  4  1  6  4  2  3  3  4  2  4  3  3  1  3  2  3  6  5  5  4  4  1  2
#> 2  1  5  5  4  4  4  5  5  2  4  1  2  4  5  4  4  5  5  3  4  4  4  5  6  2
#> 3  4  5  5  3  5  5  4  4  3  2  1  1  5  6  5  2  2  4  4  2  5  2  4  5  2
#> 4  1  5  3  6  6  4  3  1  5  5  4  2  2  5  2  3  3  2  5  6  5  2  5  6  1
#> 5  2  4  4  4  3  5  5  4  2  6  5  5  3  3  3  3  3  4  5  5  3  2  3  5  3
#> 6  3  6  4  6  5  5  5  5  2  4  5  2  4  5  5  1  1  2  2  5  4  4  5  4  4
```
