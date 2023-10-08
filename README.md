
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ackwards

<!-- badges: start -->
<!-- badges: end -->

The goal of {ackwards} is to provide functions for the traditional and
extended bass-ackwards techniques.

## Installation

You can install the development version of {ackwards} from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jmgirard/ackwards")
```

## Example

Work in progress on the extended bass-ackwards technique:

``` r
library(ackwards)
x <- wba(r = cor(skiers), nfactors = 3)
print(x, digits = 3, cut = 0.3)
#> Bass-Ackwards Analysis
#> 
#>   Method = Waller
#>   Engine = eigen
#> 
#> Cross-Level Correlations
#> 
#>          A1
#> B1    0.904
#> B2   -0.427
#> 
#>          B1      B2
#> C1    1.000   0.000
#> C2    0.000   1.000
#> C3   -0.021   0.003
#> 
#> 
#> Factor Loadings
#> 
#>              A1
#> cost     -0.500
#> lift      0.357
#> depth     0.891
#> powder    0.919
#> 
#>             B1       B2
#> cost         .    0.988
#> lift         .   -0.989
#> depth    0.997        .
#> powder   0.998        .
#> 
#>             C1       C2   C3
#> cost         .    0.987    .
#> lift         .   -0.989    .
#> depth    0.998        .    .
#> powder   0.997        .    .
#> 
#> . = Loading magnitude less than 0.3
#> Bass-Ackwards Analysis
#> 
#>   Method = Waller
#>   Engine = eigen
#> 
#> Cross-Level Correlations
#> 
#>          A1
#> B1    0.904
#> B2   -0.427
#> 
#>          B1      B2
#> C1    1.000   0.000
#> C2    0.000   1.000
#> C3   -0.021   0.003
#> 
#> 
#> Factor Loadings
#> 
#>              A1
#> cost     -0.500
#> lift      0.357
#> depth     0.891
#> powder    0.919
#> 
#>             B1       B2
#> cost         .    0.988
#> lift         .   -0.989
#> depth    0.997        .
#> powder   0.998        .
#> 
#>             C1       C2   C3
#> cost         .    0.987    .
#> lift         .   -0.989    .
#> depth    0.998        .    .
#> powder   0.997        .    .
#> 
#> . = Loading magnitude less than 0.3
```

## References

- Goldberg, L. R. (2006). Doing it all bass-ackwards: The development of
  hierarchical factor structures from the top down. *Journal of Research
  in Personality, 40*(4), 347–358. <https://doi.org/10/c9sqkd>

- Waller, N. (2007). A general method for computing hierarchical
  component structures by Goldberg’s bass-ackwards method. *Journal of
  Research in Personality, 41*(4), 745–752. <https://doi.org/10/bcz8wd>

- Forbes, M. K. (2022). Improving hierarchical models of individual
  differences: An extension of Goldberg’s bass-ackwards method.
  *Psychological Methods.* Advanced online publication.
  <https://doi.org/10.1037/met0000546>
