# Glance at an ackwards object

Returns a one-row data frame of top-level model metadata. For EFA and
ESEM objects, fit indices at the deepest converged level are included.
The same five columns (`CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`) are present
across all engines; columns unavailable for a given engine or estimator
are `NA` (e.g., `CFI` and `SRMR` are `NA` for EFA; all five are `NA` for
PCA; for ESEM, `BIC` is `NA` under `estimator = "WLSMV"`/`"ULSMV"` –
these limited-information estimators have no proper log-likelihood – and
populated under `"ML"`/`"MLR"`). Under a scaled-test estimator
(`"WLSMV"`/`"ULSMV"`/`"MLR"`) the `CFI`/`TLI`/`RMSEA` reported here are
the scaled variants (see
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
for the rationale).

## Usage

``` r
# S3 method for class 'ackwards'
glance(x, ...)
```

## Arguments

- x:

  An `ackwards` object.

- ...:

  Ignored.

## Value

A one-row `data.frame`.

## See also

[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)

## Examples

``` r
x <- ackwards(sim16, k_max = 5)
glance(x)
#>   engine rotation     cor k_max n_obs deepest_converged n_edges CFI TLI RMSEA
#> 1    pca  varimax pearson     5  1000                 5      40  NA  NA    NA
#>   SRMR BIC
#> 1   NA  NA
```
