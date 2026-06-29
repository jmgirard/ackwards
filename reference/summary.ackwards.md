# Summarise an ackwards object

Returns a structured `summary_ackwards` object that, when printed, shows
per-level variance and fit indices, a readable lineage list, and (when
present) pruning annotations. More verbose than
[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md);
designed for inspection and reporting.

## Usage

``` r
# S3 method for class 'ackwards'
summary(object, ...)
```

## Arguments

- object:

  An `ackwards` object.

- ...:

  Ignored.

## Value

An object of class `"summary_ackwards"`, printed via
[`print.summary_ackwards()`](https://jmgirard.github.io/ackwards/reference/print.summary_ackwards.md).

## See also

[`print.ackwards()`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[`glance.ackwards()`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- ackwards(psych::bfi[, 1:25], k_max = 5)
summary(x)
} # }
```
