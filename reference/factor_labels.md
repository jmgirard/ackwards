# Read the factor labels stored on an ackwards object

Read the factor labels stored on an ackwards object

## Usage

``` r
factor_labels(x)
```

## Arguments

- x:

  An `ackwards` object.

## Value

The named character vector of factor labels (names are factor IDs), or
`NULL` if none have been set. See
[`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md)
to attach them.

## See also

[`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md)

## Examples

``` r
x <- ackwards(sim16, k_max = 4, engine = "pca")
factor_labels(x) # NULL -- none set yet
#> NULL
x <- set_factor_labels(x, c(m4f1 = "Alpha"))
factor_labels(x)
#>    m4f1 
#> "Alpha" 
```
