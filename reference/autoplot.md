# autoplot generic

Create ggplot2-based visualisations from model objects.

This generic is defined in ackwards so that
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
is available without requiring
[`library(ggplot2)`](https://ggplot2.tidyverse.org). When ggplot2 is
also loaded its own `autoplot` generic takes over in the search path,
but S3 dispatch still finds `autoplot.ackwards` correctly via either
route. Defining our own generic is intentional – it avoids putting
ggplot2 in Imports while keeping the `autoplot(x)` ergonomic without a
[`library()`](https://rdrr.io/r/base/library.html) call.

## Usage

``` r
autoplot(object, ...)
```

## Arguments

- object:

  An object to visualise.

- ...:

  Additional arguments passed to methods.

## Value

A `ggplot` object.
