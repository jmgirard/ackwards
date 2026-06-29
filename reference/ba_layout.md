# Compute a layered layout for a bass-ackwards diagram

Returns tidy node coordinates and the edge table needed to render a
bass-ackwards diagram. The layout uses a two-pass barycenter algorithm:

## Usage

``` r
ba_layout(x, min_sep = 1)
```

## Arguments

- x:

  An `ackwards` object.

- min_sep:

  Minimum horizontal separation between adjacent nodes at the same
  level. Default `1.0`. Increase for wider diagrams.

## Value

A list with two data frames:

- nodes:

  One row per factor: `id`, `level`, `x`, `y`, `label`. `y = -level` so
  level 1 is at the top.

- edges:

  The tidy edge table from `x$edges$tidy`.

## Details

1.  **Top-down pass** – determines the left-to-right *order* of factors
    at each level. Each factor's ordinal rank is the \|r\|-weighted mean
    of its parents' ranks in the level above, so siblings that share a
    parent are grouped together.

2.  **Bottom-up pass** – assigns actual x coordinates. The deepest level
    (level k) is spread evenly; every upper-level factor is placed at
    the simple mean x of its **primary** children: a factor with one
    primary child lands directly above it; a factor with two primary
    children lands exactly halfway between them. Falls back to
    \|r\|-weighted mean of all children for factors with no primary
    children. Spreading resolves any remaining overlaps.

After both passes the layout is shifted so that the single level-1 node
is always at x = 0.

## See also

[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- ackwards(psych::bfi[, 1:25], k_max = 5)
lay <- ba_layout(x)
head(lay$nodes)
} # }
```
