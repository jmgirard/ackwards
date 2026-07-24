# Compute a layered layout for a bass-ackwards diagram

Returns tidy node coordinates and the edge table needed to render a
bass-ackwards diagram. The layout is a layered (Sugiyama-style)
barycenter algorithm in two stages:

## Usage

``` r
ba_layout(x, min_sep = 1, order = NULL)
```

## Arguments

- x:

  An `ackwards` object.

- min_sep:

  Minimum horizontal separation between adjacent nodes at the same
  level. Default `1.0`. Increase for wider diagrams.

- order:

  Optional manual left-to-right ordering of the **deepest** (`k_max`)
  level, overriding Stage 1. Supply either a character vector of the
  deepest level's factor IDs in the desired order, or a named list with
  an entry for the deepest level (keyed by the level number as a string,
  e.g. `list("5" = c("m5f2", "m5f1", ...))`). Fixing the leaf order
  propagates upward: every upper factor stays at the mean-x of its
  primary children, so any arrangement of the primary forest is
  reachable from the leaf order. Entries for non-deepest levels are
  ignored with a warning (an upper factor's position is derived, not
  free). The vector must be a permutation of the deepest level's IDs.
  `NULL` (default) uses the automatic ordering.

## Value

A list with two data frames:

- nodes:

  One row per factor: `id`, `level`, `x`, `y`, `label`. `y = -level` so
  level 1 is at the top.

- edges:

  The tidy edge table from `x$edges$tidy`.

## Details

1.  **Ordering** – determines the left-to-right *order* of factors at
    each level. Two candidate orderings are scored and the
    crossing-minimising one kept: a single top-down \|r\|-weighted
    barycenter sweep (the historical order), and a **primary-forest
    traversal** – each factor has exactly one primary parent, so the
    primary edges form a forest, and a depth-first, subtree-contiguous
    leaf order lays every subtree out as an unbroken run. In deep
    hierarchies (`k >= 10`) the traversal drives primary-tree crossings
    to zero (the "bent levels" a single pass leaves behind); keep-best
    scoring (lexicographic: primary crossings, then all crossings) means
    shallow layouts are never made worse.

2.  **X-assignment** – assigns actual x coordinates bottom-up. The
    deepest level is spread evenly; every upper-level factor is placed
    at the simple mean x of its **primary** children: a factor with one
    primary child lands directly above it; a factor with two primary
    children lands exactly halfway between them. Falls back to
    \|r\|-weighted mean of all children for factors with no primary
    children. Spreading resolves any remaining overlaps.

After both stages the layout is shifted so that the single level-1 node
is always at x = 0. The result is fully deterministic.

## See also

[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md),
[`tidy.ackwards()`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)

## Examples

``` r
x <- ackwards(sim16, k_max = 5)
lay <- ba_layout(x)
head(lay$nodes)
#>     id level      x  y label
#> 1 m1f1     1  0.000 -1  m1f1
#> 2 m2f1     2 -1.125 -2  m2f1
#> 3 m2f2     2  1.125 -2  m2f2
#> 4 m3f1     3  1.125 -3  m3f1
#> 5 m3f2     3 -0.625 -3  m3f2
#> 6 m3f3     3 -1.625 -3  m3f3
```
