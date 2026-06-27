# Plot a bass-ackwards diagram

Renders the layered bass-ackwards hierarchy as a ggplot2 diagram.
Factors appear as labelled boxes arranged in levels (level 1 at top,
level k at bottom). Between-level edges are drawn as arrows coloured by
sign and scaled by \|r\|. Edges below `cut_show` are hidden; edges above
`cut_strong` are solid and those between the two thresholds are dashed.

## Usage

``` r
# S3 method for class 'ackwards'
autoplot(
  object,
  cut_show = NULL,
  cut_strong = 0.5,
  color_pos = "#2166AC",
  color_neg = "#D6604D",
  node_width = 0.8,
  node_height = 0.4,
  min_sep = 1,
  show_skip = NULL,
  curvature = 0.2,
  color_pruned = "grey80",
  ...
)

# S3 method for class 'ackwards'
plot(x, ...)
```

## Arguments

- object:

  An `ackwards` object.

- cut_show:

  Edges with `|r| < cut_show` are not drawn. Defaults to the value used
  when the object was created (`x$meta$cut_show`, typically 0.3).

- cut_strong:

  Edges with `|r| >= cut_strong` are drawn solid; those between
  `cut_show` and `cut_strong` are dashed. A documented threshold
  separating "notable" from "strong" associations. Default `0.5`.

- color_pos:

  Colour for positive edges. Default `"#2166AC"` (blue).

- color_neg:

  Colour for negative edges. Default `"#D6604D"` (red).

- node_width:

  Width of factor boxes in layout units. Default `0.8`.

- node_height:

  Height of factor boxes in layout units. Default `0.4`.

- min_sep:

  Minimum horizontal separation between nodes; passed to
  [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md).
  Default `1.0`.

- show_skip:

  Whether to draw skip-level (non-adjacent) edges. `NULL` (default)
  auto-detects: `TRUE` when the object was run with `pairs = "all"`,
  `FALSE` otherwise. Skip-level edges are drawn as curved lines to
  distinguish them from the adjacent straight arrows.

- curvature:

  Curvature of skip-level edge arcs. Passed to
  [`ggplot2::geom_curve()`](https://ggplot2.tidyverse.org/reference/geom_segment.html).
  Positive values curve right; default `0.2`.

- color_pruned:

  Fill colour for nodes flagged as pruned/redundant. Default `"grey80"`.
  Only applied when the object carries pruning annotations (`x$prune` is
  non-`NULL`).

- ...:

  Ignored.

- x:

  An `ackwards` object.

## Value

A `ggplot` object.

## Details

Requires the ggplot2 package.

## See also

[`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md),
`plot.ackwards()`

## Examples

``` r
if (FALSE) { # \dontrun{
x <- ackwards(psych::bfi[, 1:25], k = 5)
autoplot(x)
autoplot(x, cut_strong = 0.6, color_pos = "steelblue")
} # }
```
