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
  show_r = FALSE,
  r_digits = 2L,
  r_label_size = 2.5,
  mono = FALSE,
  show_level_labels = TRUE,
  level_label_size = 3,
  node_labels = NULL,
  primary_only = FALSE,
  drop_pruned = FALSE,
  compress_levels = FALSE,
  show_arrows = TRUE,
  edge_linewidth = NULL,
  legend = TRUE,
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
  `cut_show` and `cut_strong` are dashed. Ignored when `mono = TRUE`.
  Default `0.5`.

- color_pos:

  Colour for positive edges. Default `"#2166AC"` (blue). Ignored when
  `mono = TRUE`.

- color_neg:

  Colour for negative edges. Default `"#D6604D"` (red). Ignored when
  `mono = TRUE`.

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
  `FALSE` otherwise. Ignored when `drop_pruned = TRUE`.

- curvature:

  Curvature of skip-level edge arcs. Passed to
  [`ggplot2::geom_curve()`](https://ggplot2.tidyverse.org/reference/geom_segment.html).
  Positive values curve right; default `0.2`. Ignored when
  `drop_pruned = TRUE`.

- color_pruned:

  Fill colour for nodes flagged as pruned/redundant. Default `"grey80"`.
  Only applied when the object carries pruning annotations (`x$prune` is
  non-`NULL`) and `drop_pruned = FALSE` (pruned nodes are omitted
  entirely when `drop_pruned = TRUE`).

- show_r:

  Whether to label each drawn edge with its correlation coefficient.
  Default `FALSE`. When `TRUE`, labels are formatted in APA style
  (leading zero stripped, e.g. `.23`, `-.30`) and placed beside each
  edge using a white-background label that clears the arrowhead.

- r_digits:

  Number of decimal places for edge labels when `show_r = TRUE`. Default
  `2L`.

- r_label_size:

  Font size for edge correlation labels when `show_r = TRUE`. Passed to
  [`ggplot2::geom_label()`](https://ggplot2.tidyverse.org/reference/geom_text.html)
  as `size`. Default `2.5`.

- mono:

  Monochrome mode. When `TRUE`, all edges are drawn in black;
  `linewidth` still encodes `|r|`; `linetype` encodes sign (`solid` =
  positive, `dashed` = negative). The `cut_strong` strong/weak linetype
  distinction is dropped in mono mode. Default `FALSE`.

- show_level_labels:

  Whether to draw level axis labels ("1 factor", "2 factors", ...) to
  the left of the diagram. Default `TRUE`.

- level_label_size:

  Font size for level axis labels. Default `3`.

- node_labels:

  A named character vector mapping factor IDs (e.g. `"m5f1"`) to custom
  display strings (e.g. `"General"`). Unspecified factors keep their
  default `m{k}f{j}` label. A warning is issued for names that match no
  factor ID in the object. Default `NULL`.

- primary_only:

  When `TRUE`, only primary-parent edges (`is_primary == TRUE`) are
  drawn. Because skip-level edges are never primary, this also
  suppresses skip arcs. Ignored when `drop_pruned = TRUE`. Default
  `FALSE`.

- drop_pruned:

  When `TRUE`, activates the Forbes (2023) pruned-view rendering path:
  pruned nodes are removed from the diagram entirely and each retained
  node is connected to its single strongest kept ancestor by a straight
  arrow (even across level gaps). Requires the object to carry pruning
  annotations (`prune != "none"` at fit time); errors if not. Overrides
  `show_skip`, `curvature`, and `primary_only`. Default `FALSE`.

- compress_levels:

  When `TRUE` under `drop_pruned = TRUE`, closes vertical gaps left by
  pruned levels so retained levels are evenly spaced; level axis
  labels (d) still show the original level numbers. Ignored when
  `drop_pruned = FALSE`. Default `FALSE`.

- show_arrows:

  When `FALSE`, edges are drawn with plain line ends instead of closed
  arrowheads (`arrow = NULL`). Applies to both straight and curved edge
  layers. Default `TRUE`. Forbes (2023) figures use plain line ends.

- edge_linewidth:

  `NULL` (default) maps `|r|` to `linewidth` via a continuous scale
  (current behaviour). A numeric value draws every edge at that constant
  width, removes the `linewidth` aesthetic mapping, and drops the `|r|`
  linewidth legend. Applies in both colour and `mono` modes and in the
  `drop_pruned` path. Forbes figures use uniform thin lines (≈ 0.5–0.6).

- legend:

  When `FALSE`, suppresses all plot legends
  (`legend.position = "none"`). Useful when `color_pos == color_neg`
  (e.g. both `"black"`) to remove an otherwise redundant Direction key.
  Default `TRUE`.

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

# Monochrome with correlation labels (for greyscale figures)
autoplot(x, mono = TRUE, show_r = TRUE)

# Custom node labels for the 5-factor level
autoplot(x, node_labels = c(m5f1 = "Neuroticism", m5f2 = "Agreeableness"))

# Primary links only — clean hierarchy tree
autoplot(x, primary_only = TRUE)

# Forbes pruned view: omit redundant nodes, straight spanning arrows
xp <- ackwards(psych::bfi[, 1:25], k = 5, prune = "redundant")
autoplot(xp, drop_pruned = TRUE)
autoplot(xp, drop_pruned = TRUE, show_r = TRUE) # with APA-style r labels
autoplot(xp, drop_pruned = TRUE, compress_levels = TRUE)

# Plain line ends without arrowheads
autoplot(x, show_arrows = FALSE)

# Uniform edge width (no |r| scaling)
autoplot(x, edge_linewidth = 0.5)

# Suppress legend
autoplot(x, legend = FALSE)

# Forbes (2023) publication style: black lines, uniform width, no arrowheads
autoplot(xp,
  drop_pruned = TRUE,
  color_pos = "black", color_neg = "black",
  edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE
)
} # }
```
