#' @title autoplot generic
#'
#' @description
#' Create ggplot2-based visualisations from model objects.
#'
#' This generic is defined in \pkg{ackwards} so that [autoplot.ackwards()] is
#' available without requiring `library(ggplot2)`. When \pkg{ggplot2} is also
#' loaded its own `autoplot` generic takes over in the search path, but S3
#' dispatch still finds `autoplot.ackwards` correctly via either route.
#' Defining our own generic is intentional -- it avoids putting ggplot2 in
#' Imports while keeping the `autoplot(x)` ergonomic without a `library()` call.
#'
#' @param object An object to visualise.
#' @param ... Additional arguments passed to methods.
#' @return A `ggplot` object.
#' @export
autoplot <- function(object, ...) UseMethod("autoplot")

#' Plot a bass-ackwards diagram
#'
#' Renders the layered bass-ackwards hierarchy as a ggplot2 diagram. Factors
#' appear as labelled boxes arranged in levels (level 1 at top, level k at
#' bottom). Between-level edges are drawn as arrows coloured by sign and scaled
#' by |r|. Edges below `cut_show` are hidden; edges above `cut_strong` are
#' solid and those between the two thresholds are dashed.
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @param object An `ackwards` object.
#' @param cut_show Edges with `|r| < cut_show` are not drawn. Defaults to the
#'   value used when the object was created (`x$meta$cut_show`, typically 0.3).
#' @param cut_strong Edges with `|r| >= cut_strong` are drawn solid; those
#'   between `cut_show` and `cut_strong` are dashed. A documented threshold
#'   separating "notable" from "strong" associations. Default `0.5`.
#' @param color_pos Colour for positive edges. Default `"#2166AC"` (blue).
#' @param color_neg Colour for negative edges. Default `"#D6604D"` (red).
#' @param node_width Width of factor boxes in layout units. Default `0.8`.
#' @param node_height Height of factor boxes in layout units. Default `0.4`.
#' @param min_sep Minimum horizontal separation between nodes; passed to
#'   [ba_layout()]. Default `1.0`.
#' @param show_skip Whether to draw skip-level (non-adjacent) edges. `NULL`
#'   (default) auto-detects: `TRUE` when the object was run with
#'   `pairs = "all"`, `FALSE` otherwise. Skip-level edges are drawn as curved
#'   lines to distinguish them from the adjacent straight arrows.
#' @param curvature Curvature of skip-level edge arcs. Passed to
#'   [ggplot2::geom_curve()]. Positive values curve right; default `0.2`.
#' @param color_pruned Fill colour for nodes flagged as pruned/redundant.
#'   Default `"grey80"`. Only applied when the object carries pruning
#'   annotations (`x$prune` is non-`NULL`).
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @seealso [ba_layout()], [plot.ackwards()]
#'
#' @examples
#' \dontrun{
#' x <- ackwards(psych::bfi[, 1:25], k = 5)
#' autoplot(x)
#' autoplot(x, cut_strong = 0.6, color_pos = "steelblue")
#' }
#'
#' @importFrom rlang .data `%||%`
#' @export
autoplot.ackwards <- function(
  object,
  cut_show = NULL,
  cut_strong = 0.5,
  color_pos = "#2166AC",
  color_neg = "#D6604D",
  node_width = 0.8,
  node_height = 0.4,
  min_sep = 1.0,
  show_skip = NULL,
  curvature = 0.2,
  color_pruned = "grey80",
  ...
) {
  rlang::check_installed("ggplot2", reason = "for autoplot.ackwards()")

  if (min_sep < node_width) {
    cli::cli_warn(
      "{.arg min_sep} ({min_sep}) is less than {.arg node_width} ({node_width}); \\
       adjacent boxes may overlap."
    )
  }

  cut_show <- cut_show %||% object$meta$cut_show %||% 0.3
  show_skip <- show_skip %||% isTRUE(object$meta$pairs == "all")

  layout <- ba_layout(object, min_sep = min_sep)
  nodes <- layout$nodes
  edges <- layout$edges

  # Annotate pruned nodes with distinct fill
  nodes$fill <- "white"
  if (!is.null(object$prune) && !is.null(object$prune$nodes)) {
    pruned_ids <- object$prune$nodes$id[object$prune$nodes$pruned]
    nodes$fill[nodes$id %in% pruned_ids] <- color_pruned
  }

  # Partition edges: adjacent (straight arrows) vs skip-level (curved arcs)
  is_adjacent <- (edges$level_to - edges$level_from) == 1L
  edges_adj <- edges[is_adjacent & abs(edges$r) >= cut_show, , drop = FALSE]
  edges_skip <- edges[!is_adjacent & abs(edges$r) >= cut_show, , drop = FALSE]

  if (nrow(edges_adj) == 0L && (!show_skip || nrow(edges_skip) == 0L)) {
    cli::cli_warn(
      "No edges with |r| {cli::symbol$geq} {cut_show} to display. \\
       Try lowering {.arg cut_show}."
    )
  }

  # Attach linetype and direction color to each edge subset
  .annotate_edges <- function(e) {
    e$linetype <- ifelse(abs(e$r) >= cut_strong, "solid", "dashed")
    e$color_group <- ifelse(e$r >= 0, "positive", "negative")
    e
  }
  edges_adj <- .annotate_edges(edges_adj)
  edges_skip <- .annotate_edges(edges_skip)

  # Attach node coordinates (named lookup avoids merge column conflicts)
  nx <- stats::setNames(nodes$x, nodes$id)
  ny <- stats::setNames(nodes$y, nodes$id)

  .attach_coords <- function(e) {
    e$x_from <- nx[e$from]
    e$y_from <- ny[e$from] - node_height / 2 # bottom of parent box
    e$x_to <- nx[e$to]
    e$y_to <- ny[e$to] + node_height / 2 # top of child box
    e
  }
  edges_adj <- .attach_coords(edges_adj)
  edges_skip <- .attach_coords(edges_skip)

  arrow_spec <- ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")

  edge_aes <- ggplot2::aes(
    x         = .data$x_from,
    y         = .data$y_from,
    xend      = .data$x_to,
    yend      = .data$y_to,
    color     = .data$color_group,
    linewidth = abs(.data$r),
    linetype  = .data$linetype
  )

  p <- ggplot2::ggplot()

  # Adjacent edges — straight arrows
  if (nrow(edges_adj) > 0L) {
    p <- p + ggplot2::geom_segment(
      data  = edges_adj,
      edge_aes,
      arrow = arrow_spec
    )
  }

  # Skip-level edges — curved arcs (drawn before nodes so nodes sit on top)
  if (show_skip && nrow(edges_skip) > 0L) {
    p <- p + ggplot2::geom_curve(
      data      = edges_skip,
      edge_aes,
      curvature = curvature,
      arrow     = arrow_spec
    )
  }

  p +
    ggplot2::geom_tile(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$fill),
      width = node_width,
      height = node_height,
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 3
    ) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::scale_color_manual(
      values = c(positive = color_pos, negative = color_neg),
      name   = "Direction",
      labels = c(positive = "Positive", negative = "Negative")
    ) +
    ggplot2::scale_linewidth_continuous(
      range  = c(0.4, 1.8),
      name   = "|r|",
      limits = c(cut_show, 1)
    ) +
    ggplot2::scale_linetype_identity(guide = "none") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin     = ggplot2::margin(10, 10, 10, 10)
    )
}

#' @rdname autoplot.ackwards
#'
#' @param x An `ackwards` object.
#'
#' @export
plot.ackwards <- function(x, ...) {
  rlang::check_installed("ggplot2", reason = "for plot.ackwards()")
  print(ggplot2::autoplot(x, ...))
  invisible(x)
}
