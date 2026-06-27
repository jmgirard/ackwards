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
#' @exportS3Method ggplot2::autoplot
autoplot.ackwards <- function(
    object,
    cut_show    = NULL,
    cut_strong  = 0.5,
    color_pos   = "#2166AC",
    color_neg   = "#D6604D",
    node_width  = 0.8,
    node_height = 0.4,
    min_sep     = 1.0,
    ...) {
  rlang::check_installed("ggplot2", reason = "for autoplot.ackwards()")

  cut_show <- cut_show %||% object$meta$cut_show %||% 0.3

  layout <- ba_layout(object, min_sep = min_sep)
  nodes  <- layout$nodes
  edges  <- layout$edges

  # Filter to displayable edges
  edges_show <- edges[abs(edges$r) >= cut_show, , drop = FALSE]
  if (nrow(edges_show) == 0L) {
    cli::cli_warn(
      "No edges with |r| {cli::symbol$geq} {cut_show} to display. \\
       Try lowering {.arg cut_show}."
    )
  }

  edges_show$linetype    <- ifelse(abs(edges_show$r) >= cut_strong, "solid", "dashed")
  edges_show$color_group <- ifelse(edges_show$r >= 0, "positive", "negative")

  # Attach node coordinates via named lookup (avoids column-name clashes in merge)
  nx <- stats::setNames(nodes$x, nodes$id)
  ny <- stats::setNames(nodes$y, nodes$id)

  # Arrow endpoints: bottom of parent box → top of child box
  edges_show$x_from <- nx[edges_show$from]
  edges_show$y_from <- ny[edges_show$from] - node_height / 2
  edges_show$x_to   <- nx[edges_show$to]
  edges_show$y_to   <- ny[edges_show$to]   + node_height / 2

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_show,
      ggplot2::aes(
        x        = .data$x_from,
        y        = .data$y_from,
        xend     = .data$x_to,
        yend     = .data$y_to,
        color    = .data$color_group,
        linewidth = abs(.data$r),
        linetype  = .data$linetype
      ),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.15, "cm"),
        type   = "closed"
      )
    ) +
    ggplot2::geom_tile(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y),
      width     = node_width,
      height    = node_height,
      fill      = "white",
      color     = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 3
    ) +
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
