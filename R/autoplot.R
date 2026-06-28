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
#'   between `cut_show` and `cut_strong` are dashed. Ignored when
#'   `mono = TRUE`. Default `0.5`.
#' @param color_pos Colour for positive edges. Default `"#2166AC"` (blue).
#'   Ignored when `mono = TRUE`.
#' @param color_neg Colour for negative edges. Default `"#D6604D"` (red).
#'   Ignored when `mono = TRUE`.
#' @param node_width Width of factor boxes in layout units. Default `0.8`.
#' @param node_height Height of factor boxes in layout units. Default `0.4`.
#' @param min_sep Minimum horizontal separation between nodes; passed to
#'   [ba_layout()]. Default `1.0`.
#' @param show_skip Whether to draw skip-level (non-adjacent) edges. `NULL`
#'   (default) auto-detects: `TRUE` when the object was run with
#'   `pairs = "all"`, `FALSE` otherwise. Ignored when `drop_pruned = TRUE`.
#' @param curvature Curvature of skip-level edge arcs. Passed to
#'   [ggplot2::geom_curve()]. Positive values curve right; default `0.2`.
#'   Ignored when `drop_pruned = TRUE`.
#' @param color_pruned Fill colour for nodes flagged as pruned/redundant.
#'   Default `"grey80"`. Only applied when the object carries pruning
#'   annotations (`x$prune` is non-`NULL`) and `drop_pruned = FALSE` (pruned
#'   nodes are omitted entirely when `drop_pruned = TRUE`).
#' @param show_r Whether to label each drawn edge with its rounded correlation
#'   coefficient. `NULL` (default) uses `FALSE` normally and `TRUE` when
#'   `drop_pruned = TRUE` (Forbes 2023 figures always label spanning arrows).
#'   Set explicitly to `TRUE` or `FALSE` to override.
#' @param r_digits Number of decimal places for edge labels when
#'   `show_r = TRUE`. Default `2L`.
#' @param mono Monochrome mode. When `TRUE`, all edges are drawn in black;
#'   `linewidth` still encodes `|r|`; `linetype` encodes sign (`solid` =
#'   positive, `dashed` = negative). The `cut_strong` strong/weak linetype
#'   distinction is dropped in mono mode. Default `FALSE`.
#' @param show_level_labels Whether to draw level axis labels
#'   ("1 factor", "2 factors", ...) to the left of the diagram. Default `TRUE`.
#' @param level_label_size Font size for level axis labels. Default `3`.
#' @param node_labels A named character vector mapping factor IDs (e.g.
#'   `"m5f1"`) to custom display strings (e.g. `"General"`). Unspecified
#'   factors keep their default `m{k}f{j}` label. A warning is issued for
#'   names that match no factor ID in the object. Default `NULL`.
#' @param primary_only When `TRUE`, only primary-parent edges (`is_primary ==
#'   TRUE`) are drawn. Because skip-level edges are never primary, this also
#'   suppresses skip arcs. Ignored when `drop_pruned = TRUE`. Default `FALSE`.
#' @param drop_pruned When `TRUE`, activates the Forbes (2023) pruned-view
#'   rendering path: pruned nodes are removed from the diagram entirely and
#'   each retained node is connected to its single strongest kept ancestor by
#'   a straight arrow (even across level gaps). Requires the object to carry
#'   pruning annotations (`prune != "none"` at fit time); errors if not.
#'   Overrides `show_skip`, `curvature`, and `primary_only`. Default `FALSE`.
#' @param compress_levels When `TRUE` under `drop_pruned = TRUE`, closes
#'   vertical gaps left by pruned levels so retained levels are evenly spaced;
#'   level axis labels (d) still show the original level numbers. Ignored when
#'   `drop_pruned = FALSE`. Default `FALSE`.
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
#'
#' # Monochrome with correlation labels (for greyscale figures)
#' autoplot(x, mono = TRUE, show_r = TRUE)
#'
#' # Custom node labels for the 5-factor level
#' autoplot(x, node_labels = c(m5f1 = "Neuroticism", m5f2 = "Agreeableness"))
#'
#' # Primary links only — clean hierarchy tree
#' autoplot(x, primary_only = TRUE)
#'
#' # Forbes pruned view: omit redundant nodes, straight spanning arrows
#' xp <- ackwards(psych::bfi[, 1:25], k = 5, prune = "redundant")
#' autoplot(xp, drop_pruned = TRUE)
#' autoplot(xp, drop_pruned = TRUE, compress_levels = TRUE)
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
  show_r = NULL,
  r_digits = 2L,
  mono = FALSE,
  show_level_labels = TRUE,
  level_label_size = 3,
  node_labels = NULL,
  primary_only = FALSE,
  drop_pruned = FALSE,
  compress_levels = FALSE,
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
  show_r <- show_r %||% drop_pruned # Forbes: label spanning arrows by default

  layout <- ba_layout(object, min_sep = min_sep)
  nodes <- layout$nodes

  # (e) Custom node labels applied before any geometry
  if (!is.null(node_labels)) {
    unknown <- setdiff(names(node_labels), nodes$id)
    if (length(unknown) > 0L) {
      cli::cli_warn(
        "Some {.arg node_labels} name{?s} match no factor ID: {.val {unknown}}"
      )
    }
    matched <- intersect(names(node_labels), nodes$id)
    nodes$label[match(matched, nodes$id)] <- node_labels[matched]
  }

  # Annotate pruned nodes with distinct fill (only visible in normal path)
  nodes$fill <- "white"
  if (!is.null(object$prune) && !is.null(object$prune$nodes)) {
    pruned_ids <- object$prune$nodes$id[object$prune$nodes$pruned]
    nodes$fill[nodes$id %in% pruned_ids] <- color_pruned
  }

  # Local helpers shared between rendering paths
  .ann <- function(e) {
    if (mono) {
      e$linetype <- ifelse(e$r >= 0, "solid", "dashed")
    } else {
      e$linetype <- ifelse(abs(e$r) >= cut_strong, "solid", "dashed")
      e$color_group <- ifelse(e$r >= 0, "positive", "negative")
    }
    e
  }
  .attach_coords <- function(e, nx, ny) {
    e$x_from <- nx[e$from]
    e$y_from <- ny[e$from] - node_height / 2 # bottom of parent box
    e$x_to <- nx[e$to]
    e$y_to <- ny[e$to] + node_height / 2 # top of child box
    e
  }

  # --------------------------------------------------------------------------
  # Build draw_edges — branching on drop_pruned
  # --------------------------------------------------------------------------

  if (drop_pruned) {
    if (is.null(object$prune)) {
      cli::cli_abort(c(
        "!" = "{.arg drop_pruned = TRUE} requires pruning annotations.",
        "i" = "Refit with {.code prune = \"redundant\"} to flag factors."
      ))
    }
    if (!any(object$prune$nodes$pruned)) {
      cli::cli_warn(c(
        "!" = "{.arg drop_pruned = TRUE}: no nodes are flagged for pruning.",
        "i" = "The pruned view will be identical to the full diagram.",
        "i" = "Use {.code prune = \"redundant\"} to flag redundant factors."
      ))
    }

    dp <- .drop_pruned_nodes(object, nodes, compress_levels = compress_levels)
    nodes <- dp$nodes # kept-only subset; y re-indexed when compress_levels
    edges <- dp$edges # reduced primary edges, no coords yet

    n_kept_levels <- length(unique(nodes$level))

    if (n_kept_levels <= 1L) {
      cli::cli_warn(
        "{n_kept_levels} level{?s} remain{?s} after pruning; \\
         no edges possible. Returning a node-only plot."
      )
      return(.ba_degenerate_plot(
        nodes, node_width, node_height,
        show_level_labels, level_label_size
      ))
    }

    edges <- edges[abs(edges$r) >= cut_show, , drop = FALSE]
    if (nrow(edges) == 0L) {
      cli::cli_warn(
        "No edges with |r| {cli::symbol$geq} {cut_show} to display. \\
         Try lowering {.arg cut_show}."
      )
    }

    nx <- stats::setNames(nodes$x, nodes$id)
    ny <- stats::setNames(nodes$y, nodes$id)
    edges <- .ann(.attach_coords(edges, nx, ny))
    edges$curved <- rep(FALSE, nrow(edges))
    draw_edges <- edges
  } else {
    # --- Normal rendering path ---
    edges <- layout$edges

    # (f) Primary-edges-only filter
    if (primary_only) {
      edges <- edges[!is.na(edges$is_primary) & edges$is_primary, , drop = FALSE]
    }

    # Partition into adjacent vs skip-level
    is_adj <- (edges$level_to - edges$level_from) == 1L
    edges_adj <- edges[is_adj & abs(edges$r) >= cut_show, , drop = FALSE]
    edges_skip <- edges[!is_adj & abs(edges$r) >= cut_show, , drop = FALSE]

    if (nrow(edges_adj) == 0L && (!show_skip || nrow(edges_skip) == 0L)) {
      cli::cli_warn(
        "No edges with |r| {cli::symbol$geq} {cut_show} to display. \\
         Try lowering {.arg cut_show}."
      )
    }

    nx <- stats::setNames(nodes$x, nodes$id)
    ny <- stats::setNames(nodes$y, nodes$id)

    edges_adj <- .ann(.attach_coords(edges_adj, nx, ny))
    edges_skip <- .ann(.attach_coords(edges_skip, nx, ny))
    edges_adj$curved <- rep(FALSE, nrow(edges_adj))
    edges_skip$curved <- rep(TRUE, nrow(edges_skip))

    parts <- list()
    if (nrow(edges_adj) > 0L) parts[[1L]] <- edges_adj
    if (show_skip && nrow(edges_skip) > 0L) parts[[2L]] <- edges_skip
    draw_edges <- if (length(parts) > 0L) {
      do.call(rbind, parts)
    } else {
      edges_adj[0L, , drop = FALSE]
    }
  }

  # --------------------------------------------------------------------------
  # Draw
  # --------------------------------------------------------------------------

  arrow_spec <- ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")

  if (mono) {
    edge_aes <- ggplot2::aes(
      x         = .data$x_from,
      y         = .data$y_from,
      xend      = .data$x_to,
      yend      = .data$y_to,
      linewidth = abs(.data$r),
      linetype  = .data$linetype
    )
  } else {
    edge_aes <- ggplot2::aes(
      x         = .data$x_from,
      y         = .data$y_from,
      xend      = .data$x_to,
      yend      = .data$y_to,
      color     = .data$color_group,
      linewidth = abs(.data$r),
      linetype  = .data$linetype
    )
  }

  p <- ggplot2::ggplot()

  # Straight edges (adjacent in normal mode; all edges in drop_pruned mode)
  ed_str <- draw_edges[!draw_edges$curved, , drop = FALSE]
  if (nrow(ed_str) > 0L) {
    if (mono) {
      p <- p + ggplot2::geom_segment(
        data = ed_str, edge_aes, arrow = arrow_spec, color = "black"
      )
    } else {
      p <- p + ggplot2::geom_segment(data = ed_str, edge_aes, arrow = arrow_spec)
    }
  }

  # Curved arcs (skip-level edges in normal path only; never in drop_pruned)
  ed_cur <- draw_edges[draw_edges$curved, , drop = FALSE]
  if (nrow(ed_cur) > 0L) {
    if (mono) {
      p <- p + ggplot2::geom_curve(
        data = ed_cur, edge_aes,
        curvature = curvature, arrow = arrow_spec, color = "black"
      )
    } else {
      p <- p + ggplot2::geom_curve(
        data = ed_cur, edge_aes, curvature = curvature, arrow = arrow_spec
      )
    }
  }

  # (a) Edge correlation labels at midpoints
  if (show_r && nrow(draw_edges) > 0L) {
    draw_edges$lx <- (draw_edges$x_from + draw_edges$x_to) / 2
    draw_edges$ly <- (draw_edges$y_from + draw_edges$y_to) / 2
    draw_edges$rl <- round(draw_edges$r, r_digits)
    p <- p + ggplot2::geom_text(
      data = draw_edges,
      ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$rl),
      size = 2.5,
      nudge_x = 0.1,
      inherit.aes = FALSE
    )
  }

  # Node tiles and labels (drawn on top of all edges)
  p <- p +
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
    ggplot2::scale_fill_identity(guide = "none")

  # Edge scales — conditional on mono
  if (mono) {
    p <- p +
      ggplot2::scale_linewidth_continuous(
        range  = c(0.4, 1.8),
        name   = "|r|",
        limits = c(cut_show, 1)
      ) +
      ggplot2::scale_linetype_manual(
        values = c("solid" = "solid", "dashed" = "dashed"),
        breaks = c("solid", "dashed"),
        labels = c("solid" = "Positive", "dashed" = "Negative"),
        name   = "Direction"
      )
  } else {
    p <- p +
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
      ggplot2::scale_linetype_identity(guide = "none")
  }

  # (d) Level axis labels in left margin
  if (show_level_labels) {
    p <- p + .ba_level_labels(nodes, node_width, level_label_size)
  }

  left_margin <- if (show_level_labels) 50 else 10

  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin     = ggplot2::margin(10, 10, 10, left_margin)
    )
}

# Build a node-only plot for the drop_pruned degenerate case (<=1 kept level).
.ba_degenerate_plot <- function(
  nodes, node_width, node_height, show_level_labels, level_label_size
) {
  p <- ggplot2::ggplot() +
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
    ggplot2::scale_fill_identity(guide = "none")

  if (show_level_labels && nrow(nodes) > 0L) {
    p <- p + .ba_level_labels(nodes, node_width, level_label_size)
  }

  left_margin <- if (show_level_labels) 50 else 10

  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin     = ggplot2::margin(10, 10, 10, left_margin)
    )
}

# Compute the level-label geom_text layer from a nodes data frame.
# The `level` column holds original level numbers (preserved even under
# compress_levels so labels read "3 factors" not "2 factors" at a re-indexed y).
.ba_level_labels <- function(nodes, node_width, level_label_size) {
  level_df <- unique(nodes[, c("level", "y"), drop = FALSE])
  pad <- node_width / 2 + 0.8
  level_df$lx <- min(nodes$x) - pad
  level_df$lt <- ifelse(
    level_df$level == 1L, "1 factor", paste(level_df$level, "factors")
  )
  ggplot2::geom_text(
    data = level_df,
    ggplot2::aes(x = .data$lx, y = .data$y, label = .data$lt),
    size = level_label_size,
    hjust = 1,
    inherit.aes = FALSE
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
