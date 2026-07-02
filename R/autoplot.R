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

#' Plot a bass-ackwards diagram or per-level fit index chart
#'
#' When `what = "hierarchy"` (default), renders the layered bass-ackwards
#' hierarchy as a ggplot2 diagram. Factors appear as labelled boxes arranged
#' in levels (level 1 at top, level k at bottom by default; set
#' `direction = "horizontal"` for a left-to-right layout). Between-level edges
#' below `cut_show` are hidden. Two aesthetics carry the edge information, each
#' explained by its own legend: **sign** (positive/negative) is shown by
#' `sign_by` (edge colour by default) and **magnitude** (`|r|`) by
#' `magnitude_by` (line thickness by default). No aesthetic is ever mapped
#' without a matching legend.
#'
#' When `what = "fit"`, renders a two-panel line plot of per-level fit indices
#' (CFI/TLI in the top panel; RMSEA/SRMR in the bottom panel) with horizontal
#' reference lines at conventional Hu & Bentler (1999) thresholds. The anchor
#' level (k = 1, saturated and always fits perfectly) is excluded. Requires an
#' EFA or ESEM engine; returns an informative empty plot for PCA (which has no
#' model-fit indices).
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @section Saving plots:
#' `autoplot()` returns a standard `ggplot` object, so save it with
#' [ggplot2::ggsave()]:
#' \preformatted{
#'   p <- autoplot(x)
#'   ggplot2::ggsave("hierarchy.png", p, width = 8, height = 6, dpi = 300)
#' }
#' `ggsave()` is not re-exported by \pkg{ackwards} (that would move
#' \pkg{ggplot2} from Suggests into Imports); call it from \pkg{ggplot2}
#' directly. For a wide slide or poster, pair it with
#' `direction = "horizontal"`.
#'
#' @param object An `ackwards` object.
#' @param what One of `"hierarchy"` (default) or `"fit"`. Controls which
#'   visualisation is produced. All other parameters are ignored when
#'   `what = "fit"`.
#' @param sign_by How edge **sign** (positive vs negative correlation) is
#'   encoded. One of `"color"` (default; positive = `color_pos`, negative =
#'   `color_neg`), `"linetype"` (positive = solid, negative = dashed),
#'   `"both"` (colour *and* linetype, with negative drawn as a distinct
#'   double-dash so it reads clearly in greyscale), or `"none"` (sign not
#'   encoded; all edges `color_edge` and solid). Whichever channel is used gets
#'   a "Direction" legend. Colour is the default because it reads sign
#'   pre-attentively and leaves linetype free; switch to `"linetype"` or
#'   `"both"` for greyscale or colour-blind-safe figures.
#' @param magnitude_by How edge **magnitude** (`|r|`) is encoded. One of
#'   `"linewidth"` (default; thicker = stronger, with a `|r|` legend) or
#'   `"none"` (constant width). A numeric `edge_linewidth` also forces constant
#'   width at that value.
#' @param cut_show Edges with `|r| < cut_show` are not drawn. Defaults to the
#'   value used when the object was created (`x$meta$cut_show`, typically 0.3).
#' @param color_pos,colour_pos Colour for positive edges when `sign_by` uses
#'   colour. Default `"#2166AC"` (blue). `colour_pos` is an accepted British
#'   alias.
#' @param color_neg,colour_neg Colour for negative edges when `sign_by` uses
#'   colour. Default `"#D6604D"` (red). `colour_neg` is an accepted British
#'   alias.
#' @param color_edge,colour_edge Single colour for all edges when `sign_by`
#'   does not use colour (`"linetype"` or `"none"`). Default `"black"`.
#'   `colour_edge` is an accepted British alias.
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
#' @param color_pruned,colour_pruned Fill colour for nodes flagged as
#'   pruned/redundant. Default `"grey80"`. Only applied when the object carries
#'   pruning annotations (`x$prune` is non-`NULL`) and `drop_pruned = FALSE`
#'   (pruned nodes are omitted entirely when `drop_pruned = TRUE`).
#'   `colour_pruned` is an accepted British alias.
#' @param show_r Whether to label each drawn edge with its correlation
#'   coefficient. Default `FALSE`. When `TRUE`, labels are formatted in APA
#'   style (leading zero stripped, e.g. `.23`, `-.30`) and placed beside each
#'   edge using a white-background label that clears the arrowhead.
#' @param r_digits Number of decimal places for edge labels when
#'   `show_r = TRUE`. Default `2L`.
#' @param r_label_size Font size for edge correlation labels when
#'   `show_r = TRUE`. Passed to `ggplot2::geom_label()` as `size`. Default
#'   `2.5`.
#' @param mono Monochrome convenience wrapper. When `TRUE`, equivalent to
#'   `sign_by = "linetype"` with black edges (it overrides `sign_by` and the
#'   colour arguments). `magnitude_by` still applies. Default `FALSE`.
#' @param show_level_labels Whether to draw level axis labels
#'   ("1 factor", "2 factors", ...) to the left of the diagram. Default `TRUE`.
#'   When the object carries pruning annotations (`x$prune` non-`NULL`) and
#'   `drop_pruned = FALSE`, a level whose factors are *all* pruned has its axis
#'   label rendered in italic to denote its status (matching the grey node
#'   fill); partially-pruned levels keep a plain label.
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
#'   pruning annotations (piped through [prune()]); errors if not.
#'   Overrides `show_skip`, `curvature`, and `primary_only`. Default `FALSE`.
#' @param compress_levels When `TRUE` under `drop_pruned = TRUE`, closes
#'   vertical gaps left by pruned levels so retained levels are evenly spaced;
#'   level axis labels (d) still show the original level numbers. Ignored when
#'   `drop_pruned = FALSE`. Default `FALSE`.
#' @param show_arrows When `FALSE`, edges are drawn with plain line ends instead
#'   of closed arrowheads (`arrow = NULL`). Applies to both straight and curved
#'   edge layers. Default `TRUE`. Forbes (2023) figures use plain line ends.
#' @param edge_linewidth `NULL` (default) uses `magnitude_by` to decide edge
#'   width. A numeric value forces every edge to that constant width (implying
#'   `magnitude_by = "none"`) and drops the `|r|` legend. Forbes figures use
#'   uniform thin lines (~= 0.5--0.6).
#' @param cut_strong **Deprecated** in M35 and
#'   ignored (with a warning). Edge magnitude is now shown by `magnitude_by`
#'   (linewidth) and sign by `sign_by`; the old strong/weak linetype split
#'   double-encoded magnitude. Retained only so existing calls do not error.
#' @param direction Layout orientation. `"vertical"` (default) stacks levels
#'   top-to-bottom (level 1 at top); `"horizontal"` lays them out left-to-right
#'   (level 1 at left), which suits wide slides or posters. Level-axis labels
#'   move to the bottom margin under `"horizontal"`.
#' @param legend When `FALSE`, suppresses all plot legends
#'   (`legend.position = "none"`). Useful when `color_pos == color_neg` (e.g.
#'   both `"black"`) to remove an otherwise redundant Direction key. Default
#'   `TRUE`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @seealso [ba_layout()], [plot.ackwards()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- ackwards(bfi25, k_max = 5)
#'   autoplot(x)
#'   autoplot(x, color_pos = "steelblue")
#'
#'   # Encode sign by linetype instead of colour, or by both
#'   autoplot(x, sign_by = "linetype")
#'   autoplot(x, sign_by = "both")
#'
#'   # Left-to-right layout (wide slides / posters)
#'   autoplot(x, direction = "horizontal")
#'
#'   # Per-level fit index chart (EFA or ESEM only)
#'   x_efa <- ackwards(bfi25, k_max = 5, engine = "efa")
#'   autoplot(x_efa, what = "fit")
#'
#'   # Monochrome with correlation labels (for greyscale figures)
#'   autoplot(x, mono = TRUE, show_r = TRUE)
#'
#'   # Custom node labels for the 5-factor level
#'   autoplot(x, node_labels = c(m5f1 = "Neuroticism", m5f2 = "Agreeableness"))
#'
#'   # Primary links only -- clean hierarchy tree
#'   autoplot(x, primary_only = TRUE)
#'
#'   # Forbes pruned view: omit redundant nodes, straight spanning arrows
#'   xp <- ackwards(bfi25, k_max = 5) |> prune("redundant")
#'   autoplot(xp, drop_pruned = TRUE)
#'   autoplot(xp, drop_pruned = TRUE, show_r = TRUE)
#'   autoplot(xp, drop_pruned = TRUE, compress_levels = TRUE)
#'
#'   # Plain line ends without arrowheads
#'   autoplot(x, show_arrows = FALSE)
#'
#'   # Uniform edge width (no |r| scaling)
#'   autoplot(x, edge_linewidth = 0.5)
#'
#'   # Suppress legend
#'   autoplot(x, legend = FALSE)
#'
#'   # Forbes (2023) publication style: black lines, uniform width, no arrowheads
#'   autoplot(xp,
#'     drop_pruned = TRUE,
#'     color_pos = "black", color_neg = "black",
#'     edge_linewidth = 0.6, show_arrows = FALSE, legend = FALSE
#'   )
#' }
#' }
#'
#' @importFrom rlang .data `%||%`
#' @export
autoplot.ackwards <- function(
  object,
  what = c("hierarchy", "fit"),
  sign_by = c("color", "linetype", "both", "none"),
  magnitude_by = c("linewidth", "none"),
  cut_show = NULL,
  color_pos = "#2166AC",
  color_neg = "#D6604D",
  color_edge = "black",
  colour_pos = NULL,
  colour_neg = NULL,
  colour_edge = NULL,
  node_width = 0.8,
  node_height = 0.4,
  min_sep = 1.0,
  show_skip = NULL,
  curvature = 0.2,
  color_pruned = "grey80",
  colour_pruned = NULL,
  show_r = FALSE,
  r_digits = 2L,
  r_label_size = 2.5,
  mono = FALSE,
  edge_linewidth = NULL,
  cut_strong = NULL,
  direction = c("vertical", "horizontal"),
  show_level_labels = TRUE,
  level_label_size = 3,
  node_labels = NULL,
  primary_only = FALSE,
  drop_pruned = FALSE,
  compress_levels = FALSE,
  show_arrows = TRUE,
  legend = TRUE,
  ...
) {
  rlang::check_installed("ggplot2", reason = "for autoplot.ackwards()")
  what <- match.arg(what)
  if (what == "fit") {
    return(.ba_fit_plot(object))
  }
  sign_by <- match.arg(sign_by)
  magnitude_by <- match.arg(magnitude_by)
  direction <- match.arg(direction)

  # British-spelling aliases: colour_* overrides color_* when supplied.
  color_pos <- colour_pos %||% color_pos
  color_neg <- colour_neg %||% color_neg
  color_edge <- colour_edge %||% color_edge
  color_pruned <- colour_pruned %||% color_pruned

  # cut_strong retired in M35: it double-encoded magnitude via linetype.
  if (!is.null(cut_strong)) {
    cli::cli_warn(c(
      "!" = "{.arg cut_strong} is deprecated (M35) and has no effect.",
      "i" = "Sign is encoded by {.arg sign_by}; magnitude by {.arg magnitude_by}."
    ))
  }

  # mono is a convenience wrapper: sign by linetype, all edges black.
  if (isTRUE(mono)) {
    sign_by <- "linetype"
    color_edge <- "black"
  }

  # Resolve which aesthetic carries which encoding.
  map_color <- sign_by %in% c("color", "both")
  map_linetype <- sign_by %in% c("linetype", "both")
  const_width <- !is.null(edge_linewidth) || magnitude_by == "none"
  map_linewidth <- !const_width
  width_val <- edge_linewidth %||% 0.7
  # sign_by = "both" draws negatives as a distinct double-dash so the redundant
  # colour+linetype pairing stays legible; plain "linetype" uses dashed.
  neg_linetype <- if (sign_by == "both") "twodash" else "dashed"

  if (min_sep < node_width) {
    cli::cli_warn(
      "{.arg min_sep} ({min_sep}) is less than {.arg node_width} ({node_width}); \\
       adjacent boxes may overlap."
    )
  }

  if (!is.null(edge_linewidth) &&
    (!is.numeric(edge_linewidth) ||
      length(edge_linewidth) != 1L ||
      edge_linewidth <= 0)) {
    cli::cli_abort(
      "{.arg edge_linewidth} must be a single positive number or {.val NULL}."
    )
  }

  if (!is.null(cut_show) &&
    (!is.numeric(cut_show) || length(cut_show) != 1L || is.na(cut_show) ||
      cut_show < 0 || cut_show > 1)) {
    cli::cli_abort(
      "{.arg cut_show} must be a single number in [0, 1] or {.val NULL}."
    )
  }
  cut_show <- cut_show %||% object$meta$cut_show %||% 0.3
  show_skip <- show_skip %||% isTRUE(object$meta$pairs == "all")

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
    e$sign_group <- ifelse(e$r >= 0, "positive", "negative")
    e
  }
  # Edge endpoints attach to the box faces that point toward the child:
  # bottom->top faces when vertical, right->left faces when horizontal.
  .attach_coords <- function(e, nx, ny) {
    if (direction == "horizontal") {
      e$x_from <- nx[e$from] + node_width / 2 # right face of parent box
      e$y_from <- ny[e$from]
      e$x_to <- nx[e$to] - node_width / 2 # left face of child box
      e$y_to <- ny[e$to]
    } else {
      e$x_from <- nx[e$from]
      e$y_from <- ny[e$from] - node_height / 2 # bottom of parent box
      e$x_to <- nx[e$to]
      e$y_to <- ny[e$to] + node_height / 2 # top of child box
    }
    e
  }
  # Transpose the layered layout to left-to-right when direction = "horizontal":
  # the level axis (encoded as -y) becomes x, the within-level spread becomes y.
  .orient_nodes <- function(nd) {
    if (direction == "horizontal") {
      old_x <- nd$x
      nd$x <- -nd$y
      nd$y <- old_x
    }
    nd
  }

  # --------------------------------------------------------------------------
  # Build draw_edges -- branching on drop_pruned
  # --------------------------------------------------------------------------

  if (drop_pruned) {
    if (is.null(object$prune)) {
      cli::cli_abort(c(
        "!" = "{.arg drop_pruned = TRUE} requires pruning annotations.",
        "i" = "Pipe through {.code prune(x, \"redundant\")} to flag factors."
      ))
    }
    if (!any(object$prune$nodes$pruned)) {
      cli::cli_warn(c(
        "!" = "{.arg drop_pruned = TRUE}: no nodes are flagged for pruning.",
        "i" = "The pruned view will be identical to the full diagram.",
        "i" = "Use {.code prune(x, \"redundant\")} to flag redundant factors."
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
        .orient_nodes(nodes), node_width, node_height,
        show_level_labels, level_label_size,
        legend = legend, direction = direction
      ))
    }

    edges <- edges[abs(edges$r) >= cut_show, , drop = FALSE]
    if (nrow(edges) == 0L) {
      cli::cli_warn(
        "No edges with |r| {cli::symbol$geq} {cut_show} to display. \\
         Try lowering {.arg cut_show}."
      )
    }

    nodes <- .orient_nodes(nodes)
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

    nodes <- .orient_nodes(nodes)
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

  arrow_spec <- if (show_arrows) {
    ggplot2::arrow(length = ggplot2::unit(0.15, "cm"), type = "closed")
  } else {
    NULL
  }

  # Build edge aesthetics. Each of sign (colour/linetype) and magnitude
  # (linewidth) enters the mapping only when the corresponding channel is
  # active; otherwise it is passed to each geom as a constant so no unmapped
  # aesthetic ever produces a legend-less encoding.
  edge_aes <- ggplot2::aes(
    x    = .data$x_from,
    y    = .data$y_from,
    xend = .data$x_to,
    yend = .data$y_to
  )
  if (map_color) {
    edge_aes <- utils::modifyList(edge_aes, ggplot2::aes(color = .data$sign_group))
  }
  if (map_linetype) {
    edge_aes <- utils::modifyList(edge_aes, ggplot2::aes(linetype = .data$sign_group))
  }
  if (map_linewidth) {
    edge_aes <- utils::modifyList(edge_aes, ggplot2::aes(linewidth = abs(.data$r)))
  }

  # Constant (non-mapped) aesthetics passed directly to each geom layer.
  const_args <- list()
  if (!map_color) const_args[["color"]] <- color_edge
  if (!map_linetype) const_args[["linetype"]] <- "solid"
  if (!map_linewidth) const_args[["linewidth"]] <- width_val

  .add_seg <- function(p, data) {
    args <- c(list(data = data, mapping = edge_aes, arrow = arrow_spec), const_args)
    p + do.call(ggplot2::geom_segment, args)
  }
  .add_cur <- function(p, data) {
    args <- c(
      list(
        data = data, mapping = edge_aes, arrow = arrow_spec,
        curvature = curvature
      ),
      const_args
    )
    p + do.call(ggplot2::geom_curve, args)
  }

  p <- ggplot2::ggplot()

  # Straight edges (adjacent in normal mode; all edges in drop_pruned mode)
  ed_str <- draw_edges[!draw_edges$curved, , drop = FALSE]
  if (nrow(ed_str) > 0L) p <- .add_seg(p, ed_str)

  # Curved arcs (skip-level edges in normal path only; never in drop_pruned)
  ed_cur <- draw_edges[draw_edges$curved, , drop = FALSE]
  if (nrow(ed_cur) > 0L) p <- .add_cur(p, ed_cur)

  # (a) Edge correlation labels -- APA style, perpendicular offset, white halo
  if (isTRUE(show_r) && nrow(draw_edges) > 0L) {
    de <- draw_edges
    de$edx <- de$x_to - de$x_from
    de$edy <- de$y_to - de$y_from
    # pmax guard: degenerate zero-length edges produce NaN coordinates without it
    de$elen <- pmax(sqrt(de$edx^2 + de$edy^2), 1e-9)
    r_nudge <- 0.15
    de$lx <- (de$x_from + de$x_to) / 2 + (-de$edy / de$elen) * r_nudge
    de$ly <- (de$y_from + de$y_to) / 2 + (de$edx / de$elen) * r_nudge
    de$rl <- .format_r(de$r, r_digits)
    p <- p + ggplot2::geom_label(
      data = de,
      ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$rl),
      size = r_label_size,
      linewidth = 0,
      label.padding = ggplot2::unit(0.1, "lines"),
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

  # Edge scales -- one per active encoding channel, each with its own legend.
  # When sign_by = "both", the colour and linetype scales share the name
  # "Direction" (and identical breaks/labels) so ggplot2 merges them into a
  # single legend key.
  if (map_color) {
    p <- p + ggplot2::scale_color_manual(
      values = c(positive = color_pos, negative = color_neg),
      breaks = c("positive", "negative"),
      labels = c(positive = "Positive", negative = "Negative"),
      name   = "Direction"
    )
  }
  if (map_linetype) {
    p <- p + ggplot2::scale_linetype_manual(
      values = c(positive = "solid", negative = neg_linetype),
      breaks = c("positive", "negative"),
      labels = c(positive = "Positive", negative = "Negative"),
      name   = "Direction"
    )
  }
  if (map_linewidth) {
    p <- p + ggplot2::scale_linewidth_continuous(
      range = c(0.4, 1.8), name = "|r|", limits = c(cut_show, 1)
    )
  }

  # (d) Level axis labels: left margin when vertical, bottom margin when horizontal.
  # Fully-pruned levels (M40) are italicised. Only the normal path retains such
  # levels' nodes; under drop_pruned they were removed, so nothing matches.
  if (show_level_labels) {
    p <- p + .ba_level_labels(
      nodes, node_width, level_label_size,
      direction = direction, node_height = node_height,
      pruned_levels = if (drop_pruned) integer(0) else .fully_pruned_levels(object)
    )
  }

  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = if (legend) "right" else "none",
      plot.margin     = .ba_label_margin(show_level_labels, direction)
    )
}

# Plot margin that leaves room for the level-axis labels: on the left when
# vertical, along the bottom when horizontal.
.ba_label_margin <- function(show_level_labels, direction) {
  pad <- if (show_level_labels) 50 else 10
  if (direction == "horizontal") {
    ggplot2::margin(10, 10, pad, 10)
  } else {
    ggplot2::margin(10, 10, 10, pad)
  }
}

# Build a node-only plot for the drop_pruned degenerate case (<=1 kept level).
.ba_degenerate_plot <- function(
  nodes, node_width, node_height, show_level_labels, level_label_size,
  legend = TRUE, direction = "vertical"
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
    p <- p + .ba_level_labels(
      nodes, node_width, level_label_size,
      direction = direction, node_height = node_height
    )
  }

  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = if (legend) "right" else "none",
      plot.margin     = .ba_label_margin(show_level_labels, direction)
    )
}

# Levels for which *every* node is flagged pruned (M40). A fully-pruned level's
# axis label is italicised in the normal render path to denote its status,
# mirroring the automatic grey node fill (a partially-pruned level keeps a plain
# label -- its retained factors are still substantive). Returns original level
# numbers (matching nodes$level). integer(0) when the object carries no pruning
# annotations or nothing is flagged.
.fully_pruned_levels <- function(object) {
  pn <- object$prune$nodes
  if (is.null(object$prune) || is.null(pn) || !any(pn$pruned)) {
    return(integer(0))
  }
  by_level <- tapply(pn$pruned, pn$level, FUN = all)
  as.integer(names(by_level)[by_level])
}

# Compute the level-label geom_text layer from a nodes data frame.
# The `level` column holds original level numbers (preserved even under
# compress_levels so labels read "3 factors" not "2 factors" at a re-indexed y).
# Labels sit to the left of the diagram when vertical, below it when horizontal.
.ba_level_labels <- function(nodes, node_width, level_label_size,
                             direction = "vertical", node_height = 0.4,
                             pruned_levels = integer(0)) {
  labels_for <- function(level) {
    ifelse(level == 1L, "1 factor", paste(level, "factors"))
  }
  # Fully-pruned levels (M40) are italicised; all others plain. fontface is a
  # per-row geom_text aesthetic, so it is carried as a column.
  face_for <- function(level) {
    ifelse(level %in% pruned_levels, "italic", "plain")
  }
  if (direction == "horizontal") {
    level_df <- unique(nodes[, c("level", "x"), drop = FALSE])
    level_df$ly <- min(nodes$y) - (node_height / 2 + 0.8)
    level_df$lt <- labels_for(level_df$level)
    level_df$lf <- face_for(level_df$level)
    ggplot2::geom_text(
      data = level_df,
      ggplot2::aes(
        x = .data$x, y = .data$ly, label = .data$lt, fontface = .data$lf
      ),
      size = level_label_size,
      vjust = 1,
      inherit.aes = FALSE
    )
  } else {
    level_df <- unique(nodes[, c("level", "y"), drop = FALSE])
    level_df$lx <- min(nodes$x) - (node_width / 2 + 0.8)
    level_df$lt <- labels_for(level_df$level)
    level_df$lf <- face_for(level_df$level)
    ggplot2::geom_text(
      data = level_df,
      ggplot2::aes(
        x = .data$lx, y = .data$y, label = .data$lt, fontface = .data$lf
      ),
      size = level_label_size,
      hjust = 1,
      inherit.aes = FALSE
    )
  }
}

# Build the per-level fit index chart (autoplot.ackwards(what = "fit")).
# Two-panel: CFI/TLI (higher = better) on top; RMSEA/SRMR (lower = better)
# on bottom. Anchor level k=1 excluded (saturated, always fits perfectly).
# Horizontal reference lines at Hu & Bentler (1999) conventional thresholds.
# Returns an empty annotated ggplot for PCA (no model-fit indices).
.ba_fit_plot <- function(object) {
  engine <- object$engine
  if (engine == "pca") {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5, y = 0.5, size = 4, color = "grey40",
          label = paste0(
            "Fit indices are not available for engine = \"pca\".\n",
            "Use engine = \"efa\" or \"esem\" to obtain model fit."
          )
        ) +
        ggplot2::theme_void()
    )
  }

  fit_long <- .tidy_fit(object)
  # Panels split by direction of "good": higher-is-better vs lower-is-better.
  # EFA reports only TLI / RMSEA; ESEM adds CFI / SRMR. Panel titles are built
  # from the indices actually present so EFA is not mislabelled "CFI / TLI".
  if (engine == "efa") {
    hi_idx <- c("TLI")
    lo_idx <- c("RMSEA")
  } else {
    hi_idx <- c("CFI", "TLI")
    lo_idx <- c("RMSEA", "SRMR")
  }
  keep_idx <- c(hi_idx, lo_idx)
  hi_lbl <- paste0(paste(hi_idx, collapse = " / "), "  (higher is better)")
  lo_lbl <- paste0(paste(lo_idx, collapse = " / "), "  (lower is better)")

  # Anchor level (k = 1) excluded: saturated baseline, always fits perfectly.
  fd <- fit_long[fit_long$level > 1L & fit_long$statistic %in% keep_idx, , drop = FALSE]

  if (nrow(fd) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5, y = 0.5, size = 4, color = "grey40",
          label = "No fit indices available."
        ) +
        ggplot2::theme_void()
    )
  }

  fd$panel <- ifelse(fd$statistic %in% hi_idx, hi_lbl, lo_lbl)
  fd$panel <- factor(fd$panel, levels = c(hi_lbl, lo_lbl))

  # Every index in keep_idx has a threshold in .fit_cutoffs(), so no NULL guard
  # is needed here.
  cuts <- .fit_cutoffs()
  ref_df <- do.call(rbind, lapply(keep_idx, function(idx) {
    panel_lbl <- if (idx %in% hi_idx) hi_lbl else lo_lbl
    data.frame(panel = panel_lbl, yintercept = cuts[[idx]]$threshold, stringsAsFactors = FALSE)
  }))
  ref_df <- unique(ref_df)
  ref_df$panel <- factor(ref_df$panel, levels = levels(fd$panel))

  ggplot2::ggplot(fd, ggplot2::aes(x = .data$level, y = .data$value, color = .data$statistic)) +
    ggplot2::geom_hline(
      data = ref_df,
      ggplot2::aes(yintercept = .data$yintercept),
      linetype = "dashed", color = "grey60", linewidth = 0.5
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~panel, ncol = 1L, scales = "free_y") +
    ggplot2::scale_x_continuous(
      name = "Number of factors (level)",
      breaks = function(x) seq(ceiling(x[1L]), floor(x[2L]), by = 1L)
    ) +
    ggplot2::scale_color_brewer(palette = "Set1", name = "Index") +
    ggplot2::labs(
      y = "Value",
      caption = "Dashed lines: Hu & Bentler (1999) thresholds (CFI/TLI >= .95, RMSEA <= .06, SRMR <= .08)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
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
