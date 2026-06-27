#' Compute a layered layout for a bass-ackwards diagram
#'
#' Returns tidy node coordinates and the edge table needed to render a
#' bass-ackwards diagram. The layout is a **layered DAG** (not a tree): level
#' 1 (one factor) sits at the top; level k (k factors) at the bottom. Horizontal
#' positions are computed with a barycenter heuristic so that each factor is
#' placed near the weighted average x-position of its parents in the level above,
#' with weights proportional to |r|. Siblings that would overlap are spread
#' symmetrically while preserving their order.
#'
#' The returned coordinates have no inherent unit. Pass the result to
#' [autoplot.ackwards()] for rendering, or use `$nodes` and `$edges` directly
#' with any graphics system.
#'
#' @param x An `ackwards` object.
#' @param min_sep Minimum horizontal separation between adjacent nodes at the
#'   same level. Default `1.0`.
#'
#' @return A list with two data frames:
#'   \item{nodes}{One row per factor: `id`, `level`, `x`, `y`, `label`.
#'     `y = -level` so level 1 is at the top.}
#'   \item{edges}{The tidy edge table from `x$edges$tidy`.}
#'
#' @seealso [autoplot.ackwards()], [tidy.ackwards()]
#'
#' @examples
#' \dontrun{
#' x <- ackwards(psych::bfi[, 1:25], k = 5)
#' lay <- ba_layout(x)
#' head(lay$nodes)
#' }
#'
#' @export
ba_layout <- function(x, min_sep = 1.0) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }

  K          <- x$k_max
  levels_lst <- x$levels
  tidy_edges <- x$edges$tidy

  # Accumulate x-positions as a named numeric per level
  node_x <- vector("list", K)
  names(node_x) <- as.character(seq_len(K))

  # Level 1: single node anchored at x = 0
  node_x[["1"]] <- stats::setNames(0, levels_lst[["1"]]$labels)

  for (k in seq(2L, K)) {
    labs_k   <- levels_lst[[as.character(k)]]$labels
    x_km1    <- node_x[[as.character(k - 1L)]]

    # Edges from level k-1 to level k
    ep <- tidy_edges[
      tidy_edges$level_from == k - 1L & tidy_edges$level_to == k, ,
      drop = FALSE
    ]

    # Barycenter: weighted-mean x of parents, weights = |r|
    bary <- vapply(labs_k, function(lab) {
      pe <- ep[ep$to == lab, , drop = FALSE]
      if (nrow(pe) == 0L) return(mean(x_km1))
      w <- abs(pe$r)
      if (sum(w) == 0) return(mean(x_km1))
      sum(w * x_km1[pe$from]) / sum(w)
    }, numeric(1L))

    node_x[[as.character(k)]] <- stats::setNames(
      .spread_positions(bary, min_sep),
      labs_k
    )
  }

  # Build nodes data frame
  nodes <- do.call(rbind, lapply(seq_len(K), function(k) {
    labs <- levels_lst[[as.character(k)]]$labels
    xs   <- node_x[[as.character(k)]]
    data.frame(
      id             = labs,
      level          = k,
      x              = xs[labs],
      y              = -k,
      label          = labs,
      stringsAsFactors = FALSE,
      row.names      = NULL
    )
  }))

  list(nodes = nodes, edges = tidy_edges)
}

# Enforce minimum separation between positions while preserving order.
# Re-centres the spread positions around the original barycenter mean.
.spread_positions <- function(bary, min_sep) {
  n <- length(bary)
  if (n == 1L) return(bary)

  bary_mean <- mean(bary)
  ord <- order(bary)
  p   <- bary[ord]

  # Forward pass: push right if too close
  for (i in seq(2L, n)) {
    if (p[i] - p[i - 1L] < min_sep) p[i] <- p[i - 1L] + min_sep
  }

  # Re-centre around the original weighted mean so parents stay central
  p <- p - mean(p) + bary_mean

  result       <- numeric(n)
  result[ord]  <- p
  result
}
