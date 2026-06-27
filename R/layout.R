#' Compute a layered layout for a bass-ackwards diagram
#'
#' Returns tidy node coordinates and the edge table needed to render a
#' bass-ackwards diagram. The layout uses a two-pass barycenter algorithm:
#'
#' 1. **Top-down pass** -- determines the left-to-right *order* of factors at
#'    each level. Each factor's ordinal rank is the |r|-weighted mean of its
#'    parents' ranks in the level above, so siblings that share a parent are
#'    grouped together.
#'
#' 2. **Bottom-up pass** -- assigns actual x coordinates. The deepest level
#'    (level k) is spread evenly; every upper-level factor is placed at the
#'    |r|-weighted mean x of its *children* in the level below. This ensures
#'    parents sit directly above (or between) their children rather than merely
#'    near them.
#'
#' After both passes the layout is shifted so that the single level-1 node is
#' always at x = 0.
#'
#' @param x An `ackwards` object.
#' @param min_sep Minimum horizontal separation between adjacent nodes at the
#'   same level. Default `1.0`. Increase for wider diagrams.
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

  # --- Pass 1: Top-down -- determine ordinal order at each level --------------
  # For each level k, rank factors by the |r|-weighted mean of their parents'
  # ranks. This groups siblings together and minimises crossings.
  ordinals        <- vector("list", K)
  names(ordinals) <- as.character(seq_len(K))
  ordinals[["1"]] <- stats::setNames(1L, levels_lst[["1"]]$labels)

  for (k in seq(2L, K)) {
    labs_k      <- levels_lst[[as.character(k)]]$labels
    parent_ords <- ordinals[[as.character(k - 1L)]]

    ep <- tidy_edges[
      tidy_edges$level_from == k - 1L & tidy_edges$level_to == k, ,
      drop = FALSE
    ]

    bary <- vapply(labs_k, function(lab) {
      pe <- ep[ep$to == lab, , drop = FALSE]
      w  <- abs(pe$r)
      if (nrow(pe) == 0L || sum(w) == 0) return(mean(as.numeric(parent_ords)))
      sum(w * parent_ords[pe$from]) / sum(w)
    }, numeric(1L))

    ordinals[[as.character(k)]] <- stats::setNames(
      rank(bary, ties.method = "first"),
      labs_k
    )
  }

  # --- Pass 2: Bottom-up -- assign x coordinates ----------------------------
  # The deepest level is evenly spread; every upper-level factor is placed at
  # the |r|-weighted mean of its children's x positions so parents sit directly
  # above (or centred between) their children.
  node_x        <- vector("list", K)
  names(node_x) <- as.character(seq_len(K))

  # Level K: evenly spaced by ordinal rank, centred at 0
  labs_K  <- levels_lst[[as.character(K)]]$labels
  ords_K  <- ordinals[[as.character(K)]]
  center  <- (K + 1) / 2
  node_x[[as.character(K)]] <- stats::setNames(
    (ords_K[labs_K] - center) * min_sep,
    labs_K
  )

  # Levels K-1 down to 1: barycenter of children, then spread to prevent overlap
  for (k in seq(K - 1L, 1L)) {
    labs_k  <- levels_lst[[as.character(k)]]$labels
    x_below <- node_x[[as.character(k + 1L)]]

    ep <- tidy_edges[
      tidy_edges$level_from == k & tidy_edges$level_to == k + 1L, ,
      drop = FALSE
    ]

    bary <- vapply(labs_k, function(lab) {
      ce <- ep[ep$from == lab, , drop = FALSE]   # edges to children
      w  <- abs(ce$r)
      if (nrow(ce) == 0L || sum(w) == 0) return(mean(x_below))
      sum(w * x_below[ce$to]) / sum(w)
    }, numeric(1L))

    node_x[[as.character(k)]] <- stats::setNames(
      .spread_positions(bary, min_sep),
      labs_k
    )
  }

  # --- Global shift: level-1 node is always at x = 0 -------------------------
  offset <- node_x[["1"]][[1L]]
  node_x <- lapply(node_x, function(xs) xs - offset)

  # --- Build nodes data frame -------------------------------------------------
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

  # Forward pass: push right to enforce min_sep
  for (i in seq(2L, n)) {
    if (p[i] - p[i - 1L] < min_sep) p[i] <- p[i - 1L] + min_sep
  }

  # Re-centre so the group doesn't drift from its natural barycenter
  p <- p - mean(p) + bary_mean

  result      <- numeric(n)
  result[ord] <- p
  result
}
