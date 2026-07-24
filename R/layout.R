#' Compute a layered layout for a bass-ackwards diagram
#'
#' Returns tidy node coordinates and the edge table needed to render a
#' bass-ackwards diagram. The layout is a layered (Sugiyama-style) barycenter
#' algorithm in two stages:
#'
#' 1. **Ordering** -- determines the left-to-right *order* of factors at each
#'    level. Two candidate orderings are scored and the crossing-minimising one
#'    kept: a single top-down |r|-weighted barycenter sweep (the historical
#'    order), and a **primary-forest traversal** -- each factor has exactly one
#'    primary parent, so the primary edges form a forest, and a depth-first,
#'    subtree-contiguous leaf order lays every subtree out as an unbroken run.
#'    In deep hierarchies (`k >= 10`) the traversal drives primary-tree crossings
#'    to zero (the "bent levels" a single pass leaves behind); keep-best scoring
#'    (lexicographic: primary crossings, then all crossings) means shallow
#'    layouts are never made worse.
#'
#' 2. **X-assignment** -- assigns actual x coordinates bottom-up. The deepest
#'    level is spread evenly; every upper-level factor is placed at the simple
#'    mean x of its **primary** children: a factor with one primary child lands
#'    directly above it; a factor with two primary children lands exactly halfway
#'    between them. Falls back to |r|-weighted mean of all children for factors
#'    with no primary children. Spreading resolves any remaining overlaps.
#'
#' After both stages the layout is shifted so that the single level-1 node is
#' always at x = 0. The result is fully deterministic.
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
#' x <- ackwards(sim16, k_max = 5)
#' lay <- ba_layout(x)
#' head(lay$nodes)
#'
#' @export
ba_layout <- function(x, min_sep = 1.0) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }

  K <- x$k_max
  levels_lst <- x$levels
  tidy_edges <- x$edges$tidy

  # --- Stage 1: Ordering ------------------------------------------------------
  # Two candidate orderings are scored and the crossing-minimising one kept:
  #   (a) the historical single top-down barycenter sweep, and
  #   (b) a primary-forest traversal that lays sibling subtrees out contiguously.
  # Each node has exactly one primary parent, so the primary edges form a forest;
  # a depth-first leaf order makes every subtree an unbroken run, which drives
  # primary-tree crossings to zero -- the deep-hierarchy "bent levels" the single
  # pass leaves behind. Scoring is lexicographic (primary crossings, then all
  # crossings) and keep-best guarantees the result is never worse than the seed.
  seed_ord <- .seed_order(levels_lst, tidy_edges, K)
  dfs_ord <- .primary_forest_order(seed_ord, levels_lst, tidy_edges, K)

  prim_edges <- tidy_edges[
    !is.na(tidy_edges$is_primary) & tidy_edges$is_primary, ,
    drop = FALSE
  ]
  score <- function(nx) {
    xmap <- do.call(c, unname(nx))
    c(
      .count_crossings_xmap(xmap, prim_edges),
      .count_crossings_xmap(xmap, tidy_edges)
    )
  }
  best_x <- .assign_x(seed_ord, levels_lst, tidy_edges, K, min_sep)
  best_cross <- score(best_x)
  cand_x <- .assign_x(dfs_ord, levels_lst, tidy_edges, K, min_sep)
  cand_cross <- score(cand_x)
  # Lexicographic improvement: fewer primary crossings, or equal primary and
  # fewer total crossings.
  if (cand_cross[1L] < best_cross[1L] ||
    (cand_cross[1L] == best_cross[1L] && cand_cross[2L] < best_cross[2L])) {
    best_x <- cand_x
  }
  node_x <- best_x

  # --- Global shift: level-1 node is always at x = 0 -------------------------
  offset <- node_x[["1"]][[1L]]
  node_x <- lapply(node_x, function(xs) xs - offset)

  # --- Build nodes data frame -------------------------------------------------
  nodes <- do.call(rbind, lapply(seq_len(K), function(k) {
    labs <- levels_lst[[as.character(k)]]$labels
    xs <- node_x[[as.character(k)]]
    data.frame(
      id = labs,
      level = k,
      x = xs[labs],
      y = -k,
      label = labs,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  list(nodes = nodes, edges = tidy_edges)
}

# Given an ackwards object with pruning annotations and a ba_layout() nodes
# data frame, returns the kept-only node set, a reduced primary-edge table, and
# the secondary-edge set for the drop_pruned rendering path.
#
# Edge selection: for each kept node, picks the single edge with the largest
# |r| to any kept node at a shallower level -- a primary-parent recomputation
# on the reduced graph. The original is_primary must NOT be reused (it was
# computed on the full adjacent lineage).
#
# Secondary edges (M79): every kept cross-level pair that is *not* a primary
# edge -- the between-level correlations the single-strongest-ancestor primary
# view hides. Returned unfiltered by |r| (autoplot applies cut_show); the set
# includes both cross-branch second parents (e.g. f1 -> d2) and same-lineage
# skip arcs (e.g. f1 -> b1 when the primary path is f1 -> d1 -> b1), consistent
# with D-032/D-017 that a direct skip-level r is a distinct, non-transitive fact.
#
# All-pairs edges are recomputed fresh from the stored levels/R (M42, fixing
# an M34 regression): x$edges$tidy holds only adjacent pairs under the default
# pairs = "adjacent", so a kept node whose adjacent ancestors were all pruned
# would find no skip-level candidate there and render edge-less. Recomputing
# mirrors prune.ackwards() (Invariant 1: one edge path via compute_edges();
# Invariant 3: recomputable from the light core) and is cheap (W'RW algebra).
# The stored weights are already sign-aligned, so recomputed adjacent edges
# are identical to x$edges$tidy's.
.drop_pruned_nodes <- function(x, nodes, compress_levels = FALSE) {
  prune_tbl <- x$prune$nodes
  tidy_edges <- compute_edges(
    levels = x$levels, R = x$r, edge_method = "auto",
    pairs = "all",
    cut_show = x$meta$cut_show %||% 0.3
  )$tidy

  kept_ids <- prune_tbl$id[!prune_tbl$pruned]
  nodes_kept <- nodes[nodes$id %in% kept_ids, , drop = FALSE]

  if (compress_levels && nrow(nodes_kept) > 0L) {
    kept_levels <- sort(unique(nodes_kept$level))
    nodes_kept$y <- -match(nodes_kept$level, kept_levels)
  }

  edge_list <- lapply(kept_ids, function(nid) {
    node_level <- nodes_kept$level[nodes_kept$id == nid][[1L]]
    candidates <- tidy_edges[
      tidy_edges$to == nid &
        tidy_edges$from %in% kept_ids &
        tidy_edges$level_from < node_level, ,
      drop = FALSE
    ]
    if (nrow(candidates) == 0L) {
      return(NULL)
    }
    candidates[which.max(abs(candidates$r)), , drop = FALSE]
  })
  edge_list <- Filter(Negate(is.null), edge_list)
  edges_kept <- if (length(edge_list) > 0L) {
    do.call(rbind, edge_list)
  } else {
    tidy_edges[0L, , drop = FALSE]
  }

  # Secondary edges: all kept cross-level pairs minus the per-node primary
  # edges just selected. tidy_edges is directed shallow -> deep, so filtering
  # to kept endpoints already restricts to cross-level pairs; the anti-join on
  # (from, to) removes exactly the primary rows (each pair is unique).
  kept_cross <- tidy_edges[
    tidy_edges$from %in% kept_ids & tidy_edges$to %in% kept_ids, ,
    drop = FALSE
  ]
  primary_key <- paste(edges_kept$from, edges_kept$to)
  is_primary_row <- paste(kept_cross$from, kept_cross$to) %in% primary_key
  secondary_kept <- kept_cross[!is_primary_row, , drop = FALSE]

  list(nodes = nodes_kept, edges = edges_kept, secondary = secondary_kept)
}

# Order every level by a depth-first traversal of the primary forest. Each node
# (level >= 2) has exactly one primary parent, so primary edges form a forest
# rooted at level 1; a pre-order walk emits each subtree as a contiguous run, so
# siblings never interleave and primary-tree crossings vanish. Children are
# visited in seed-ordinal order, keeping the result deterministic and close to
# the barycenter seed where the tree leaves order free. Returns a per-level list
# of ordinal ranks (the same shape .barycenter_sweep produces).
.primary_forest_order <- function(seed_ord, levels_lst, tidy_edges, K) {
  prim <- tidy_edges[
    !is.na(tidy_edges$is_primary) & tidy_edges$is_primary, ,
    drop = FALSE
  ]
  kids <- split(prim$to, prim$from)

  node_level <- integer(0)
  for (k in seq_len(K)) {
    labs <- levels_lst[[as.character(k)]]$labels
    node_level[labs] <- k
  }

  # Order a set of same-level nodes by their seed ordinal.
  by_seed <- function(nodes) {
    lev <- node_level[[nodes[[1L]]]]
    nodes[order(seed_ord[[as.character(lev)]][nodes])]
  }
  visit <- function(node) {
    ch <- kids[[node]]
    if (length(ch) == 0L) {
      return(node)
    }
    c(node, unlist(lapply(by_seed(ch), visit), use.names = FALSE))
  }

  roots <- levels_lst[["1"]]$labels
  seq_all <- unlist(lapply(roots, visit), use.names = FALSE)
  # Defensive: append any node the forest walk missed (deterministically), so
  # every level is fully ranked even if a primary parent were absent.
  missing <- setdiff(unlist(lapply(levels_lst, `[[`, "labels")), seq_all)
  if (length(missing) > 0L) seq_all <- c(seq_all, sort(missing)) # nocov

  pos <- stats::setNames(seq_along(seq_all), seq_all)
  ordinals <- vector("list", K)
  names(ordinals) <- as.character(seq_len(K))
  for (k in seq_len(K)) {
    labs <- levels_lst[[as.character(k)]]$labels
    ordinals[[as.character(k)]] <- stats::setNames(
      rank(pos[labs], ties.method = "first"),
      labs
    )
  }
  ordinals
}

# Seed ordering: a single top-down pass ranking each level by the |r|-weighted
# mean of its parents' ranks (the historical single-pass order). This is the
# baseline the primary-forest traversal is scored against. Deterministic -- rank
# ties break "first".
.seed_order <- function(levels_lst, tidy_edges, K) {
  ordinals <- vector("list", K)
  names(ordinals) <- as.character(seq_len(K))
  ordinals[["1"]] <- stats::setNames(1L, levels_lst[["1"]]$labels)
  for (k in seq(2L, K)) {
    labs_k <- levels_lst[[as.character(k)]]$labels
    parent_ords <- ordinals[[as.character(k - 1L)]]
    ep <- tidy_edges[
      tidy_edges$level_from == k - 1L & tidy_edges$level_to == k, ,
      drop = FALSE
    ]
    bary <- vapply(labs_k, function(lab) {
      pe <- ep[ep$to == lab, , drop = FALSE]
      w <- abs(pe$r)
      if (nrow(pe) == 0L || sum(w) == 0) { # nocov start
        return(mean(as.numeric(parent_ords)))
      } # nocov end
      sum(w * parent_ords[pe$from]) / sum(w)
    }, numeric(1L))
    ordinals[[as.character(k)]] <- stats::setNames(
      rank(bary, ties.method = "first"),
      labs_k
    )
  }
  ordinals
}

# Assign x coordinates bottom-up from a fixed per-level ordering. The deepest
# level is spread evenly by ordinal rank; every upper-level factor is placed at
# the simple mean x of its *primary* children so that a factor with one primary
# child sits directly above it and one with several sits halfway between them.
# Falls back to |r|-weighted mean of all children for orphaned parents (matching
# can assign all level-k+1 factors elsewhere). Spreading resolves overlaps.
# Returns a per-level list of named x vectors (pre global-shift).
.assign_x <- function(ordinals, levels_lst, tidy_edges, K, min_sep) {
  node_x <- vector("list", K)
  names(node_x) <- as.character(seq_len(K))

  labs_K <- levels_lst[[as.character(K)]]$labels
  ords_K <- ordinals[[as.character(K)]]
  center <- (K + 1) / 2
  node_x[[as.character(K)]] <- stats::setNames(
    (ords_K[labs_K] - center) * min_sep,
    labs_K
  )

  for (k in seq(K - 1L, 1L)) {
    labs_k <- levels_lst[[as.character(k)]]$labels
    x_below <- node_x[[as.character(k + 1L)]]

    primary_ep <- tidy_edges[
      tidy_edges$level_from == k & tidy_edges$level_to == k + 1L &
        !is.na(tidy_edges$is_primary) & tidy_edges$is_primary, ,
      drop = FALSE
    ]
    all_ep <- tidy_edges[
      tidy_edges$level_from == k & tidy_edges$level_to == k + 1L, ,
      drop = FALSE
    ]

    bary <- vapply(labs_k, function(lab) {
      pc <- primary_ep[primary_ep$from == lab, , drop = FALSE]
      if (nrow(pc) > 0L) {
        return(mean(x_below[pc$to]))
      }
      ce <- all_ep[all_ep$from == lab, , drop = FALSE] # nocov start
      w <- abs(ce$r)
      if (nrow(ce) == 0L || sum(w) == 0) {
        return(mean(x_below))
      }
      sum(w * x_below[ce$to]) / sum(w) # nocov end
    }, numeric(1L))

    node_x[[as.character(k)]] <- stats::setNames(
      .spread_positions(bary, min_sep),
      labs_k
    )
  }
  node_x
}

# Count edge crossings given an id -> x map and the tidy edge table (the Sugiyama
# crossing number the barycenter ordering minimises). Considers only edges
# spanning a single adjacent level band (level_to == level_from + 1) -- the bands
# the ordering controls; skip-level arcs (present under pairs = "all") are
# excluded because their crossings are a rendering property of the arc, not of
# the row ordering. Two same-band edges (a1->b1), (a2->b2) cross iff their
# endpoint x-orders invert: sign(x[a1]-x[a2]) != sign(x[b1]-x[b2]), neither zero.
.count_crossings_xmap <- function(xmap, edges) {
  band <- edges[edges$level_to == edges$level_from + 1L, , drop = FALSE]
  total <- 0L
  for (lf in sort(unique(band$level_from))) {
    be <- band[band$level_from == lf, , drop = FALSE]
    n <- nrow(be)
    if (n < 2L) next
    xf <- xmap[be$from]
    xt <- xmap[be$to]
    for (i in seq_len(n - 1L)) {
      for (j in seq(i + 1L, n)) {
        df <- xf[i] - xf[j]
        dt <- xt[i] - xt[j]
        if (df != 0 && dt != 0 && sign(df) != sign(dt)) {
          total <- total + 1L
        }
      }
    }
  }
  total
}

# Count crossings from a ba_layout() result (nodes + edges). Thin wrapper over
# .count_crossings_xmap for callers/tests holding a full layout.
.count_crossings <- function(layout) {
  .count_crossings_xmap(
    stats::setNames(layout$nodes$x, layout$nodes$id),
    layout$edges
  )
}

# Enforce minimum separation between positions while preserving order.
# Re-centres the spread positions around the original barycenter mean.
# Rationale: barycenter placement can pack siblings closer than a node is wide,
# so their tiles would overprint; spreading to min_sep then re-centring keeps the
# nodes legible while the group stays visually centred under its parents.
.spread_positions <- function(bary, min_sep) {
  n <- length(bary)
  if (n == 1L) {
    return(bary)
  }

  bary_mean <- mean(bary)
  ord <- order(bary)
  p <- bary[ord]

  # Forward pass: push right to enforce min_sep
  for (i in seq(2L, n)) {
    if (p[i] - p[i - 1L] < min_sep) p[i] <- p[i - 1L] + min_sep
  }

  # Re-centre so the group doesn't drift from its natural barycenter
  p <- p - mean(p) + bary_mean

  result <- numeric(n)
  result[ord] <- p
  result
}

# Dodge overlapping edge-label anchors. Edge correlation labels are drawn at
# edge midpoints; in dense (deep) hierarchies several midpoints collide (the
# .99/1.00/.64 pile-ups). Any pair closer than `threshold` is repelled apart
# along the line joining them (coincident points separate deterministically
# along x, by index) in simultaneous force-directed steps, iterating until no
# pair collides or a cap is reached. Vectorised (matrix force step, no
# interpreted pairwise loop), deterministic (no RNG) and pure -- takes/returns
# plain numeric vectors, so it is testable without ggplot2. Returns a
# list(lx, ly) of adjusted anchors.
.dodge_edge_labels <- function(lx, ly, threshold = 0.4, max_iter = 200L) {
  n <- length(lx)
  if (n < 2L) {
    return(list(lx = lx, ly = ly))
  }
  idx_sign <- outer(seq_len(n), seq_len(n), function(a, b) sign(a - b))
  for (iter in seq_len(max_iter)) {
    dx <- outer(lx, lx, "-") # dx[i, j] = lx[i] - lx[j]
    dy <- outer(ly, ly, "-")
    dist <- sqrt(dx^2 + dy^2)
    close <- dist < threshold
    diag(close) <- FALSE
    if (!any(close)) break

    # Unit vectors pointing from j toward i; coincident pairs (dist ~ 0) push
    # along x with a deterministic sign so duplicates always separate.
    coincident <- close & dist < 1e-9
    ux <- ifelse(coincident, idx_sign, dx / dist)
    uy <- ifelse(coincident, 0, dy / dist)
    shift <- (threshold - dist) / 2 + 1e-6
    ux[!close] <- 0
    uy[!close] <- 0
    shift[!close] <- 0

    lx <- lx + rowSums(ux * shift)
    ly <- ly + rowSums(uy * shift)
  }
  list(lx = lx, ly = ly)
}
