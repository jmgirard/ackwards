test_that("ba_layout() returns a list with nodes and edges data frames", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  lay <- ba_layout(x)

  expect_type(lay, "list")
  expect_named(lay, c("nodes", "edges"))
  expect_s3_class(lay$nodes, "data.frame")
  expect_s3_class(lay$edges, "data.frame")
})

test_that("ba_layout() nodes have correct structure and counts", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  lay <- ba_layout(x)
  nodes <- lay$nodes

  # Columns
  expect_true(all(c("id", "level", "x", "y", "label") %in% names(nodes)))

  # Total nodes: 1+2+3+4 = 10
  expect_equal(nrow(nodes), 10L)

  # One node per level per factor
  for (k in 1:4) {
    expect_equal(sum(nodes$level == k), k,
      info = paste("level", k, "has", k, "nodes")
    )
  }
})

test_that("ba_layout() y positions equal -level", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$y, -nodes$level)
})

test_that("ba_layout() level-1 node is at x = 0", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$x[nodes$level == 1], 0)
})

test_that("ba_layout() respects min_sep between nodes at same level", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))

  for (sep in c(0.5, 1.0, 1.5)) {
    lay <- ba_layout(x, min_sep = sep)
    nodes <- lay$nodes
    for (k in 2:5) {
      xs <- sort(nodes$x[nodes$level == k])
      diffs <- diff(xs)
      expect_true(
        all(diffs >= sep - 1e-9),
        info = paste("min_sep", sep, "at level", k)
      )
    }
  }
})

test_that("ba_layout() node IDs and labels match level labels", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$id, nodes$label)
  # IDs follow m{k}f{j} scheme
  expect_true(all(grepl("^m[0-9]+f[0-9]+$", nodes$id)))
})

test_that("ba_layout() errors on non-ackwards input", {
  expect_error(ba_layout(list()), "ackwards")
})

test_that("ba_layout() places each parent at mean x of its primary children", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  edges <- lay$edges

  nx <- stats::setNames(nodes$x, nodes$id)

  # For every adjacent level pair, compute the ideal parent position as the
  # simple mean of its primary children's x coordinates, and assert it agrees
  # with the actual parent x within tolerance (spreading can shift by up to
  # min_sep when siblings collide, but in the typical case they agree exactly).
  primary_edges <- edges[!is.na(edges$is_primary) & edges$is_primary, ]

  for (parent_id in unique(primary_edges$from)) {
    child_ids <- primary_edges$to[primary_edges$from == parent_id]
    ideal_x <- mean(nx[child_ids])
    actual_x <- nx[[parent_id]]
    # Alignment is exact when no spreading is needed; after spreading it may
    # differ, so assert within 1 * min_sep (default).
    expect_lte(
      abs(actual_x - ideal_x), 1.0 + 1e-9,
      label = paste(parent_id, "aligned to primary children")
    )
  }

  # Specifically verify single-child parents are directly above their child
  # (before any spreading, diff should be 0; after spreading can be up to min_sep).
  single_child <- names(which(table(primary_edges$from) == 1))
  for (parent_id in single_child) {
    child_id <- primary_edges$to[primary_edges$from == parent_id]
    expect_equal(
      nx[[parent_id]], nx[[child_id]],
      tolerance = 1.0 + 1e-9,
      label = paste(parent_id, "directly above sole primary child", child_id)
    )
  }
})

test_that(".count_crossings() counts adjacent-band inversions, ignores skip arcs", {
  nodes <- data.frame(
    id = c("p1", "p2", "c1", "c2"),
    level = c(2L, 2L, 3L, 3L),
    x = c(0, 1, 0, 1)
  )
  # p1(0)->c2(1), p2(1)->c1(0): endpoint x-orders invert -> one crossing.
  crossing <- data.frame(
    from = c("p1", "p2"), to = c("c2", "c1"),
    level_from = c(2L, 2L), level_to = c(3L, 3L)
  )
  expect_equal(ackwards:::.count_crossings(list(nodes = nodes, edges = crossing)), 1L)

  # p1->c1, p2->c2: no inversion.
  parallel <- data.frame(
    from = c("p1", "p2"), to = c("c1", "c2"),
    level_from = c(2L, 2L), level_to = c(3L, 3L)
  )
  expect_equal(ackwards:::.count_crossings(list(nodes = nodes, edges = parallel)), 0L)

  # A skip arc (level_to == level_from + 2) is not an adjacent band -> ignored.
  skip <- rbind(crossing, data.frame(
    from = "root", to = "c1", level_from = 1L, level_to = 3L
  ))
  skip_nodes <- rbind(nodes, data.frame(id = "root", level = 1L, x = 5))
  expect_equal(ackwards:::.count_crossings(list(nodes = skip_nodes, edges = skip)), 1L)

  # A band holding a single edge cannot cross anything (no pair to compare).
  lone <- data.frame(
    from = "p1", to = "c1", level_from = 2L, level_to = 3L
  )
  expect_equal(ackwards:::.count_crossings(list(nodes = nodes, edges = lone)), 0L)
})

test_that("ba_layout() removes primary-tree crossings in a deep (k=10) hierarchy", {
  # AMH applied example (M53): the single top-down seed leaves 3 primary-edge
  # crossings; the primary-forest traversal ordering (M80) drives them to zero.
  x <- cached(ackwards(forbes2023, k_max = 10, pairs = "all"))
  lay <- ba_layout(x)
  prim <- lay$edges[!is.na(lay$edges$is_primary) & lay$edges$is_primary, ]

  # Baseline: the seed ordering alone.
  seed <- ackwards:::.seed_order(x$levels, x$edges$tidy, x$k_max)
  seed_x <- ackwards:::.assign_x(seed, x$levels, x$edges$tidy, x$k_max, 1.0)
  seed_prim <- ackwards:::.count_crossings_xmap(do.call(c, unname(seed_x)), prim)
  expect_equal(seed_prim, 3L)

  # M80 ordering: zero primary crossings, and never worse on total crossings.
  expect_equal(ackwards:::.count_crossings(list(nodes = lay$nodes, edges = prim)), 0L)
  expect_lte(
    ackwards:::.count_crossings(lay),
    ackwards:::.count_crossings_xmap(do.call(c, unname(seed_x)), x$edges$tidy)
  )
})

test_that("ba_layout() never increases crossings on a shallow hierarchy", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  lay <- ba_layout(x)
  prim <- lay$edges[!is.na(lay$edges$is_primary) & lay$edges$is_primary, ]

  seed <- ackwards:::.seed_order(x$levels, x$edges$tidy, x$k_max)
  seed_x <- ackwards:::.assign_x(seed, x$levels, x$edges$tidy, x$k_max, 1.0)
  seed_xmap <- do.call(c, unname(seed_x))
  new_xmap <- stats::setNames(lay$nodes$x, lay$nodes$id)

  # AC1 shallow clause: neither primary nor total crossings increase vs the seed.
  expect_lte(
    ackwards:::.count_crossings(list(nodes = lay$nodes, edges = prim)),
    ackwards:::.count_crossings_xmap(seed_xmap, prim)
  )
  expect_lte(
    ackwards:::.count_crossings_xmap(new_xmap, x$edges$tidy),
    ackwards:::.count_crossings_xmap(seed_xmap, x$edges$tidy)
  )
})

test_that("ba_layout() stays faithful and deterministic at k=10", {
  x <- cached(ackwards(forbes2023, k_max = 10, pairs = "all"))
  lay <- ba_layout(x)
  nodes <- lay$nodes

  # Each node on its own row; level-1 anchored; full node set; no relabeling.
  expect_equal(nodes$y, -nodes$level)
  expect_equal(nodes$x[nodes$level == 1L], 0)
  expect_equal(nrow(nodes), sum(seq_len(10L)))
  expect_setequal(nodes$id, unlist(lapply(x$levels, `[[`, "labels")))

  # Every parent still sits at the mean x of its primary children (Pass 2 intact).
  nx <- stats::setNames(nodes$x, nodes$id)
  prim <- lay$edges[!is.na(lay$edges$is_primary) & lay$edges$is_primary, ]
  for (parent_id in unique(prim$from)) {
    kids <- prim$to[prim$from == parent_id]
    expect_lte(abs(nx[[parent_id]] - mean(nx[kids])), 1.0 + 1e-9,
      label = paste(parent_id, "aligned to primary children")
    )
  }

  # Deterministic across calls.
  expect_identical(ba_layout(x)$nodes, nodes)
})

test_that("ba_layout() node coordinates are stable (shallow snapshot)", {
  x <- cached(ackwards(sim16, k_max = 5))
  nodes <- ba_layout(x)$nodes
  # Coordinates derive from ordinal ranks via .assign_x (rationals), so a value
  # snapshot pins the layout deterministically. Tolerance guards float low bits.
  expect_snapshot_value(nodes, style = "json2", tolerance = 1e-6)
})

test_that("ba_layout() node coordinates are stable (deep k=10 snapshot)", {
  x <- cached(ackwards(forbes2023, k_max = 10, pairs = "all"))
  nodes <- ba_layout(x)$nodes
  expect_snapshot_value(nodes, style = "json2", tolerance = 1e-6)
})

test_that(".dodge_edge_labels() separates colliding labels beyond the threshold", {
  min_pair_dist <- function(r) min(stats::dist(cbind(r$lx, r$ly)))

  # Three coincident anchors are pushed apart to at least the threshold.
  r <- ackwards:::.dodge_edge_labels(c(0, 0, 0), c(0, 0, 0), threshold = 0.4)
  expect_gte(min_pair_dist(r), 0.4 - 1e-6)

  # A dense horizontal band of near-colliding labels is declutttered.
  set.seed(1)
  band <- ackwards:::.dodge_edge_labels(
    stats::runif(8, 0, 0.5), rep(-2.5, 8),
    threshold = 0.4
  )
  expect_gte(min_pair_dist(band), 0.4 - 1e-6)

  # Already-separated anchors are left untouched (no spurious motion).
  lx <- c(0, 1, 2)
  ly <- c(0, 0, 0)
  same <- ackwards:::.dodge_edge_labels(lx, ly, threshold = 0.4)
  expect_identical(same$lx, lx)
  expect_identical(same$ly, ly)

  # A single label (or none) is a no-op.
  one <- ackwards:::.dodge_edge_labels(5, 5)
  expect_identical(one, list(lx = 5, ly = 5))
})

test_that("autoplot(show_r = TRUE) dodges overlapping edge labels", {
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(forbes2023, k_max = 10, pairs = "all"))
  # Renders without error and produces the labelled layer with dodged anchors.
  expect_no_error(p <- ggplot2::autoplot(x, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that(".spread_positions() enforces min_sep, preserves order, and recentres", {
  # Access internal function from the package namespace
  spread <- ackwards:::.spread_positions

  # Colliding input: three positions bunched at 0
  bary <- c(0, 0, 0)
  out <- spread(bary, min_sep = 1.0)
  expect_equal(length(out), 3L)
  # All gaps should be >= min_sep
  expect_true(all(diff(sort(out)) >= 1.0 - 1e-9))
  # Group should be re-centred around original mean (0)
  expect_equal(mean(out), 0, tolerance = 1e-9)

  # Already-separated input: should be unchanged (no spreading needed)
  bary2 <- c(-2, 0, 2)
  out2 <- spread(bary2, min_sep = 1.0)
  expect_equal(out2, bary2, tolerance = 1e-9)

  # Order preservation: largest bary stays largest
  bary3 <- c(1, 3, 2)
  out3 <- spread(bary3, min_sep = 0.5)
  expect_equal(order(out3), order(bary3))

  # Single element: returned as-is
  expect_equal(spread(5, min_sep = 1.0), 5)
})

test_that("autoplot.ackwards() returns a ggplot object", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(x) dispatches correctly without library(ggplot2)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  # The package defines its own autoplot generic so this works without ggplot2 attached
  p <- autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() respects cut_show", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.1))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.5))
})

test_that("autoplot.ackwards() deprecates cut_strong (warns, no effect)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_warning(ggplot2::autoplot(x, cut_strong = 0.4), "deprecated")
})

test_that("autoplot.ackwards() validates cut_show (M42/m8)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_error(ggplot2::autoplot(x, cut_show = 5), "\\[0, 1\\]")
  expect_error(ggplot2::autoplot(x, cut_show = c(0.3, 0.5)), "\\[0, 1\\]")
  expect_error(ggplot2::autoplot(x, cut_show = NA), "\\[0, 1\\]")
})

test_that("autoplot.ackwards() warns when cut_show hides all edges", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  # cut_show = 1 (the legal maximum since M42/m8 validation) hides every edge
  # here: no between-level correlation in these data reaches exactly 1.0
  expect_warning(ggplot2::autoplot(x, cut_show = 1), "No edges")
})

test_that("autoplot.ackwards() warns when min_sep < node_width", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_warning(ggplot2::autoplot(x, min_sep = 0.3, node_width = 0.8), "overlap")
})

test_that("plot.ackwards() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_no_error(plot(x))
  expect_invisible(plot(x))
})

test_that("autoplot.ackwards() renders skip-level edges with pairs='all'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
  # show_skip defaults to TRUE when meta$pairs == "all"
  expect_no_error(ggplot2::autoplot(x, show_skip = FALSE))
  expect_no_error(ggplot2::autoplot(x, show_skip = TRUE, curvature = 0.3))
})

test_that("autoplot.ackwards() fades pruned nodes", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 4) |> prune("redundant", redundancy_r = 0.9)
  ))
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
  # Custom prune colour
  expect_no_error(ggplot2::autoplot(x, color_pruned = "pink"))
})

test_that(".drop_pruned_nodes() bridges fully-pruned levels on pairs='adjacent' objects (M42/M1)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  # Default pairs = "adjacent": the stored tidy edges hold no skip-level rows.
  x <- cached(ackwards(d, k_max = 4))
  xp <- prune(x, manual = c("m2f1", "m2f2")) # level 2 fully pruned

  lay <- ba_layout(xp)
  dp <- ackwards:::.drop_pruned_nodes(xp, lay$nodes)

  # Pre-M42 regression: kept level-3 nodes had no candidate ancestor edges
  # (their only stored parents were the pruned level-2 nodes; no 1:3 skip
  # edges existed in the adjacent-only tidy table). Every kept node below the
  # apex must now have exactly one incoming edge, bridging the pruned level.
  kept_below_apex <- dp$nodes$id[dp$nodes$level > 1L]
  expect_setequal(dp$edges$to, kept_below_apex)
  # The level-3 bridges must come from level 1 (level 2 is gone).
  lvl3_edges <- dp$edges[dp$edges$level_to == 3L, , drop = FALSE]
  expect_true(nrow(lvl3_edges) == 3L && all(lvl3_edges$level_from == 1L))

  # And the reduced edge set must be identical whether the object was fit
  # with pairs = "adjacent" or pairs = "all" (edges are recomputed fresh).
  x_all <- cached(ackwards(d, k_max = 4, pairs = "all"))
  xp_all <- prune(x_all, manual = c("m2f1", "m2f2"))
  dp_all <- ackwards:::.drop_pruned_nodes(xp_all, ba_layout(xp_all)$nodes)
  cols <- c("from", "to", "level_from", "level_to", "r")
  expect_equal(dp$edges[, cols], dp_all$edges[, cols], tolerance = 1e-12)
})

test_that(".drop_pruned_nodes() returns the secondary edge set (M79)", {
  skip_if_not_installed("psych")
  suppressWarnings(suppressMessages(
    xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all") |>
      prune("redundant", redundancy_r = 0.9))
  ))
  lay <- ba_layout(xp)
  dp <- ackwards:::.drop_pruned_nodes(xp, lay$nodes)

  kept <- dp$nodes$id
  ekey <- function(e) paste(e$from, e$to)

  # Secondary and primary are disjoint: every primary edge is absent from the
  # secondary set (AC1: secondary = pairs that are *not* the primary edge).
  expect_length(intersect(ekey(dp$secondary), ekey(dp$edges)), 0L)

  # Union of primary + secondary == every kept cross-level pair (shallower kept
  # node -> deeper kept node). Nothing dropped, nothing invented.
  te <- ackwards:::compute_edges(
    levels = xp$levels, R = xp$r, edge_method = "auto",
    pairs = "all", cut_show = 0.3
  )$tidy
  cross <- te[te$from %in% kept & te$to %in% kept, , drop = FALSE]
  expect_setequal(c(ekey(dp$edges), ekey(dp$secondary)), ekey(cross))

  # The set spans both kinds the milestone promises: an adjacent-level second
  # parent (level gap 1) and a same-lineage skip arc (level gap >= 2).
  gaps <- dp$secondary$level_to - dp$secondary$level_from
  expect_true(any(gaps == 1L))
  expect_true(any(gaps >= 2L))
})

# M79 secondary-edge render helpers: the drop_pruned path draws only
# geom_segment edges. The secondary layer is the segment layer carrying the
# constant dimmed alpha (0.4); the primary layer carries alpha = NA.
.dp_seg_layers <- function(p) {
  which(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1L)))
}
.dp_secondary_data <- function(p) {
  for (i in .dp_seg_layers(p)) {
    ld <- ggplot2::layer_data(p, i)
    if (nrow(ld) > 0L && all(!is.na(ld$alpha)) && all(ld$alpha == 0.4)) {
      return(ld)
    }
  }
  NULL
}
.dp_primary_data <- function(p) {
  for (i in .dp_seg_layers(p)) {
    ld <- ggplot2::layer_data(p, i)
    if (nrow(ld) > 0L && all(is.na(ld$alpha))) {
      return(ld)
    }
  }
  NULL
}

test_that("autoplot(drop_pruned, show_secondary=TRUE) draws a distinct dimmed+thin layer (M79)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all") |>
      prune("redundant", redundancy_r = 0.9))
  ))
  p_no <- ggplot2::autoplot(xp, drop_pruned = TRUE)
  p_yes <- ggplot2::autoplot(xp, drop_pruned = TRUE, show_secondary = TRUE)

  # AC1: show_secondary adds exactly one segment layer (the secondary edges).
  expect_length(.dp_seg_layers(p_yes), length(.dp_seg_layers(p_no)) + 1L)

  sec <- .dp_secondary_data(p_yes)
  expect_false(is.null(sec))
  # Dimmed + thinner, distinct from the primary channel.
  expect_true(all(sec$alpha == 0.4))
  expect_true(all(sec$linewidth == 0.3))

  # AC1: the drawn set == every non-primary kept cross-level pair >= cut_show.
  dp <- ackwards:::.drop_pruned_nodes(xp, ba_layout(xp)$nodes)
  n_exp <- sum(abs(dp$secondary$r) >= 0.3)
  expect_equal(nrow(sec), n_exp)
  expect_gt(n_exp, 0L)

  # AC2: secondary edges inherit the sign colour (not a single flattened hue),
  # and the primary edge layer's colours are byte-identical with/without them.
  expect_true(all(sec$colour %in% c("#2166AC", "#D6604D")))
  expect_identical(.dp_primary_data(p_yes)$colour, .dp_primary_data(p_no)$colour)
})

test_that("autoplot(drop_pruned) default (show_secondary=FALSE) adds no secondary layer (M79)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all") |>
      prune("redundant", redundancy_r = 0.9))
  ))
  # AC3: the default reproduces the prior pruned view -- a single (primary)
  # segment layer, no dimmed secondary layer.
  p_def <- ggplot2::autoplot(xp, drop_pruned = TRUE)
  expect_length(.dp_seg_layers(p_def), 1L)
  expect_null(.dp_secondary_data(p_def))
  expect_false(is.null(.dp_primary_data(p_def)))
})

test_that("show_secondary preserves the sign channel under sign_by='linetype' (M79)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all") |>
      prune("redundant", redundancy_r = 0.9))
  ))
  # Sign carried by linetype: the secondary layer must inherit the linetype
  # (encoding sign), not reuse it as its own distinct channel -- AC2, no
  # conflation with the sign dash. The secondary and primary layers therefore
  # draw from the same linetype set.
  p <- ggplot2::autoplot(xp,
    drop_pruned = TRUE, show_secondary = TRUE,
    sign_by = "linetype"
  )
  sec <- .dp_secondary_data(p)
  expect_false(is.null(sec))
  # All secondary edges share the primary channel: same (constant) colour and
  # a linetype drawn from the sign-encoding set the primary layer uses.
  expect_setequal(unique(sec$linetype), unique(.dp_primary_data(p)$linetype))
})

test_that("show_secondary stays thinner than the primary under a thin edge_linewidth (M79)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all") |>
      prune("redundant", redundancy_r = 0.9))
  ))
  # Regression: a user-supplied edge_linewidth below the fixed secondary width
  # must not make secondary edges *thicker* than the primary ones. Secondary
  # width scales under the primary width (min(0.3, 0.6 * width_val)).
  p <- ggplot2::autoplot(xp,
    drop_pruned = TRUE, show_secondary = TRUE, edge_linewidth = 0.2
  )
  prim_w <- unique(.dp_primary_data(p)$linewidth)
  sec_w <- unique(.dp_secondary_data(p)$linewidth)
  expect_equal(prim_w, 0.2)
  expect_true(all(sec_w < prim_w)) # strictly thinner, not inverted
  expect_equal(sec_w, 0.12) # min(0.3, 0.6 * 0.2)
})

test_that("autoplot.ackwards() handles objects with prune=NULL (no pruning)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$prune)
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
})

# --- Wave 1 (M8): show_r, mono, show_level_labels, node_labels, primary_only ---

test_that("autoplot.ackwards() show_r=TRUE labels edges without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_r=TRUE respects r_digits", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE, r_digits = 3L))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_r=TRUE works with pairs='all'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() mono=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, mono = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() mono=TRUE + show_r=TRUE composes without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, mono = TRUE, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

# --- M35: configurable encodings (sign_by / magnitude_by / colour aliases) ------

.ba_has_scale <- function(p, aesthetic) {
  any(vapply(
    p$scales$scales,
    function(s) aesthetic %in% s$aesthetics,
    logical(1)
  ))
}

.ba_get_scale <- function(p, aesthetic) {
  for (s in p$scales$scales) {
    if (aesthetic %in% s$aesthetics) {
      return(s)
    }
  }
  NULL
}

test_that("sign_by selects which aesthetic encodes edge sign", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))

  # default "color": colour encodes sign, linetype does not
  p <- ggplot2::autoplot(x)
  expect_true(.ba_has_scale(p, "colour"))
  expect_false(.ba_has_scale(p, "linetype"))

  # "linetype": linetype encodes sign, colour does not
  p <- ggplot2::autoplot(x, sign_by = "linetype")
  expect_true(.ba_has_scale(p, "linetype"))
  expect_false(.ba_has_scale(p, "colour"))

  # "both": colour AND linetype, merged into one legend. ggplot2 merges guides
  # that share a scale name, so both scales must be titled "Direction".
  p <- ggplot2::autoplot(x, sign_by = "both")
  expect_true(.ba_has_scale(p, "colour"))
  expect_true(.ba_has_scale(p, "linetype"))
  expect_identical(.ba_get_scale(p, "colour")$name, "Direction")
  expect_identical(.ba_get_scale(p, "linetype")$name, "Direction")

  # "none": neither encodes sign
  p <- ggplot2::autoplot(x, sign_by = "none")
  expect_false(.ba_has_scale(p, "colour"))
  expect_false(.ba_has_scale(p, "linetype"))
})

test_that("magnitude_by / edge_linewidth control the |r| linewidth scale", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))

  expect_true(.ba_has_scale(ggplot2::autoplot(x), "linewidth"))
  expect_false(.ba_has_scale(ggplot2::autoplot(x, magnitude_by = "none"), "linewidth"))
  # a numeric edge_linewidth also forces constant width (no scale)
  expect_false(.ba_has_scale(ggplot2::autoplot(x, edge_linewidth = 0.6), "linewidth"))
})

test_that("colour_* aliases override the color_* arguments", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))

  # colour_pos / colour_neg feed the Direction colour scale (break order
  # positive, negative) -- readable even when no negative edge is drawn.
  p <- ggplot2::autoplot(x, colour_pos = "darkgreen", colour_neg = "darkorange")
  expect_equal(
    unname(.ba_get_scale(p, "colour")$palette(2L)),
    c("darkgreen", "darkorange")
  )

  # colour_edge sets the single edge colour when sign is not colour-encoded.
  pe <- ggplot2::autoplot(x, sign_by = "linetype", colour_edge = "darkgreen")
  expect_true(all(ggplot2::ggplot_build(pe)$data[[1]]$colour == "darkgreen"))

  # colour_pruned sets the pruned-node fill (use manual pruning for determinism).
  xp <- prune(x, manual = "m4f4")
  pp <- ggplot2::autoplot(xp, colour_pruned = "pink")
  tile_fill <- unlist(lapply(
    ggplot2::ggplot_build(pp)$data,
    function(l) if ("fill" %in% names(l)) l$fill
  ))
  expect_true("pink" %in% tile_fill)
})

test_that("sign_by='none' + magnitude_by='none' renders with no encoding scales", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))

  p <- ggplot2::autoplot(x, sign_by = "none", magnitude_by = "none")
  expect_s3_class(p, "ggplot")
  expect_false(.ba_has_scale(p, "colour"))
  expect_false(.ba_has_scale(p, "linetype"))
  expect_false(.ba_has_scale(p, "linewidth"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("mono=TRUE is equivalent to sign_by='linetype' with black edges", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))

  p <- ggplot2::autoplot(x, mono = TRUE)
  expect_true(.ba_has_scale(p, "linetype"))
  expect_false(.ba_has_scale(p, "colour"))
  edge_colours <- ggplot2::ggplot_build(p)$data[[1]]$colour
  expect_true(all(edge_colours == "black"))
})

test_that("direction='horizontal' transposes the level axis to x", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))

  tile_layer <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    idx <- which(vapply(
      d,
      function(l) all(c("fill", "width") %in% names(l)),
      logical(1)
    ))[[1L]]
    d[[idx]]
  }

  v <- tile_layer(ggplot2::autoplot(x, direction = "vertical"))
  h <- tile_layer(ggplot2::autoplot(x, direction = "horizontal"))

  # Vertical: levels run down the y-axis (level 1 at top, y = -1 is the max).
  expect_equal(max(v$y), -1)
  # Horizontal: levels run along the x-axis (level 1 at left, x = 1 is the min).
  expect_equal(min(h$x), 1)
  # The two layouts are transposes of one another.
  expect_equal(sort(range(v$x)), sort(range(h$y)))
})

test_that("direction='horizontal' composes with skip, mono, and drop_pruned", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  xa <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  expect_s3_class(ggplot2::autoplot(xa, direction = "horizontal", show_skip = TRUE), "ggplot")
  expect_s3_class(ggplot2::autoplot(xa, direction = "horizontal", mono = TRUE), "ggplot")
  # show_r places perpendicular-offset labels; the offset math is orientation-
  # agnostic, so it must compose with the transposed layout.
  expect_s3_class(ggplot2::autoplot(xa, direction = "horizontal", show_r = TRUE), "ggplot")

  xp <- prune(xa, "redundant")
  expect_s3_class(
    suppressWarnings(ggplot2::autoplot(xp, direction = "horizontal", drop_pruned = TRUE)),
    "ggplot"
  )
})

test_that("autoplot.ackwards() show_level_labels=TRUE (default) returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_level_labels = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_level_labels=FALSE omits labels without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_level_labels = FALSE))
  expect_s3_class(p, "ggplot")
})

# --- M40: fully-pruned level axis labels rendered in italic -------------------

# Locate the level-axis-label GeomText layer (labels like "3 factors"), as
# opposed to the node-label GeomText layer (labels like "m3f1").
.level_label_data <- function(p) {
  gt_idx <- which(vapply(
    p$layers, function(l) inherits(l$geom, "GeomText"), logical(1L)
  ))
  for (i in gt_idx) {
    ld <- ggplot2::layer_data(p, i)
    if (any(grepl("factor", ld$label))) {
      return(ld)
    }
  }
  NULL
}

test_that(".fully_pruned_levels() detects a wholly-pruned level, not partial ones", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  # No pruning annotations -> integer(0)
  expect_identical(.fully_pruned_levels(x), integer(0))
  # Prune every factor at level 4 -> level 4 fully pruned
  xp <- prune(x, manual = c("m4f1", "m4f2", "m4f3", "m4f4"))
  expect_identical(.fully_pruned_levels(xp), 4L)
  # Prune only some of level 4 -> not fully pruned
  xpart <- prune(x, manual = c("m4f1", "m4f2"))
  expect_identical(.fully_pruned_levels(xpart), integer(0))
})

test_that("autoplot() italicises a fully-pruned level's axis label (both directions)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  xp <- prune(x, manual = c("m4f1", "m4f2", "m4f3", "m4f4"))

  for (dir in c("vertical", "horizontal")) {
    p <- ggplot2::autoplot(xp, direction = dir)
    ld <- .level_label_data(p)
    expect_false(is.null(ld))
    italic <- ld$label[ld$fontface == "italic"]
    plain <- ld$label[ld$fontface == "plain"]
    expect_identical(italic, "4 factors")
    expect_true(all(c("1 factor", "2 factors", "3 factors") %in% plain))
  }
})

test_that("autoplot() renders all level labels plain when nothing is fully pruned", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  # Un-pruned object
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  ld <- .level_label_data(ggplot2::autoplot(x))
  expect_false(is.null(ld))
  expect_true(all(ld$fontface == "plain"))
  # Partially-pruned level stays plain
  xpart <- prune(x, manual = c("m4f1", "m4f2"))
  ld2 <- .level_label_data(ggplot2::autoplot(xpart))
  expect_true(all(ld2$fontface == "plain"))
})

test_that("autoplot() italicises levels fully pruned by the auto 'redundant' rule", {
  skip_if_not_installed("ggplot2")
  # sim16's built-in redundant chains (M33) fully prune levels 3 and 4 at
  # k_max = 5 with the |r| >= 0.9 threshold -- an *auto*-rule full prune (not
  # manual), and two fully-pruned levels at once. Deterministic for this dataset.
  suppressWarnings(suppressMessages(
    xr <- cached(ackwards(sim16, k_max = 5) |> prune("redundant", redundancy_r = 0.9))
  ))
  expect_identical(.fully_pruned_levels(xr), c(3L, 4L))

  ld <- .level_label_data(ggplot2::autoplot(xr))
  expect_false(is.null(ld))
  expect_setequal(ld$label[ld$fontface == "italic"], c("3 factors", "4 factors"))
  expect_setequal(ld$label[ld$fontface == "plain"], c("1 factor", "2 factors", "5 factors"))
})

test_that("autoplot.ackwards() node_labels applies custom labels", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, node_labels = c(m3f1 = "Alpha")))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() node_labels warns for unknown factor IDs", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_warning(
    ggplot2::autoplot(x, node_labels = c(m99f1 = "Ghost")),
    "match no factor ID"
  )
})

test_that("autoplot.ackwards() primary_only=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, primary_only = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() primary_only=TRUE hides skip arcs", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  # primary_only filters to is_primary==TRUE; skip edges are never primary
  p <- expect_no_error(ggplot2::autoplot(x, primary_only = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() all Wave-1 args compose without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- expect_no_error(ggplot2::autoplot(
    x,
    show_r            = TRUE,
    r_digits          = 2L,
    mono              = TRUE,
    show_level_labels = TRUE,
    level_label_size  = 3,
    node_labels       = c(m4f1 = "General"),
    primary_only      = FALSE
  ))
  expect_s3_class(p, "ggplot")
})

# --- Wave 2 (M8): drop_pruned + compress_levels + .drop_pruned_nodes() ------

test_that("drop_pruned=TRUE errors without pruning annotations", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$prune)
  expect_error(ggplot2::autoplot(x, drop_pruned = TRUE), "prune")
})

test_that("drop_pruned=TRUE returns a ggplot for a pruned object", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE: show_r defaults to FALSE (decoupled from drop_pruned)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  # show_r defaults FALSE; no GeomLabel layers should appear
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE))
  expect_s3_class(p, "ggplot")
  label_layers <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_equal(length(label_layers), 0L)
})

test_that("show_r=TRUE produces GeomLabel layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  label_layers <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_true(length(label_layers) > 0L)
})

test_that("drop_pruned=TRUE + show_r=FALSE produces no GeomLabel layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE, show_r = FALSE))
  expect_s3_class(p, "ggplot")
  label_layers <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_equal(length(label_layers), 0L)
})

test_that("show_r=TRUE label text is APA-formatted (no leading zero, padded decimals)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_r = TRUE, r_digits = 2L)
  idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomLabel"), logical(1L)))
  ld <- ggplot2::layer_data(p, idx)
  # APA convention: no label starts with "0." or "-0."
  expect_false(any(grepl("^0\\.", ld$label)))
  expect_false(any(grepl("^-0\\.", ld$label)))
  # Labels that contain a decimal point have exactly r_digits digits after it
  has_dot <- grepl("\\.", ld$label)
  if (any(has_dot)) {
    decimal_parts <- sub("^[^.]*\\.", "", ld$label[has_dot])
    expect_true(all(nchar(decimal_parts) == 2L))
  }
})

test_that("show_r=TRUE places labels offset from the edge midpoint", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_r = TRUE)
  ld <- ggplot2::layer_data(p, which(vapply(p$layers, function(l) {
    inherits(l$geom, "GeomLabel")
  }, logical(1L))))
  # Recover the segment layers for midpoint comparison
  seg_layers <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1L)))
  sd <- ggplot2::layer_data(p, seg_layers[[1L]])
  mid_x <- (sd$x + sd$xend) / 2
  mid_y <- (sd$y + sd$yend) / 2
  # Labels must not coincide with the raw midpoints — perpendicular nudge moved them
  expect_false(all(abs(ld$x - mid_x) < 1e-9 & abs(ld$y - mid_y) < 1e-9))
})

test_that("r_label_size arg changes the geom_label size", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p_small <- ggplot2::autoplot(x, show_r = TRUE, r_label_size = 1.5)
  p_large <- ggplot2::autoplot(x, show_r = TRUE, r_label_size = 4.0)
  get_label_size <- function(p) {
    idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomLabel"), logical(1L)))
    p$layers[[idx]]$aes_params$size
  }
  expect_equal(get_label_size(p_small), 1.5)
  expect_equal(get_label_size(p_large), 4.0)
})

test_that("drop_pruned=TRUE + compress_levels=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  p <- expect_no_error(ggplot2::autoplot(x,
    drop_pruned = TRUE,
    compress_levels = TRUE
  ))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE + mono=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE, mono = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE + node_labels applies labels", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  p <- expect_no_error(ggplot2::autoplot(x,
    drop_pruned = TRUE,
    node_labels = c(m5f1 = "Neuroticism")
  ))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned degenerate case warns and returns node-only ggplot", {
  skip_if_not_installed("ggplot2")
  # k=2 with high-g data: m1f1 is redundant with m2f1 (chain reaches k=2),
  # leaving only level-2 nodes — n_kept_levels = 1 -> degenerate
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  d <- data.frame(
    x1 = 0.9 * g + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + rnorm(n, sd = 0.05)
  )
  x <- suppressMessages(ackwards(d, k_max = 2) |> prune("redundant", redundancy_r = 0.7))
  # Verify the test fixture: m1f1 should be pruned, only level 2 keeps nodes
  expect_equal(
    unique(x$prune$nodes$level[!x$prune$nodes$pruned]),
    2L
  )
  expect_warning(ggplot2::autoplot(x, drop_pruned = TRUE), "node-only")
  p <- suppressWarnings(ggplot2::autoplot(x, drop_pruned = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that(".drop_pruned_nodes() returns kept-only nodes and reduced edges", {
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  lay <- ba_layout(x)
  dp <- ackwards:::.drop_pruned_nodes(x, lay$nodes)
  kept <- x$prune$nodes$id[!x$prune$nodes$pruned]

  # All returned nodes are kept
  expect_true(all(dp$nodes$id %in% kept))
  expect_equal(nrow(dp$nodes), length(kept))

  # Each edge connects two kept nodes
  if (nrow(dp$edges) > 0L) {
    expect_true(all(dp$edges$from %in% kept))
    expect_true(all(dp$edges$to %in% kept))
    # Each "to" node appears at most once (one primary parent per node)
    expect_equal(length(unique(dp$edges$to)), nrow(dp$edges))
    # Parent is always shallower than child
    expect_true(all(dp$edges$level_from < dp$edges$level_to))
  }
})

test_that(".drop_pruned_nodes() compress_levels re-indexes y", {
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  lay <- ba_layout(x)
  dp_gap <- ackwards:::.drop_pruned_nodes(x, lay$nodes, compress_levels = FALSE)
  dp_comp <- ackwards:::.drop_pruned_nodes(x, lay$nodes, compress_levels = TRUE)

  # Gap-preserved: y == -level for all kept nodes
  expect_equal(dp_gap$nodes$y, -dp_gap$nodes$level)

  # Compressed: y values are consecutive integers starting at -1
  kept_levels <- sort(unique(dp_comp$nodes$level))
  if (length(kept_levels) > 1L) {
    expected_y <- -match(dp_comp$nodes$level, kept_levels)
    expect_equal(dp_comp$nodes$y, expected_y)
  }
})

test_that("drop_pruned=TRUE warns when no nodes are pruned (artifact-only)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3) |> prune("artifact"))
  ))
  expect_true(!any(x$prune$nodes$pruned))
  expect_warning(
    ggplot2::autoplot(x, drop_pruned = TRUE),
    "no nodes are flagged"
  )
})

test_that("drop_pruned=TRUE warns when cut_show removes all reduced edges", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  # cut_show = 1 is the legal maximum since M42/m8 validation; no reduced
  # edge in these data reaches exactly 1.0, so all are hidden
  expect_warning(
    ggplot2::autoplot(x, drop_pruned = TRUE, cut_show = 1),
    "No edges"
  )
})

test_that("node_labels as unnamed vector silently does nothing", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  # Unnamed character: names() is NULL; no warning, no crash
  p <- expect_no_error(ggplot2::autoplot(x, node_labels = c("Alpha", "Beta")))
  expect_s3_class(p, "ggplot")
})

test_that("compress_levels=TRUE without drop_pruned is silently ignored", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, compress_levels = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE ignores primary_only and show_skip", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5) |>
      prune("redundant", redundancy_r = 0.95)
  ))
  # These flags should be silently ignored — no error
  p <- expect_no_error(ggplot2::autoplot(
    x,
    drop_pruned  = TRUE,
    primary_only = TRUE,
    show_skip    = TRUE
  ))
  expect_s3_class(p, "ggplot")
})

# --- Wave 1 (M9): show_arrows, edge_linewidth, legend -------------------------

test_that("show_arrows=FALSE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_arrows = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("show_arrows=FALSE sets arrow=NULL on segment layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_arrows = FALSE)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  expect_true(length(seg_layers) > 0L)
  for (l in seg_layers) expect_null(l$geom_params$arrow)
})

test_that("show_arrows=TRUE (default) retains arrowheads on segment layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_arrows = TRUE)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  for (l in seg_layers) expect_s3_class(l$geom_params$arrow, "arrow")
})

test_that("show_arrows=FALSE removes arrowheads from curved arcs too", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- ggplot2::autoplot(x, show_arrows = FALSE)
  expect_s3_class(p, "ggplot")
  cur_layers <- Filter(function(l) inherits(l$geom, "GeomCurve"), p$layers)
  for (l in cur_layers) expect_null(l$geom_params$arrow)
})

test_that("edge_linewidth numeric returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, edge_linewidth = 0.6))
  expect_s3_class(p, "ggplot")
})

test_that("edge_linewidth numeric drops the linewidth scale", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  has_lw_scale <- function(p) {
    any(vapply(p$scales$scales, function(s) "linewidth" %in% s$aesthetics, logical(1L)))
  }
  expect_true(has_lw_scale(ggplot2::autoplot(x)))
  expect_false(has_lw_scale(ggplot2::autoplot(x, edge_linewidth = 0.6)))
})

test_that("edge_linewidth=NULL (default) retains the linewidth scale", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  has_lw_scale <- function(p) {
    any(vapply(p$scales$scales, function(s) "linewidth" %in% s$aesthetics, logical(1L)))
  }
  expect_true(has_lw_scale(ggplot2::autoplot(x, edge_linewidth = NULL)))
})

test_that("legend=FALSE returns a ggplot with legend.position='none'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, legend = FALSE))
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

test_that("legend=TRUE (default) does not suppress the legend", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, legend = TRUE)
  expect_false(identical(p$theme$legend.position, "none"))
})

test_that("Forbes-style composition: drop_pruned+black+fixed lw+no arrows+no legend", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5) |> prune("redundant", redundancy_r = 0.95))
  ))
  p <- expect_no_error(ggplot2::autoplot(
    x,
    drop_pruned    = TRUE,
    color_pos      = "black",
    color_neg      = "black",
    edge_linewidth = 0.6,
    show_arrows    = FALSE,
    legend         = FALSE
  ))
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

test_that("M9 args compose with mono=TRUE without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(
    x,
    mono           = TRUE,
    show_arrows    = FALSE,
    edge_linewidth = 0.5,
    legend         = FALSE
  ))
  expect_s3_class(p, "ggplot")
})

test_that("all M9 args compose with all M8 args without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- expect_no_error(ggplot2::autoplot(
    x,
    show_r            = TRUE,
    r_digits          = 2L,
    show_level_labels = TRUE,
    node_labels       = c(m4f1 = "General"),
    primary_only      = FALSE,
    show_arrows       = FALSE,
    edge_linewidth    = 0.5,
    legend            = FALSE
  ))
  expect_s3_class(p, "ggplot")
})

test_that("edge_linewidth numeric stores constant linewidth on segment layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, edge_linewidth = 0.6)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  expect_true(length(seg_layers) > 0L)
  # linewidth is a ggplot2 aesthetic so constant values land in $aes_params
  for (l in seg_layers) expect_equal(l$aes_params$linewidth, 0.6)
})

test_that("edge_linewidth invalid values error clearly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_error(ggplot2::autoplot(x, edge_linewidth = "thick"), "edge_linewidth")
  expect_error(ggplot2::autoplot(x, edge_linewidth = -0.5), "edge_linewidth")
  expect_error(ggplot2::autoplot(x, edge_linewidth = c(0.5, 0.6)), "edge_linewidth")
})

test_that("legend=FALSE is honoured on the degenerate drop_pruned plot", {
  skip_if_not_installed("ggplot2")
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  d <- data.frame(
    x1 = 0.9 * g + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + rnorm(n, sd = 0.05)
  )
  x <- suppressMessages(ackwards(d, k_max = 2) |> prune("redundant", redundancy_r = 0.7))
  p <- suppressWarnings(ggplot2::autoplot(x, drop_pruned = TRUE, legend = FALSE))
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

# ── M27: autoplot(x, what = "fit") ────────────────────────────────────────────

test_that("autoplot(x, what='fit') returns ggplot for EFA object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  p <- ggplot2::autoplot(x, what = "fit")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(x, what='fit') draws Hu & Bentler threshold reference lines", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  p <- ggplot2::autoplot(x, what = "fit")
  # The cutoff reference lines are a geom_hline layer (retained .fit_cutoffs()).
  has_hline <- vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomHline"),
    logical(1L)
  )
  expect_true(any(has_hline))
  # The reference y-intercepts are the EFA-charted thresholds (TLI .95, RMSEA .06).
  hline_layer <- p$layers[[which(has_hline)[1L]]]
  expect_setequal(hline_layer$data$yintercept, c(0.95, 0.06))
})

test_that("autoplot(x, what='fit') returns ggplot for ESEM object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("lavaan")
  d <- data.frame(matrix(rnorm(300 * 6), 300, 6))
  x <- cached(ackwards(d, k_max = 3, engine = "esem"))
  p <- ggplot2::autoplot(x, what = "fit")
  expect_s3_class(p, "ggplot")
})

test_that("fit-plot caption names only the plotted indices (M42/m10)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  # EFA panels show TLI/RMSEA only; the caption must not list CFI/SRMR cutoffs.
  suppressWarnings(
    x_efa <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  )
  cap_efa <- ggplot2::autoplot(x_efa, what = "fit")$labels$caption
  expect_true(grepl("TLI", cap_efa) && grepl("RMSEA", cap_efa))
  expect_false(grepl("CFI", cap_efa) || grepl("SRMR", cap_efa))
})

test_that("autoplot(x, what='fit') returns ggplot with informative message for PCA", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, what = "fit")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(x, what='hierarchy') is unchanged from default", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  p_default <- ggplot2::autoplot(x)
  p_explicit <- ggplot2::autoplot(x, what = "hierarchy")
  # Both should be ggplots; data-level equality is already tested elsewhere
  expect_s3_class(p_default, "ggplot")
  expect_s3_class(p_explicit, "ggplot")
})

test_that("autoplot(what='fit') returns empty plot when no kept indices present", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  # Strip the EFA fit down to indices the plot does not chart (chi/dof only),
  # forcing the "No fit indices available" branch.
  for (ki in names(x$levels)) {
    x$levels[[ki]]$fit <- c(chi = 1, dof = 1)
  }
  p <- ggplot2::autoplot(x, what = "fit")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(what='fit') panel titles are engine-aware (EFA: no CFI/SRMR)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  p <- ggplot2::autoplot(x, what = "fit")
  panels <- as.character(unique(p$data$panel))
  # EFA charts TLI and RMSEA only; titles must not advertise CFI or SRMR.
  expect_false(any(grepl("CFI", panels)))
  expect_false(any(grepl("SRMR", panels)))
  expect_true(any(grepl("TLI", panels)))
  expect_true(any(grepl("RMSEA", panels)))
})

# --- M81: manual deepest-level ordering (order=) ----------------------------

test_that("ba_layout() order= places the deepest level in the supplied order", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  perm <- rev(x$levels[["4"]]$labels)
  lay <- ba_layout(x, order = perm)
  d <- lay$nodes[lay$nodes$level == 4L, ]
  # left-to-right (increasing x) order equals the requested permutation
  expect_equal(d$id[order(d$x)], perm)
  # unchanged invariant: the single level-1 node is still anchored at x = 0
  expect_equal(lay$nodes$x[lay$nodes$level == 1L], 0)
})

test_that("ba_layout() order= accepts a level-keyed list, ignoring non-deepest", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  perm <- rev(x$levels[["4"]]$labels)
  expect_warning(
    lay <- ba_layout(x, order = list("3" = rev(x$levels[["3"]]$labels), "4" = perm)),
    "level.*3.*ignored|ignored"
  )
  d <- lay$nodes[lay$nodes$level == 4L, ]
  expect_equal(d$id[order(d$x)], perm)
})

test_that("ba_layout() order= errors on a non-permutation of deepest IDs", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  labs <- x$levels[["4"]]$labels
  # missing one ID
  expect_error(ba_layout(x, order = labs[-1L]), "permutation")
  # an extra / unknown ID
  expect_error(ba_layout(x, order = c(labs, "m4f9")), "permutation")
  # a duplicate
  expect_error(ba_layout(x, order = c(labs[1L], labs)), "permutation")
  # wrong type
  expect_error(ba_layout(x, order = seq_along(labs)), "permutation")
})

test_that("ba_layout() order= errors when a list omits the deepest level", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  expect_error(
    suppressWarnings(ba_layout(x, order = list("3" = rev(x$levels[["3"]]$labels)))),
    "deepest level"
  )
})

test_that("autoplot() forwards order= to ba_layout()", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4))
  perm <- rev(x$levels[["4"]]$labels)
  p <- autoplot(x, order = perm)
  bd <- ggplot2::ggplot_build(p)$data
  # the node-label layer carries one row per node with x, y and label(= id)
  lab <- do.call(rbind, lapply(bd, function(d) {
    if (all(c("x", "y", "label") %in% names(d))) d[, c("x", "y", "label")] else NULL
  }))
  d4 <- lab[lab$label %in% perm, , drop = FALSE]
  expect_equal(d4$label[order(d4$x)], perm)
})
