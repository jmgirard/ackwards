test_that("ba_layout() returns a list with nodes and edges data frames", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4))
  lay <- ba_layout(x)

  expect_type(lay, "list")
  expect_named(lay, c("nodes", "edges"))
  expect_s3_class(lay$nodes, "data.frame")
  expect_s3_class(lay$edges, "data.frame")
})

test_that("ba_layout() nodes have correct structure and counts", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$y, -nodes$level)
})

test_that("ba_layout() level-1 node is at x = 0", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$x[nodes$level == 1], 0)
})

test_that("ba_layout() respects min_sep between nodes at same level", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 5))

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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 5))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(x) dispatches correctly without library(ggplot2)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  # The package defines its own autoplot generic so this works without ggplot2 attached
  p <- autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() respects cut_show and cut_strong", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.1, cut_strong = 0.4))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.5, cut_strong = 0.7))
})

test_that("autoplot.ackwards() warns when cut_show hides all edges", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  # cut_show > 1.0 always hides all edges (correlations are bounded by [-1, 1])
  expect_warning(ggplot2::autoplot(x, cut_show = 1.01), "No edges")
})

test_that("autoplot.ackwards() warns when min_sep < node_width", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_warning(ggplot2::autoplot(x, min_sep = 0.3, node_width = 0.8), "overlap")
})

test_that("plot.ackwards() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_no_error(plot(x))
  expect_invisible(plot(x))
})

test_that("autoplot.ackwards() renders skip-level edges with pairs='all'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
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
    ackwards(data, k_max = 4, prune = "redundant", redundancy_r = 0.9)
  ))
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
  # Custom prune colour
  expect_no_error(ggplot2::autoplot(x, color_pruned = "pink"))
})

test_that("autoplot.ackwards() handles objects with prune=NULL (no pruning)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$prune)
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
})

# --- Wave 1 (M8): show_r, mono, show_level_labels, node_labels, primary_only ---

test_that("autoplot.ackwards() show_r=TRUE labels edges without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_r=TRUE respects r_digits", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE, r_digits = 3L))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_r=TRUE works with pairs='all'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() mono=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, mono = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() mono=TRUE + show_r=TRUE composes without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, mono = TRUE, show_r = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_level_labels=TRUE (default) returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_level_labels = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() show_level_labels=FALSE omits labels without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_level_labels = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() node_labels applies custom labels", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, node_labels = c(m3f1 = "Alpha")))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() node_labels warns for unknown factor IDs", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_warning(
    ggplot2::autoplot(x, node_labels = c(m99f1 = "Ghost")),
    "match no factor ID"
  )
})

test_that("autoplot.ackwards() primary_only=TRUE returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, primary_only = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() primary_only=TRUE hides skip arcs", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  # primary_only filters to is_primary==TRUE; skip edges are never primary
  p <- expect_no_error(ggplot2::autoplot(x, primary_only = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() all Wave-1 args compose without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$prune)
  expect_error(ggplot2::autoplot(x, drop_pruned = TRUE), "prune")
})

test_that("drop_pruned=TRUE returns a ggplot for a pruned object", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE: show_r defaults to FALSE (decoupled from drop_pruned)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_r = TRUE))
  label_layers <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_true(length(label_layers) > 0L)
})

test_that("drop_pruned=TRUE + show_r=FALSE produces no GeomLabel layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE, show_r = FALSE))
  expect_s3_class(p, "ggplot")
  label_layers <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_equal(length(label_layers), 0L)
})

test_that("show_r=TRUE label text is APA-formatted (no leading zero, padded decimals)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
  ))
  p <- expect_no_error(ggplot2::autoplot(x, drop_pruned = TRUE, mono = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE + node_labels applies labels", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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
  x <- suppressMessages(ackwards(d,
    k_max = 2, prune = "redundant",
    redundancy_r = 0.7
  ))
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
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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

test_that("drop_pruned=TRUE warns when no nodes are pruned (artefact-only)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 3, prune = "artefact")
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
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
  ))
  expect_warning(
    ggplot2::autoplot(x, drop_pruned = TRUE, cut_show = 1.01),
    "No edges"
  )
})

test_that("node_labels as unnamed vector silently does nothing", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  # Unnamed character: names() is NULL; no warning, no crash
  p <- expect_no_error(ggplot2::autoplot(x, node_labels = c("Alpha", "Beta")))
  expect_s3_class(p, "ggplot")
})

test_that("compress_levels=TRUE without drop_pruned is silently ignored", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, compress_levels = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("drop_pruned=TRUE ignores primary_only and show_skip", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25],
      k_max = 5, prune = "redundant",
      redundancy_r = 0.95
    )
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, show_arrows = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("show_arrows=FALSE sets arrow=NULL on segment layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_arrows = FALSE)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  expect_true(length(seg_layers) > 0L)
  for (l in seg_layers) expect_null(l$geom_params$arrow)
})

test_that("show_arrows=TRUE (default) retains arrowheads on segment layers", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, show_arrows = TRUE)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  for (l in seg_layers) expect_s3_class(l$geom_params$arrow, "arrow")
})

test_that("show_arrows=FALSE removes arrowheads from curved arcs too", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
  p <- ggplot2::autoplot(x, show_arrows = FALSE)
  expect_s3_class(p, "ggplot")
  cur_layers <- Filter(function(l) inherits(l$geom, "GeomCurve"), p$layers)
  for (l in cur_layers) expect_null(l$geom_params$arrow)
})

test_that("edge_linewidth numeric returns a ggplot", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, edge_linewidth = 0.6))
  expect_s3_class(p, "ggplot")
})

test_that("edge_linewidth numeric drops the linewidth scale", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  has_lw_scale <- function(p) {
    any(vapply(p$scales$scales, function(s) "linewidth" %in% s$aesthetics, logical(1L)))
  }
  expect_true(has_lw_scale(ggplot2::autoplot(x)))
  expect_false(has_lw_scale(ggplot2::autoplot(x, edge_linewidth = 0.6)))
})

test_that("edge_linewidth=NULL (default) retains the linewidth scale", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  has_lw_scale <- function(p) {
    any(vapply(p$scales$scales, function(s) "linewidth" %in% s$aesthetics, logical(1L)))
  }
  expect_true(has_lw_scale(ggplot2::autoplot(x, edge_linewidth = NULL)))
})

test_that("legend=FALSE returns a ggplot with legend.position='none'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- expect_no_error(ggplot2::autoplot(x, legend = FALSE))
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})

test_that("legend=TRUE (default) does not suppress the legend", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, legend = TRUE)
  expect_false(identical(p$theme$legend.position, "none"))
})

test_that("Forbes-style composition: drop_pruned+black+fixed lw+no arrows+no legend", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(suppressMessages(
    x <- ackwards(psych::bfi[, 1:25], k_max = 5, prune = "redundant", redundancy_r = 0.95)
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 4, pairs = "all"))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  p <- ggplot2::autoplot(x, edge_linewidth = 0.6)
  seg_layers <- Filter(function(l) inherits(l$geom, "GeomSegment"), p$layers)
  expect_true(length(seg_layers) > 0L)
  # linewidth is a ggplot2 aesthetic so constant values land in $aes_params
  for (l in seg_layers) expect_equal(l$aes_params$linewidth, 0.6)
})

test_that("edge_linewidth invalid values error clearly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
  x <- suppressMessages(ackwards(d, k_max = 2, prune = "redundant", redundancy_r = 0.7))
  p <- suppressWarnings(ggplot2::autoplot(x, drop_pruned = TRUE, legend = FALSE))
  expect_s3_class(p, "ggplot")
  expect_equal(p$theme$legend.position, "none")
})
