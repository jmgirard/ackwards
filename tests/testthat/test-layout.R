test_that("ba_layout() returns a list with nodes and edges data frames", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 4))
  lay <- ba_layout(x)

  expect_type(lay, "list")
  expect_named(lay, c("nodes", "edges"))
  expect_s3_class(lay$nodes, "data.frame")
  expect_s3_class(lay$edges, "data.frame")
})

test_that("ba_layout() nodes have correct structure and counts", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 4))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$y, -nodes$level)
})

test_that("ba_layout() level-1 node is at x = 0", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  lay <- ba_layout(x)
  nodes <- lay$nodes
  expect_equal(nodes$x[nodes$level == 1], 0)
})

test_that("ba_layout() respects min_sep between nodes at same level", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 5))

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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 5))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  p <- ggplot2::autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot(x) dispatches correctly without library(ggplot2)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  # The package defines its own autoplot generic so this works without ggplot2 attached
  p <- autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() respects cut_show and cut_strong", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.1, cut_strong = 0.4))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.5, cut_strong = 0.7))
})

test_that("autoplot.ackwards() warns when cut_show hides all edges", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  # cut_show > 1.0 always hides all edges (correlations are bounded by [-1, 1])
  expect_warning(ggplot2::autoplot(x, cut_show = 1.01), "No edges")
})

test_that("autoplot.ackwards() warns when min_sep < node_width", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_warning(ggplot2::autoplot(x, min_sep = 0.3, node_width = 0.8), "overlap")
})

test_that("plot.ackwards() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_no_error(plot(x))
  expect_invisible(plot(x))
})

test_that("autoplot.ackwards() renders skip-level edges with pairs='all'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 4, pairs = "all"))
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
    ackwards(data, k = 4, prune = "redundant", redundancy_r = 0.9)
  ))
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
  # Custom prune colour
  expect_no_error(ggplot2::autoplot(x, color_pruned = "pink"))
})

test_that("autoplot.ackwards() handles objects with prune=NULL (no pruning)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_null(x$prune)
  p <- expect_no_error(ggplot2::autoplot(x))
  expect_s3_class(p, "ggplot")
})
