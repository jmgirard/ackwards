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
                 info = paste("level", k, "has", k, "nodes"))
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
    lay   <- ba_layout(x, min_sep = sep)
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

test_that("autoplot.ackwards() returns a ggplot object", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  p <- ggplot2::autoplot(x)
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.ackwards() respects cut_show and cut_strong", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  # Should run without error for various threshold combos
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.1, cut_strong = 0.4))
  expect_no_error(ggplot2::autoplot(x, cut_show = 0.5, cut_strong = 0.7))
})

test_that("plot.ackwards() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_no_error(plot(x))
  expect_invisible(plot(x))
})
