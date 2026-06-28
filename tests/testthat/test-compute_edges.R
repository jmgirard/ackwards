test_that("algebra and scores paths agree for PCA engine on all adjacent pairs", {
  skip_if_not_installed("psych")

  set.seed(42)
  n <- 150
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  data <- data.frame(
    x1 = f1 + 0.3 * rnorm(n),
    x2 = f1 + 0.3 * rnorm(n),
    x3 = f1 + 0.3 * rnorm(n),
    x4 = f2 + 0.3 * rnorm(n),
    x5 = f2 + 0.3 * rnorm(n),
    x6 = f2 + 0.3 * rnorm(n)
  )
  R <- cor(data)

  suppressWarnings(x <- ackwards(data, k_max = 4))

  E_scores_all <- compute_edges(
    levels = x$levels,
    R = R,
    edge_method = "scores",
    pairs = "adjacent",
    data = data,
    align = FALSE
  )$matrices

  for (key in names(x$edges$matrices)) {
    E_alg <- x$edges$matrices[[key]]
    E_sc <- E_scores_all[[key]]
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 1e-6,
      label = paste("algebra vs scores for pair", key)
    )
  }
})

test_that("compute_edges errors when algebra forced but R missing", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(
    compute_edges(x$levels, R = NULL, edge_method = "algebra"),
    "conditions not met"
  )
})

test_that("compute_edges errors when scores needed but data absent", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(
    compute_edges(x$levels, R = NULL, edge_method = "scores", data = NULL),
    "data"
  )
})

test_that("pairs='all' produces skip-level edge matrices in ackwards()", {
  skip_if_not_installed("psych")

  set.seed(1)
  n <- 200
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = g + s1 + rnorm(n, sd = 0.3),
    x2 = g + s1 + rnorm(n, sd = 0.3),
    x3 = g + s1 + rnorm(n, sd = 0.3),
    x4 = g + s2 + rnorm(n, sd = 0.3),
    x5 = g + s2 + rnorm(n, sd = 0.3),
    x6 = g + s2 + rnorm(n, sd = 0.3)
  )

  x_adj <- suppressWarnings(ackwards(data, k_max = 4, pairs = "adjacent"))
  x_all <- suppressWarnings(ackwards(data, k_max = 4, pairs = "all"))

  # Default is adjacent
  expect_equal(x_adj$meta$pairs, "adjacent")
  expect_equal(x_all$meta$pairs, "all")

  # adjacent: only consecutive level keys
  adj_keys <- names(x_adj$edges$matrices)
  expect_setequal(adj_keys, c("1:2", "2:3", "3:4"))

  # all: includes skip-level keys
  all_keys <- names(x_all$edges$matrices)
  expect_true(all(c("1:2", "2:3", "3:4", "1:3", "1:4", "2:4") %in% all_keys))

  # Adjacent edges are consistent between the two modes (same aligned weights)
  for (key in c("1:2", "2:3", "3:4")) {
    expect_equal(
      x_adj$edges$matrices[[key]],
      x_all$edges$matrices[[key]],
      tolerance = 1e-10,
      label = paste("adjacent edge", key, "consistent between pairs modes")
    )
  }
})

test_that("skip-level edges in pairs='all' have is_primary=FALSE and correct level_from/to", {
  skip_if_not_installed("psych")

  set.seed(2)
  n <- 200
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = g + s1 + rnorm(n, sd = 0.3),
    x2 = g + s1 + rnorm(n, sd = 0.3),
    x3 = g + s1 + rnorm(n, sd = 0.3),
    x4 = g + s2 + rnorm(n, sd = 0.3),
    x5 = g + s2 + rnorm(n, sd = 0.3),
    x6 = g + s2 + rnorm(n, sd = 0.3)
  )

  x <- suppressWarnings(ackwards(data, k_max = 3, pairs = "all"))
  tidy <- x$edges$tidy

  # Skip-level pair 1:3 should exist
  skip <- tidy[tidy$level_from == 1 & tidy$level_to == 3, ]
  expect_gt(nrow(skip), 0)

  # All skip-level edges must be non-primary
  expect_true(all(!skip$is_primary))

  # level_from < level_to throughout
  expect_true(all(tidy$level_from < tidy$level_to))
})

test_that("algebra and scores agree for all-levels edges (PCA)", {
  skip_if_not_installed("psych")

  set.seed(3)
  n <- 200
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = g + s1 + rnorm(n, sd = 0.3),
    x2 = g + s1 + rnorm(n, sd = 0.3),
    x3 = g + s1 + rnorm(n, sd = 0.3),
    x4 = g + s2 + rnorm(n, sd = 0.3),
    x5 = g + s2 + rnorm(n, sd = 0.3),
    x6 = g + s2 + rnorm(n, sd = 0.3)
  )
  R <- cor(data)

  x <- suppressWarnings(ackwards(data, k_max = 4, pairs = "all"))

  E_scores_all <- compute_edges(
    levels = x$levels,
    R = R,
    edge_method = "scores",
    pairs = "all",
    data = data,
    align = FALSE
  )$matrices

  for (key in names(x$edges$matrices)) {
    E_alg <- x$edges$matrices[[key]]
    E_sc <- E_scores_all[[key]]
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 1e-6,
      label = paste("algebra vs scores for all-pairs edge", key)
    )
  }
})

# --- match_parents unit tests ---------------------------------------------------
# Adjacent levels always have n_b = n_a + 1 (non-square). LSAP (bijection) would
# require a padding row that can return index > nrow(E) → subscript OOB.
# The greedy argmax is correct: each child independently picks its closest parent;
# multiple children may share a parent, which is normal in the hierarchy.

test_that("match_parents: non-square (n_b > n_a) assigns correct greedy parents", {
  # Level k-1 has 2 factors (rows), level k has 3 factors (cols).
  # Col 1 and 2 both peak at row 1; col 3 peaks at row 2.
  E <- matrix(
    c(
      0.9, 0.8, 0.1,
      0.1, 0.2, 0.7
    ),
    nrow = 2
  )
  result <- ackwards:::match_parents(E)
  expect_equal(result, c(1L, 1L, 2L))
  # All returned indices must be within 1..nrow(E)
  expect_true(all(result >= 1L & result <= nrow(E)))
})

test_that("match_parents: square case assigns correct parents", {
  E <- matrix(
    c(
      0.9, 0.1, 0.1,
      0.2, 0.8, 0.2,
      0.1, 0.2, 0.7
    ),
    nrow = 3, byrow = TRUE
  )
  result <- ackwards:::match_parents(E)
  expect_equal(result, c(1L, 2L, 3L))
})

test_that("ackwards() parent indices are always within bounds (regression: LSAP padding)", {
  skip_if_not_installed("psych")
  # This call previously triggered subscript OOB when clue was installed,
  # because LSAP padding rows returned indices > nrow(E) for non-square levels.
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 5))
  expect_s3_class(x, "ackwards")
  expect_equal(x$k_max, 5L)
  # Verify every lineage index is in bounds for its level
  for (ki in seq(2L, x$k_max)) {
    parents <- x$lineage[[as.character(ki)]]
    n_parents <- ki - 1L # nrow of edge matrix for level ki-1 → ki
    expect_true(
      all(parents >= 1L & parents <= n_parents),
      label = paste("lineage indices in bounds at level", ki)
    )
  }
})
