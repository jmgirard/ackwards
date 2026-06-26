test_that("algebra and scores paths agree for PCA engine", {
  skip_if_not_installed("psych")

  # Small dataset with known structure (no external dependency needed)
  set.seed(42)
  n <- 150
  f1 <- rnorm(n); f2 <- rnorm(n)
  data <- data.frame(
    x1 = f1 + 0.3 * rnorm(n),
    x2 = f1 + 0.3 * rnorm(n),
    x3 = f1 + 0.3 * rnorm(n),
    x4 = f2 + 0.3 * rnorm(n),
    x5 = f2 + 0.3 * rnorm(n),
    x6 = f2 + 0.3 * rnorm(n)
  )
  R <- cor(data)

  suppressWarnings(x <- ackwards(data, k = 3))

  E_algebra <- x$edges$matrices[["1:2"]]

  E_scores <- compute_edges(
    levels  = x$levels,
    R       = R,
    method  = "scores",
    pairs   = "adjacent",
    data    = data,
    align   = FALSE
  )$matrices[["1:2"]]

  expect_lt(max(abs(abs(E_algebra) - abs(E_scores))), 1e-6)
})

test_that("compute_edges errors when algebra forced but R missing", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 2))
  expect_error(
    compute_edges(x$levels, R = NULL, method = "algebra"),
    "conditions not met"
  )
})

test_that("compute_edges errors when scores needed but data absent", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 2))
  expect_error(
    compute_edges(x$levels, R = NULL, method = "scores", data = NULL),
    "data"
  )
})
