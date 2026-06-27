test_that("ackwards() returns a valid ackwards object for the smoke-test case", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 5))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)  # internal: checks all required fields

  expect_equal(x$method,   "pca")
  expect_equal(x$k_max,    5L)
  expect_equal(x$n_obs,    2800L)
  expect_equal(x$cor_type, "pearson")

  # All five levels present and converged
  expect_equal(length(x$levels), 5L)
  expect_true(all(vapply(x$levels, `[[`, logical(1), "converged")))

  # Four adjacent edge matrices (1:2, 2:3, 3:4, 4:5)
  expect_equal(length(x$edges$matrices), 4L)
  expect_named(x$edges$matrices, c("1:2", "2:3", "3:4", "4:5"))

  # Edge matrix dimensions
  for (i in 1:4) {
    E <- x$edges$matrices[[paste0(i, ":", i + 1)]]
    expect_equal(nrow(E), i,     info = paste("rows of", i, ":", i + 1))
    expect_equal(ncol(E), i + 1, info = paste("cols of", i, ":", i + 1))
    expect_true(all(abs(E) <= 1 + 1e-9), info = "correlations in [-1, 1]")
  }
})

test_that("PCA edge correlations match psych::bassAckward within tolerance", {
  skip_if_not_installed("psych")

  ba  <- psych::bassAckward(
    psych::bfi[, 1:25], nfactors = 5, fm = "pca", rotate = "varimax", plot = FALSE
  )
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 5))

  for (i in 1:4) {
    # psych stores (k+1) × k; ours is k × (k+1) — transpose to align
    psych_mat <- t(ba$bass.ack[[i + 1L]])
    our_mat   <- x$edges$matrices[[paste0(i, ":", i + 1)]]

    max_diff <- max(abs(abs(psych_mat) - abs(our_mat)))
    expect_lt(max_diff, 1e-4, label = paste("Level", i, "->", i + 1))
  }
})

test_that("levels have correct structure and label formats", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))

  for (ki in 1:3) {
    lev <- x$levels[[as.character(ki)]]
    expect_equal(ncol(lev$loadings), ki)
    expect_equal(nrow(lev$loadings), 25L)
    expect_equal(lev$labels, paste0("m", ki, "f", seq_len(ki)))
    expect_equal(colnames(lev$loadings), lev$labels)
    expect_true(isTRUE(lev$converged))
    expect_equal(dim(lev$scoring$weights), c(25L, ki))
  }
})

test_that("ackwards() errors informatively on bad inputs", {
  skip_if_not_installed("psych")
  d <- psych::bfi[, 1:5]
  expect_error(ackwards(d, k = 0),     "integer >= 2")
  expect_error(ackwards(d, k = 1),     "integer >= 2")
  expect_error(ackwards(d, k = 100),   "cannot exceed")
  expect_error(ackwards(d, k = 1.5),   "integer >= 2")
  expect_error(ackwards(list(), k = 2), "data frame")
})

test_that("rotation = 'cfQ' errors with a clear message for PCA engine", {
  skip_if_not_installed("psych")
  d <- psych::bfi[, 1:5]
  expect_error(ackwards(d, k = 2, rotation = "cfQ"), "not yet implemented")
})

test_that("scores = TRUE warns that storage is not yet implemented", {
  skip_if_not_installed("psych")
  expect_warning(
    ackwards(psych::bfi[, 1:25], k = 2, scores = TRUE),
    "not yet implemented"
  )
})

test_that("keep_fits = TRUE warns that storage is not yet implemented", {
  skip_if_not_installed("psych")
  expect_warning(
    ackwards(psych::bfi[, 1:25], k = 2, keep_fits = TRUE),
    "not yet implemented"
  )
})

test_that("meta$cut_show stores the cut_show value", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 2, cut_show = 0.4))
  expect_equal(x$meta$cut_show, 0.4)
})

test_that("detect_ordinal() flags bfi columns", {
  skip_if_not_installed("psych")
  expect_true(ackwards:::detect_ordinal(psych::bfi[, 1:25]))
})

test_that("detect_ordinal() returns FALSE for continuous data", {
  set.seed(1)
  d <- data.frame(matrix(rnorm(500), 100, 5))
  expect_false(ackwards:::detect_ordinal(d))
})

test_that("meta$ordinal_warned is TRUE for bfi data", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 2))
  expect_true(x$meta$ordinal_warned)
})
