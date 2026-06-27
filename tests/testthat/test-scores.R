# tests/testthat/test-scores.R
# M6: scores storage, keep_fits, augment.ackwards(), tidy(what="scores"), cfQ

# ── cfQ unsupported (all engines) ─────────────────────────────────────────────

test_that("cfQ errors before any engine runs (no engine-specific message)", {
  skip_if_not_installed("psych")
  d <- as.data.frame(matrix(rnorm(300), 100, 6))
  # Error fires at ackwards() validation, not inside the engine
  expect_error(ackwards(d, k = 2, rotation = "cfQ"), "not supported")
  expect_error(ackwards(d, k = 2, method = "efa", rotation = "cfQ"), "not supported")
})

# ── scores = FALSE (default) ──────────────────────────────────────────────────

test_that("scores = FALSE (default) leaves x$scores as NULL", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 3))
  expect_null(x$scores)
})

# ── scores = TRUE: storage ────────────────────────────────────────────────────

test_that("scores = TRUE stores list of n × k_j matrices with correct dims", {
  skip_if_not_installed("psych")
  n <- nrow(psych::bfi)
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 3, scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  expect_equal(nrow(x$scores[["1"]]), n)
  expect_equal(nrow(x$scores[["2"]]), n)
  expect_equal(nrow(x$scores[["3"]]), n)
  expect_equal(ncol(x$scores[["1"]]), 1L)
  expect_equal(ncol(x$scores[["2"]]), 2L)
  expect_equal(ncol(x$scores[["3"]]), 3L)
})

test_that("stored score column names match factor labels", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 3, scores = TRUE))
  expect_equal(colnames(x$scores[["1"]]), "m1f1")
  expect_equal(colnames(x$scores[["2"]]), c("m2f1", "m2f2"))
  expect_equal(colnames(x$scores[["3"]]), c("m3f1", "m3f2", "m3f3"))
})

# ── Invariant 1: standardization ──────────────────────────────────────────────

test_that("PCA tenBerge scores are approximately unit variance (Inv. 1)", {
  skip_if_not_installed("psych")
  set.seed(42)
  data <- as.data.frame(matrix(rnorm(600), 200, 6))
  x <- suppressWarnings(ackwards(data, k = 3, scores = TRUE))
  # For orthogonal cfT scores, SD ≈ 1 (tenBerge: W'RW = I exactly)
  for (ki in names(x$scores)) {
    sds <- apply(x$scores[[ki]], 2, sd, na.rm = TRUE)
    expect_true(all(abs(sds - 1) < 0.05),
      label = paste("SDs ≈ 1 at level", ki)
    )
  }
})

# ── Algebra vs. materialized scores cross-check (Inv. 2) ─────────────────────

test_that("between-level correlations from stored scores agree with algebra edges", {
  skip_if_not_installed("psych")
  set.seed(7)
  n <- 300
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = g + s1 + rnorm(n, sd = 0.2),
    x2 = g + s1 + rnorm(n, sd = 0.2),
    x3 = g + s1 + rnorm(n, sd = 0.2),
    x4 = g + s2 + rnorm(n, sd = 0.2),
    x5 = g + s2 + rnorm(n, sd = 0.2),
    x6 = g + s2 + rnorm(n, sd = 0.2)
  )
  x <- suppressWarnings(ackwards(data, k = 3, scores = TRUE))

  # For each adjacent pair, correlation of materialized scores ≈ algebra edge
  for (ki in seq(2L, x$k_max)) {
    key <- paste0(ki - 1L, ":", ki)
    E_alg <- x$edges$matrices[[key]]
    Sa <- x$scores[[as.character(ki - 1L)]]
    Sb <- x$scores[[as.character(ki)]]
    E_sc <- cor(Sa, Sb)
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 0.05,
      label = paste("algebra vs. score correlations at level pair", key)
    )
  }
})

# ── keep_fits = FALSE (default) ───────────────────────────────────────────────

test_that("keep_fits = FALSE (default) leaves x$fits as NULL", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 2))
  expect_null(x$fits)
})

# ── keep_fits = TRUE: storage ─────────────────────────────────────────────────

test_that("keep_fits = TRUE stores list of raw psych objects for PCA", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 3, keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2", "3"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "psych"),
      label = paste("psych object at level", ki)
    )
  }
})

test_that("keep_fits = TRUE and scores = TRUE can be combined", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(
    ackwards(psych::bfi[, 1:25], k = 2, scores = TRUE, keep_fits = TRUE)
  )
  expect_false(is.null(x$scores))
  expect_false(is.null(x$fits))
})

# ── augment.ackwards() ────────────────────────────────────────────────────────

test_that("augment(x, data) returns data frame with score columns appended", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k = 3))
  out <- augment(x, data = bfi_items)
  expect_s3_class(out, "data.frame")
  # Same number of rows as input
  expect_equal(nrow(out), nrow(bfi_items))
  # Original columns present
  expect_true(all(names(bfi_items) %in% names(out)))
  # Score columns present: .m1f1, .m2f1, .m2f2, .m3f1, .m3f2, .m3f3
  score_cols <- grep("^\\.m", names(out), value = TRUE)
  expect_length(score_cols, 1L + 2L + 3L)
})

test_that("augment(x) without data uses stored scores", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 2, scores = TRUE))
  out <- augment(x)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), x$n_obs)
  score_cols <- grep("^\\.m", names(out), value = TRUE)
  expect_length(score_cols, 1L + 2L)
})

test_that("augment(x) without data or stored scores gives informative error", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 2))
  expect_error(augment(x), "not stored")
})

test_that("augment(x, data) and augment(x) [with stored scores] agree", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k = 2, scores = TRUE))
  out_stored <- augment(x)
  out_recomp <- augment(x, data = bfi_items)
  # Score columns should be numerically identical
  score_cols_stored <- grep("^\\.m", names(out_stored), value = TRUE)
  score_cols_recomp <- grep("^\\.m", names(out_recomp), value = TRUE)
  for (col in score_cols_stored) {
    expect_equal(out_stored[[col]], out_recomp[[col]],
      tolerance = 1e-10,
      label = paste("stored vs. recomputed scores:", col)
    )
  }
})

# ── tidy(what = "scores") ─────────────────────────────────────────────────────

test_that("tidy(x, what='scores') returns long data frame with correct columns", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 3, scores = TRUE))
  out <- tidy(x, what = "scores")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("obs", "level", "factor", "score") %in% names(out)))
  # n_obs × (1 + 2 + 3) rows
  expect_equal(nrow(out), x$n_obs * (1L + 2L + 3L))
  expect_setequal(unique(out$level), 1:3)
})

test_that("tidy(x, what='scores') errors informatively when no scores stored", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k = 2))
  expect_error(tidy(x, what = "scores"), "not stored")
})
