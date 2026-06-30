# Tests for correlation-matrix input (M22)

# Helper: build a valid small correlation matrix from the first 6 bfi25 columns
.bfi6_R <- function() {
  dat <- bfi25[, 1:6]
  cor(dat, use = "pairwise.complete.obs")
}

# ---- ackwards(): detection ---------------------------------------------------

test_that("ackwards() detects R-matrix input and sets input_type", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_equal(x$meta$input_type, "cor_matrix")
})

test_that("ackwards() sets cor = NA for R-matrix input", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_true(is.na(x$cor))
})

test_that("ackwards() data path still sets input_type = 'data'", {
  x <- suppressWarnings(ackwards(bfi25[, 1:6], k_max = 3))
  expect_equal(x$meta$input_type, "data")
})

# ---- ackwards(): edge correctness (algebra parity) --------------------------

test_that("R-matrix PCA edges match raw-data PCA edges within tolerance", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  x_R <- suppressMessages(ackwards(R, k_max = 5))
  x_d <- suppressWarnings(ackwards(bfi25, k_max = 5))
  expect_equal(tidy(x_R)$r, tidy(x_d)$r, tolerance = 1e-10)
})

test_that("R-matrix EFA edges match raw-data EFA edges within tolerance", {
  skip_if_not_installed("psych")
  set.seed(42)
  R <- cor(bfi25, use = "pairwise.complete.obs")
  x_R <- suppressMessages(ackwards(R, k_max = 5, engine = "efa", n_obs = 875L))
  set.seed(42)
  x_d <- suppressWarnings(ackwards(bfi25, k_max = 5, engine = "efa"))
  expect_equal(tidy(x_R)$r, tidy(x_d)$r, tolerance = 1e-8)
})

# ---- ackwards(): n_obs handling ----------------------------------------------

test_that("EFA + R-matrix without n_obs errors with helpful message", {
  R <- .bfi6_R()
  expect_error(
    suppressMessages(ackwards(R, k_max = 3, engine = "efa")),
    "n_obs"
  )
})

test_that("PCA + R-matrix without n_obs is allowed and stores NA", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_true(is.na(x$n_obs))
})

test_that("PCA + R-matrix with n_obs stores the supplied value", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3, n_obs = 100L))
  expect_equal(x$n_obs, 100L)
})

test_that("n_obs supplied with raw data warns and is ignored", {
  expect_warning(
    ackwards(bfi25[, 1:6], k_max = 3, n_obs = 999L),
    "n_obs.*ignored"
  )
})

# ---- ackwards(): engine gating -----------------------------------------------

test_that("engine = 'esem' + R-matrix errors clearly", {
  R <- .bfi6_R()
  expect_error(
    suppressMessages(ackwards(R, k_max = 3, engine = "esem")),
    "esem.*requires raw"
  )
})

# ---- ackwards(): cor/missing warnings ----------------------------------------

test_that("cor != default + R-matrix warns it is ignored", {
  R <- .bfi6_R()
  expect_warning(
    suppressMessages(ackwards(R, k_max = 3, cor = "polychoric")),
    "cor.*ignored"
  )
})

test_that("missing != default + R-matrix warns it is ignored", {
  R <- .bfi6_R()
  expect_warning(
    suppressMessages(ackwards(R, k_max = 3, missing = "listwise")),
    "missing.*ignored"
  )
})

test_that("no ordinal warning fires for R-matrix input", {
  R <- .bfi6_R()
  # Should not warn about ordinal data even though bfi6 items are Likert
  expect_no_warning(
    withCallingHandlers(
      suppressMessages(ackwards(R, k_max = 3)),
      simpleWarning = function(w) {
        if (grepl("ordinal|Likert", conditionMessage(w))) stop(w)
        invokeRestart("muffleWarning")
      }
    )
  )
})

# ---- ackwards(): scores blocked for R-matrix ---------------------------------

test_that("keep_scores = TRUE + R-matrix errors clearly", {
  R <- .bfi6_R()
  expect_error(
    suppressMessages(ackwards(R, k_max = 3, keep_scores = TRUE)),
    "keep_scores.*requires raw"
  )
})

test_that("augment() on R-matrix fit with no data errors informatively", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_error(augment(x), "not stored")
})

test_that("tidy(what = 'scores') on R-matrix fit errors informatively", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_error(tidy(x, what = "scores"), "not stored")
})

# ---- ackwards(): print / summary / glance / autoplot work -------------------

test_that("print() works on R-matrix fit object without error", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  # cli writes to stderr; verify no error + returns invisibly
  expect_no_error(print(x))
  expect_invisible(print(x))
  # Content verified via the NA cor field that drives the "(user-supplied matrix)" label
  expect_true(is.na(x$cor))
})

test_that("summary() works on R-matrix fit object", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  s <- summary(x)
  expect_s3_class(s, "summary_ackwards")
  expect_true(is.na(s$cor))
})

test_that("glance() works on R-matrix fit object", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  g <- glance(x)
  expect_s3_class(g, "data.frame")
  expect_true("n_obs" %in% names(g))
})

test_that("tidy() works on R-matrix fit object", {
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  expect_s3_class(tidy(x), "data.frame")
})

test_that("autoplot() works on R-matrix fit object", {
  skip_if_not_installed("ggplot2")
  R <- .bfi6_R()
  x <- suppressMessages(ackwards(R, k_max = 3))
  p <- autoplot(x)
  expect_s3_class(p, "gg")
})

# ---- ackwards(): R-matrix validation errors ----------------------------------
# Note: non-square, non-symmetric, and non-unit-diagonal matrices fail
# .is_cor_matrix() and are routed to the raw-data path (Branch B), not Branch A.
# Their validation is tested directly in test-utils.R via .validate_cor_matrix().
# Here we test the one case that passes .is_cor_matrix() but fails validation:
# a symmetric matrix with unit diagonal but off-diagonal |r| > 1.

test_that("R-like matrix with |r| > 1 off-diagonal errors clearly", {
  # Passes .is_cor_matrix() (symmetric, unit diagonal) but fails validation
  R <- matrix(c(1, 1.5, 1.5, 1), 2, 2)
  expect_error(suppressMessages(ackwards(R, k_max = 2)), "off-diagonal")
})

test_that("covariance matrix gives a targeted error not a confusing one", {
  S <- cov(bfi25, use = "pairwise.complete.obs")
  expect_error(ackwards(S, k_max = 3), "covariance matrix")
})

# ---- ackwards(): prune works with R-matrix input ----------------------------

test_that("prune='redundant' works with R-matrix input", {
  R <- .bfi6_R()
  x <- suppressMessages(
    ackwards(R, k_max = 3, prune = "redundant")
  )
  expect_s3_class(x, "ackwards")
  expect_false(is.null(x$prune))
})

# ---- suggest_k(): R-matrix input --------------------------------------------

test_that("suggest_k() accepts R-matrix input with n_obs", {
  skip_if_not_installed("psych")
  R <- .bfi6_R()
  sk <- suppressMessages(suggest_k(R, n_obs = 875L, n_iter = 5L))
  expect_s3_class(sk, "suggest_k")
  expect_equal(sk$input_type, "cor_matrix")
  expect_equal(sk$n_obs, 875L)
  expect_true(is.na(sk$cor))
})

test_that("suggest_k() R-matrix sets cd_available = FALSE", {
  skip_if_not_installed("psych")
  R <- .bfi6_R()
  sk <- suppressMessages(suggest_k(R, n_obs = 875L, n_iter = 5L))
  expect_false(sk$cd_available)
})

test_that("suggest_k() R-matrix without n_obs errors", {
  R <- .bfi6_R()
  expect_error(suggest_k(R), "n_obs.*required")
})

test_that("suggest_k() print works without error for R input", {
  skip_if_not_installed("psych")
  R <- .bfi6_R()
  sk <- suppressMessages(suggest_k(R, n_obs = 875L, n_iter = 5L))
  # cli writes to stderr; verify no error + NA cor field drives the label
  expect_no_error(print(sk))
  expect_invisible(print(sk))
  expect_true(is.na(sk$cor))
})

test_that("suggest_k() n_obs ignored for raw data with warning", {
  skip_if_not_installed("psych")
  expect_warning(
    suggest_k(bfi25[, 1:6], n_obs = 999L, n_iter = 5L),
    "n_obs.*ignored"
  )
})

test_that("suggest_k() errors on non-positive n_obs", {
  R <- .bfi6_R()
  expect_error(suggest_k(R, n_obs = -1L), "positive integer")
})

test_that("suggest_k() errors on non-integer n_obs", {
  R <- .bfi6_R()
  expect_error(suggest_k(R, n_obs = 100.5), "positive integer")
})

test_that("suggest_k() covariance matrix gives a targeted error", {
  S <- cov(bfi25, use = "pairwise.complete.obs")
  expect_error(suggest_k(S), "covariance matrix")
})

# ---- ackwards(): additional validation for R-matrix branch ------------------

test_that("ackwards() R-matrix: invalid n_obs (negative) errors", {
  R <- .bfi6_R()
  expect_error(
    suppressMessages(ackwards(R, k_max = 3, n_obs = -1L)),
    "positive integer"
  )
})

test_that("ackwards() R-matrix: invalid n_obs (non-numeric) errors", {
  R <- .bfi6_R()
  expect_error(
    suppressMessages(ackwards(R, k_max = 3, n_obs = "abc")),
    "positive integer"
  )
})

test_that("ackwards() R-matrix: k_max > p errors clearly", {
  R <- .bfi6_R() # 6 x 6 matrix
  expect_error(
    suppressMessages(ackwards(R, k_max = 100L)),
    "k_max.*exceed"
  )
})

test_that("ackwards() R-matrix with seed runs without error", {
  R <- .bfi6_R()
  expect_no_error(suppressMessages(ackwards(R, k_max = 3, seed = 42L)))
})

# ---- ackwards(): raw data validation branches --------------------------------

test_that("ackwards() errors on non-numeric raw data (character columns)", {
  d <- data.frame(
    x1 = c("a", "b", "c"),
    x2 = c("d", "e", "f"),
    x3 = c("g", "h", "i"),
    stringsAsFactors = FALSE
  )
  expect_error(ackwards(d, k_max = 2), "numeric")
})

test_that("ackwards() raw data PCA with seed runs without error", {
  set.seed(1L)
  d <- as.data.frame(matrix(rnorm(120L), 20L, 6L))
  expect_no_error(suppressWarnings(ackwards(d, k_max = 3, engine = "pca", seed = 7L)))
})

test_that("autoplot.suggest_k() works on R-matrix suggest_k object", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  R <- .bfi6_R()
  sk <- suppressMessages(suggest_k(R, n_obs = 875L, n_iter = 5L))
  p <- suppressMessages(autoplot(sk))
  expect_s3_class(p, "gg")
  expect_false("CD (RMSE, minimize)" %in% levels(p$data$panel))
})
