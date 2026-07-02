test_that("ackwards() returns a valid ackwards object for the smoke-test case", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x) # internal: checks all required fields

  expect_equal(x$engine, "pca")
  expect_equal(x$k_max, 5L)
  expect_equal(x$n_obs, 2800L)
  expect_equal(x$cor, "pearson")

  # All five levels present and converged
  expect_equal(length(x$levels), 5L)
  expect_true(all(vapply(x$levels, `[[`, logical(1), "converged")))

  # Four adjacent edge matrices (1:2, 2:3, 3:4, 4:5)
  expect_equal(length(x$edges$matrices), 4L)
  expect_named(x$edges$matrices, c("1:2", "2:3", "3:4", "4:5"))

  # Edge matrix dimensions
  for (i in 1:4) {
    E <- x$edges$matrices[[paste0(i, ":", i + 1)]]
    expect_equal(nrow(E), i, info = paste("rows of", i, ":", i + 1))
    expect_equal(ncol(E), i + 1, info = paste("cols of", i, ":", i + 1))
    expect_true(all(abs(E) <= 1 + 1e-9), info = "correlations in [-1, 1]")
  }
})

test_that("PCA edge correlations match psych::bassAckward within tolerance", {
  skip_if_not_installed("psych")

  ba <- psych::bassAckward(
    psych::bfi[, 1:25],
    nfactors = 5, fm = "pca", rotate = "varimax", plot = FALSE
  )
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))

  for (i in 1:4) {
    # psych stores (k+1) × k; ours is k × (k+1) — transpose to align
    psych_mat <- t(ba$bass.ack[[i + 1L]])
    our_mat <- x$edges$matrices[[paste0(i, ":", i + 1)]]

    max_diff <- max(abs(abs(psych_mat) - abs(our_mat)))
    expect_lt(max_diff, 1e-4, label = paste("Level", i, "->", i + 1))
  }
})

test_that("levels have correct structure and label formats", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))

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
  expect_error(ackwards(d, k_max = 0), "integer >= 2")
  expect_error(ackwards(d, k_max = 1), "integer >= 2")
  expect_error(ackwards(d, k_max = 100), "cannot exceed")
  expect_error(ackwards(d, k_max = 1.5), "integer >= 2")
  expect_error(ackwards(list(), k_max = 2), "data frame")
  # cut_show thresholds |r|, so it must be a single number in [0, 1] (M42/m8)
  expect_error(ackwards(d, k_max = 2, cut_show = 5), "\\[0, 1\\]")
  expect_error(ackwards(d, k_max = 2, cut_show = -0.1), "\\[0, 1\\]")
  expect_error(ackwards(d, k_max = 2, cut_show = c(0.3, 0.5)), "\\[0, 1\\]")
  expect_error(ackwards(d, k_max = 2, cut_show = NA), "\\[0, 1\\]")
  expect_error(ackwards(d, k_max = 2, cut_show = "0.3"), "\\[0, 1\\]")
})

test_that("keep_scores = TRUE stores a named list of score matrices", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2"))
  expect_equal(nrow(x$scores[["1"]]), nrow(psych::bfi))
  expect_equal(ncol(x$scores[["1"]]), 1L)
  expect_equal(ncol(x$scores[["2"]]), 2L)
})

test_that("keep_fits = TRUE stores a named list of raw fit objects", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2, keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2"))
  expect_true(inherits(x$fits[["1"]], "psych"))
})

test_that("meta$cut_show stores the cut_show value", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2, cut_show = 0.4))
  expect_equal(x$meta$cut_show, 0.4)
})

test_that("detect_ordinal() flags bfi columns", {
  skip_if_not_installed("psych")
  expect_length(ackwards:::detect_ordinal(psych::bfi[, 1:25]), 25L)
})

test_that("ordinal-detection warning names the flagged columns (M42/e3)", {
  skip_if_not_installed("psych")
  d <- data.frame(
    L1 = sample(1:5, 100, replace = TRUE),
    L2 = sample(1:5, 100, replace = TRUE),
    c1 = rnorm(100), c2 = rnorm(100), c3 = rnorm(100)
  )
  rlang::reset_warning_verbosity("ackwards_ordinal_warning")
  expect_warning(ackwards(d, k_max = 2), "L1")
})

test_that("detect_ordinal() flags nothing for continuous data", {
  set.seed(1)
  d <- data.frame(matrix(rnorm(500), 100, 5))
  expect_identical(ackwards:::detect_ordinal(d), character(0))
})

test_that("meta$ordinal_warned is TRUE for bfi data", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_true(x$meta$ordinal_warned)
})

test_that("kappa is not a parameter of ackwards() (M13: removed dead arg)", {
  expect_false("kappa" %in% names(formals(ackwards)))
})

test_that("ackwards() meta does not store kappa (M13: removed)", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_false("kappa" %in% names(x$meta))
})

# ── tidy(sort = "strength") ────────────────────────────────────────────────────

test_that("tidy(x, sort = 'strength') returns edges in descending |r| order", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  out <- tidy(x, sort = "strength")
  expect_s3_class(out, "data.frame")
  expect_true(all(diff(abs(out$r)) <= 0))
})

test_that("tidy(x, sort = 'none') is byte-identical to the default", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  expect_identical(tidy(x, sort = "none"), tidy(x))
})

test_that("tidy(x, sort = 'strength', what = 'loadings') errors informatively", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "loadings", sort = "strength"), "edges")
})

# ── tidy(primary_only = TRUE) ──────────────────────────────────────────────────

test_that("tidy(x, primary_only = TRUE) returns only is_primary edges", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  full <- tidy(x, what = "edges")
  prim <- tidy(x, what = "edges", primary_only = TRUE)
  expect_true(all(prim$is_primary))
  expect_identical(nrow(prim), sum(full$is_primary))
})

test_that("tidy(x, primary_only = TRUE, sort = 'strength') filters and sorts", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5))
  out <- tidy(x, what = "edges", primary_only = TRUE, sort = "strength")
  expect_true(all(out$is_primary))
  expect_true(all(diff(abs(out$r)) <= 0))
})

test_that("tidy(x, primary_only = TRUE, what = 'loadings') errors informatively", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(
    tidy(x, what = "loadings", primary_only = TRUE),
    "edges"
  )
})

# ── tidy(what = "loadings") has SE/CI columns, NA for PCA ──────────────────────

test_that("tidy(x, what = 'loadings') has se/ci cols present but NA for PCA", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 2))
  ld <- tidy(x, what = "loadings")
  expect_true(all(c("se", "ci_lower", "ci_upper") %in% names(ld)))
  expect_true(all(is.na(ld$se)))
})

# ── glance() fit columns ───────────────────────────────────────────────────────

test_that("glance() for PCA has fit columns present but all NA", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  g <- generics::glance(x)
  expect_true(all(c("CFI", "TLI", "RMSEA", "SRMR", "BIC") %in% names(g)))
  expect_true(is.na(g$CFI))
  expect_true(is.na(g$TLI))
  expect_true(is.na(g$RMSEA))
  expect_true(is.na(g$SRMR))
  expect_true(is.na(g$BIC))
})

# ── .glance_fit defensive branches (Invariant 7: truncation) ──────────────────

test_that("glance() returns all-NA fit when deepest_converged < 2", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  # Simulate a run that truncated to the 1-factor anchor only.
  x$meta$deepest_converged <- 1L
  g <- generics::glance(x)
  expect_true(is.na(g$CFI))
  expect_true(is.na(g$RMSEA))
  expect_true(is.na(g$BIC))
})

test_that("glance() returns all-NA fit when the deepest level carries no fit", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3))
  # Force the empty-fit guard: deepest level has no stored fit vector.
  x$levels[[as.character(x$meta$deepest_converged)]]$fit <- numeric(0)
  g <- generics::glance(x)
  expect_true(is.na(g$CFI))
  expect_true(is.na(g$RMSEA))
  expect_true(is.na(g$BIC))
})
