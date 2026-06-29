# tests/testthat/test-scores.R
# M6: scores storage, keep_fits, augment.ackwards(), tidy(what="scores")

# ── keep_scores = FALSE (default) ──────────────────────────────────────────────────

test_that("keep_scores = FALSE (default) leaves x$scores as NULL", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$scores)
})

# ── keep_scores = TRUE: storage ────────────────────────────────────────────────────

test_that("keep_scores = TRUE stores list of n × k_j matrices with correct dims", {
  skip_if_not_installed("psych")
  n <- nrow(psych::bfi)
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
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
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
  expect_equal(colnames(x$scores[["1"]]), "m1f1")
  expect_equal(colnames(x$scores[["2"]]), c("m2f1", "m2f2"))
  expect_equal(colnames(x$scores[["3"]]), c("m3f1", "m3f2", "m3f3"))
})

# ── Invariant 1: standardization ──────────────────────────────────────────────

test_that("PCA tenBerge scores are approximately unit variance (Inv. 1)", {
  skip_if_not_installed("psych")
  set.seed(42)
  data <- as.data.frame(matrix(rnorm(600), 200, 6))
  x <- suppressWarnings(ackwards(data, k_max = 3, keep_scores = TRUE))
  # For orthogonal (varimax) scores, SD ≈ 1 (tenBerge: W'RW = I exactly)
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
  x <- suppressWarnings(ackwards(data, k_max = 3, keep_scores = TRUE))

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
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_null(x$fits)
})

# ── keep_fits = TRUE: storage ─────────────────────────────────────────────────

test_that("keep_fits = TRUE stores list of raw psych objects for PCA", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2", "3"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "psych"),
      label = paste("psych object at level", ki)
    )
  }
})

test_that("keep_fits = TRUE and keep_scores = TRUE can be combined", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(
    ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE, keep_fits = TRUE)
  )
  expect_false(is.null(x$scores))
  expect_false(is.null(x$fits))
})

# ── augment.ackwards() ────────────────────────────────────────────────────────

test_that("augment(x, data) returns data frame with score columns appended", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  out <- suppressWarnings(augment(x, data = bfi_items))
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
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE))
  out <- augment(x)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), x$n_obs)
  score_cols <- grep("^\\.m", names(out), value = TRUE)
  expect_length(score_cols, 1L + 2L)
})

test_that("augment(x) without data or stored scores gives informative error", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(augment(x), "not stored")
})

test_that("augment(x, data) and augment(x) [with stored scores] agree", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2, keep_scores = TRUE))
  out_stored <- suppressWarnings(augment(x))
  out_recomp <- suppressWarnings(augment(x, data = bfi_items))
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
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
  out <- tidy(x, what = "scores")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("obs", "level", "factor", "score") %in% names(out)))
  # n_obs × (1 + 2 + 3) rows
  expect_equal(nrow(out), x$n_obs * (1L + 2L + 3L))
  expect_setequal(unique(out$level), 1:3)
})

test_that("tidy(x, what='scores') errors informatively when no scores stored", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "scores"), "not stored")
})

# ── B1: Engine coverage for scores and keep_fits ──────────────────────────────

test_that("EFA keep_scores = TRUE stores correctly shaped matrices", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  x <- suppressWarnings(ackwards(d, k_max = 3, engine = "efa", keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  for (ki in 1:3) {
    expect_equal(nrow(x$scores[[as.character(ki)]]), 300L,
      label = paste("EFA score nrow at level", ki)
    )
    expect_equal(ncol(x$scores[[as.character(ki)]]), ki,
      label = paste("EFA score ncol at level", ki)
    )
  }
})

test_that("EFA keep_fits = TRUE stores psych objects for all levels", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  x <- suppressWarnings(ackwards(d, k_max = 2, engine = "efa", keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "psych"),
      label = paste("EFA psych object at level", ki)
    )
  }
})

test_that("ESEM keep_scores = TRUE stores correctly shaped matrices", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem", keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  for (ki in 1:3) {
    expect_equal(nrow(x$scores[[as.character(ki)]]), nrow(d),
      label = paste("ESEM score nrow at level", ki)
    )
    expect_equal(ncol(x$scores[[as.character(ki)]]), ki,
      label = paste("ESEM score ncol at level", ki)
    )
  }
})

test_that("ESEM keep_fits = TRUE stores lavaan objects for all levels", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 2, engine = "esem", keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "lavaan"),
      label = paste("ESEM lavaan object at level", ki)
    )
  }
})

# ── B3: Truncation — scores only cover converged levels ───────────────────────

test_that("keep_scores = TRUE only covers converged levels when model is truncated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # 6 variables → lavaan::efa() can only fit k <= 3; k = 5 triggers truncation
  suppressWarnings(x <- ackwards(d, k_max = 5, engine = "esem", keep_scores = TRUE))
  expect_equal(x$k_max, 3L)
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  expect_false("4" %in% names(x$scores))
  expect_false("5" %in% names(x$scores))
})

# ── B5: augment() column validation ──────────────────────────────────────────

test_that("augment(x, data) errors when data has wrong column count (unnamed)", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  # Use a plain matrix (no colnames) with wrong column count → dimension path
  d_wrong <- matrix(rnorm(100 * 5), 100, 5)
  expect_error(augment(x, data = d_wrong), "5.*column|column.*5")
})

test_that("augment(x, data) errors when named data is missing expected columns", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  # Supply only 20 of the 25 named BFI columns
  d_missing <- psych::bfi[1:50, 1:20]
  expect_error(augment(x, data = d_missing), "[Mm]issing")
})

test_that("augment(x, data) works when data has extra named columns (supersets)", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[1:50, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2))
  # Add an extra column not in the model
  d_extra <- cbind(bfi_items, extra = rnorm(50))
  out <- suppressWarnings(augment(x, data = d_extra))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 50L)
  expect_true("extra" %in% names(out))
  expect_true(".m1f1" %in% names(out))
})

test_that("augment(x, data) errors on non-numeric data", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[1:50, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2))
  d_chr <- as.data.frame(lapply(bfi_items, as.character))
  expect_error(augment(x, data = d_chr), "[Nn]umeric")
})

# ── NA data: warning + correct propagation ────────────────────────────────────

test_that("augment(x, data) warns when data has missing rows", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  x <- ackwards(d, k_max = 2L)
  d_na <- d
  d_na[5L, 1L] <- NA_real_
  expect_warning(augment(x, data = d_na), "missing")
})

test_that("augment(x, data) produces NA scores for rows with missing values", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  x <- ackwards(d, k_max = 2L)
  d_na <- d
  d_na[5L, 1L] <- NA_real_
  out <- suppressWarnings(augment(x, data = d_na))
  expect_true(is.na(out$.m1f1[5L]))
})

test_that("keep_scores = TRUE warns about NA score propagation", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  d[5L, 1L] <- NA_real_
  # Collect all warnings (pairwise advisory + score NA) and verify the score
  # NA warning is present; suppress so no warnings leak to testthat.
  warns <- character(0L)
  withCallingHandlers(
    ackwards(d, k_max = 2L, keep_scores = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("NA scores", warns)))
})

test_that("keep_scores = TRUE produces NA scores for rows with missing values", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  d[5L, 1L] <- NA_real_
  x <- suppressWarnings(ackwards(d, k_max = 2L, keep_scores = TRUE))
  expect_true(any(is.na(x$scores[["1"]])))
})

# ── Non-Pearson basis: warning on scoring ────────────────────────────────────

test_that("scoring with non-Pearson basis warns about basis mismatch", {
  skip_if_not_installed("psych")
  set.seed(42)
  # Build a small ordinal data set so polychoric is valid
  d <- as.data.frame(matrix(
    sample(1L:5L, 150L * 6L, replace = TRUE), 150L, 6L
  ))
  x <- suppressWarnings(ackwards(d, k_max = 2, cor = "polychoric"))
  expect_warning(augment(x, data = d), "polychoric")
})

# ── B4: keep_fits truncation ─────────────────────────────────────────────────

test_that("keep_fits = TRUE only stores fits for converged levels when truncated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # 6 variables → lavaan::efa() truncates at k = 3; k = 5 requested
  suppressWarnings(x <- ackwards(d, k_max = 5, engine = "esem", keep_fits = TRUE))
  expect_equal(x$k_max, 3L)
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2", "3"))
  expect_false("4" %in% names(x$fits))
  expect_false("5" %in% names(x$fits))
})
