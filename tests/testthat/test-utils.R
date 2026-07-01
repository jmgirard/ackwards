test_that(".format_r() strips leading zero and pads trailing zeros", {
  fmt <- ackwards:::.format_r

  expect_equal(fmt(0.23, 2L), ".23")
  expect_equal(fmt(-0.23, 2L), "-.23")
  expect_equal(fmt(0.3, 2L), ".30")
  expect_equal(fmt(-0.3, 2L), "-.30")
  expect_equal(fmt(0, 2L), ".00")
  expect_equal(fmt(1, 2L), "1.00")
  expect_equal(fmt(-1, 2L), "-1.00")
})

test_that(".format_r() does not produce '-.00' when magnitude rounds to zero", {
  fmt <- ackwards:::.format_r

  # -0.003 rounds to .00 at 2 digits; sign should be suppressed
  expect_equal(fmt(-0.003, 2L), ".00")
  # Exact zero is also clean
  expect_equal(fmt(-0, 2L), ".00")
})

test_that(".format_r() respects digits argument", {
  fmt <- ackwards:::.format_r

  expect_equal(fmt(0.234, 3L), ".234")
  expect_equal(fmt(0.2, 1L), ".2")
})

test_that(".format_r() is vectorised", {
  fmt <- ackwards:::.format_r

  result <- fmt(c(0.23, -0.45, 1, 0), 2L)
  expect_equal(result, c(".23", "-.45", "1.00", ".00"))
})

test_that("detect_ordinal returns FALSE for an all-NA column", {
  # all-NA column alongside a clearly continuous column
  d <- data.frame(x = rep(NA_real_, 5L), y = c(1.1, 2.2, 3.3, 4.4, 5.5))
  expect_false(ackwards:::detect_ordinal(d))
})

test_that("detect_ordinal returns FALSE for a data frame with no ordinal columns", {
  d <- data.frame(x = c(1.1, 2.2, 3.3), y = c(4.4, 5.5, 6.6))
  expect_false(ackwards:::detect_ordinal(d))
})

test_that("detect_ordinal returns TRUE for integer columns with few distinct values", {
  d <- data.frame(x = sample(1:5, 20, replace = TRUE))
  expect_true(ackwards:::detect_ordinal(d))
})

test_that(".standardize produces NA only in cells where input was NA", {
  m <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 3, ncol = 2)
  z <- ackwards:::.standardize(m)
  # Cell [2,1] was NA -> stays NA
  expect_true(is.na(z[2L, 1L]))
  # Cell [2,2] was not NA -> should be finite
  expect_false(is.na(z[2L, 2L]))
  # Rows 1 and 3 have no NAs -> all finite
  expect_false(anyNA(z[c(1L, 3L), ]))
})

test_that(".is_cor_matrix() returns TRUE for a valid correlation matrix", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  expect_true(ackwards:::.is_cor_matrix(R))
})

test_that(".is_cor_matrix() returns FALSE for raw data", {
  expect_false(ackwards:::.is_cor_matrix(as.matrix(bfi25)))
})

test_that(".is_cor_matrix() returns FALSE for non-square matrix", {
  expect_false(ackwards:::.is_cor_matrix(matrix(1:6, 2, 3)))
})

test_that(".is_cor_matrix() returns FALSE for non-numeric input", {
  expect_false(ackwards:::.is_cor_matrix(list(a = 1, b = 2)))
})

test_that(".is_cor_matrix() returns FALSE for asymmetric matrix", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  R[1, 2] <- 0.999
  expect_false(ackwards:::.is_cor_matrix(R))
})

test_that(".validate_cor_matrix() errors on non-matrix input", {
  expect_error(
    ackwards:::.validate_cor_matrix(list(a = 1)),
    "numeric matrix"
  )
})

test_that(".validate_cor_matrix() errors on non-square matrix", {
  expect_error(
    ackwards:::.validate_cor_matrix(matrix(1:6, 2, 3)),
    "square"
  )
})

test_that(".validate_cor_matrix() errors when NA present", {
  R <- diag(3)
  R[1, 2] <- R[2, 1] <- NA_real_
  expect_error(
    ackwards:::.validate_cor_matrix(R),
    "NA value"
  )
})

test_that(".validate_cor_matrix() errors on non-symmetric matrix", {
  R <- diag(3)
  R[1, 2] <- 0.5
  expect_error(
    ackwards:::.validate_cor_matrix(R),
    "symmetric"
  )
})

test_that(".validate_cor_matrix() errors when diagonal != 1", {
  R <- matrix(c(1, 0.5, 0.5, 2), 2, 2) # diagonal element is 2
  expect_error(
    ackwards:::.validate_cor_matrix(R),
    "diagonal must be all 1s"
  )
})

test_that(".validate_cor_matrix() errors when |r| > 1", {
  R <- matrix(c(1, 1.5, 1.5, 1), 2, 2)
  expect_error(
    ackwards:::.validate_cor_matrix(R),
    "off-diagonal"
  )
})

test_that(".validate_cor_matrix() synthesises dimnames when absent", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  dimnames(R) <- NULL
  R2 <- ackwards:::.validate_cor_matrix(R)
  expect_equal(rownames(R2), paste0("V", seq_len(ncol(R2))))
})

test_that(".validate_cor_matrix() passes valid matrix through unchanged", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  R2 <- ackwards:::.validate_cor_matrix(R)
  expect_equal(R2, R)
})

test_that(".validate_cor_matrix() warns on non-positive-definite matrix", {
  # Inconsistent triangle: Var1-Var2 and Var1-Var3 are highly positive,
  # but Var2-Var3 are highly negative -- impossible under transitivity,
  # so det < 0 and the matrix has a negative eigenvalue.
  R <- matrix(
    c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3,
    dimnames = list(paste0("V", 1:3), paste0("V", 1:3))
  )
  expect_warning(
    ackwards:::.validate_cor_matrix(R),
    "not positive definite"
  )
})

test_that(".check_maybe_cov_matrix() errors on covariance-like matrix", {
  S <- cov(bfi25, use = "pairwise.complete.obs")
  expect_error(
    ackwards:::.check_maybe_cov_matrix(S),
    "covariance matrix"
  )
})

test_that(".check_maybe_cov_matrix() passes silently for a valid correlation matrix", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  expect_no_error(ackwards:::.check_maybe_cov_matrix(R))
  expect_no_error(ackwards:::.check_maybe_cov_matrix(as.data.frame(bfi25)))
})

test_that(".check_maybe_cov_matrix() returns invisibly for a non-square numeric matrix", {
  # Non-square: rows != cols, so not a square symmetric matrix → early return
  m <- matrix(1:6, nrow = 2, ncol = 3)
  expect_no_error(ackwards:::.check_maybe_cov_matrix(m))
})

test_that(".validate_cor_matrix() errors on matrix containing Inf", {
  R <- diag(3L)
  dimnames(R) <- list(c("a", "b", "c"), c("a", "b", "c"))
  R[1L, 2L] <- R[2L, 1L] <- Inf
  expect_error(
    ackwards:::.validate_cor_matrix(R),
    "non-finite"
  )
})

test_that("validate_ackwards() errors when required fields are missing", {
  x <- structure(list(call = NULL), class = "ackwards")
  expect_error(validate_ackwards(x), "missing fields")
})

test_that("validate_ackwards() errors when object does not inherit ackwards class", {
  required_fields <- c(
    "call", "engine", "rotation", "cor", "n_obs", "k_max",
    "seed", "pkg_version", "levels", "edges", "lineage",
    "scores", "fits", "r", "data", "meta", "prune"
  )
  x <- structure(
    stats::setNames(vector("list", length(required_fields)), required_fields),
    class = "not_ackwards"
  )
  expect_error(validate_ackwards(x), "class")
})

test_that("inst/CITATION has a single Girard software entry, no Goldberg", {
  skip_if_not_installed("ackwards")
  cites <- citation("ackwards")
  expect_length(cites, 1)
  # Package entry: authored by Girard, not Goldberg -- the method paper is
  # cited in DESCRIPTION and roxygen @references, not the software citation.
  expect_match(paste(format(cites[[1]]), collapse = "\n"), "Girard")
  expect_false(any(grepl("Goldberg", format(cites))))
  # Year is set, version note is present, no "????" fallback
  pkg_year <- cites[[1]]$year
  expect_false(is.null(pkg_year))
  expect_false(is.na(pkg_year))
  expect_match(pkg_year, "^[0-9]{4}$")
  expect_match(cites[[1]]$note, "R package version")
})
