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

test_that(".fmt_pct() keeps a fixed single decimal with trailing zeros", {
  pct <- ackwards:::.fmt_pct

  # The drift the shared helper prevents: a whole/half number keeps its ".0"
  # rather than collapsing to "30" (the old round()-based print path).
  expect_equal(pct(0.30), "30.0")
  expect_equal(pct(0.302), "30.2")
  expect_equal(pct(0.1355), "13.6") # rounds to 1 dp
  expect_equal(pct(1), "100.0")
  expect_equal(pct(0), "0.0")
})

test_that(".ok_glyph() returns a tick when ok and a cross otherwise", {
  glyph <- ackwards:::.ok_glyph

  # Maps TRUE -> tick, FALSE -> cross using cli::symbol, so it is terminal-
  # adaptive (UTF-8 ✔/✖, ASCII v/x) rather than a hard-coded literal.
  expect_equal(glyph(TRUE), cli::symbol$tick)
  expect_equal(glyph(FALSE), cli::symbol$cross)
  expect_false(identical(glyph(TRUE), glyph(FALSE)))
})

test_that("detect_ordinal flags nothing for an all-NA column", {
  # all-NA column alongside a clearly continuous column
  d <- data.frame(x = rep(NA_real_, 5L), y = c(1.1, 2.2, 3.3, 4.4, 5.5))
  expect_identical(ackwards:::detect_ordinal(d), character(0))
})

test_that("detect_ordinal flags nothing for a data frame with no ordinal columns", {
  d <- data.frame(x = c(1.1, 2.2, 3.3), y = c(4.4, 5.5, 6.6))
  expect_identical(ackwards:::detect_ordinal(d), character(0))
})

test_that("detect_ordinal names integer columns with few distinct values (M42/e3)", {
  d <- data.frame(
    x = sample(1:5, 20, replace = TRUE),
    y = rnorm(20)
  )
  expect_identical(ackwards:::detect_ordinal(d), "x")
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

test_that(".is_cor_matrix() returns FALSE for a square, non-unit-diagonal matrix", {
  # Square and symmetric but diagonal != 1 (e.g. a covariance-like matrix):
  # the unit-diagonal check short-circuits before the symmetry check.
  M <- matrix(c(2, 0.5, 0.5, 2), 2, 2)
  expect_false(ackwards:::.is_cor_matrix(M))
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
  R2 <- ackwards:::.validate_cor_matrix(R)$R
  expect_equal(rownames(R2), paste0("V", seq_len(ncol(R2))))
})

test_that(".validate_cor_matrix() passes valid matrix through unchanged", {
  R <- cor(bfi25, use = "pairwise.complete.obs")
  out <- ackwards:::.validate_cor_matrix(R)
  expect_equal(out$R, R)
  # min_eig rides along with the matrix it describes (M60)
  expect_equal(out$min_eig, min(eigen(R, symmetric = TRUE, only.values = TRUE)$values))
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

test_that(".check_maybe_cov_matrix() returns invisibly for a square, non-symmetric matrix", {
  # Square but not symmetric: not covariance-like, so no error and no abort
  # (the symmetry check short-circuits before the diagonal check).
  m <- matrix(c(1, 0.5, 0.2, 1), nrow = 2, ncol = 2)
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

test_that(".capture_item_labels() returns named labels only for labelled columns", {
  cap <- ackwards:::.capture_item_labels
  d <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  attr(d$a, "label") <- "First item"
  attr(d$c, "label") <- "Third item"
  labs <- cap(d)
  expect_equal(labs, c(a = "First item", c = "Third item"))
})

test_that(".capture_item_labels() returns NULL when nothing is labelled or input is not a data.frame", {
  cap <- ackwards:::.capture_item_labels
  # Unlabelled data.frame
  expect_null(cap(data.frame(a = 1:3, b = 4:6)))
  # Matrix input (labels can't survive a matrix)
  expect_null(cap(matrix(1:6, nrow = 3)))
  # Non-scalar / empty / NA labels are ignored
  d <- data.frame(a = 1:3, b = 4:6, e = 7:9)
  attr(d$a, "label") <- ""
  attr(d$b, "label") <- NA_character_
  attr(d$e, "label") <- c("two", "labels")
  expect_null(cap(d))
})

test_that("ackwards() captures data.frame variable labels into meta$item_labels", {
  d <- as.data.frame(na.omit(bfi25)[, 1:6])
  attr(d$A1, "label") <- "Am indifferent to the feelings of others"
  attr(d$A3, "label") <- "Know how to comfort others"
  x <- cached(ackwards(d, k_max = 3))
  expect_equal(
    x$meta$item_labels,
    c(
      A1 = "Am indifferent to the feelings of others",
      A3 = "Know how to comfort others"
    )
  )
  # No labels -> NULL
  x2 <- cached(ackwards(as.data.frame(na.omit(bfi25)[, 1:6]), k_max = 3))
  expect_null(x2$meta$item_labels)
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

test_that(".factor_model_df() matches ((p-k)^2 - (p+k))/2", {
  df <- ackwards:::.factor_model_df
  expect_equal(df(6L, 3L), 0) # just-identified
  expect_equal(df(6L, 4L), -3) # under-identified
  expect_equal(df(16L, 10L), 5)
  expect_equal(df(16L, 11L), -1)
  # vectorised over k (used by .ledermann_bound)
  expect_equal(df(6L, 1:3), c(9, 4, 0))
})

test_that(".ledermann_bound() returns the largest identifiable factor count", {
  lb <- ackwards:::.ledermann_bound
  expect_equal(lb(6L), 3L) # df(6,3)=0, df(6,4)<0
  expect_equal(lb(16L), 10L) # df(16,10)=5, df(16,11)<0
  expect_equal(lb(25L), 18L) # df(25,18)=3, df(25,19)=-4
  # Degenerate small p: no latent model is identified
  expect_equal(lb(3L), 1L) # df(3,1)=0
  expect_equal(lb(2L), 0L)
  expect_equal(lb(1L), 0L)
})

test_that(".as_numeric_matrix() coerces or aborts with the shared messages", {
  anm <- ackwards:::.as_numeric_matrix
  # data frame -> numeric matrix
  expect_identical(anm(data.frame(a = 1:2, b = 3:4)), matrix(c(1L, 2L, 3L, 4L),
    nrow = 2, dimnames = list(NULL, c("a", "b"))
  ))
  # numeric matrix passes through
  m <- matrix(1:6, 2)
  expect_equal(anm(m), m)
  # non-data/non-matrix input
  expect_error(anm(list(a = 1)), "data frame or numeric matrix")
  expect_error(anm(1:5), "data frame or numeric matrix")
  # non-numeric columns
  expect_error(anm(data.frame(a = 1:2, b = c("x", "y"))), "only numeric columns")
  # arg name is threaded into the message
  expect_error(anm(list(), arg = "foo"), "foo")
})

test_that(".check_count() validates a single positive integer and returns int", {
  cc <- ackwards:::.check_count
  expect_identical(cc(5, "n"), 5L)
  expect_identical(cc(5L, "n"), 5L)
  # the drift-bug case: numeric NA must give the intended message, not the
  # base-R "missing value where TRUE/FALSE needed"
  expect_error(cc(NA_real_, "n_obs"), "positive integer")
  expect_error(cc(NA, "n_obs"), "positive integer")
  expect_error(cc(0, "n"), "positive integer")
  expect_error(cc(-1, "n"), "positive integer")
  expect_error(cc(2.5, "n"), "positive integer")
  expect_error(cc(c(1, 2), "n"), "positive integer")
  expect_error(cc("3", "n"), "positive integer")
  # min > 1 switches to the ">= min" message (n_boot uses 2)
  expect_identical(cc(2, "n_boot", min = 2L), 2L)
  expect_error(cc(1, "n_boot", min = 2L), "integer >= 2")
  expect_error(cc(NA_real_, "n_boot", min = 2L), "integer >= 2")
})
