test_that("print.ackwards runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  # cli writes to stderr in non-interactive mode; just verify no error + invisible
  expect_no_error(print(x))
  expect_invisible(print(x))
})

test_that("tidy(x, what = 'edges') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  ed <- generics::tidy(x, what = "edges")
  expect_s3_class(ed, "data.frame")
  expect_true(all(c("from", "to", "level_from", "level_to", "r", "is_primary", "above_cut") %in% names(ed)))
  # 3 levels → 2 adjacent pairs → (1*2) + (2*3) = 8 edges
  expect_equal(nrow(ed), 8L)
  expect_true(all(abs(ed$r) <= 1 + 1e-9))
})

test_that("tidy(x, what = 'loadings') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 2))
  ld <- generics::tidy(x, what = "loadings")
  expect_s3_class(ld, "data.frame")
  expect_true(all(c("level", "factor", "item", "loading") %in% names(ld)))
  # k=1 has 25 rows, k=2 has 50 rows → 75 total
  expect_equal(nrow(ld), 75L)
})

test_that("tidy(x, what = 'variance') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  vr <- generics::tidy(x, what = "variance")
  expect_s3_class(vr, "data.frame")
  expect_true(all(c("level", "factor", "variance_pct", "cumulative_pct") %in% names(vr)))
  # k=1 has 1 row, k=2 has 2, k=3 has 3 → 6 total
  expect_equal(nrow(vr), 6L)
})

test_that("glance.ackwards returns a one-row data frame with expected columns", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  gl <- generics::glance(x)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c(
    "method", "rotation", "cor_type", "k_max", "n_obs",
    "deepest_converged", "n_edges"
  ) %in% names(gl)))
  expect_equal(gl$k_max, 3L)
  expect_equal(gl$n_obs, 2800L)
  expect_equal(gl$n_edges, 8L)
})
