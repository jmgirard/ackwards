test_that("suggest_k() returns a suggest_k object with expected fields", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5)

  expect_s3_class(sk, "suggest_k")
  expect_named(
    sk,
    c("k_parallel", "k_map", "criteria", "k_max", "n_obs", "n_vars", "cor")
  )
  expect_equal(sk$k_max,   4L)
  expect_equal(sk$n_obs,  2800L)
  expect_equal(sk$n_vars,   25L)
  expect_equal(sk$cor,   "pearson")
})

test_that("suggest_k() criteria table has correct structure and row count", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 5, n_iter = 5)
  cr <- sk$criteria

  expect_s3_class(cr, "data.frame")
  expect_equal(nrow(cr), 5L)
  expect_true(all(c("k", "map", "pa_suggested") %in% names(cr)))
  expect_equal(cr$k, 1:5)
  expect_type(cr$map, "double")
  expect_type(cr$pa_suggested, "logical")
})

test_that("suggest_k() MAP values decrease then increase (optimal in middle)", {
  skip_if_not_installed("psych")
  # MAP should have a minimum somewhere in the range for bfi data
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 7, n_iter = 5)
  # k_map should be between 1 and k_max
  expect_gte(sk$k_map, 1L)
  expect_lte(sk$k_map, 7L)
  # The MAP at k_map should be <= all other MAP values
  expect_true(all(sk$criteria$map[sk$k_map] <= sk$criteria$map))
})

test_that("suggest_k() k_parallel is capped at k_max", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 3, n_iter = 5)
  expect_lte(sk$k_parallel, 3L)
})

test_that("suggest_k() errors on bad inputs", {
  skip_if_not_installed("psych")
  expect_error(suggest_k(list()),                         "data frame")
  expect_error(suggest_k(psych::bfi[, 1:5], k_max = 100), "k_max")
  expect_error(suggest_k(psych::bfi[, 1:5], k_max = 0),   "k_max")
})

test_that("print.suggest_k() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5)
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})
