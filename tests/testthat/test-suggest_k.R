test_that("suggest_k() returns a suggest_k object with expected fields", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)

  expect_s3_class(sk, "suggest_k")
  expect_named(
    sk,
    c(
      "k_parallel_pc", "k_parallel_fa",
      "k_map", "k_vss1", "k_vss2",
      "k_cd", "cd_available",
      "criteria",
      "k_max", "n_obs", "n_vars", "cor"
    )
  )
  expect_equal(sk$k_max, 4L)
  expect_equal(sk$n_obs, 2800L)
  expect_equal(sk$n_vars, 25L)
  expect_equal(sk$cor, "pearson")
  expect_type(sk$cd_available, "logical")
})

test_that("suggest_k() criteria table has correct structure and row count", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 5, n_iter = 5, seed = 1L)
  cr <- sk$criteria

  expect_s3_class(cr, "data.frame")
  expect_equal(nrow(cr), 5L)
  expect_named(
    cr,
    c(
      "k",
      "ev_obs",
      "pa_pc_quant", "pa_pc_suggested",
      "pa_fa_quant", "pa_fa_suggested",
      "map", "vss1", "vss2"
    )
  )
  expect_equal(cr$k, 1:5)
  expect_type(cr$ev_obs, "double")
  expect_type(cr$pa_pc_quant, "double")
  expect_type(cr$pa_fa_quant, "double")
  expect_type(cr$pa_pc_suggested, "logical")
  expect_type(cr$pa_fa_suggested, "logical")
  expect_type(cr$map, "double")
  expect_type(cr$vss1, "double")
  expect_type(cr$vss2, "double")
})

test_that("suggest_k() MAP optimal is the row minimising map", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 7, n_iter = 5, seed = 1L)
  expect_gte(sk$k_map, 1L)
  expect_lte(sk$k_map, 7L)
  expect_equal(sk$k_map, which.min(sk$criteria$map))
})

test_that("suggest_k() VSS optima match which.max of vss columns", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 7, n_iter = 5, seed = 1L)
  expect_equal(sk$k_vss1, which.max(sk$criteria$vss1))
  expect_equal(sk$k_vss2, which.max(sk$criteria$vss2))
})

test_that("suggest_k() PA-PC and PA-FA suggestions are capped at k_max", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 3, n_iter = 5, seed = 1L)
  expect_lte(sk$k_parallel_pc, 3L)
  expect_lte(sk$k_parallel_fa, 3L)
  expect_gte(sk$k_parallel_pc, 1L)
  # PA-FA may legitimately be 0 on small k_max; allow >= 0
  expect_gte(sk$k_parallel_fa, 0L)
})

test_that("suggest_k() pa_pc_suggested matches k_parallel_pc boundary", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 5, n_iter = 5, seed = 1L)
  cr <- sk$criteria
  expect_equal(cr$pa_pc_suggested, cr$k <= sk$k_parallel_pc)
})

test_that("suggest_k() seed makes deterministic outputs identical", {
  skip_if_not_installed("psych")
  sk1 <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 99L)
  sk2 <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 99L)
  # Observed eigenvalues, MAP, and VSS are computed from the correlation
  # matrix and are fully deterministic regardless of seed.
  expect_equal(sk1$criteria$ev_obs, sk2$criteria$ev_obs)
  expect_equal(sk1$criteria$map, sk2$criteria$map)
  expect_equal(sk1$criteria$vss1, sk2$criteria$vss1)
  expect_equal(sk1$criteria$vss2, sk2$criteria$vss2)
  # Suggested k from PA should also be stable across repeated calls.
  expect_equal(sk1$k_parallel_pc, sk2$k_parallel_pc)
  expect_equal(sk1$k_parallel_fa, sk2$k_parallel_fa)
})

test_that("suggest_k() errors on bad inputs", {
  skip_if_not_installed("psych")
  expect_error(suggest_k(list()), "data frame")
  expect_error(suggest_k(psych::bfi[, 1:5], k_max = 100), "k_max")
  expect_error(suggest_k(psych::bfi[, 1:5], k_max = 0), "k_max")
})

test_that("suggest_k() k_cd is NA when EFAtools is not installed", {
  skip_if_not_installed("psych")
  skip_if(rlang::is_installed("EFAtools"), "EFAtools is installed; skipping absence test")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  expect_false(sk$cd_available)
  expect_identical(sk$k_cd, NA_integer_)
})

test_that("suggest_k() k_cd is within range when EFAtools is installed", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 6, n_iter = 5, seed = 1L)
  expect_true(sk$cd_available)
  expect_true(!is.na(sk$k_cd))
  expect_gte(sk$k_cd, 1L)
  expect_lte(sk$k_cd, 6L)
})

test_that("print.suggest_k() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

test_that("autoplot.suggest_k() returns a ggplot with three facet panels", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  panel_vals <- levels(p$data$panel)
  expect_equal(
    panel_vals,
    c("Scree / Parallel Analysis", "MAP (minimize)", "VSS (maximize)")
  )
})

test_that("autoplot.suggest_k() marks optimal k with star points", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  p <- autoplot(sk)
  # There should be exactly one row marked is_opt per criterion
  opt_rows <- p$data[p$data$is_opt, ]
  # At minimum: one MAP + at least one VSS + one scree
  expect_gte(nrow(opt_rows), 3L)
})
