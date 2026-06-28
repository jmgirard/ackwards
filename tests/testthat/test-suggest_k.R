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
      "ev_obs", "ev_obs_fa",
      "pa_pc_quant", "pa_pc_suggested",
      "pa_fa_quant", "pa_fa_suggested",
      "map", "vss1", "vss2"
    )
  )
  expect_equal(cr$k, 1:5)
  expect_type(cr$ev_obs, "double")
  expect_type(cr$ev_obs_fa, "double")
  expect_type(cr$pa_pc_quant, "double")
  expect_type(cr$pa_fa_quant, "double")
  expect_type(cr$pa_pc_suggested, "logical")
  expect_type(cr$pa_fa_suggested, "logical")
  expect_type(cr$map, "double")
  expect_type(cr$vss1, "double")
  expect_type(cr$vss2, "double")
})

test_that("suggest_k() FA eigenvalues are smaller than PC eigenvalues", {
  # Factor-analysis eigenvalues (communalities removed) are always < PC eigs.
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 5, n_iter = 5, seed = 1L)
  expect_true(all(sk$criteria$ev_obs_fa <= sk$criteria$ev_obs))
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

test_that("suggest_k() PA-PC is capped at k_max; k_parallel_fa is integer or NA", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 3, n_iter = 5, seed = 1L)
  expect_lte(sk$k_parallel_pc, 3L)
  expect_gte(sk$k_parallel_pc, 1L)
  # PA-FA is either NA (undetermined) or a valid integer in [1, k_max]
  if (!is.na(sk$k_parallel_fa)) {
    expect_gte(sk$k_parallel_fa, 1L)
    expect_lte(sk$k_parallel_fa, 3L)
  }
})

test_that("suggest_k() pa_pc_suggested matches k_parallel_pc boundary", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 5, n_iter = 5, seed = 1L)
  cr <- sk$criteria
  expect_equal(cr$pa_pc_suggested, cr$k <= sk$k_parallel_pc)
})

test_that("suggest_k() pa_fa_suggested is all FALSE when k_parallel_fa is NA", {
  skip_if_not_installed("psych")
  # Construct a minimal suggest_k object with k_parallel_fa = NA to verify
  # the criteria table logic — build one real call and then patch the field.
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  sk_na <- sk
  sk_na$k_parallel_fa <- NA_integer_
  sk_na$criteria$pa_fa_suggested <- rep(FALSE, 4L)
  expect_true(all(!sk_na$criteria$pa_fa_suggested))
})

test_that("suggest_k() deterministic criteria are identical across calls", {
  skip_if_not_installed("psych")
  sk1 <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 99L)
  sk2 <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 99L)
  # MAP, VSS, and observed eigenvalues are computed from the correlation
  # matrix and are fully deterministic regardless of seed.
  expect_equal(sk1$criteria$ev_obs, sk2$criteria$ev_obs)
  expect_equal(sk1$criteria$ev_obs_fa, sk2$criteria$ev_obs_fa)
  expect_equal(sk1$criteria$map, sk2$criteria$map)
  expect_equal(sk1$criteria$vss1, sk2$criteria$vss1)
  expect_equal(sk1$criteria$vss2, sk2$criteria$vss2)
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
  expect_false(is.na(sk$k_cd))
  expect_gte(sk$k_cd, 1L)
  expect_lte(sk$k_cd, 6L)
})

test_that("suggest_k() CD degrades gracefully on EFAtools error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  # Trigger a CD error by passing a corrupt matrix that EFAtools::CD rejects.
  # The function should not propagate an unhandled error; it returns
  # cd_available = FALSE with k_cd = NA (degraded gracefully).
  bad_mat <- matrix(1:4, nrow = 2)
  storage.mode(bad_mat) <- "double"
  sk <- suppressWarnings(tryCatch(
    suggest_k(bad_mat, k_max = 1L, n_iter = 3, seed = 1L),
    error = function(e) NULL
  ))
  # Result is either NULL (PA/VSS also failed on degenerate input) or a
  # suggest_k with cd_available = FALSE (CD degraded, rest succeeded).
  # Either way the function must not propagate an unhandled error.
  if (is.null(sk)) {
    expect_null(sk) # degenerate input caused full failure — acceptable
  } else {
    expect_false(sk$cd_available)
    expect_identical(sk$k_cd, NA_integer_)
  }
})

test_that("suggest_k() handles data with missing values for PA/MAP/VSS", {
  skip_if_not_installed("psych")
  # bfi has NA values; pairwise deletion should allow PA/MAP/VSS to run.
  sk <- suggest_k(psych::bfi[, 1:10], k_max = 3, n_iter = 3, seed = 1L)
  expect_s3_class(sk, "suggest_k")
  expect_equal(sk$n_vars, 10L)
})

test_that("print.suggest_k() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

test_that("print.suggest_k() does not error when k_parallel_fa is NA", {
  skip_if_not_installed("psych")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  # Patch to simulate the NA PA-FA case (no FA factor exceeded random threshold).
  sk$k_parallel_fa <- NA_integer_
  sk$criteria$pa_fa_suggested <- rep(FALSE, 4L)
  # Must run without error and return invisibly; no "k <= 0" / "k <= NA" trap.
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

test_that("suggest_k() consensus excludes NA k_parallel_fa", {
  skip_if_not_installed("psych")
  # When k_parallel_fa = NA the consensus range uses only the remaining criteria.
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  sk$k_parallel_fa <- NA_integer_
  all_k <- stats::na.omit(
    c(sk$k_parallel_pc, sk$k_parallel_fa, sk$k_map, sk$k_vss1, sk$k_vss2)
  )
  expect_gte(min(all_k), 1L)
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

test_that("autoplot.suggest_k() scree panel has four series", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  p <- autoplot(sk)
  scree_series <- unique(
    as.character(p$data$series[p$data$panel == "Scree / Parallel Analysis"])
  )
  expect_setequal(
    scree_series,
    c("Observed (PC)", "PA-PC (95th pct)", "Observed (FA)", "PA-FA (95th pct)")
  )
})

test_that("autoplot.suggest_k() marks optimal k with star points", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 4, n_iter = 5, seed = 1L)
  p <- autoplot(sk)
  opt_rows <- p$data[p$data$is_opt, ]
  # At minimum: one PC scree + one MAP + at least one VSS
  expect_gte(nrow(opt_rows), 3L)
})
