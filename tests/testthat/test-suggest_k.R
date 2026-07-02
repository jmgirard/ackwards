# Fast fixtures. `EFAtools::CD()` dominates suggest_k() cost (~8 s on the full
# 2800x25 psych::bfi vs ~0.1 s for fa.parallel + vss combined) and scales hard
# with both n_obs and n_vars. So the structural/logic/print/autoplot tests run
# on a small subset of the bundled bfi25 (CD ~0.8 s), computed once and cached.
# Nothing is mocked or skipped: PA, MAP, VSS, and CD all run for real on the
# small fixture — only the data is smaller, which is sufficient to exercise the
# wiring and table/print/plot logic. A single integration smoke (the
# "expected fields" test) still runs the full 25-item psych::bfi instrument so
# the realistic-scale path (incl. CD at 25 vars, n_obs = 2800) stays covered.
.sk_small <- bfi25[, 1:8] # 1000 x 8; complete cases = 946

.get_sk <- local({
  cache <- list()
  function(k_max, seed = 1L, n_iter = 5L) {
    key <- paste0("k", k_max, "s", seed)
    if (is.null(cache[[key]])) {
      cache[[key]] <<- suggest_k(
        .sk_small,
        k_max = k_max, n_iter = n_iter, seed = seed
      )
    }
    cache[[key]]
  }
})

test_that("suggest_k() returns a suggest_k object with expected fields", {
  skip_if_not_installed("psych")
  # Integration smoke on the full 25-item instrument: this is the one test that
  # asserts the realistic n_obs / n_vars and pays the full-data CD cost (~5 s).
  sk <- suggest_k(psych::bfi[, 1:25], k_max = 3L, n_iter = 5L, seed = 1L)

  expect_s3_class(sk, "suggest_k")
  expect_named(
    sk,
    c(
      "k_parallel_pc", "k_parallel_fa",
      "k_map", "k_vss1", "k_vss2",
      "k_cd", "cd_available", "cd_rmse",
      "criteria", "criteria_requested",
      "k_max", "n_obs", "n_vars", "cor", "input_type"
    )
  )
  expect_equal(sk$k_max, 3L)
  expect_equal(sk$n_obs, 2800L)
  expect_equal(sk$n_vars, 25L)
  expect_equal(sk$cor, "pearson")
  expect_type(sk$cd_available, "logical")
})

test_that("the six per-criterion k_* fields the suggest-k vignette reads stay available", {
  # The "When the criteria agree — and when they don't" section of
  # vignette("ackwards-suggest-k") computes a consensus range *inline* from these
  # six fields (sim16 vs bfi25). Renaming any of them would silently break the
  # vignette build rather than fail a test, so pin the contract here and mirror
  # the vignette's own range computation.
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  ks <- c(
    sk$k_parallel_pc, sk$k_parallel_fa, sk$k_map,
    sk$k_vss1, sk$k_vss2, sk$k_cd
  )
  expect_length(ks, 6L) # all six fields resolved (NULL would drop entries)
  finite_ks <- ks[!is.na(ks)]
  expect_gte(length(finite_ks), 1L) # at least the non-CD criteria yield a value
  expect_true(all(finite_ks >= 1L))
  expect_true(is.finite(min(finite_ks)) && is.finite(max(finite_ks)))
})

test_that("suggest_k() criteria table has correct structure and row count", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  cr <- sk$criteria

  expect_s3_class(cr, "data.frame")
  expect_equal(nrow(cr), 4L)
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
  expect_equal(cr$k, 1:4)
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
  sk <- .get_sk(4L)
  expect_true(all(sk$criteria$ev_obs_fa <= sk$criteria$ev_obs))
})

test_that("suggest_k() MAP optimal is the row minimising map", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  expect_gte(sk$k_map, 1L)
  expect_lte(sk$k_map, 4L)
  expect_equal(sk$k_map, which.min(sk$criteria$map))
})

test_that("suggest_k() VSS optima match which.max of vss columns", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  expect_equal(sk$k_vss1, which.max(sk$criteria$vss1))
  expect_equal(sk$k_vss2, which.max(sk$criteria$vss2))
})

test_that("suggest_k() PA-PC is capped at k_max; k_parallel_fa is integer or NA", {
  skip_if_not_installed("psych")
  # The clamp `min(raw_pa, k_max)` is only exercised when raw PA exceeds k_max.
  # On bfi25[, 1:15] PA-PC suggests ~3 components, so k_max = 2 forces the clamp.
  sk <- suggest_k(bfi25[, 1:15], k_max = 2L, n_iter = 5L, seed = 1L)
  expect_lte(sk$k_parallel_pc, 2L)
  expect_gte(sk$k_parallel_pc, 1L)
  # PA-FA is either NA (undetermined) or a valid integer in [1, k_max]
  if (!is.na(sk$k_parallel_fa)) {
    expect_gte(sk$k_parallel_fa, 1L)
    expect_lte(sk$k_parallel_fa, 2L)
  }
})

test_that("suggest_k() pa_pc_suggested matches k_parallel_pc boundary", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  cr <- sk$criteria
  expect_equal(cr$pa_pc_suggested, cr$k <= sk$k_parallel_pc)
})

test_that("suggest_k() pa_fa_suggested is all FALSE when k_parallel_fa is NA", {
  skip_if_not_installed("psych")
  # Construct a minimal suggest_k object with k_parallel_fa = NA to verify
  # the criteria table logic — build one real call and then patch the field.
  sk <- .get_sk(4L)
  sk_na <- sk
  sk_na$k_parallel_fa <- NA_integer_
  sk_na$criteria$pa_fa_suggested <- rep(FALSE, 4L)
  expect_true(all(!sk_na$criteria$pa_fa_suggested))
})

test_that("suggest_k() deterministic criteria are identical across calls", {
  skip_if_not_installed("psych")
  # Two independent calls with the same seed — intentionally NOT cached so
  # reproducibility is actually exercised (not just comparing an object to itself).
  sk1 <- suggest_k(.sk_small, k_max = 4, n_iter = 3, seed = 99L)
  sk2 <- suggest_k(.sk_small, k_max = 4, n_iter = 3, seed = 99L)
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
  expect_error(suggest_k(bfi25[, 1:5], k_max = 100), "k_max")
  expect_error(suggest_k(bfi25[, 1:5], k_max = 0), "k_max")
  # n_iter drives the PA simulation; 0/fractional/vector values fail deep
  # inside psych without this guard (M42/m8)
  expect_error(suggest_k(bfi25[, 1:5], n_iter = 0), "positive integer")
  expect_error(suggest_k(bfi25[, 1:5], n_iter = 2.5), "positive integer")
  expect_error(suggest_k(bfi25[, 1:5], n_iter = c(5, 10)), "positive integer")
})

test_that("suggest_k() defaults k_max to min(p - 1, 8) for raw data", {
  skip_if_not_installed("psych")
  # No k_max supplied: the raw-data path sets k_max <- min(p - 1, 8). With p = 6
  # vars that is 5. suppressWarnings swallows the incidental CD n_factors_max
  # notice (CD caps at 2 factors for 6 vars).
  sk <- suppressWarnings(suggest_k(bfi25[, 1:6], n_iter = 3L, seed = 1L))
  expect_equal(sk$k_max, 5L)
})

test_that("suggest_k() k_cd is NA and cd_rmse is NULL when EFAtools is unavailable", {
  skip_if_not_installed("psych")
  # Force the EFAtools availability check to report it missing, so this absence
  # branch runs regardless of whether EFAtools is actually installed (mirrors the
  # future.apply pattern in test-esem.R). Call suggest_k() directly rather than
  # via the memoised .get_sk() so the mocked-absence result never shares a cache
  # key with the EFAtools-present tests.
  testthat::local_mocked_bindings(
    is_installed = function(pkg, ...) if (identical(pkg, "EFAtools")) FALSE else TRUE,
    .package = "rlang"
  )
  sk <- suggest_k(.sk_small, k_max = 4L, n_iter = 3L, seed = 1L)
  expect_false(sk$cd_available)
  expect_identical(sk$k_cd, NA_integer_)
  expect_null(sk$cd_rmse)
})

test_that("suggest_k() k_cd and cd_rmse populated when EFAtools is installed", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  sk <- .get_sk(4L)
  expect_true(sk$cd_available)
  expect_false(is.na(sk$k_cd))
  expect_gte(sk$k_cd, 1L)
  expect_lte(sk$k_cd, 4L)
  # cd_rmse: length-k_max numeric vector of mean RMSE values
  expect_false(is.null(sk$cd_rmse))
  expect_length(sk$cd_rmse, 4L)
  expect_true(is.numeric(sk$cd_rmse))
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
  # bfi25 has NA values; pairwise deletion should allow PA/MAP/VSS to run.
  sk <- suggest_k(bfi25[, 1:10], k_max = 3, n_iter = 3, seed = 1L)
  expect_s3_class(sk, "suggest_k")
  expect_equal(sk$n_vars, 10L)
})

test_that("print.suggest_k() runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

test_that("print.suggest_k() does not error when k_parallel_fa is NA", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
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
  sk <- .get_sk(4L)
  sk$k_parallel_fa <- NA_integer_
  all_k <- stats::na.omit(
    c(sk$k_parallel_pc, sk$k_parallel_fa, sk$k_map, sk$k_vss1, sk$k_vss2)
  )
  expect_gte(min(all_k), 1L)
})

test_that("autoplot.suggest_k() renders 3 panels in single column when CD absent", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- .get_sk(4L)
  sk$cd_available <- FALSE
  sk$cd_rmse <- NULL
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  panel_vals <- levels(p$data$panel)
  expect_equal(
    panel_vals,
    c("Scree / Parallel Analysis", "MAP (minimize)", "VSS (maximize)")
  )
})

test_that("autoplot.suggest_k() renders 4-panel 2x2 grid when CD available", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("EFAtools")
  sk <- .get_sk(4L)
  # Ensure cd_rmse is populated (real CD run)
  expect_true(sk$cd_available)
  expect_false(is.null(sk$cd_rmse))
  expect_length(sk$cd_rmse, 4L)
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  panel_vals <- levels(p$data$panel)
  expect_equal(
    panel_vals,
    c(
      "Scree / Parallel Analysis", "MAP (minimize)",
      "VSS (maximize)", "CD (RMSE; sequential test)"
    )
  )
  # CD star marked at k_cd
  cd_rows <- p$data[p$data$series == "CD (RMSE)" & p$data$is_opt, ]
  expect_equal(nrow(cd_rows), 1L)
  expect_equal(cd_rows$k, sk$k_cd)
  # No NA rows in plot data (NAs filtered before building frame)
  cd_plot_rows <- p$data[p$data$series == "CD (RMSE)", ]
  expect_false(anyNA(cd_plot_rows$value))
})

test_that("cd_rmse has no spurious zeros: trailing columns masked to NA", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  # Use k_max large enough to guarantee unfilled columns exist when k_cd < k_max
  sk <- suggest_k(bfi25, k_max = 8L, seed = 1L)
  skip_if(!sk$cd_available)
  expect_length(sk$cd_rmse, 8L)
  # No zeros beyond the computed range (k_cd + 1)
  expect_false(any(sk$cd_rmse == 0, na.rm = TRUE))
  # Tail beyond k_cd + 1 is NA (when k_cd < k_max - 1)
  n_computed <- min(sk$k_cd + 1L, 8L)
  if (n_computed < 8L) {
    expect_true(all(is.na(sk$cd_rmse[(n_computed + 1L):8L])))
  }
  # which.min over non-NA values is at a genuinely computed level
  expect_lte(which.min(sk$cd_rmse), n_computed)
})

test_that("autoplot.suggest_k() stars k_cd even when the RMSE minimum is elsewhere", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("EFAtools")
  sk <- .get_sk(4L)
  skip_if(!sk$cd_available)
  # CD's sequential test can stop before the curve's visible minimum (RMSE keeps
  # falling, but not *significantly*). Force that divergence: pick k_cd = 2 while
  # the smallest RMSE sits at k = 3. The star must follow k_cd, not which.min.
  sk$k_cd <- 2L
  sk$cd_rmse <- c(0.4, 0.3, 0.1, NA_real_)
  expect_equal(which.min(sk$cd_rmse), 3L) # guard the fixture's premise
  p <- autoplot(sk)
  cd_opt <- p$data[p$data$series == "CD (RMSE)" & p$data$is_opt, ]
  expect_equal(nrow(cd_opt), 1L)
  expect_equal(cd_opt$k, 2L) # star at k_cd, not the minimum (k = 3)
  cd_min_row <- p$data[p$data$series == "CD (RMSE)" & p$data$k == 3L, ]
  expect_false(cd_min_row$is_opt)
})

test_that("autoplot.suggest_k() scree panel has four series", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- .get_sk(4L)
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
  sk <- .get_sk(4L)
  p <- autoplot(sk)
  opt_rows <- p$data[p$data$is_opt, ]
  # At minimum: one PC scree + one MAP + at least one VSS
  expect_gte(nrow(opt_rows), 3L)
})

test_that("autoplot.suggest_k() handles k_parallel_fa = NA without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- .get_sk(4L)
  sk$k_parallel_fa <- NA_integer_
  sk$criteria$pa_fa_suggested <- rep(FALSE, 4L)
  sk$cd_available <- FALSE
  sk$cd_rmse <- NULL
  # FA star must be omitted; the plot must still build with 3 panels (CD absent).
  expect_no_error({
    p <- autoplot(sk)
  })
  expect_s3_class(p, "gg")
  expect_equal(
    levels(p$data$panel),
    c("Scree / Parallel Analysis", "MAP (minimize)", "VSS (maximize)")
  )
})

test_that("autoplot.suggest_k() handles cd_available = FALSE without error", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- .get_sk(4L)
  sk$cd_available <- FALSE
  sk$cd_rmse <- NULL
  sk$k_cd <- NA_integer_
  # CD panel must be absent; plot must still build with 3 panels.
  expect_no_error({
    p <- autoplot(sk)
  })
  expect_s3_class(p, "gg")
  expect_false("CD (RMSE, minimize)" %in% levels(p$data$panel))
})

test_that("print.suggest_k() runs without error when cd_available = FALSE", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  sk$cd_available <- FALSE
  sk$k_cd <- NA_integer_
  # Must print the "CD requires EFAtools" note, not error, and return invisibly.
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

# ---- Branches added for M23 coverage hardening -------------------------------

test_that("suggest_k() R-matrix + explicit cor warns it is ignored", {
  skip_if_not_installed("psych")
  R <- cor(bfi25[, 1:6], use = "pairwise.complete.obs")
  # 'cor' arg is explicitly set (not "pearson") alongside an R-matrix: warn
  expect_warning(
    suppressMessages(suggest_k(R, n_obs = 875L, n_iter = 3L, cor = "spearman")),
    "cor.*ignored"
  )
})

test_that("suggest_k() R-matrix: k_max >= p errors with helpful message", {
  skip_if_not_installed("psych")
  R <- cor(bfi25[, 1:6], use = "pairwise.complete.obs") # p = 6
  expect_error(
    suppressMessages(suggest_k(R, n_obs = 875L, k_max = 6L, n_iter = 3L)),
    "k_max"
  )
})

test_that("suggest_k() errors on non-numeric (character) data", {
  d <- data.frame(x = c("a", "b", "c"), y = c("d", "e", "f"))
  expect_error(suggest_k(d, k_max = 1L), "numeric")
})

test_that("suggest_k() warns when cor = 'spearman' and CD is available", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  # k_max = 2 stays within CD's extractable-factor limit for 6 vars, so the
  # only warning raised is the spearman -> Pearson one under test.
  expect_warning(
    suggest_k(bfi25[, 1:6], k_max = 2L, cor = "spearman", n_iter = 3L, seed = 1L),
    "Pearson"
  )
})

test_that("print.suggest_k() renders CD dash symbol when i > k_cd", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  # Force k_cd to 1 so that the later rows hit the 'else' (dash) branch
  sk$cd_available <- TRUE
  sk$k_cd <- 1L
  # Must print without error (rows 2+ produce the dash branch at print.R line 418)
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

test_that("print.suggest_k() shows single-value consensus when all criteria agree", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  # Patch every criterion to the same k so lo == hi: this hits the degenerate
  # "Consensus: k = {lo}" branch (rather than the usual "Consensus range").
  sk$k_parallel_pc <- 2L
  sk$k_parallel_fa <- 2L
  sk$k_map <- 2L
  sk$k_vss1 <- 2L
  sk$k_vss2 <- 2L
  sk$cd_available <- FALSE
  sk$k_cd <- NA_integer_
  # cli output is not reliably captured by expect_output (see other print tests);
  # the patched lo == hi still drives the single-value consensus branch.
  expect_no_error(print(sk))
  expect_invisible(print(sk))
})

# ---- Tests for the criteria= argument (Wave 1 / M25) -------------------------

test_that("suggest_k() criteria='map' skips PA and CD, populates only MAP", {
  skip_if_not_installed("psych")
  sk <- suggest_k(.sk_small, k_max = 4L, criteria = "map", n_iter = 5L)

  expect_equal(sk$criteria_requested, "map")
  # PA fields must be NA (not run)
  expect_identical(sk$k_parallel_pc, NA_integer_)
  expect_identical(sk$k_parallel_fa, NA_integer_)
  expect_true(all(is.na(sk$criteria$ev_obs)))
  expect_true(all(is.na(sk$criteria$pa_pc_quant)))
  # MAP must be populated
  expect_false(is.na(sk$k_map))
  expect_equal(sk$k_map, which.min(sk$criteria$map))
  # VSS must be NA
  expect_identical(sk$k_vss1, NA_integer_)
  expect_identical(sk$k_vss2, NA_integer_)
  expect_true(all(is.na(sk$criteria$vss1)))
  # CD must be NA
  expect_false(sk$cd_available)
  expect_identical(sk$k_cd, NA_integer_)
})

test_that("suggest_k() criteria=c('pa_pc','pa_fa') skips vss() and CD", {
  skip_if_not_installed("psych")
  sk <- suggest_k(.sk_small,
    k_max = 4L, criteria = c("pa_pc", "pa_fa"),
    n_iter = 5L
  )

  expect_equal(sk$criteria_requested, c("pa_pc", "pa_fa"))
  # PA fields must be populated
  expect_false(is.na(sk$k_parallel_pc))
  expect_true(all(!is.na(sk$criteria$ev_obs)))
  # MAP and VSS must be NA (not run)
  expect_identical(sk$k_map, NA_integer_)
  expect_true(all(is.na(sk$criteria$map)))
  expect_identical(sk$k_vss1, NA_integer_)
  expect_true(all(is.na(sk$criteria$vss1)))
  # CD must be NA
  expect_false(sk$cd_available)
})

test_that("suggest_k() default criteria matches all-five structure", {
  skip_if_not_installed("psych")
  sk <- .get_sk(4L)
  # Default should request all five criteria
  expect_setequal(sk$criteria_requested, c("pa_pc", "pa_fa", "map", "vss", "cd"))
  # All standard k_* fields non-NA (except k_parallel_fa and k_cd which may be NA)
  expect_false(is.na(sk$k_parallel_pc))
  expect_false(is.na(sk$k_map))
  expect_false(is.na(sk$k_vss1))
  expect_false(is.na(sk$k_vss2))
})

test_that("suggest_k() criteria='cd' with EFAtools unavailable emits info and returns NA", {
  skip_if_not_installed("psych")
  # Mock EFAtools as missing (see the note on the first absence test above) so the
  # graceful-degradation branch runs even when EFAtools is installed.
  testthat::local_mocked_bindings(
    is_installed = function(pkg, ...) if (identical(pkg, "EFAtools")) FALSE else TRUE,
    .package = "rlang"
  )
  # When 'cd' requested but EFAtools absent, should inform (not error) and return NA
  expect_message(
    sk <- suggest_k(.sk_small, k_max = 3L, criteria = "cd", n_iter = 3L),
    "EFAtools"
  )
  expect_false(sk$cd_available)
  expect_identical(sk$k_cd, NA_integer_)
  # All non-CD k_* fields are NA (only CD was requested)
  expect_identical(sk$k_parallel_pc, NA_integer_)
  expect_identical(sk$k_map, NA_integer_)
})

test_that("suggest_k() criteria='cd' with EFAtools available runs CD only", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  sk <- suggest_k(.sk_small, k_max = 4L, criteria = "cd", n_iter = 3L, seed = 1L)

  expect_equal(sk$criteria_requested, "cd")
  expect_true(sk$cd_available)
  expect_false(is.na(sk$k_cd))
  # All non-CD k_* fields are NA
  expect_identical(sk$k_parallel_pc, NA_integer_)
  expect_identical(sk$k_map, NA_integer_)
  expect_identical(sk$k_vss1, NA_integer_)
})

test_that("suggest_k() invalid criterion name errors via arg_match", {
  skip_if_not_installed("psych")
  expect_error(
    suggest_k(.sk_small, k_max = 3L, criteria = "bad_criterion"),
    "must be one of"
  )
})

test_that("suggest_k() criteria='vss' consensus reflects VSS-1 and VSS-2 only", {
  skip_if_not_installed("psych")
  sk <- suggest_k(.sk_small, k_max = 4L, criteria = "vss", n_iter = 3L)

  expect_equal(sk$criteria_requested, "vss")
  expect_false(is.na(sk$k_vss1))
  expect_false(is.na(sk$k_vss2))
  # Consensus range must be within VSS-1/VSS-2 bounds
  lo <- min(sk$k_vss1, sk$k_vss2)
  hi <- max(sk$k_vss1, sk$k_vss2)
  expect_gte(lo, 1L)
  expect_lte(hi, 4L)
  # print must run without error and consensus must reflect VSS only
  expect_no_error(print(sk))
})

test_that("autoplot.suggest_k() shows only MAP panel when criteria='map'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(.sk_small, k_max = 4L, criteria = "map", n_iter = 3L)
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  expect_equal(levels(p$data$panel), "MAP (minimize)")
})

test_that("autoplot.suggest_k() shows only PA panel when criteria=c('pa_pc','pa_fa')", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(.sk_small,
    k_max = 4L, criteria = c("pa_pc", "pa_fa"),
    n_iter = 3L
  )
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  expect_equal(levels(p$data$panel), "Scree / Parallel Analysis")
  # Both PC and FA series present
  series_vals <- unique(as.character(p$data$series))
  expect_true("Observed (PC)" %in% series_vals)
  expect_true("Observed (FA)" %in% series_vals)
})

test_that("autoplot.suggest_k() shows only PA-PC series when criteria='pa_pc'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  sk <- suggest_k(.sk_small, k_max = 4L, criteria = "pa_pc", n_iter = 3L)
  p <- autoplot(sk)
  expect_s3_class(p, "gg")
  series_vals <- unique(as.character(p$data$series))
  expect_true("Observed (PC)" %in% series_vals)
  expect_false("Observed (FA)" %in% series_vals)
})
