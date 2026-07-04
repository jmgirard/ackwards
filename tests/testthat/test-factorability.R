# factorability(): pre-analysis KMO / Bartlett / N:p / Ledermann screen, plus
# the shared internal screen ackwards() runs (Ledermann + adequacy warnings).

test_that("factorability() on raw data returns the expected shape", {
  f <- factorability(sim16)
  expect_s3_class(f, "factorability")
  expect_named(
    f,
    c(
      "cor", "n_obs", "n_vars", "np_ratio", "kmo_overall",
      "kmo_items", "bartlett", "ledermann"
    ),
    ignore.order = TRUE
  )
  expect_identical(f$cor, "pearson")
  expect_identical(f$n_obs, nrow(sim16))
  expect_identical(f$n_vars, ncol(sim16))
  expect_equal(f$np_ratio, nrow(sim16) / ncol(sim16))
  expect_identical(f$ledermann, 10L) # p = 16
  expect_equal(nrow(f$kmo_items), ncol(sim16))
  expect_named(f$kmo_items, c("item", "msa"))
})

test_that("factorability() KMO and Bartlett match psych on a fixed matrix", {
  R <- cor(sim16)
  n <- nrow(sim16)
  f <- factorability(R, n_obs = n)

  k <- psych::KMO(R)
  expect_equal(f$kmo_overall, unname(k$MSA))
  expect_equal(f$kmo_items$msa, unname(k$MSAi))

  b <- psych::cortest.bartlett(R, n = n)
  expect_equal(f$bartlett$chisq, unname(b$chisq))
  expect_equal(f$bartlett$df, unname(b$df))
  expect_equal(f$bartlett$p_value, unname(b$p.value))
})

test_that("N-based diagnostics are NA/NULL for a matrix without n_obs", {
  R <- cor(sim16)
  f <- factorability(R)
  expect_true(is.na(f$n_obs))
  expect_true(is.na(f$np_ratio))
  expect_null(f$bartlett)
  # KMO and Ledermann do not need N and are still reported
  expect_false(is.na(f$kmo_overall))
  expect_identical(f$ledermann, 10L)
})

test_that("factorability() honours the correlation basis", {
  f_p <- factorability(sim16, cor = "pearson")
  f_s <- factorability(sim16, cor = "spearman")
  expect_identical(f_s$cor, "spearman")
  # Spearman and Pearson KMO differ (different R)
  expect_false(isTRUE(all.equal(f_p$kmo_overall, f_s$kmo_overall)))
})

test_that("factorability() errors on a constant item and warns on ignored n_obs", {
  d <- sim16
  d$dead <- 3
  expect_error(factorability(d), "no variance")

  expect_warning(factorability(sim16, n_obs = 999), "ignored for raw data")
})

test_that("factorability() validates n_obs and rejects too-few columns", {
  R <- cor(sim16)
  expect_error(factorability(R, n_obs = -5), "positive integer")
  expect_error(factorability(R, n_obs = 2.5), "positive integer")
  expect_error(factorability(sim16[, 1, drop = FALSE]), "at least two columns")
})

test_that("print.factorability() renders without error and returns invisibly", {
  f <- factorability(cor(sim16), n_obs = nrow(sim16))
  expect_no_error(print(f))
  expect_invisible(print(f))
  # cli output is captured on the message connection (see test-print.R)
  txt <- cli::ansi_strip(capture.output(print(f), type = "message"))
  expect_true(any(grepl("Factorability screen", txt)))
  expect_true(any(grepl("Ledermann bound", txt)))
  expect_true(any(grepl("Bartlett", txt)))
  # No-N variant prints the 'not supplied' path
  txt2 <- cli::ansi_strip(capture.output(print(factorability(cor(sim16))), type = "message"))
  expect_true(any(grepl("not supplied", txt2)))
})

test_that(".kmo_band() maps values to Kaiser bands", {
  band <- ackwards:::.kmo_band
  expect_identical(band(0.95), "marvellous")
  expect_identical(band(0.85), "meritorious")
  expect_identical(band(0.75), "middling")
  expect_identical(band(0.65), "mediocre")
  expect_identical(band(0.55), "miserable")
  expect_identical(band(0.45), "unacceptable")
  expect_identical(band(NA_real_), "uncomputable")
})

# --- Internal screen wired into ackwards() -----------------------------------

test_that("ackwards() warns when k_max exceeds the Ledermann bound (EFA only, not PCA)", {
  R <- cor(sim16) # p = 16, Ledermann bound = 10

  rlang::reset_warning_verbosity("ackwards_ledermann")
  w <- suppressMessages(testthat::capture_warnings(
    ackwards(R, k_max = 13L, engine = "efa", n_obs = 1000L)
  ))
  expect_true(any(grepl("Ledermann bound", w)))
  expect_true(any(grepl("at most 10 common factor", w)))

  # PCA at the same k_max: components have no identification limit -> no warning
  rlang::reset_warning_verbosity("ackwards_ledermann")
  w_pca <- suppressMessages(testthat::capture_warnings(
    ackwards(R, k_max = 13L, engine = "pca", n_obs = 1000L)
  ))
  expect_false(any(grepl("Ledermann", w_pca)))
})

test_that("ackwards() does not warn about Ledermann when k_max is within the bound", {
  R <- cor(sim16)
  rlang::reset_warning_verbosity("ackwards_ledermann")
  w <- suppressMessages(testthat::capture_warnings(
    ackwards(R, k_max = 5L, engine = "efa", n_obs = 1000L)
  ))
  expect_false(any(grepl("Ledermann", w)))
})

test_that("ackwards() emits a consolidated adequacy warning below thresholds, silent when adequate", {
  set.seed(42)
  d <- as.data.frame(matrix(stats::rnorm(40 * 10), 40, 10)) # N:p = 4, N < 100
  names(d) <- paste0("v", seq_len(10))

  rlang::reset_warning_verbosity("ackwards_factorability")
  w <- suppressMessages(testthat::capture_warnings(
    ackwards(d, k_max = 3L, engine = "pca")
  ))
  expect_true(any(grepl("poorly suited to factor analysis", w)))
  expect_true(any(grepl("N:p", w)))

  # Well-powered, well-structured data: no adequacy warning
  rlang::reset_warning_verbosity("ackwards_factorability")
  w_ok <- suppressMessages(testthat::capture_warnings(
    ackwards(sim16, k_max = 4L, engine = "pca")
  ))
  expect_false(any(grepl("poorly suited", w_ok)))
})

test_that("the internal screen skips N-based checks for a correlation matrix without n_obs", {
  R <- cor(sim16)
  rlang::reset_warning_verbosity("ackwards_factorability")
  # PCA cor-matrix input, no n_obs -> n_obs_eff is NA; N:p / Bartlett skipped,
  # no error from the screen.
  expect_no_error(suppressMessages(ackwards(R, k_max = 3L, engine = "pca")))
})

# --- Coverage of input variants and print branches ---------------------------

test_that("factorability() accepts a raw numeric matrix and rejects bad input", {
  expect_s3_class(factorability(as.matrix(sim16)), "factorability")
  expect_error(factorability("not data"), "data frame or numeric matrix")
  d <- sim16
  d$chr <- "a" # character column -> non-numeric matrix
  expect_error(factorability(d), "only numeric columns")
})

test_that("factorability() computes on the polychoric basis", {
  skip_if_not_installed("psych")
  f <- factorability(bfi25[, 1:6], cor = "polychoric")
  expect_identical(f$cor, "polychoric")
  expect_false(is.na(f$kmo_overall))
  # matches a direct polychoric KMO
  R <- suppressWarnings(psych::polychoric(bfi25[, 1:6])$rho)
  expect_equal(f$kmo_overall, unname(psych::KMO(R)$MSA))
})

test_that("print.factorability() shows low-MSA items and a non-significant Bartlett", {
  set.seed(7)
  R <- cor(matrix(stats::rnorm(45 * 6), 45, 6)) # weak, small N
  colnames(R) <- rownames(R) <- paste0("v", 1:6)
  f <- factorability(R, n_obs = 45)
  # collapse to one string: cli may wrap a phrase across output lines
  txt <- paste(cli::ansi_strip(capture.output(print(f), type = "message")), collapse = " ")
  expect_match(txt, "MSA < .50", fixed = TRUE)
  expect_match(txt, "NOT significant", fixed = TRUE)
})
