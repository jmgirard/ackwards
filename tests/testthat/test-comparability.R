# Tests for comparability() -- split-half factor comparability (M46).

# Shared fixture: one seeded PCA run on sim16, reused across structural,
# discriminating, and print/plot tests (memoised to keep the suite fast).
.cmp_cache <- new.env(parent = emptyenv())
.get_cmp <- function() {
  if (is.null(.cmp_cache$cmp)) {
    .cmp_cache$cmp <- suppressMessages(
      comparability(sim16, k_max = 5, n_splits = 4, seed = 1)
    )
  }
  .cmp_cache$cmp
}

# ---- Object structure --------------------------------------------------------

test_that("comparability() returns a well-formed object", {
  cmp <- .get_cmp()

  expect_s3_class(cmp, "comparability")
  expect_named(
    cmp,
    c(
      "coefficients", "summary", "k_max", "k_requested", "n_splits", "n_half",
      "engine", "cor", "fm", "n_obs", "n_vars", "seed", "call"
    )
  )
  expect_equal(cmp$k_max, 5L)
  expect_equal(cmp$k_requested, 5L)
  expect_equal(cmp$n_splits, 4L)
  expect_equal(cmp$n_half, 500L)
  expect_equal(cmp$engine, "pca")
  expect_equal(cmp$n_obs, 1000L)
  expect_equal(cmp$n_vars, 16L)

  # One row per split x level x factor; one summary row per level x factor.
  n_factors <- sum(seq_len(5L))
  expect_equal(nrow(cmp$coefficients), 4L * n_factors)
  expect_equal(nrow(cmp$summary), n_factors)
  expect_named(cmp$coefficients, c("split", "level", "factor", "r", "phi"))
  expect_named(
    cmp$summary,
    c("level", "factor", "r_median", "r_min", "phi_median", "phi_min", "n_splits_ok")
  )

  # Coefficients are correlations / congruences: bounded by 1 in magnitude.
  expect_true(all(abs(cmp$coefficients$r) <= 1 + 1e-8))
  expect_true(all(abs(cmp$coefficients$phi) <= 1 + 1e-8))
  expect_true(all(cmp$summary$n_splits_ok == 4L))

  # Factor labels match the full-sample ackwards() labelling scheme.
  expect_true(all(grepl("^m\\d+f\\d+$", cmp$summary$factor)))
})

# ---- Discriminating behaviour (the acceptance-criterion test) ----------------

test_that("comparability() separates replicable from non-replicable structure", {
  cmp <- .get_cmp()
  sm <- cmp$summary

  # sim16 has a known 1 -> 2 -> 4 hierarchy: those levels replicate cleanly.
  expect_true(all(sm$r_median[sm$level %in% c(1L, 2L, 4L)] > 0.95))
  expect_true(all(sm$phi_median[sm$level %in% c(1L, 2L, 4L)] > 0.95))

  # Level 3 sits between the true 2- and 4-factor structures (an arbitrary
  # blend), so at least one of its factors fails to replicate.
  expect_lt(min(sm$r_median[sm$level == 3L]), 0.9)

  # Level 5 is overextracted by design: the first four factors are the true
  # 4-structure (replicate), the fifth is a sample-idiosyncratic fragment.
  lev5 <- sm[sm$level == 5L, , drop = FALSE]
  expect_true(all(lev5$r_median[lev5$factor != "m5f5"] > 0.95))
  expect_lt(lev5$r_median[lev5$factor == "m5f5"], 0.5)
  expect_lt(lev5$phi_median[lev5$factor == "m5f5"], 0.5)
})

test_that("comparability() flags a pure-noise factor as non-replicable", {
  set.seed(99)
  dat <- cbind(
    sim16[, 1:8],
    noise1 = rnorm(1000), noise2 = rnorm(1000), noise3 = rnorm(1000)
  )
  cmp <- suppressMessages(comparability(dat, k_max = 3, n_splits = 4, seed = 1))
  sm <- cmp$summary

  # The two real factors at level 3 replicate; the noise-anchored one does not.
  lev3 <- sm[sm$level == 3L, , drop = FALSE]
  expect_lt(min(lev3$r_median), 0.8)
  expect_equal(sum(lev3$r_median > 0.95), 2L)
})

# ---- Internal algebra: identity, invariance, oracle ---------------------------

test_that(".cross_cor() of a solution with itself has a unit diagonal", {
  R <- stats::cor(sim16)
  lev <- pca_levels(R, k_max = 3L)$levels[["3"]]

  E <- .cross_cor(lev, lev, R)
  expect_equal(diag(E), rep(1, 3L), tolerance = 1e-12, ignore_attr = TRUE)

  # Self-comparability: every factor matches itself with r = phi = 1.
  out <- .level_comparability(lev, lev, lev, R)
  expect_equal(out$r, rep(1, 3L), tolerance = 1e-12)
  expect_equal(out$phi, rep(1, 3L), tolerance = 1e-12)
})

test_that(".cross_cor() agrees with materialised half-solution scores", {
  # Invariant 2 extended to cross-solution correlations: the W'RW algebra on
  # the pooled R must match correlating actual scores computed by applying
  # both halves' weights to the full (complete-data) sample.
  dat <- as.matrix(sim16)
  R <- stats::cor(dat)
  half_a <- dat[1:500, , drop = FALSE]
  half_b <- dat[501:1000, , drop = FALSE]
  lev_a <- pca_levels(stats::cor(half_a), k_max = 3L)$levels[["3"]]
  lev_b <- pca_levels(stats::cor(half_b), k_max = 3L)$levels[["3"]]

  E_algebra <- .cross_cor(lev_a, lev_b, R)
  Z <- scale(dat)
  E_scores <- stats::cor(Z %*% lev_a$scoring$weights, Z %*% lev_b$scoring$weights)
  expect_equal(unname(E_algebra), unname(E_scores), tolerance = 1e-12)
})

test_that(".level_comparability() is invariant to column order and sign", {
  dat <- as.matrix(sim16)
  R <- stats::cor(dat)
  lev_f <- pca_levels(R, k_max = 3L)$levels[["3"]]
  lev_a <- pca_levels(stats::cor(dat[1:500, ]), k_max = 3L)$levels[["3"]]
  lev_b <- pca_levels(stats::cor(dat[501:1000, ]), k_max = 3L)$levels[["3"]]

  base <- .level_comparability(lev_f, lev_a, lev_b, R)

  # Permute and sign-flip half B's factors; the anchored coefficients must
  # not change (matching + sign alignment absorb both transformations).
  perm <- c(3L, 1L, 2L)
  flips <- c(-1, 1, -1)
  lev_b2 <- lev_b
  lev_b2$loadings <- sweep(lev_b$loadings[, perm, drop = FALSE], 2L, flips, "*")
  lev_b2$scoring$weights <-
    sweep(lev_b$scoring$weights[, perm, drop = FALSE], 2L, flips, "*")

  permuted <- .level_comparability(lev_f, lev_a, lev_b2, R)
  expect_equal(permuted$r, base$r, tolerance = 1e-12)
  expect_equal(permuted$phi, base$phi, tolerance = 1e-12)
})

test_that(".match_square() picks the global-max bijection", {
  E <- rbind(
    c(0.2, 0.9, 0.1),
    c(0.85, 0.3, 0.1),
    c(0.4, 0.5, -0.95) # negative: matching is on |r|
  )
  expect_identical(.match_square(E), c(2L, 1L, 3L))
  expect_identical(.match_square(diag(3L)), 1:3)
})

test_that(".level_comparability() returns NA rows for a missing half-level", {
  R <- stats::cor(sim16)
  lev <- pca_levels(R, k_max = 2L)$levels[["2"]]
  out <- .level_comparability(lev, NULL, lev, R)
  expect_equal(out$factor, lev$labels)
  expect_true(all(is.na(out$r)))
  expect_true(all(is.na(out$phi)))
})

# ---- Reproducibility ----------------------------------------------------------

test_that("comparability() is reproducible under a seed", {
  a <- suppressMessages(comparability(sim16, k_max = 3, n_splits = 2, seed = 42))
  b <- suppressMessages(comparability(sim16, k_max = 3, n_splits = 2, seed = 42))
  expect_identical(a$coefficients, b$coefficients)
  expect_identical(a$summary, b$summary)

  c <- suppressMessages(comparability(sim16, k_max = 3, n_splits = 2, seed = 43))
  expect_false(identical(a$coefficients$r, c$coefficients$r))
})

# ---- EFA engine ---------------------------------------------------------------

test_that("comparability() works with the EFA engine", {
  cmp <- suppressMessages(
    comparability(sim16[, 1:8], k_max = 3, engine = "efa", n_splits = 2, seed = 1)
  )
  expect_s3_class(cmp, "comparability")
  expect_equal(cmp$engine, "efa")
  expect_equal(cmp$fm, "minres")
  # The true 2-factor structure of these 8 items replicates.
  expect_true(all(cmp$summary$r_median[cmp$summary$level == 2L] > 0.95))
})

# ---- Validation ----------------------------------------------------------------

test_that("comparability() rejects unsupported engines and bases with pointers", {
  expect_error(
    comparability(sim16, k_max = 3, engine = "esem"),
    "not supported"
  )
  expect_error(
    comparability(sim16, k_max = 3, cor = "polychoric"),
    "not supported"
  )
  expect_error(comparability(sim16, k_max = 3, engine = "midi"), "must be one of")
  expect_error(comparability(sim16, k_max = 3, cor = "kendall"), "must be one of")
})

test_that("comparability() rejects non-raw-data input", {
  R <- stats::cor(sim16)
  expect_error(comparability(R, k_max = 3), "raw item data")
  expect_error(comparability(list(a = 1), k_max = 3), "data frame or numeric matrix")
  ch <- matrix(letters[1:20], nrow = 5)
  expect_error(comparability(ch, k_max = 3), "numeric columns")
})

test_that("comparability() validates n_splits and half-sample size", {
  expect_error(comparability(sim16, k_max = 3, n_splits = 0), "positive integer")
  expect_error(comparability(sim16, k_max = 3, n_splits = 1.5), "positive integer")
  expect_error(comparability(sim16, k_max = 3, n_splits = c(2, 3)), "positive integer")

  # n/2 must exceed p, or half-sample solutions cannot be scored.
  expect_error(
    comparability(sim16[1:30, ], k_max = 3),
    "more rows than variables"
  )
})

# ---- Convergence shortfalls (mocked) -------------------------------------------

test_that("a level missing from some half-fits yields NA + a shortfall message", {
  real_fit_half <- .fit_half
  calls <- 0L
  testthat::local_mocked_bindings(
    .fit_half = function(data_half, k_max, engine, cor, fm) {
      calls <<- calls + 1L
      out <- real_fit_half(data_half, k_max, engine, cor, fm)
      if (calls == 1L) out[["3"]] <- NULL # split 1, half A: level 3 "fails"
      out
    },
    .package = "ackwards"
  )
  expect_message(
    cmp <- comparability(sim16, k_max = 3, n_splits = 2, seed = 1),
    "did not converge"
  )
  lev3 <- cmp$summary[cmp$summary$level == 3L, , drop = FALSE]
  expect_true(all(lev3$n_splits_ok == 1L))
  expect_false(anyNA(lev3$r_median)) # one usable split still aggregates
  expect_equal(sum(is.na(cmp$coefficients$r)), 3L) # split 1's level-3 rows
  expect_no_error(print(cmp)) # "[1/2 splits usable]" suffix branch
})

test_that("a level with zero usable splits aggregates to NA and prints", {
  real_fit_half <- .fit_half
  testthat::local_mocked_bindings(
    .fit_half = function(data_half, k_max, engine, cor, fm) {
      out <- real_fit_half(data_half, k_max, engine, cor, fm)
      out[["3"]] <- NULL # level 3 "fails" in every half
      out
    },
    .package = "ackwards"
  )
  expect_message(
    cmp <- comparability(sim16, k_max = 3, n_splits = 2, seed = 1),
    "did not converge"
  )
  lev3 <- cmp$summary[cmp$summary$level == 3L, , drop = FALSE]
  expect_true(all(is.na(lev3$r_median)))
  expect_true(all(lev3$n_splits_ok == 0L))
  expect_no_error(print(cmp)) # "no usable splits" branch
})

# ---- print + autoplot -----------------------------------------------------------

test_that("print.comparability() runs and returns invisibly", {
  cmp <- .get_cmp()
  # cli output is not reliably captured by expect_output (see test-print.R);
  # assert the method runs and honours the print contract.
  expect_no_error(print(cmp))
  expect_invisible(print(cmp))
})

test_that("autoplot.comparability() builds a ggplot", {
  skip_if_not_installed("ggplot2")
  cmp <- .get_cmp()
  p <- autoplot(cmp)
  expect_s3_class(p, "ggplot")
  # Both panels present: score comparability and loading congruence.
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(levels(built$plot$data$panel)), 2L)
})

test_that("benchmark citations name the verified lineage, not Goldberg (1990)", {
  # Regression test (2026-07-16): the .90/.95 benchmark lines used to cite
  # "Goldberg, 1990", which contains no split-half comparability analyses.
  # Verified lineage (cairn/references/): Everett 1983 (procedure + .90
  # rationale) and Saucier et al. 2005 (the Goldberg-lab .90 split-half gate).
  cmp <- .get_cmp()

  msgs <- character()
  withCallingHandlers(
    txt <- utils::capture.output(print(cmp)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  out <- cli::ansi_strip(paste(c(txt, msgs), collapse = " "))
  expect_match(out, "Everett", fixed = TRUE)
  expect_match(out, "Saucier", fixed = TRUE)
  expect_no_match(out, "Goldberg, 1990", fixed = TRUE)

  skip_if_not_installed("ggplot2")
  cap <- autoplot(cmp)$labels$caption
  expect_match(cap, "Saucier", fixed = TRUE)
  expect_no_match(cap, "Goldberg, 1990", fixed = TRUE)
})

# ---- Post-review follow-up coverage ---------------------------------------------

test_that(".level_comparability() returns NA rows when cross matrices carry NA", {
  R <- stats::cor(sim16)
  lev <- pca_levels(R, k_max = 2L)$levels[["2"]]
  # An NA cell in the pooled R (pathological pairwise missingness) propagates
  # into the cross-solution matrices; matching must degrade to NA, not crash.
  R_bad <- R
  R_bad[1L, 2L] <- R_bad[2L, 1L] <- NA_real_
  out <- .level_comparability(lev, lev, lev, R_bad)
  expect_equal(out$factor, lev$labels)
  expect_true(all(is.na(out$r)))
  expect_true(all(is.na(out$phi)))
})

test_that(".fit_half() returns an empty levels list when a half cannot be factored", {
  d <- as.matrix(sim16[1:100, 1:6])
  d[, 6] <- 1 # zero-variance column -> NA row/col in the half R
  expect_identical(.fit_half(d, k_max = 3L, engine = "pca", cor = "pearson", fm = "minres"), list())
  expect_identical(.fit_half(d, k_max = 3L, engine = "efa", cor = "pearson", fm = "minres"), list())
})

test_that("unknown arguments in ... error loudly (typo guard)", {
  expect_error(comparability(sim16, k_max = 3, nsplits = 20), "Unknown argument")
  expect_error(suggest_k(sim16, kmax = 4), "Unknown argument")
  expect_error(ackwards(sim16, k_max = 3, alignsigns = FALSE), "Unknown argument")
  # Unnamed extras are rejected too.
  expect_error(comparability(sim16, 3, "pca", "pearson", "minres", 2, NULL, 42), "unnamed extra")
  # The M34 moved-args guard keeps its more specific message.
  expect_error(ackwards(sim16, k_max = 3, prune = "redundant"), "no longer")
})

test_that("comparability() works with a single split and the spearman basis", {
  cmp1 <- suppressMessages(comparability(sim16, k_max = 2, n_splits = 1, seed = 5))
  expect_equal(cmp1$n_splits, 1L)
  # With one split the median and min across splits coincide.
  expect_identical(cmp1$summary$r_median, cmp1$summary$r_min)

  cmp_sp <- suppressMessages(
    comparability(sim16, k_max = 2, cor = "spearman", n_splits = 2, seed = 5)
  )
  expect_equal(cmp_sp$cor, "spearman")
  expect_true(all(cmp_sp$summary$r_median > 0.9)) # sim16's 2-level structure
})

test_that("a truncated full-sample anchor is recorded and printed", {
  real_ackwards <- ackwards
  testthat::local_mocked_bindings(
    # Simulate the anchor truncating (non-convergence below the request):
    # whatever k_max comparability() asks for, the full-sample fit stops at 2.
    ackwards = function(data, k_max, ...) real_ackwards(data, k_max = 2L, ...),
    .package = "ackwards"
  )
  cmp <- suppressMessages(comparability(sim16, k_max = 3, n_splits = 2, seed = 1))
  expect_equal(cmp$k_max, 2L)
  expect_equal(cmp$k_requested, 3L)
  expect_true(all(cmp$coefficients$level <= 2L))
  expect_no_error(print(cmp)) # "(requested 1-3; ...)" suffix branch
})

test_that("autoplot.comparability() handles all-NA levels from shortfalls", {
  skip_if_not_installed("ggplot2")
  real_fit_half <- .fit_half
  testthat::local_mocked_bindings(
    .fit_half = function(data_half, k_max, engine, cor, fm) {
      out <- real_fit_half(data_half, k_max, engine, cor, fm)
      out[["3"]] <- NULL # level 3 "fails" in every half
      out
    },
    .package = "ackwards"
  )
  cmp <- suppressMessages(comparability(sim16, k_max = 3, n_splits = 2, seed = 1))
  p <- autoplot(cmp)
  expect_s3_class(p, "ggplot")
  # NA coefficients are dropped from the plot data, not drawn.
  expect_false(anyNA(p$data$value))
  expect_true(all(p$data$level <= 2L))
})
