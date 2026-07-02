# Tests for boot_edges() -- bootstrap confidence intervals on edges (M47).

# Shared fixture: one seeded PCA bootstrap on sim16, reused across the
# structure / seed / integration tests (memoised to keep the suite fast).
.boot_cache <- new.env(parent = emptyenv())
.get_boot <- function() {
  if (is.null(.boot_cache$x)) {
    x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))
    .boot_cache$x <- x
    .boot_cache$xb <- suppressMessages(
      boot_edges(x, sim16, n_boot = 200, seed = 1)
    )
  }
  .boot_cache$xb
}

# ---- Object structure --------------------------------------------------------

test_that("boot_edges() returns the object with a well-formed $boot", {
  xb <- .get_boot()

  expect_s3_class(xb, "ackwards")
  expect_type(xb$boot, "list")
  expect_named(xb$boot, c("edges", "n_boot", "conf", "seed"))

  be <- xb$boot$edges
  expect_s3_class(be, "data.frame")
  expect_named(be, c(
    "from", "to", "level_from", "level_to", "r", "se", "lo", "hi", "n_boot_ok"
  ))
  # One row per edge, aligned with the object's tidy edge table.
  expect_identical(nrow(be), nrow(xb$edges$tidy))
  expect_identical(be$from, xb$edges$tidy$from)
  expect_identical(be$to, xb$edges$tidy$to)
  expect_identical(be$r, xb$edges$tidy$r)
})

test_that("boot_edges() intervals bracket the point estimate and endpoints order", {
  be <- .get_boot()$boot$edges
  expect_true(all(be$lo <= be$hi))
  # Percentile endpoints stay inside the correlation bound.
  expect_true(all(be$lo >= -1 & be$hi <= 1))
  expect_true(all(be$se >= 0))
  # The point estimate lies within the (wider) interval for stable edges;
  # allow a small slack because percentile endpoints are order statistics.
  expect_true(all(be$r >= be$lo - 0.05 & be$r <= be$hi + 0.05))
})

test_that("boot_edges() is pipeable off ackwards()", {
  xb <- suppressMessages(suppressWarnings(
    ackwards(sim16, k_max = 3) |> boot_edges(sim16, n_boot = 30, seed = 4)
  ))
  expect_false(is.null(xb$boot))
})

# ---- Reproducibility + serial == parallel ------------------------------------

test_that("boot_edges() is reproducible under a seed", {
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))
  a <- suppressMessages(boot_edges(x, sim16, n_boot = 80, seed = 11))
  b <- suppressMessages(boot_edges(x, sim16, n_boot = 80, seed = 11))
  expect_identical(a$boot$edges, b$boot$edges)

  # A different seed gives a different (but same-shaped) result.
  c <- suppressMessages(boot_edges(x, sim16, n_boot = 80, seed = 12))
  expect_false(isTRUE(all.equal(a$boot$edges$se, c$boot$edges$se)))
})

test_that("serial and parallel dispatch agree exactly (upfront indices)", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))

  skip_if(!future::supportsMulticore(), "multicore plan unavailable here")

  future::plan(future::sequential)
  serial <- suppressMessages(boot_edges(x, sim16, n_boot = 60, seed = 21))

  future::plan(future::multicore, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  parallel <- suppressMessages(boot_edges(x, sim16, n_boot = 60, seed = 21))

  expect_equal(serial$boot$edges, parallel$boot$edges)
})

# ---- Anchoring (label switching + sign flipping) -----------------------------

test_that(".anchor_levels() restores full-sample factor order and sign", {
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))
  lev <- x$levels

  # Build a "replicate" that is the full-sample solution with, at each level,
  # columns permuted and some signs flipped -- exactly the label-switching /
  # sign-flipping an unanchored replicate would exhibit. Anchoring must undo it.
  rep_levels <- lapply(lev, function(l) {
    k <- l$k
    perm <- rev(seq_len(k))
    signs <- rep(c(1, -1), length.out = k)
    l$loadings <- sweep(l$loadings[, perm, drop = FALSE], 2L, signs, "*")
    l$scoring$weights <- sweep(l$scoring$weights[, perm, drop = FALSE], 2L, signs, "*")
    l$scoring$score_var <- l$scoring$score_var[perm]
    l
  })

  anchored <- ackwards:::.anchor_levels(lev, rep_levels, x$r)
  expect_false(any(vapply(anchored, is.null, logical(1L))))

  # After anchoring, each replicate factor aligns (positively) with its
  # full-sample counterpart, so the cross-solution matrix is ~ identity.
  for (ki in names(lev)) {
    cc <- ackwards:::.cross_cor(lev[[ki]], anchored[[ki]], x$r)
    expect_equal(cc, diag(nrow(cc)), tolerance = 1e-8, ignore_attr = TRUE)
  }

  # And edges recomputed from the anchored replicate reproduce the object's
  # own edges (the permutation/sign scrambling is fully undone).
  E_anchored <- compute_edges(
    levels = anchored, R = x$r, edge_method = "algebra",
    pairs = "adjacent", align = FALSE
  )$matrices
  for (key in names(x$edges$matrices)) {
    expect_equal(E_anchored[[key]], x$edges$matrices[[key]],
      tolerance = 1e-8, ignore_attr = TRUE
    )
  }
})

test_that("bootstrap primary-parent edges center on the full-sample estimate", {
  xb <- .get_boot()
  be <- xb$boot$edges
  prim <- xb$edges$tidy$is_primary
  # Primary edges are the stable spine; their point estimate should sit inside
  # the bootstrap interval (no sign-flip/label-switch corruption pulling the
  # distribution off the estimate).
  dprim <- be[prim, ]
  expect_true(all(dprim$r >= dprim$lo - 1e-6 & dprim$r <= dprim$hi + 1e-6))
})

# ---- Statistical oracle ------------------------------------------------------

test_that("CI arithmetic matches Fisher-z on fixed weights; full pipeline adds refit variance", {
  # The full boot_edges pipeline re-extracts factors each replicate, so its SE
  # legitimately EXCEEDS the fixed-weights (Fisher-z) SE that treats the scores
  # as observed variables. This test asserts both halves of that relationship:
  #   (a) a fixed-weights percentile bootstrap of the materialised scores
  #       reproduces the Fisher-z analytic CI (validates the percentile + SE
  #       arithmetic), and
  #   (b) the full pipeline centers on the estimate and its SE is >= the
  #       fixed-weights SE and of the same order (validates the machinery and
  #       the added refit variance).
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 2)))
  sc <- augment(x, data = sim16, append = FALSE)
  s1 <- sc$.m1f1
  s2 <- sc$.m2f1
  n <- length(s1)
  r0 <- cor(s1, s2)

  z <- atanh(r0)
  se_z <- 1 / sqrt(n - 3)
  fisher_lo <- tanh(z - stats::qnorm(0.975) * se_z)
  fisher_hi <- tanh(z + stats::qnorm(0.975) * se_z)
  fisher_se_r <- se_z * (1 - r0^2) # delta-method r-space SE

  # (a) fixed-weights percentile bootstrap ~ Fisher-z
  set.seed(2024)
  rb <- replicate(3000, {
    i <- sample.int(n, n, replace = TRUE)
    cor(s1[i], s2[i])
  })
  expect_equal(quantile(rb, 0.025, names = FALSE), fisher_lo, tolerance = 0.02)
  expect_equal(quantile(rb, 0.975, names = FALSE), fisher_hi, tolerance = 0.02)
  expect_equal(sd(rb), fisher_se_r, tolerance = 0.1) # Monte Carlo slack

  # (b) full pipeline: centered on r0, SE >= fixed-weights SE, same order
  xb <- suppressMessages(boot_edges(x, sim16, n_boot = 400, seed = 7))
  row <- xb$boot$edges[
    xb$boot$edges$from == "m1f1" & xb$boot$edges$to == "m2f1",
  ]
  expect_true(row$r >= row$lo && row$r <= row$hi)
  expect_gt(row$se, fisher_se_r) # refit variance strictly adds
  expect_lt(row$se, 4 * fisher_se_r) # but stays the same order of magnitude
})

test_that("strong edges exclude 0, a near-zero edge covers 0 (discriminating)", {
  xb <- .get_boot()
  be <- xb$boot$edges
  prim <- xb$edges$tidy$is_primary

  # Every primary-parent edge on sim16 is strong; its 95% CI excludes 0.
  dprim <- be[prim, ]
  expect_true(all(dprim$lo > 0))

  # There is at least one near-zero off-lineage edge whose CI covers 0.
  weakest <- be[which.min(abs(be$r)), ]
  expect_lt(abs(weakest$r), 0.1)
  expect_true(weakest$lo < 0 && weakest$hi > 0)
})

# ---- Invariant 7: a failing replicate is dropped, not fatal ------------------

test_that("a failing replicate is dropped and counted; the run completes", {
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))

  # Rig ~1/4 of the replicate refits to fail. .boot_replicate wraps the engine
  # in tryCatch, so a failure must produce an all-NA replicate (dropped), never
  # an abort (Invariant 7). Deterministic given the counter.
  real_pca <- ackwards:::pca_levels
  call_i <- 0L
  fail_pca <- function(R, k_max, cor = "pearson", keep_fits = FALSE) {
    call_i <<- call_i + 1L
    if (call_i %% 4L == 0L) stop("rigged replicate failure")
    real_pca(R, k_max = k_max, cor = cor, keep_fits = keep_fits)
  }

  testthat::local_mocked_bindings(pca_levels = fail_pca)
  xb <- suppressMessages(boot_edges(x, sim16, n_boot = 40, seed = 5))

  expect_s3_class(xb, "ackwards")
  be <- xb$boot$edges
  # Some replicates were dropped, so usable counts fall below n_boot but stay
  # high enough to form intervals.
  expect_true(all(be$n_boot_ok < 40))
  expect_true(all(be$n_boot_ok > 1))
  expect_false(anyNA(be$lo))
  expect_false(anyNA(be$hi))
})

# ---- Guards ------------------------------------------------------------------

test_that("boot_edges() rejects correlation-matrix-input objects", {
  R <- cor(sim16)
  x <- suppressMessages(suppressWarnings(ackwards(R, k_max = 3, n_obs = 1000)))
  expect_error(boot_edges(x, sim16, n_boot = 10), "raw item data")
})

test_that("boot_edges() rejects ESEM and polychoric objects", {
  skip_if_not_installed("lavaan")
  ord <- .make_ordinal_data()
  xe <- suppressMessages(suppressWarnings(
    ackwards(ord, k_max = 2, engine = "esem")
  ))
  expect_error(boot_edges(xe, ord, n_boot = 10), "esem")

  xp <- suppressMessages(suppressWarnings(
    ackwards(ord, k_max = 2, engine = "efa", cor = "polychoric")
  ))
  expect_error(boot_edges(xp, ord, n_boot = 10), "olychoric")
})

test_that("boot_edges() validates data and its arguments", {
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))

  expect_error(boot_edges(x, n_boot = 10), "required")
  expect_error(boot_edges(x, sim16[, 1:5], n_boot = 10), "missing")
  expect_error(boot_edges(x, sim16, n_boot = 1), "n_boot")
  expect_error(boot_edges(x, sim16, n_boot = 10, conf = 1.5), "conf")
  expect_error(
    boot_edges(x, sim16, n_boot = 10, notanarg = TRUE),
    "unused|unknown|notanarg"
  )
})

test_that("boot_edges() warns when data do not look like the fit data", {
  x <- suppressMessages(suppressWarnings(ackwards(sim16, k_max = 3)))
  # A different sample of the columns: right variables, wrong rows/moments.
  wrong <- sim16
  wrong[] <- lapply(wrong, function(v) v + 5)
  expect_warning(
    suppressMessages(boot_edges(x, wrong, n_boot = 10, seed = 1)),
    "does not look like"
  )
})

test_that("boot_edges() works with the EFA engine", {
  x <- suppressMessages(suppressWarnings(
    ackwards(sim16, k_max = 3, engine = "efa")
  ))
  xb <- suppressMessages(suppressWarnings(
    boot_edges(x, sim16, n_boot = 30, seed = 9)
  ))
  expect_false(is.null(xb$boot))
  expect_true(all(xb$boot$edges$lo <= xb$boot$edges$hi))
})

# ---- Downstream integration --------------------------------------------------

test_that("tidy(what = 'edges') gains bootstrap columns after boot_edges()", {
  xb <- .get_boot()
  te <- tidy(xb)
  expect_true(all(c("se", "lo", "hi", "n_boot_ok") %in% names(te)))
  # A plain object's edge table is unchanged.
  expect_false(any(c("se", "lo", "hi") %in% names(tidy(.boot_cache$x))))

  # sort = "strength" preserves the joined CI columns.
  ts <- tidy(xb, sort = "strength")
  expect_true(all(c("se", "lo", "hi") %in% names(ts)))
  expect_identical(nrow(ts), nrow(te))
})

test_that("print() and summary() note the bootstrap coverage", {
  # cli output is not reliably captured by expect_output (see test-print.R);
  # capture the cli stream directly and assert the boot line is present.
  xb <- .get_boot()
  print_out <- cli::cli_fmt(print(xb))
  expect_true(any(grepl("bootstrap", print_out, ignore.case = TRUE)))
  summ_out <- cli::cli_fmt(print(summary(xb)))
  expect_true(any(grepl("bootstrap", summ_out, ignore.case = TRUE)))
})
