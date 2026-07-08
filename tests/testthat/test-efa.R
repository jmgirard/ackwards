test_that("ackwards() with method = 'efa' returns a valid ackwards object", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 5, engine = "efa"))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$engine, "efa")
  expect_equal(x$k_max, 5L)
  expect_equal(x$n_obs, 2800L)
  expect_equal(x$cor, "pearson")

  # All five levels present and converged
  expect_equal(length(x$levels), 5L)
  expect_true(all(vapply(x$levels, `[[`, logical(1), "converged")))

  # Four adjacent edge matrices
  expect_equal(length(x$edges$matrices), 4L)
  expect_named(x$edges$matrices, c("1:2", "2:3", "3:4", "4:5"))

  for (i in 1:4) {
    E <- x$edges$matrices[[paste0(i, ":", i + 1)]]
    expect_equal(nrow(E), i, info = paste("rows of", i, ":", i + 1))
    expect_equal(ncol(E), i + 1, info = paste("cols of", i, ":", i + 1))
    expect_true(all(abs(E) <= 1 + 1e-9), info = "correlations in [-1, 1]")
  }
})

test_that("EFA levels use tenBerge scoring (linear, method = 'tenBerge')", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))

  for (ki in seq_len(x$k_max)) {
    sc <- x$levels[[as.character(ki)]]$scoring
    expect_true(isTRUE(sc$linear), info = paste("level", ki, "linear"))
    expect_equal(sc$method, "tenBerge", info = paste("level", ki, "method"))
    expect_equal(dim(sc$weights), c(25L, ki), info = paste("level", ki, "weight dims"))
    # tenBerge weights yield unit-variance scores: diag(W'RW) ≈ 1
    sv <- sc$score_var
    expect_true(all(abs(sv - 1) < 1e-8), info = paste("level", ki, "score_var ≈ 1"))
  }
})

test_that(".tenBerge_weights warns (once) on a rank-deficient loading matrix", {
  # Degenerate input no shipped engine emits, but the guard must be loud, not
  # silent (Invariant 6). R is PD; L has two identical columns, so B = L'R^-1 L
  # is singular and the relative-tolerance clamp fires.
  R <- matrix(0.3, 4L, 4L)
  diag(R) <- 1
  L_full <- matrix(c(0.8, 0.7, 0.1, 0.1, 0.1, 0.1, 0.8, 0.7), ncol = 2L)
  L_dup <- cbind(L_full[, 1L], L_full[, 1L])

  # Full-rank L: no warning, and W'RW = I (unit-variance scores).
  expect_no_warning(W_full <- .tenBerge_weights(R, L_full))
  expect_equal(diag(crossprod(W_full, R %*% W_full)), c(1, 1), tolerance = 1e-8)

  # Rank-deficient L: warns, returns finite weights, and W'RW departs from I --
  # exactly the case the comment now documents (edges stay valid downstream
  # because compute_edges() standardizes by the actual score SDs).
  expect_warning(W_dup <- .tenBerge_weights(R, L_dup), "rank-deficient")
  expect_true(all(is.finite(W_dup)))
  expect_false(isTRUE(all.equal(
    diag(crossprod(W_dup, R %*% W_dup)), c(1, 1)
  )))
})

test_that("EFA algebra and scores paths agree (algebra-vs-scores cross-check)", {
  skip_if_not_installed("psych")

  set.seed(42)
  n <- 200
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  data <- data.frame(
    x1 = f1 + 0.3 * rnorm(n),
    x2 = f1 + 0.3 * rnorm(n),
    x3 = f1 + 0.3 * rnorm(n),
    x4 = f2 + 0.3 * rnorm(n),
    x5 = f2 + 0.3 * rnorm(n),
    x6 = f2 + 0.3 * rnorm(n)
  )
  R <- cor(data)

  x <- cached(ackwards(data, k_max = 4, engine = "efa"))

  E_scores <- compute_edges(
    levels      = x$levels,
    R           = R,
    edge_method = "scores",
    pairs       = "adjacent",
    data        = data
  )$matrices

  for (key in names(x$edges$matrices)) {
    E_alg <- x$edges$matrices[[key]]
    E_sc <- E_scores[[key]]
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 1e-6,
      label = paste("EFA algebra vs scores for pair", key)
    )
  }
})

test_that("EFA factor_cor is identity (orthogonal rotation)", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))

  for (ki in seq_len(x$k_max)) {
    Phi <- x$levels[[as.character(ki)]]$factor_cor
    expect_equal(Phi, diag(ki),
      info = paste("factor_cor is I at level", ki)
    )
  }
})

test_that("EFA levels have correct structure and label formats", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))

  for (ki in 1:3) {
    lev <- x$levels[[as.character(ki)]]
    expect_equal(ncol(lev$loadings), ki)
    expect_equal(nrow(lev$loadings), 25L)
    expect_equal(lev$labels, paste0("m", ki, "f", seq_len(ki)))
    expect_equal(colnames(lev$loadings), lev$labels)
    expect_true(isTRUE(lev$converged))
    expect_equal(dim(lev$scoring$weights), c(25L, ki))
    # Fit indices named correctly
    expect_named(lev$fit, c("chi", "dof", "p_value", "RMSEA", "TLI", "BIC"))
  }
})

test_that("EFA chi and p_value are one consistent pair (M42/C1)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  x <- cached(ackwards(d, k_max = 3, engine = "efa"))

  for (ki in 2:3) {
    fv <- x$levels[[as.character(ki)]]$fit
    # p_value must be the tail probability of the reported chi at the reported
    # dof -- the pre-M42 pairing (empirical chi + likelihood p) violated this.
    expect_equal(
      unname(stats::pchisq(fv[["chi"]], fv[["dof"]], lower.tail = FALSE)),
      unname(fv[["p_value"]]),
      tolerance = 1e-8
    )
  }

  # chi is psych::fa()'s likelihood-ratio STATISTIC, not its empirical $chi.
  R <- stats::cor(d)
  ref <- psych::fa(R,
    nfactors = 3, rotate = "varimax", fm = "minres",
    n.obs = nrow(d)
  )
  fv3 <- x$levels[["3"]]$fit
  expect_equal(unname(fv3[["chi"]]), unname(ref$STATISTIC)[[1L]],
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(
    unname(fv3[["chi"]]), unname(ref$chi)[[1L]],
    tolerance = 1e-6
  )))
})

test_that("print, tidy, glance work for EFA objects", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))

  expect_no_error(print(x))
  expect_s3_class(generics::tidy(x), "data.frame")
  expect_s3_class(generics::glance(x), "data.frame")
  expect_equal(nrow(generics::glance(x)), 1L)
})

test_that("EFA tidy() edges contain is_primary column", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  td <- generics::tidy(x)
  expect_true("is_primary" %in% names(td))
  expect_type(td$is_primary, "logical")
  # Exactly one primary parent per factor (one TRUE per to-label per adjacent pair)
  for (k in 2:3) {
    td_k <- td[td$level_to == k, ]
    for (lab in unique(td_k$to)) {
      n_primary <- sum(td_k$is_primary[td_k$to == lab], na.rm = TRUE)
      expect_equal(n_primary, 1L, info = paste("one primary parent for", lab))
    }
  }
})

test_that("meta$k_requested stores original k even when not truncated", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 4, engine = "efa"))
  expect_equal(x$meta$k_requested, 4L)
  expect_equal(x$k_max, 4L) # no truncation on well-formed data
})

test_that("ackwards errors informatively when all EFA levels fail (k_eff < 2)", {
  skip_if_not_installed("psych")
  # Perfectly collinear data: x4 = x1 + x2 + x3 → singular R.
  # fm = "ml" errors on this with "L-BFGS-B needs finite values of 'fn'",
  # which causes efa_levels to return 0 levels → k_eff < 2 guard fires.
  set.seed(42)
  x1 <- rnorm(100)
  x2 <- rnorm(100)
  x3 <- rnorm(100)
  d_sing <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x1 + x2 + x3)

  expect_error(
    suppressWarnings(ackwards(d_sing, k_max = 3, engine = "efa", fm = "ml")),
    regexp = "at least 2" # from our k_eff < 2 guard
  )
})

test_that("fm argument is validated", {
  skip_if_not_installed("psych")
  d <- psych::bfi[, 1:5]
  expect_error(
    ackwards(d, k_max = 2, engine = "efa", fm = "bad_fm"),
    "fm"
  )
})

test_that("EFA edge correlations match psych::bassAckward(fm='minres') within tolerance", {
  skip_if_not_installed("psych")
  # Use complete cases so both implementations work from the same correlation matrix
  d <- as.data.frame(na.omit(psych::bfi[, 1:25]))
  ba <- suppressWarnings(
    psych::bassAckward(d, nfactors = 4L, fm = "minres", rotate = "varimax", plot = FALSE)
  )
  x <- cached(ackwards(d, k_max = 4L, engine = "efa", fm = "minres"))

  # psych::bassAckward mis-standardized any k=1 level: diag(1/sqrt(rs)) with a
  # scalar rs builds a 1×1 *identity* (dropping the scale factor), not the
  # intended 1×1 diagonal, so its level 1:2 edge was left unstandardized. (The
  # PCA oracle in test-pca.R is unaffected — its k=1 score variance is 1, so the
  # missing scaling is a no-op.) We reported this; Revelle fixed it in psych
  # 2.6.6 (diag(..., nrow = length(rs))). We verified that under 2.6.6 level 1:2
  # then agrees to ~8e-16, like 2:3 and 3:4. So compare 1:2 only when a fixed
  # psych is installed: this stays green on the still-buggy CRAN 2.6.5 and
  # auto-activates full k=1 coverage once psych >= 2.6.6 reaches the user/CRAN.
  lo <- if (utils::packageVersion("psych") >= "2.6.6") 1L else 2L
  for (i in lo:3) {
    psych_mat <- t(ba$bass.ack[[i + 1L]])
    our_mat <- x$edges$matrices[[paste0(i, ":", i + 1L)]]
    max_diff <- max(abs(abs(psych_mat) - abs(our_mat)))
    expect_lt(max_diff, 1e-4, label = paste("EFA oracle level", i, "->", i + 1))
  }
})

test_that("fm = 'ml' produces a valid object and passes algebra-vs-scores", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 200
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  d <- data.frame(
    x1 = f1 + 0.3 * rnorm(n), x2 = f1 + 0.3 * rnorm(n),
    x3 = f2 + 0.3 * rnorm(n), x4 = f2 + 0.3 * rnorm(n)
  )
  x <- cached(ackwards(d, k_max = 3L, engine = "efa", fm = "ml"))
  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  # Cross-check: algebra vs scores
  E_sc <- compute_edges(x$levels,
    R = cor(d), edge_method = "scores",
    pairs = "adjacent", data = d
  )$matrices
  for (key in names(x$edges$matrices)) {
    expect_lt(max(abs(abs(x$edges$matrices[[key]]) - abs(E_sc[[key]]))), 1e-4,
      label = paste("ml cross-check", key)
    )
  }
})

test_that("fm = 'pa' produces a valid object and passes algebra-vs-scores", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 200
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  d <- data.frame(
    x1 = f1 + 0.3 * rnorm(n), x2 = f1 + 0.3 * rnorm(n),
    x3 = f2 + 0.3 * rnorm(n), x4 = f2 + 0.3 * rnorm(n)
  )
  x <- cached(ackwards(d, k_max = 3L, engine = "efa", fm = "pa"))
  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  E_sc <- compute_edges(x$levels,
    R = cor(d), edge_method = "scores",
    pairs = "adjacent", data = d
  )$matrices
  for (key in names(x$edges$matrices)) {
    expect_lt(max(abs(abs(x$edges$matrices[[key]]) - abs(E_sc[[key]]))), 1e-4,
      label = paste("pa cross-check", key)
    )
  }
})

test_that("scoring$basis reflects the actual cor= argument", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 6), 200, 6))
  # Pearson
  xp <- cached(ackwards(d, k_max = 2L, engine = "efa"))
  expect_equal(xp$levels[["1"]]$scoring$basis, "pearson")
  expect_equal(xp$levels[["2"]]$scoring$basis, "pearson")
  # Spearman
  xs <- cached(ackwards(d, k_max = 2L, engine = "efa", cor = "spearman"))
  expect_equal(xs$levels[["1"]]$scoring$basis, "spearman")
  expect_equal(xs$levels[["2"]]$scoring$basis, "spearman")
  # PCA engine too
  xpca_s <- cached(ackwards(d, k_max = 2L, engine = "pca", cor = "spearman"))
  expect_equal(xpca_s$levels[["1"]]$scoring$basis, "spearman")
})

test_that("tidy(x, what = 'fit') returns one row per level with named indices", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3L, engine = "efa"))
  td <- generics::tidy(x, what = "fit")
  expect_s3_class(td, "data.frame")
  expect_true(all(c("level", "statistic", "value") %in% names(td)))
  # 3 levels × 6 indices each
  expect_equal(nrow(td), 18L)
  expect_equal(sort(unique(td$level)), 1:3)
  # PCA fit returns different indices (eigenvalues) — should also work
  xp <- cached(ackwards(psych::bfi[, 1:25], k_max = 3L, engine = "pca"))
  td_pca <- generics::tidy(xp, what = "fit")
  expect_s3_class(td_pca, "data.frame")
  expect_equal(nrow(td_pca), 1L + 2L + 3L) # 1+2+3 eigenvalues across 3 levels
})

# ── glance() fit columns for EFA ──────────────────────────────────────────────

test_that("glance() for EFA has TLI/RMSEA/BIC at deepest level; CFI/SRMR NA", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  g <- generics::glance(x)
  expect_true(all(c("CFI", "TLI", "RMSEA", "SRMR", "BIC") %in% names(g)))
  # EFA carries TLI, RMSEA, BIC but not CFI or SRMR
  expect_true(is.na(g$CFI))
  expect_true(is.na(g$SRMR))
  expect_false(is.na(g$TLI))
  expect_false(is.na(g$RMSEA))
  expect_false(is.na(g$BIC))
  # fit values come from the deepest converged level
  expect_equal(g$deepest_converged, 3L)
})

# --- Polychoric continuity correction (`correct`) ---------------------------
# Real-world ordinal data with a singleton-category item makes
# psych::polychoric() fail under its default continuity correction (0.5); psych
# itself suggests correct = 0. Build such a dataset and check the escape hatch.
.singleton_ordinal <- function(n = 300, seed = 1) {
  set.seed(seed)
  mk <- function(p) sample(seq_along(p), n, replace = TRUE, prob = p)
  d <- as.data.frame(lapply(1:5, function(i) mk(rep(0.2, 5))))
  names(d) <- paste0("i", 1:5)
  d$i6 <- c(rep(2L, n - 2L), 3L, 4L) # 298, 1, 1 -> singleton categories
  d
}

test_that("correct = 0 rescues polychoric EFA where the default (0.5) fails", {
  skip_if_not_installed("psych")
  d <- .singleton_ordinal()

  # Default correct = 0.5: psych::polychoric() fails, and ackwards() surfaces an
  # actionable error naming the correct = 0 remedy.
  expect_error(
    suppressWarnings(ackwards(d, k_max = 3, engine = "efa", cor = "polychoric")),
    "correct = 0"
  )

  # correct = 0 succeeds and keeps the polychoric basis (not the cor-matrix NA).
  x <- suppressWarnings(
    ackwards(d, k_max = 3, engine = "efa", cor = "polychoric", correct = 0)
  )
  expect_s3_class(x, "ackwards")
  expect_equal(x$cor, "polychoric")
  expect_equal(x$k_max, 3L)
})

test_that("correct is validated as a single non-negative number", {
  skip_if_not_installed("psych")
  d <- .singleton_ordinal()
  expect_error(
    ackwards(d, k_max = 3, cor = "polychoric", correct = -1),
    "non-negative"
  )
  expect_error(
    ackwards(d, k_max = 3, cor = "polychoric", correct = c(0, 0.5)),
    "single non-negative"
  )
})

# --- Near-singular polychoric matrix: one clear warning, no psych flood -------
.triplicate_ordinal <- function(n = 250, seed = 3) {
  set.seed(seed)
  mk <- function(p) sample(seq_along(p), n, replace = TRUE, prob = p)
  base <- as.data.frame(lapply(1:6, function(i) mk(c(.5, .25, .15, .07, .03))))
  d <- cbind(base, base, base) # triplicate -> perfectly collinear -> rank-deficient
  names(d) <- paste0("v", 1:18)
  d
}

test_that("ackwards() warns once on a near-singular polychoric matrix", {
  skip_if_not_installed("psych")
  d <- .triplicate_ordinal()
  rlang::reset_warning_verbosity("ackwards_near_singular")
  w <- testthat::capture_warnings(
    x <- suppressMessages(
      ackwards(d, k_max = 2, engine = "efa", cor = "polychoric", correct = 0)
    )
  )
  expect_true(any(grepl("near-singular", w)))
  expect_s3_class(x, "ackwards")
})

test_that("ackwards() does not leak psych's per-level chatter to the console", {
  skip_if_not_installed("psych")
  d <- .triplicate_ordinal()
  rlang::reset_warning_verbosity("ackwards_near_singular")
  msgs <- capture.output(
    suppressWarnings(
      ackwards(d, k_max = 2, engine = "efa", cor = "polychoric", correct = 0)
    ),
    type = "message"
  )
  expect_false(any(grepl("determinant|smcs|objective function", msgs)))
})

test_that("near-singular fit is recorded in meta and re-surfaced by print/summary", {
  skip_if_not_installed("psych")
  d <- .triplicate_ordinal()
  rlang::reset_warning_verbosity("ackwards_near_singular")
  x <- suppressWarnings(suppressMessages(
    ackwards(d, k_max = 2, engine = "efa", cor = "polychoric", correct = 0)
  ))
  expect_true(x$meta$near_singular)
  expect_lt(x$meta$min_eigenvalue, 1e-4)
  # Durable caution appears in both print() and summary(), not just at fit time.
  pr <- cli::ansi_strip(capture.output(print(x), type = "message"))
  expect_true(any(grepl("Near-singular", pr)))
  sm <- cli::ansi_strip(capture.output(print(summary(x)), type = "message"))
  expect_true(any(grepl("Near-singular", sm)))
})

test_that("a well-conditioned fit is not flagged near-singular", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(psych::bfi[, 1:25], k_max = 3, engine = "efa"))
  expect_false(isTRUE(x$meta$near_singular))
  expect_false(any(grepl(
    "Near-singular",
    cli::ansi_strip(capture.output(print(x), type = "message"))
  )))
})
