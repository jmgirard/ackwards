test_that("ackwards() with method = 'efa' returns a valid ackwards object", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 5, method = "efa"))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$method, "efa")
  expect_equal(x$k_max, 5L)
  expect_equal(x$n_obs, 2800L)
  expect_equal(x$cor_type, "pearson")

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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3, method = "efa"))

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

  suppressWarnings(x <- ackwards(data, k = 4, method = "efa"))

  E_scores <- compute_edges(
    levels  = x$levels,
    R       = R,
    method  = "scores",
    pairs   = "adjacent",
    data    = data,
    align   = FALSE
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3, method = "efa"))

  for (ki in seq_len(x$k_max)) {
    Phi <- x$levels[[as.character(ki)]]$factor_cor
    expect_equal(Phi, diag(ki),
      info = paste("factor_cor is I at level", ki)
    )
  }
})

test_that("EFA levels have correct structure and label formats", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3, method = "efa"))

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

test_that("print, tidy, glance work for EFA objects", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3, method = "efa"))

  expect_no_error(print(x))
  expect_s3_class(generics::tidy(x), "data.frame")
  expect_s3_class(generics::glance(x), "data.frame")
  expect_equal(nrow(generics::glance(x)), 1L)
})

test_that("EFA tidy() edges contain is_primary column", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3, method = "efa"))
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
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 4, method = "efa"))
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
    suppressWarnings(ackwards(d_sing, k = 3, method = "efa", fm = "ml")),
    regexp = "at least 2" # from our k_eff < 2 guard
  )
})

test_that("fm argument is validated", {
  skip_if_not_installed("psych")
  d <- psych::bfi[, 1:5]
  expect_error(
    ackwards(d, k = 2, method = "efa", fm = "bad_fm"),
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
  suppressWarnings(x <- ackwards(d, k = 4L, method = "efa", fm = "minres"))

  # Skip i=1 (level 1:2): psych::bassAckward has a diag(scalar) bug that silently
  # skips score-variance standardization when a k=1 level is involved — the call
  # diag(1/sqrt(rs)) with scalar rs creates a 1×1 identity rather than the intended
  # 1×1 diagonal. Levels 2:3 and 3:4 agree to machine precision; 1:2 does not.
  for (i in 2:3) {
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
  suppressWarnings(x <- ackwards(d, k = 3L, method = "efa", fm = "ml"))
  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  # Cross-check: algebra vs scores
  E_sc <- compute_edges(x$levels,
    R = cor(d), method = "scores",
    pairs = "adjacent", data = d, align = FALSE
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
  suppressWarnings(x <- ackwards(d, k = 3L, method = "efa", fm = "pa"))
  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  E_sc <- compute_edges(x$levels,
    R = cor(d), method = "scores",
    pairs = "adjacent", data = d, align = FALSE
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
  xp <- ackwards(d, k = 2L, method = "efa")
  expect_equal(xp$levels[["1"]]$scoring$basis, "pearson")
  expect_equal(xp$levels[["2"]]$scoring$basis, "pearson")
  # Spearman
  xs <- ackwards(d, k = 2L, method = "efa", cor = "spearman")
  expect_equal(xs$levels[["1"]]$scoring$basis, "spearman")
  expect_equal(xs$levels[["2"]]$scoring$basis, "spearman")
  # PCA engine too
  xpca_s <- ackwards(d, k = 2L, method = "pca", cor = "spearman")
  expect_equal(xpca_s$levels[["1"]]$scoring$basis, "spearman")
})

test_that("tidy(x, what = 'fit') returns one row per level with named indices", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3L, method = "efa"))
  td <- generics::tidy(x, what = "fit")
  expect_s3_class(td, "data.frame")
  expect_true(all(c("level", "index", "value") %in% names(td)))
  # 3 levels × 6 indices each
  expect_equal(nrow(td), 18L)
  expect_equal(sort(unique(td$level)), 1:3)
  # PCA fit returns different indices (eigenvalues) — should also work
  suppressWarnings(xp <- ackwards(psych::bfi[, 1:25], k = 3L, method = "pca"))
  td_pca <- generics::tidy(xp, what = "fit")
  expect_s3_class(td_pca, "data.frame")
  expect_equal(nrow(td_pca), 1L + 2L + 3L) # 1+2+3 eigenvalues across 3 levels
})
