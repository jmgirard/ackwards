# Tests for M16: estimator-aware missing-data handling
# Covers: .resolve_missing() validation, pairwise/listwise/fiml paths,
# n_obs/n_complete threading, pairwise warning, score NA propagation.

# ── .resolve_missing() validation ─────────────────────────────────────────────

test_that(".resolve_missing() is silent for valid pairwise/listwise combinations", {
  expect_silent(.resolve_missing("pairwise", "pca", NULL))
  expect_silent(.resolve_missing("listwise", "pca", NULL))
  expect_silent(.resolve_missing("pairwise", "efa", NULL))
  expect_silent(.resolve_missing("listwise", "efa", NULL))
  expect_silent(.resolve_missing("pairwise", "esem", "ML"))
  expect_silent(.resolve_missing("listwise", "esem", "ML"))
  expect_silent(.resolve_missing("pairwise", "esem", "WLSMV"))
  expect_silent(.resolve_missing("listwise", "esem", "WLSMV"))
})

test_that(".resolve_missing() errors for fiml + pca", {
  expect_error(.resolve_missing("fiml", "pca", NULL), "pca")
})

test_that(".resolve_missing() errors for fiml + efa", {
  expect_error(.resolve_missing("fiml", "efa", NULL), "efa")
})

test_that(".resolve_missing() errors for fiml + WLSMV", {
  expect_error(.resolve_missing("fiml", "esem", "WLSMV"), "WLSMV")
})

test_that(".resolve_missing() errors for fiml + ULSMV", {
  expect_error(.resolve_missing("fiml", "esem", "ULSMV"), "ULSMV")
})

test_that(".resolve_missing() is silent for fiml + ML", {
  expect_silent(.resolve_missing("fiml", "esem", "ML"))
})

test_that(".resolve_missing() is silent for fiml + MLR", {
  expect_silent(.resolve_missing("fiml", "esem", "MLR"))
})

# ── ackwards() missing= argument validation via ackwards() ────────────────────

test_that("ackwards() accepts missing = 'pairwise' (default)", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_no_error(suppressWarnings(ackwards(d, k_max = 2L, missing = "pairwise")))
})

test_that("ackwards() errors for missing = 'fiml' with engine = 'pca'", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_error(ackwards(d, k_max = 2L, missing = "fiml", engine = "pca"), "pca")
})

test_that("ackwards() errors for missing = 'fiml' with engine = 'efa'", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_error(ackwards(d, k_max = 2L, missing = "fiml", engine = "efa"), "efa")
})

test_that("ackwards() errors for missing = 'fiml' with esem + WLSMV", {
  skip_if_not_installed("psych")
  skip_if_not_installed("lavaan")
  set.seed(1)
  d <- .make_ordinal_data()
  expect_error(
    suppressWarnings(ackwards(d,
      k_max = 2L, engine = "esem",
      cor = "polychoric", missing = "fiml"
    )),
    "WLSMV"
  )
})

test_that("ackwards() errors for unrecognised missing value", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_error(ackwards(d, k_max = 2L, missing = "em"), "arg_match|must be one of")
})

# ── Regression: default output unchanged ──────────────────────────────────────

test_that("missing = 'pairwise' (default) reproduces prior PCA output exactly", {
  skip_if_not_installed("psych")
  suppressWarnings({
    x_default <- ackwards(psych::bfi[, 1:25], k_max = 3)
    x_explicit <- ackwards(psych::bfi[, 1:25], k_max = 3, missing = "pairwise")
  })
  expect_equal(x_default$edges$tidy, x_explicit$edges$tidy)
  expect_equal(x_default$r, x_explicit$r)
})

# ── PCA listwise: n_obs and complete-case R ────────────────────────────────────

test_that("PCA listwise uses complete-case n_obs", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[1:10, 1] <- NA_real_ # 10 incomplete rows
  x <- suppressWarnings(ackwards(d, k_max = 2L, missing = "listwise"))
  expect_equal(x$n_obs, 190L)
  expect_equal(x$meta$n_complete, 190L)
})

test_that("PCA pairwise n_obs is total row count", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[1:10, 1] <- NA_real_
  x <- suppressWarnings(ackwards(d, k_max = 2L, missing = "pairwise"))
  expect_equal(x$n_obs, 200L)
  expect_equal(x$meta$n_complete, 190L)
})

test_that("PCA listwise correlation matrix has no NAs", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[1:10, 1] <- NA_real_
  x <- suppressWarnings(ackwards(d, k_max = 2L, missing = "listwise"))
  expect_false(any(is.na(x$r)))
})

# ── EFA listwise: n_obs threading ─────────────────────────────────────────────

test_that("EFA listwise uses complete-case n_obs", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(300 * 8), 300, 8))
  d[1:20, 2] <- NA_real_ # 20 incomplete rows
  x <- suppressWarnings(ackwards(d, k_max = 2L, engine = "efa", missing = "listwise"))
  expect_equal(x$n_obs, 280L)
  expect_equal(x$meta$n_complete, 280L)
})

# ── Pairwise warning fires when NAs present ───────────────────────────────────

test_that("pairwise with NAs emits a warning about missing values", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[5, 1] <- NA_real_
  expect_warning(
    suppressMessages(ackwards(d, k_max = 2L, missing = "pairwise")),
    "missing"
  )
})

test_that("pairwise with complete data emits no pairwise-missing warning", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  # Should get no pairwise-missing warning; may get ordinal warning — that's ok
  expect_no_warning(
    suppressMessages(ackwards(d, k_max = 2L, missing = "pairwise"))
  )
})

# ── meta$missing and meta$n_complete populated ────────────────────────────────

test_that("meta$missing records the missing argument", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  x_pw <- suppressWarnings(ackwards(d, k_max = 2L, missing = "pairwise"))
  x_lw <- suppressWarnings(ackwards(d, k_max = 2L, missing = "listwise"))
  expect_equal(x_pw$meta$missing, "pairwise")
  expect_equal(x_lw$meta$missing, "listwise")
})

test_that("meta$n_complete records complete case count", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[1:5, 1] <- NA_real_
  x <- suppressWarnings(ackwards(d, k_max = 2L))
  expect_equal(x$meta$n_complete, 195L)
})

# ── Scores still NA-propagate under every missing= ───────────────────────────

test_that("keep_scores = TRUE with listwise: scores for reduced data have no NAs", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[1:5, 1] <- NA_real_
  x <- suppressWarnings(ackwards(d,
    k_max = 2L, missing = "listwise",
    keep_scores = TRUE
  ))
  # listwise removes incomplete rows, so stored scores should be complete
  expect_false(any(is.na(x$scores[["1"]])))
  expect_equal(nrow(x$scores[["1"]]), 195L)
})

test_that("keep_scores = TRUE with pairwise: NA rows produce NA scores", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  d[5, 1] <- NA_real_
  x <- suppressWarnings(ackwards(d,
    k_max = 2L, missing = "pairwise",
    keep_scores = TRUE
  ))
  expect_true(any(is.na(x$scores[["1"]])))
})

# ── ESEM FIML path (ML estimator) ────────────────────────────────────────────

test_that("ESEM missing = 'fiml' succeeds with ML estimator", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:5, 1] <- NA_real_
  expect_no_error(
    suppressWarnings(
      ackwards(d,
        k_max = 2L, engine = "esem", missing = "fiml",
        estimator = "ML"
      )
    )
  )
})

test_that("ESEM FIML result has valid edges (no NAs in r)", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:5, 1] <- NA_real_
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "ML"
    )
  )
  expect_false(any(is.na(x$edges$tidy$r)))
})

test_that("ESEM FIML n_obs equals total rows (all rows used)", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:5, 1] <- NA_real_
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "ML"
    )
  )
  expect_equal(x$n_obs, nrow(d))
})

# ── ESEM listwise consistent fit-vs-edges ─────────────────────────────────────

test_that("ESEM listwise: n_obs equals complete-case count", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:10, 1] <- NA_real_ # 10 incomplete rows
  x <- suppressWarnings(
    ackwards(d, k_max = 2L, engine = "esem", missing = "listwise")
  )
  expect_equal(x$n_obs, 190L)
  expect_equal(x$meta$n_complete, 190L)
})

test_that("ESEM listwise: correlation matrix has no NAs", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:10, 1] <- NA_real_
  x <- suppressWarnings(
    ackwards(d, k_max = 2L, engine = "esem", missing = "listwise")
  )
  expect_false(any(is.na(x$r)))
})

# ── ESEM WLSMV pairwise uses available.cases (full N) ────────────────────────

test_that("ESEM WLSMV pairwise uses full N (available.cases), not complete-case N", {
  skip_if_not_installed("lavaan")
  d <- .make_ordinal_data()
  d[1:20, 1] <- NA_integer_ # 20 incomplete rows; 280 complete
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", cor = "polychoric",
      missing = "pairwise"
    )
  )
  # available.cases uses all rows; n_obs should reflect total, not complete-case N
  expect_equal(x$n_obs, nrow(d))
  expect_equal(x$meta$n_complete, nrow(d) - 20L)
})

test_that("ESEM WLSMV listwise uses complete-case N, pairwise uses full N", {
  skip_if_not_installed("lavaan")
  d <- .make_ordinal_data()
  d[1:20, 1] <- NA_integer_
  x_pw <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", cor = "polychoric",
      missing = "pairwise"
    )
  )
  x_lw <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", cor = "polychoric",
      missing = "listwise"
    )
  )
  expect_equal(x_pw$n_obs, nrow(d)) # pairwise: full N
  expect_equal(x_lw$n_obs, nrow(d) - 20L) # listwise: complete-case N
})
