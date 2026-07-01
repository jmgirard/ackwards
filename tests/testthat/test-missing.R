# Tests for M16: estimator-aware missing-data handling
# Covers: .resolve_missing() validation, pairwise/listwise/fiml paths,
# n_obs/n_complete threading, pairwise warning, score NA propagation.

# ── .resolve_missing() validation ─────────────────────────────────────────────

test_that(".resolve_missing() is silent for valid pairwise/listwise combinations", {
  expect_silent(.resolve_missing("pairwise", "pca", NULL, "pearson"))
  expect_silent(.resolve_missing("listwise", "pca", NULL, "pearson"))
  expect_silent(.resolve_missing("pairwise", "efa", NULL, "pearson"))
  expect_silent(.resolve_missing("listwise", "efa", NULL, "pearson"))
  expect_silent(.resolve_missing("pairwise", "esem", "ML", "pearson"))
  expect_silent(.resolve_missing("listwise", "esem", "ML", "pearson"))
  expect_silent(.resolve_missing("pairwise", "esem", "WLSMV", "polychoric"))
  expect_silent(.resolve_missing("listwise", "esem", "WLSMV", "polychoric"))
})

# M38: fiml + pca/efa is now VALID on the Pearson basis (corFiml route),
# but still errors for a non-Pearson basis.
test_that(".resolve_missing() is silent for fiml + pca/efa + pearson", {
  expect_silent(.resolve_missing("fiml", "pca", NULL, "pearson"))
  expect_silent(.resolve_missing("fiml", "efa", NULL, "pearson"))
})

test_that(".resolve_missing() errors for fiml + pca/efa + non-pearson basis", {
  expect_error(.resolve_missing("fiml", "pca", NULL, "spearman"), "pearson")
  expect_error(.resolve_missing("fiml", "efa", NULL, "polychoric"), "pearson")
})

test_that(".resolve_missing() errors for fiml + WLSMV", {
  expect_error(.resolve_missing("fiml", "esem", "WLSMV", "polychoric"), "WLSMV")
})

test_that(".resolve_missing() errors for fiml + ULSMV", {
  expect_error(.resolve_missing("fiml", "esem", "ULSMV", "polychoric"), "ULSMV")
})

test_that(".resolve_missing() is silent for fiml + ML", {
  expect_silent(.resolve_missing("fiml", "esem", "ML", "pearson"))
})

test_that(".resolve_missing() is silent for fiml + MLR", {
  expect_silent(.resolve_missing("fiml", "esem", "MLR", "pearson"))
})

# ── ackwards() missing= argument validation via ackwards() ────────────────────

test_that("ackwards() accepts missing = 'pairwise' (default)", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_no_error(suppressWarnings(ackwards(d, k_max = 2L, missing = "pairwise")))
})

test_that("ackwards() errors for missing = 'fiml' with pca/efa + non-pearson basis", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(200 * 8), 200, 8))
  expect_error(
    ackwards(d, k_max = 2L, missing = "fiml", engine = "pca", cor = "spearman"),
    "pearson"
  )
  expect_error(
    ackwards(d, k_max = 2L, missing = "fiml", engine = "efa", cor = "polychoric"),
    "pearson"
  )
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

# ── ESEM polychoric + listwise: complete-case consistency ─────────────────────

test_that("ESEM polychoric listwise: n_obs equals complete-case count", {
  skip_if_not_installed("lavaan")
  d <- .make_ordinal_data()
  d[1:15, 2] <- NA_integer_
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", cor = "polychoric",
      missing = "listwise"
    )
  )
  expect_equal(x$n_obs, nrow(d) - 15L)
  expect_equal(x$meta$n_complete, nrow(d) - 15L)
})

test_that("ESEM polychoric listwise: correlation matrix has no NAs", {
  skip_if_not_installed("lavaan")
  d <- .make_ordinal_data()
  d[1:15, 2] <- NA_integer_
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", cor = "polychoric",
      missing = "listwise"
    )
  )
  expect_false(any(is.na(x$r)))
})

# ── ESEM FIML + MLR end-to-end ────────────────────────────────────────────────

test_that("ESEM FIML with MLR estimator succeeds and produces valid edges", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:5, 1] <- NA_real_
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "MLR"
    )
  )
  expect_false(any(is.na(x$edges$tidy$r)))
  expect_equal(x$meta$missing, "fiml")
})

# ── FIML criterion 4: edge R derived from lavaan h1 saturated model ──────────

test_that("ESEM FIML r matches lavaan h1 saturated-model correlation", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[1:10, 1] <- NA_real_
  # keep_fits = TRUE so we can inspect the k=1 lavaan fit directly
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "ML", keep_fits = TRUE
    )
  )
  # r_lv in esem_levels() is populated from the k=1 fit's h1 model;
  # compare x$r against the same h1 extraction to verify the source
  fit_k1 <- x$fits[["1"]]
  h1 <- lavaan::lavInspect(fit_k1, "h1")
  cov_h1 <- if (is.list(h1[[1L]])) h1[[1L]]$cov else h1$cov
  r_from_h1 <- stats::cov2cor(cov_h1)
  expect_equal(x$r, r_from_h1, tolerance = 1e-8)
})

test_that("ESEM FIML edge R differs from pairwise R under substantial missingness", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data(n = 400)
  # 20% MCAR: two different variables, non-overlapping rows
  set.seed(99)
  rows <- sample(400, 160)
  d[rows[1:80], 1] <- NA_real_
  d[rows[81:160], 2] <- NA_real_
  x_fiml <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "ML"
    )
  )
  x_pw <- suppressWarnings(
    ackwards(d, k_max = 2L, engine = "esem", missing = "pairwise")
  )
  # FIML uses all rows for estimation; pairwise uses a separate stats::cor().
  # With 20% MCAR on different variables the correlation matrices will differ.
  expect_false(isTRUE(all.equal(x_fiml$r, x_pw$r, tolerance = 1e-3)))
})

# ── Pairwise warning mentions "fiml" for ESEM ML path ────────────────────────

test_that("pairwise warning for ESEM ML mentions fiml as an alternative", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[5, 1] <- NA_real_
  warns <- character(0L)
  withCallingHandlers(
    suppressMessages(ackwards(d,
      k_max = 2L, engine = "esem",
      estimator = "ML", missing = "pairwise"
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # The ESEM ML branch (supports_fiml = TRUE) should mention fiml in the advisory
  expect_true(any(grepl("fiml", warns, ignore.case = TRUE)))
})

# ── FIML score NA propagation ────────────────────────────────────────────────

test_that("ESEM FIML keep_scores: NA rows still produce NA scores", {
  skip_if_not_installed("lavaan")
  set.seed(42)
  d <- .make_esem_data()
  d[5, 1] <- NA_real_ # one incomplete row
  x <- suppressWarnings(
    ackwards(d,
      k_max = 2L, engine = "esem", missing = "fiml",
      estimator = "ML", keep_scores = TRUE
    )
  )
  # FIML improves estimation but does not impute item responses;
  # score projection still propagates NAs row-wise
  expect_true(any(is.na(x$scores[["1"]])))
  expect_true(any(is.na(x$scores[["2"]])))
})

# ── available.cases vs listwise: R matrices differ ────────────────────────────

test_that("ESEM WLSMV pairwise (available.cases) R differs from listwise R with NAs", {
  skip_if_not_installed("lavaan")
  d <- .make_ordinal_data()
  d[1:30, 1] <- NA_integer_ # 10% missingness concentrated on one variable
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
  # available.cases uses all rows for each pair; listwise uses only the 270
  # complete rows. With 10% missingness on one variable the polychoric
  # correlation matrices will differ.
  expect_false(isTRUE(all.equal(x_pw$r, x_lw$r, tolerance = 1e-3)))
})

# ── M38: FIML via psych::corFiml() for PCA/EFA ───────────────────────────────

# Continuous data with MCAR missingness on non-overlapping rows/vars, so the
# FIML matrix differs from the pairwise one and no ordinal warning fires.
.make_fiml_data <- function(n = 300, p = 8, miss = 40) {
  set.seed(7)
  d <- as.data.frame(matrix(stats::rnorm(n * p), n, p))
  for (j in seq_len(p)) d[sample(n, miss), j] <- NA_real_
  d
}

test_that(".corfiml_R() returns a valid correlation matrix (unit diagonal, PD)", {
  skip_if_not_installed("psych")
  d <- as.matrix(.make_fiml_data())
  R <- .corfiml_R(d)
  expect_equal(dim(R), c(ncol(d), ncol(d)))
  expect_equal(unname(diag(R)), rep(1, ncol(d)), tolerance = 1e-6)
  expect_gt(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("EFA missing='fiml' builds and routes R through psych::corFiml()", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  x <- suppressMessages(ackwards(d, k_max = 3L, engine = "efa", missing = "fiml"))
  expect_s3_class(x, "ackwards")
  expect_equal(x$meta$missing, "fiml")
  # x$r should match a direct corFiml() call, not the pairwise cor()
  expect_equal(x$r, .corfiml_R(as.matrix(d)), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(
    x$r, stats::cor(as.matrix(d), use = "pairwise.complete.obs"),
    tolerance = 1e-3
  )))
})

test_that("PCA missing='fiml' builds and uses the corFiml matrix", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  x <- suppressMessages(ackwards(d, k_max = 3L, engine = "pca", missing = "fiml"))
  expect_s3_class(x, "ackwards")
  expect_equal(x$r, .corfiml_R(as.matrix(d)), tolerance = 1e-8)
})

test_that("FIML n_obs = 'total' (default) uses all rows; 'complete' uses complete cases", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  n_complete <- sum(stats::complete.cases(d))
  x_total <- suppressMessages(ackwards(d, k_max = 3L, engine = "efa", missing = "fiml"))
  x_comp <- suppressMessages(
    ackwards(d, k_max = 3L, engine = "efa", missing = "fiml", n_obs = "complete")
  )
  expect_equal(x_total$n_obs, nrow(d))
  expect_equal(x_comp$n_obs, n_complete)
  expect_lt(n_complete, nrow(d)) # sanity: missingness reduced the complete count
})

test_that("PCA FIML n_obs = 'complete' is honored (engine-agnostic switch)", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  n_complete <- sum(stats::complete.cases(d))
  x_total <- suppressMessages(ackwards(d, k_max = 3L, engine = "pca", missing = "fiml"))
  x_comp <- suppressMessages(
    ackwards(d, k_max = 3L, engine = "pca", missing = "fiml", n_obs = "complete")
  )
  expect_equal(x_total$n_obs, nrow(d))
  expect_equal(x_comp$n_obs, n_complete)
  # PCA fit does not use N, so the R matrix (and thus edges) is identical either way
  expect_equal(x_total$r, x_comp$r, tolerance = 1e-10)
})

test_that("FIML point estimates (edges) are unchanged by the n_obs choice", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  x_total <- suppressMessages(ackwards(d, k_max = 3L, engine = "efa", missing = "fiml"))
  x_comp <- suppressMessages(
    ackwards(d, k_max = 3L, engine = "efa", missing = "fiml", n_obs = "complete")
  )
  # Only fit-index N differs; loadings/edges come from the same R.
  expect_equal(x_total$edges$tidy$r, x_comp$edges$tidy$r, tolerance = 1e-10)
  expect_equal(x_total$r, x_comp$r, tolerance = 1e-10)
})

test_that("FIML route announces itself and the approximate-fit caveat via cli", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  expect_message(
    suppressWarnings(ackwards(d, k_max = 3L, engine = "efa", missing = "fiml")),
    "corFiml"
  )
  expect_message(
    suppressWarnings(ackwards(d, k_max = 3L, engine = "efa", missing = "fiml")),
    "approximate"
  )
})

test_that("string n_obs errors when not on the raw-data FIML pca/efa path", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  # non-FIML raw data
  expect_error(
    suppressMessages(ackwards(d, k_max = 3L, engine = "efa", n_obs = "total")),
    "only valid"
  )
  # correlation-matrix input rejects a string n_obs (numeric only)
  expect_error(
    ackwards(stats::cor(as.matrix(d), use = "pairwise.complete.obs"),
      k_max = 3L, engine = "efa", n_obs = "complete"
    ),
    "positive integer"
  )
})

test_that("numeric n_obs with raw-data FIML warns and defaults to total", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  # The advisory uses .frequency = "once"; force it to fire even if an earlier
  # test in the suite already triggered the same frequency id.
  rlang::local_options(rlib_warning_verbosity = "verbose")
  expect_warning(
    x <- suppressMessages(
      ackwards(d, k_max = 3L, engine = "efa", missing = "fiml", n_obs = 999)
    ),
    "ignored"
  )
  expect_equal(x$n_obs, nrow(d))
})

test_that("PCA/EFA FIML keep_scores: incomplete rows still produce NA scores", {
  skip_if_not_installed("psych")
  d <- .make_fiml_data()
  x <- suppressMessages(suppressWarnings(
    ackwards(d, k_max = 2L, engine = "efa", missing = "fiml", keep_scores = TRUE)
  ))
  # corFiml estimates R but does not impute item responses -> NA score rows
  expect_true(any(is.na(x$scores[["1"]])))
})
