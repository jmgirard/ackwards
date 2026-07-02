# tests/testthat/test-scores.R
# M6: scores storage, keep_fits, augment.ackwards(), tidy(what="scores")

# ── M45: fit-time moments + externally supplied standardization ───────────────────

test_that("fit-time item moments are stored in meta (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 2))
  expect_equal(x$meta$item_means, colMeans(as.matrix(d)))
  expect_equal(x$meta$item_sds, apply(as.matrix(d), 2, stats::sd))

  # listwise: moments reflect the reduced data actually fit
  suppressWarnings(xl <- ackwards(ackwards::bfi25, k_max = 2, missing = "listwise"))
  cc <- as.matrix(na.omit(ackwards::bfi25))
  expect_equal(unname(xl$meta$item_means), unname(colMeans(cc)))

  # correlation-matrix input: raw-data moments do not exist
  suppressWarnings(suppressMessages(xr <- ackwards(cor(d), k_max = 2)))
  expect_null(xr$meta$item_means)
  expect_null(xr$meta$item_sds)
})

test_that(".standardize() honours supplied center/scale moments (M45)", {
  set.seed(1)
  x <- matrix(rnorm(40, mean = 5, sd = 3), 20, 2, dimnames = list(NULL, c("a", "b")))
  ctr <- c(a = 1, b = 2)
  scl <- c(a = 2, b = 4)
  Z <- ackwards:::.standardize(x, center = ctr, scale = scl)
  expect_equal(Z, sweep(sweep(x, 2, ctr, "-"), 2, scl, "/"))
  # NULL default unchanged: sample moments
  Z0 <- ackwards:::.standardize(x)
  expect_equal(unname(colMeans(Z0)), c(0, 0))
})

# ── M45: out-of-sample scoring via augment(scaling =) ─────────────────────────────

test_that("scaling='fit' scores a test split by the training moments (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  train <- d[1:500, ]
  test <- d[501:nrow(d), ]
  suppressWarnings(x <- ackwards(train, k_max = 3))

  sc <- augment(x, data = test, append = FALSE)

  # Hand computation: standardize the TEST data by the TRAINING moments,
  # apply the stored weights, divide by the model-implied score SDs.
  mu <- colMeans(as.matrix(train))
  sg <- apply(as.matrix(train), 2, stats::sd)
  Z <- sweep(sweep(as.matrix(test), 2, mu, "-"), 2, sg, "/")
  for (ki in names(x$levels)) {
    W <- x$levels[[ki]]$scoring$weights
    S <- sweep(Z %*% W, 2, sqrt(x$levels[[ki]]$scoring$score_var), "/")
    for (j in seq_len(ncol(S))) {
      expect_equal(
        sc[[paste0(".", colnames(W)[j])]], unname(S[, j]),
        tolerance = 1e-12
      )
    }
  }
})

test_that("scaling='fit' is metric-consistent: subsets score like the full set (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 3))

  full <- augment(x, data = d, append = FALSE)
  sub <- augment(x, data = d[1:50, ], append = FALSE)
  expect_equal(sub, full[1:50, ], ignore_attr = TRUE)

  # Re-scoring the training data itself reproduces the fit-time stored scores.
  suppressWarnings(xs <- ackwards(d, k_max = 3, keep_scores = TRUE))
  stored <- do.call(cbind, lapply(xs$scores, unclass))
  again <- augment(xs, data = d, append = FALSE)
  expect_equal(unname(as.matrix(again)), unname(stored), tolerance = 1e-12)

  # Under scaling='sample', a subset deliberately scores differently (its own
  # moments) -- the pre-M45 behavior, retained as an explicit opt-in.
  sub_sample <- augment(x, data = d[1:50, ], append = FALSE, scaling = "sample")
  expect_false(isTRUE(all.equal(sub, sub_sample, ignore_attr = TRUE)))
})

test_that("scaling='sample' reproduces sample-moment scoring (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 2))
  sc <- augment(x, data = d, append = FALSE, scaling = "sample")
  Z <- ackwards:::.standardize(as.matrix(d))
  W <- x$levels[["2"]]$scoring$weights
  S <- sweep(Z %*% W, 2, sqrt(x$levels[["2"]]$scoring$score_var), "/")
  expect_equal(sc[[".m2f1"]], unname(S[, 1]), tolerance = 1e-12)
})

test_that("scaling='fit' errors informatively for correlation-matrix objects (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(suppressMessages(xr <- ackwards(cor(d), k_max = 2)))
  expect_error(augment(xr, data = d), "fit-time item means")
  # scaling='sample' remains available there
  expect_no_error(sc <- augment(xr, data = d, append = FALSE, scaling = "sample"))
  expect_equal(ncol(sc), 3L) # 1 + 2 score columns
})

test_that("scaling='fit' handles unnamed matrix input positionally (M45)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 2))
  m <- unname(as.matrix(d[1:20, ]))
  sc_unnamed <- augment(x, data = m, append = FALSE)
  sc_named <- augment(x, data = d[1:20, ], append = FALSE)
  expect_equal(sc_unnamed, sc_named, ignore_attr = TRUE)
})

# ── keep_scores = FALSE (default) ──────────────────────────────────────────────────

test_that("keep_scores = FALSE (default) leaves x$scores as NULL", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3))
  expect_null(x$scores)
})

# ── keep_scores = TRUE: storage ────────────────────────────────────────────────────

test_that("keep_scores = TRUE stores list of n × k_j matrices with correct dims", {
  skip_if_not_installed("psych")
  n <- nrow(psych::bfi)
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  expect_equal(nrow(x$scores[["1"]]), n)
  expect_equal(nrow(x$scores[["2"]]), n)
  expect_equal(nrow(x$scores[["3"]]), n)
  expect_equal(ncol(x$scores[["1"]]), 1L)
  expect_equal(ncol(x$scores[["2"]]), 2L)
  expect_equal(ncol(x$scores[["3"]]), 3L)
})

test_that("stored score column names match factor labels", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
  expect_equal(colnames(x$scores[["1"]]), "m1f1")
  expect_equal(colnames(x$scores[["2"]]), c("m2f1", "m2f2"))
  expect_equal(colnames(x$scores[["3"]]), c("m3f1", "m3f2", "m3f3"))
})

# ── Invariant 1: standardization ──────────────────────────────────────────────

test_that("PCA tenBerge scores are approximately unit variance (Inv. 1)", {
  skip_if_not_installed("psych")
  set.seed(42)
  data <- as.data.frame(matrix(rnorm(600), 200, 6))
  x <- suppressWarnings(ackwards(data, k_max = 3, keep_scores = TRUE))
  # For orthogonal (varimax) scores, SD ≈ 1 (tenBerge: W'RW = I exactly)
  for (ki in names(x$scores)) {
    sds <- apply(x$scores[[ki]], 2, sd, na.rm = TRUE)
    expect_true(all(abs(sds - 1) < 0.05),
      label = paste("SDs ≈ 1 at level", ki)
    )
  }
})

# ── Algebra vs. materialized scores cross-check (Inv. 2) ─────────────────────

test_that("between-level correlations from stored scores agree with algebra edges", {
  skip_if_not_installed("psych")
  set.seed(7)
  n <- 300
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = g + s1 + rnorm(n, sd = 0.2),
    x2 = g + s1 + rnorm(n, sd = 0.2),
    x3 = g + s1 + rnorm(n, sd = 0.2),
    x4 = g + s2 + rnorm(n, sd = 0.2),
    x5 = g + s2 + rnorm(n, sd = 0.2),
    x6 = g + s2 + rnorm(n, sd = 0.2)
  )
  x <- suppressWarnings(ackwards(data, k_max = 3, keep_scores = TRUE))

  # For each adjacent pair, correlation of materialized scores ≈ algebra edge
  for (ki in seq(2L, x$k_max)) {
    key <- paste0(ki - 1L, ":", ki)
    E_alg <- x$edges$matrices[[key]]
    Sa <- x$scores[[as.character(ki - 1L)]]
    Sb <- x$scores[[as.character(ki)]]
    E_sc <- cor(Sa, Sb)
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 0.05,
      label = paste("algebra vs. score correlations at level pair", key)
    )
  }
})

# ── keep_fits = FALSE (default) ───────────────────────────────────────────────

test_that("keep_fits = FALSE (default) leaves x$fits as NULL", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_null(x$fits)
})

# ── keep_fits = TRUE: storage ─────────────────────────────────────────────────

test_that("keep_fits = TRUE stores list of raw psych objects for PCA", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2", "3"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "psych"),
      label = paste("psych object at level", ki)
    )
  }
})

test_that("keep_fits = TRUE and keep_scores = TRUE can be combined", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(
    ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE, keep_fits = TRUE)
  )
  expect_false(is.null(x$scores))
  expect_false(is.null(x$fits))
})

# ── augment.ackwards() ────────────────────────────────────────────────────────

test_that("augment(x, data) returns data frame with score columns appended", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  out <- suppressWarnings(augment(x, data = bfi_items))
  expect_s3_class(out, "data.frame")
  # Same number of rows as input
  expect_equal(nrow(out), nrow(bfi_items))
  # Original columns present
  expect_true(all(names(bfi_items) %in% names(out)))
  # Score columns present: .m1f1, .m2f1, .m2f2, .m3f1, .m3f2, .m3f3
  score_cols <- grep("^\\.m", names(out), value = TRUE)
  expect_length(score_cols, 1L + 2L + 3L)
})

test_that("augment(x) without data uses stored scores", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE))
  out <- augment(x)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), x$n_obs)
  score_cols <- grep("^\\.m", names(out), value = TRUE)
  expect_length(score_cols, 1L + 2L)
})

test_that("augment(x) without data or stored scores gives informative error", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(augment(x), "not stored")
})

test_that("augment(x, data) and augment(x) [with stored scores] agree", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2, keep_scores = TRUE))
  out_stored <- suppressWarnings(augment(x))
  out_recomp <- suppressWarnings(augment(x, data = bfi_items))
  # Score columns should be numerically identical
  score_cols_stored <- grep("^\\.m", names(out_stored), value = TRUE)
  score_cols_recomp <- grep("^\\.m", names(out_recomp), value = TRUE)
  for (col in score_cols_stored) {
    expect_equal(out_stored[[col]], out_recomp[[col]],
      tolerance = 1e-10,
      label = paste("stored vs. recomputed scores:", col)
    )
  }
})

# ── augment(append=, id_cols=) : scores-only output (M36) ─────────────────────

test_that("augment(append = FALSE) returns only score columns, row order preserved", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  full <- suppressWarnings(augment(x, data = bfi_items))
  only <- suppressWarnings(augment(x, data = bfi_items, append = FALSE))
  # Only the .m score columns, none of the items
  expect_true(all(grepl("^\\.m", names(only))))
  expect_length(names(only), 1L + 2L + 3L)
  expect_equal(nrow(only), nrow(bfi_items))
  # cbind reproduces the appended output exactly (positional identity)
  expect_equal(cbind(bfi_items, only), full, ignore_attr = TRUE)
})

test_that("augment(append = FALSE, id_cols = ) carries identifier columns through", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  df <- data.frame(id = seq_len(nrow(bfi_items)), bfi_items)
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  out <- suppressWarnings(augment(x, data = df, append = FALSE, id_cols = "id"))
  expect_identical(names(out)[1L], "id")
  expect_identical(out$id, df$id)
  expect_true(all(grepl("^\\.m", names(out)[-1L])))
})

test_that("augment(append = FALSE) on stored scores drops the .obs index", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2, keep_scores = TRUE))
  full <- augment(x)
  only <- augment(x, append = FALSE)
  expect_true(".obs" %in% names(full))
  expect_false(".obs" %in% names(only))
  expect_true(all(grepl("^\\.m", names(only))))
})

test_that("augment() argument guards fire", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2))
  x2 <- suppressWarnings(ackwards(bfi_items, k_max = 2, keep_scores = TRUE))
  # id_cols with append = TRUE is a conflict
  expect_error(
    suppressWarnings(augment(x, data = bfi_items, id_cols = "A1")),
    "append = FALSE"
  )
  # id_cols with no data has no source columns
  expect_error(
    augment(x2, append = FALSE, id_cols = "A1"),
    "needs.*data|no source"
  )
  # id_cols naming an absent column
  expect_error(
    suppressWarnings(augment(x, data = bfi_items, append = FALSE, id_cols = "nope")),
    "not found"
  )
  # append must be a scalar logical
  expect_error(
    suppressWarnings(augment(x, data = bfi_items, append = NA)),
    "append"
  )
  # id_cols must be character
  expect_error(
    suppressWarnings(augment(x, data = bfi_items, append = FALSE, id_cols = 1L)),
    "character vector"
  )
})

test_that("augment(id_cols) works with matrix input carrying an extra id column", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  # Numeric matrix: the model item columns plus an extra numeric "id" column.
  mat <- cbind(id = seq_len(nrow(bfi_items)), as.matrix(bfi_items))
  out <- suppressWarnings(augment(x, data = mat, append = FALSE, id_cols = "id"))
  expect_identical(names(out)[1L], "id")
  expect_equal(out$id, seq_len(nrow(bfi_items)), ignore_attr = TRUE)
  expect_true(all(grepl("^\\.m", names(out)[-1L])))
})

test_that("augment(id_cols) may name a column that is also a model item", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  out <- suppressWarnings(augment(x, data = bfi_items, append = FALSE, id_cols = "A1"))
  # A1 is carried through verbatim AND still contributes to the scores.
  expect_identical(names(out)[1L], "A1")
  expect_identical(out$A1, bfi_items$A1)
  expect_true(all(grepl("^\\.m", names(out)[-1L])))
})

test_that("augment(append = FALSE) preserves NA-score rows and their id_cols", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[, 1:25]
  df <- data.frame(id = seq_len(nrow(bfi_items)), bfi_items)
  df$A1[1L] <- NA # force row 1 to score NA (listwise NA propagation)
  x <- suppressWarnings(ackwards(bfi_items, k_max = 3))
  out <- suppressWarnings(augment(x, data = df, append = FALSE, id_cols = "id"))
  # Row count preserved; the NA-scored row is retained, not dropped...
  expect_identical(nrow(out), nrow(df))
  expect_identical(out$id[1L], 1L)
  # ...with NA scores, so it still aligns positionally for a rejoin.
  expect_true(is.na(out$.m1f1[1L]))
})

# ── tidy(what = "scores") ─────────────────────────────────────────────────────

test_that("tidy(x, what='scores') returns long data frame with correct columns", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 3, keep_scores = TRUE))
  out <- tidy(x, what = "scores")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("obs", "level", "factor", "score") %in% names(out)))
  # n_obs × (1 + 2 + 3) rows
  expect_equal(nrow(out), x$n_obs * (1L + 2L + 3L))
  expect_setequal(unique(out$level), 1:3)
})

test_that("tidy(x, what='scores') errors informatively when no scores stored", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "scores"), "not stored")
})

# ── B1: Engine coverage for scores and keep_fits ──────────────────────────────

test_that("EFA keep_scores = TRUE stores correctly shaped matrices", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  x <- suppressWarnings(ackwards(d, k_max = 3, engine = "efa", keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  for (ki in 1:3) {
    expect_equal(nrow(x$scores[[as.character(ki)]]), 300L,
      label = paste("EFA score nrow at level", ki)
    )
    expect_equal(ncol(x$scores[[as.character(ki)]]), ki,
      label = paste("EFA score ncol at level", ki)
    )
  }
})

test_that("EFA keep_fits = TRUE stores psych objects for all levels", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300, 6))
  x <- suppressWarnings(ackwards(d, k_max = 2, engine = "efa", keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "psych"),
      label = paste("EFA psych object at level", ki)
    )
  }
})

test_that("ESEM keep_scores = TRUE stores correctly shaped matrices", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem", keep_scores = TRUE))
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  for (ki in 1:3) {
    expect_equal(nrow(x$scores[[as.character(ki)]]), nrow(d),
      label = paste("ESEM score nrow at level", ki)
    )
    expect_equal(ncol(x$scores[[as.character(ki)]]), ki,
      label = paste("ESEM score ncol at level", ki)
    )
  }
})

test_that("ESEM keep_fits = TRUE stores lavaan objects for all levels", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 2, engine = "esem", keep_fits = TRUE))
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2"))
  for (ki in names(x$fits)) {
    expect_true(inherits(x$fits[[ki]], "lavaan"),
      label = paste("ESEM lavaan object at level", ki)
    )
  }
})

# ── B3: Truncation — scores only cover converged levels ───────────────────────

test_that("keep_scores = TRUE only covers converged levels when model is truncated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # 6 variables → lavaan::efa() can only fit k <= 3; k = 5 triggers truncation
  suppressWarnings(x <- ackwards(d, k_max = 5, engine = "esem", keep_scores = TRUE))
  expect_equal(x$k_max, 3L)
  expect_false(is.null(x$scores))
  expect_named(x$scores, c("1", "2", "3"))
  expect_false("4" %in% names(x$scores))
  expect_false("5" %in% names(x$scores))
})

# ── B5: augment() column validation ──────────────────────────────────────────

test_that("augment(x, data) errors when data has wrong column count (unnamed)", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  # Use a plain matrix (no colnames) with wrong column count → dimension path
  d_wrong <- matrix(rnorm(100 * 5), 100, 5)
  expect_error(augment(x, data = d_wrong), "5.*column|column.*5")
})

test_that("augment(x, data) errors when named data is missing expected columns", {
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(psych::bfi[, 1:25], k_max = 2))
  # Supply only 20 of the 25 named BFI columns
  d_missing <- psych::bfi[1:50, 1:20]
  expect_error(augment(x, data = d_missing), "[Mm]issing")
})

test_that("augment(x, data) works when data has extra named columns (supersets)", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[1:50, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2))
  # Add an extra column not in the model
  d_extra <- cbind(bfi_items, extra = rnorm(50))
  out <- suppressWarnings(augment(x, data = d_extra))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 50L)
  expect_true("extra" %in% names(out))
  expect_true(".m1f1" %in% names(out))
})

test_that("augment(x, data) errors on non-numeric data", {
  skip_if_not_installed("psych")
  bfi_items <- psych::bfi[1:50, 1:25]
  x <- suppressWarnings(ackwards(bfi_items, k_max = 2))
  d_chr <- as.data.frame(lapply(bfi_items, as.character))
  expect_error(augment(x, data = d_chr), "[Nn]umeric")
})

# ── NA data: warning + correct propagation ────────────────────────────────────

test_that("augment(x, data) warns when data has missing rows", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  x <- ackwards(d, k_max = 2L)
  d_na <- d
  d_na[5L, 1L] <- NA_real_
  expect_warning(augment(x, data = d_na), "missing")
})

test_that("augment(x, data) produces NA scores for rows with missing values", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  x <- ackwards(d, k_max = 2L)
  d_na <- d
  d_na[5L, 1L] <- NA_real_
  out <- suppressWarnings(augment(x, data = d_na))
  expect_true(is.na(out$.m1f1[5L]))
  # Complete rows must NOT be NA (distinguishes row-NA from column-NA propagation)
  expect_false(anyNA(out$.m1f1[-5L]))
})

test_that("keep_scores = TRUE warns about NA score propagation", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  d[5L, 1L] <- NA_real_
  # Collect all warnings (pairwise advisory + score NA) and verify the score
  # NA warning is present; suppress so no warnings leak to testthat.
  warns <- character(0L)
  withCallingHandlers(
    ackwards(d, k_max = 2L, keep_scores = TRUE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("NA scores", warns)))
})

test_that("keep_scores = TRUE produces NA scores only for incomplete rows", {
  skip_if_not_installed("psych")
  set.seed(42)
  d <- as.data.frame(matrix(rnorm(200L * 10L), 200L, 10L))
  d[5L, 1L] <- NA_real_
  x <- suppressWarnings(ackwards(d, k_max = 2L, keep_scores = TRUE))
  # Row 5 (incomplete) must be NA at all factors
  expect_true(all(is.na(x$scores[["1"]][5L, ])))
  # All other rows (complete) must NOT be NA
  expect_false(anyNA(x$scores[["1"]][-5L, ]))
})

# ── Non-Pearson basis: warning on scoring ────────────────────────────────────

test_that("scoring with non-Pearson basis warns about basis mismatch", {
  skip_if_not_installed("psych")
  set.seed(42)
  # Build a small ordinal data set so polychoric is valid
  d <- as.data.frame(matrix(
    sample(1L:5L, 150L * 6L, replace = TRUE), 150L, 6L
  ))
  x <- suppressWarnings(ackwards(d, k_max = 2, cor = "polychoric"))
  expect_warning(augment(x, data = d), "polychoric")
})

# ── B4: keep_fits truncation ─────────────────────────────────────────────────

test_that("keep_fits = TRUE only stores fits for converged levels when truncated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # 6 variables → lavaan::efa() truncates at k = 3; k = 5 requested
  suppressWarnings(x <- ackwards(d, k_max = 5, engine = "esem", keep_fits = TRUE))
  expect_equal(x$k_max, 3L)
  expect_false(is.null(x$fits))
  expect_named(x$fits, c("1", "2", "3"))
  expect_false("4" %in% names(x$fits))
  expect_false("5" %in% names(x$fits))
})
