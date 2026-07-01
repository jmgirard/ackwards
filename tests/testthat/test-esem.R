# (Shared helpers .make_esem_data() and .make_ordinal_data() live in helper-data.R)

# ── Basic ESEM validity ────────────────────────────────────────────────────────

test_that("ackwards() with method = 'esem' returns a valid ackwards object", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$engine, "esem")
  expect_equal(x$cor, "pearson")
  expect_equal(length(x$levels), 3L)
  expect_true(all(vapply(x$levels, `[[`, logical(1), "converged")))
})

test_that("ESEM levels have correct structure and label formats", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  for (ki in 1:3) {
    lev <- x$levels[[as.character(ki)]]
    expect_equal(ncol(lev$loadings), ki, info = paste("ncol loadings level", ki))
    expect_equal(nrow(lev$loadings), 6L, info = paste("nrow loadings level", ki))
    expect_equal(lev$labels, paste0("m", ki, "f", seq_len(ki)))
    expect_equal(colnames(lev$loadings), lev$labels)
    expect_true(isTRUE(lev$converged))
  }

  # Edge matrices
  expect_equal(length(x$edges$matrices), 2L)
  expect_named(x$edges$matrices, c("1:2", "2:3"))
  for (i in 1:2) {
    E <- x$edges$matrices[[paste0(i, ":", i + 1)]]
    expect_equal(nrow(E), i, info = paste("rows of", i, ":", i + 1))
    expect_equal(ncol(E), i + 1, info = paste("cols of", i, ":", i + 1))
    expect_true(all(abs(E) <= 1 + 1e-9), info = "correlations in [-1, 1]")
  }
})

# ── loadings_se ───────────────────────────────────────────────────────────────

test_that("ESEM levels have a loadings_se matrix with correct dimensions", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  for (ki in seq_len(x$k_max)) {
    lev <- x$levels[[as.character(ki)]]
    expect_false(is.null(lev$loadings_se),
      info = paste("loadings_se non-NULL at level", ki)
    )
    expect_true(is.matrix(lev$loadings_se),
      info = paste("loadings_se is matrix at level", ki)
    )
    expect_equal(dim(lev$loadings_se), c(6L, ki),
      info = paste("loadings_se dims at level", ki)
    )
  }
})

test_that("PCA and EFA engines have loadings_se = NULL", {
  skip_if_not_installed("psych")
  d <- .make_esem_data()
  xp <- suppressWarnings(ackwards(d, k_max = 2, engine = "pca"))
  xe <- suppressWarnings(ackwards(d, k_max = 2, engine = "efa"))

  for (ki in 1:2) {
    expect_null(xp$levels[[as.character(ki)]]$loadings_se,
      info = paste("PCA loadings_se NULL level", ki)
    )
    expect_null(xe$levels[[as.character(ki)]]$loadings_se,
      info = paste("EFA loadings_se NULL level", ki)
    )
  }
})

# ── Fit indices ───────────────────────────────────────────────────────────────

test_that("ESEM fit indices are named correctly (chi, dof, p_value, CFI, TLI, RMSEA, SRMR, BIC)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  expected_names <- c("chi", "dof", "p_value", "CFI", "TLI", "RMSEA", "SRMR", "BIC")
  for (ki in seq_len(x$k_max)) {
    expect_named(x$levels[[as.character(ki)]]$fit, expected_names,
      info = paste("fit names at level", ki)
    )
  }
})

# ── Scoring ───────────────────────────────────────────────────────────────────

test_that("ESEM levels use tenBerge scoring (linear, method = 'tenBerge')", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  for (ki in seq_len(x$k_max)) {
    sc <- x$levels[[as.character(ki)]]$scoring
    expect_true(isTRUE(sc$linear), info = paste("level", ki, "linear"))
    expect_equal(sc$method, "tenBerge", info = paste("level", ki, "method"))
    expect_equal(sc$basis, "pearson", info = paste("level", ki, "basis"))
    expect_equal(dim(sc$weights), c(6L, ki), info = paste("level", ki, "weight dims"))
  }
})

test_that("ESEM factor_cor is identity (orthogonal rotation)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  for (ki in seq_len(x$k_max)) {
    Phi <- unname(x$levels[[as.character(ki)]]$factor_cor)
    expect_equal(Phi, diag(ki),
      tolerance = 1e-6,
      info = paste("factor_cor is I at level", ki)
    )
  }
})

# ── Algebra-vs-scores cross-check ─────────────────────────────────────────────

test_that("ESEM algebra and scores paths agree (algebra-vs-scores cross-check)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  # Scope: continuous (Pearson) paths only. For polychoric ESEM the algebra side
  # uses lavaan's polychoric R while the scores route applies Pearson
  # standardization to the raw ordinal data — they diverge by design and the
  # cross-check oracle (DESIGN.md §5.4) does not apply there.
  E_scores <- compute_edges(
    levels      = x$levels,
    R           = x$r,
    edge_method = "scores",
    pairs       = "adjacent",
    data        = d,
    align       = FALSE
  )$matrices

  for (key in names(x$edges$matrices)) {
    E_alg <- x$edges$matrices[[key]]
    E_sc <- E_scores[[key]]
    expect_lt(
      max(abs(abs(E_alg) - abs(E_sc))), 1e-4,
      label = paste("ESEM algebra vs scores for pair", key)
    )
  }
})

# ── Convergence truncation ────────────────────────────────────────────────────

test_that("ESEM warns and truncates when deep levels fail to converge", {
  skip_if_not_installed("lavaan")
  # With 6 variables, lavaan::efa() errors for k >= 4 ("maximum number of
  # factors is 3"). Requesting k = 5 exercises the warn-and-truncate path.
  d <- .make_esem_data()
  warn_msgs <- character(0L)
  withCallingHandlers(
    x <- ackwards(d, k_max = 5, engine = "esem"),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Engine emits a truncation warning and ackwards() emits a hierarchy warning
  expect_true(any(grepl("ESEM failed", warn_msgs)),
    info = "engine-level warning for failed levels"
  )
  expect_true(any(grepl("[Tt]runcated|truncated|did not converge", warn_msgs)),
    info = "ackwards-level hierarchy truncation warning"
  )

  # Object still builds to deepest converged level (k_eff = 3)
  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$k_max, 3L)
  expect_true(all(vapply(x$levels, `[[`, logical(1), "converged")))
})

# ── tidy / glance / print ─────────────────────────────────────────────────────

test_that("print, tidy, glance work for ESEM objects", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))

  expect_no_error(print(x))
  expect_s3_class(generics::tidy(x), "data.frame")
  expect_s3_class(generics::glance(x), "data.frame")
  expect_equal(nrow(generics::glance(x)), 1L)
})

test_that("tidy(x, what = 'fit') returns 8 indices per level for ESEM", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  td <- generics::tidy(x, what = "fit")
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 3L * 8L) # 3 levels × 8 indices (incl. BIC)
  expect_equal(sort(unique(td$level)), 1:3)
})

# ── M27: loading SEs/CIs folded into tidy(what = "loadings") ─────────────────

test_that("tidy(what='loadings') for ESEM has se/ci_lower/ci_upper populated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  ld <- generics::tidy(x, what = "loadings")
  expect_true(all(c("se", "ci_lower", "ci_upper") %in% names(ld)))
  expect_true(all(is.finite(ld$se)))
  expect_true(all(ld$se >= 0))
  # CIs bracket the loading
  expect_true(all(ld$ci_lower <= ld$loading))
  expect_true(all(ld$ci_upper >= ld$loading))
})

test_that("tidy(what='loadings', conf_level=0.99) widens intervals vs 0.95", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 2, engine = "esem"))
  ld95 <- generics::tidy(x, what = "loadings", conf_level = 0.95)
  ld99 <- generics::tidy(x, what = "loadings", conf_level = 0.99)
  width95 <- ld95$ci_upper - ld95$ci_lower
  width99 <- ld99$ci_upper - ld99$ci_lower
  expect_true(all(width99 >= width95, na.rm = TRUE))
})

test_that("tidy(what='loadings') for PCA/EFA has se/ci cols present but all NA", {
  skip_if_not_installed("psych")
  suppressWarnings(x_pca <- ackwards(psych::bfi[, 1:25], k_max = 2))
  suppressWarnings(x_efa <- ackwards(psych::bfi[, 1:25], k_max = 2, engine = "efa"))
  for (x in list(x_pca, x_efa)) {
    ld <- generics::tidy(x, what = "loadings")
    expect_true(all(c("se", "ci_lower", "ci_upper") %in% names(ld)))
    expect_true(all(is.na(ld$se)))
    expect_true(all(is.na(ld$ci_lower)))
    expect_true(all(is.na(ld$ci_upper)))
  }
})

test_that("conf_level passed with wrong what= errors", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "edges", conf_level = 0.90), "conf_level")
})

# ── M27: wide fit table ───────────────────────────────────────────────────────

test_that("tidy(what='fit', format='wide') gives one row per non-anchor level", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  wide <- generics::tidy(x, what = "fit", format = "wide")
  expect_s3_class(wide, "data.frame")
  # Anchor (k = 1) dropped: levels 2 and 3 only.
  expect_equal(nrow(wide), 2L)
  expect_false(1L %in% wide$level)
  expect_equal(sort(wide$level), c(2L, 3L))
  expect_true("level" %in% names(wide))
  expect_true(all(c("CFI", "TLI", "RMSEA", "SRMR") %in% names(wide)))
})

test_that("tidy(what='fit', format='long') is byte-identical to the default", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  expect_identical(
    generics::tidy(x, what = "fit"),
    generics::tidy(x, what = "fit", format = "long")
  )
})

test_that("format passed with wrong what= errors", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "edges", format = "wide"), "format")
})

# ── M27: cutoff flags ─────────────────────────────────────────────────────────

test_that("tidy(what='fit', cutoffs=TRUE) adds meets column for ESEM", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  td <- generics::tidy(x, what = "fit", cutoffs = TRUE)
  expect_true("meets" %in% names(td))
  # CFI/TLI/RMSEA/SRMR rows should have non-NA meets
  idx_with_cutoff <- c("CFI", "TLI", "RMSEA", "SRMR")
  meets_for_cutoff <- td$meets[td$index %in% idx_with_cutoff]
  expect_true(all(!is.na(meets_for_cutoff)))
  # chi/dof/p_value get NA
  meets_no_cutoff <- td$meets[td$index %in% c("chi", "dof", "p_value")]
  expect_true(all(is.na(meets_no_cutoff)))
})

test_that("tidy(what='fit', cutoffs=FALSE) has no meets column", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  td <- generics::tidy(x, what = "fit")
  expect_false("meets" %in% names(td))
})

test_that("cutoffs=TRUE with format='wide' produces *_meets columns only for thresholded indices", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  wide <- generics::tidy(x, what = "fit", format = "wide", cutoffs = TRUE)
  expect_true("CFI_meets" %in% names(wide))
  expect_true("TLI_meets" %in% names(wide))
  expect_true("RMSEA_meets" %in% names(wide))
  expect_true("SRMR_meets" %in% names(wide))
  # chi/dof/p_value/BIC have no defined threshold; the pivot must not
  # generate an always-NA *_meets column for them (M31 cleanup).
  expect_false("chi_meets" %in% names(wide))
  expect_false("dof_meets" %in% names(wide))
  expect_false("p_value_meets" %in% names(wide))
  expect_false("BIC_meets" %in% names(wide))
})

test_that("cutoffs passed with wrong what= errors", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  expect_error(tidy(x, what = "edges", cutoffs = TRUE), "cutoffs")
})


# ── estimator argument ────────────────────────────────────────────────────────

test_that("estimator argument is validated", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  expect_error(
    ackwards(d, k_max = 2, engine = "esem", estimator = "bad"),
    "estimator"
  )
})

# ── Polychoric basis: PCA and EFA ─────────────────────────────────────────────

test_that("cor = 'polychoric' works for PCA engine", {
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "pca", cor = "polychoric"))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$cor, "polychoric")
  expect_equal(length(x$levels), 3L)

  # Scoring basis recorded correctly
  for (ki in seq_len(x$k_max)) {
    expect_equal(x$levels[[as.character(ki)]]$scoring$basis, "polychoric",
      info = paste("level", ki, "basis")
    )
  }

  # R stored in result is a p×p correlation matrix
  expect_true(is.matrix(x$r))
  expect_equal(dim(x$r), c(6L, 6L))
  expect_true(all(abs(x$r) <= 1 + 1e-9))
})

test_that("cor = 'polychoric' works for EFA engine", {
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "efa", cor = "polychoric"))

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$cor, "polychoric")
  expect_equal(length(x$levels), 3L)

  for (ki in seq_len(x$k_max)) {
    expect_equal(x$levels[[as.character(ki)]]$scoring$basis, "polychoric",
      info = paste("level", ki, "basis")
    )
  }
})

# ── Polychoric basis: ESEM ────────────────────────────────────────────────────

test_that("cor = 'polychoric' with method = 'esem' uses WLSMV and returns valid object", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("psych") # needed for detect_ordinal helper path
  d <- .make_ordinal_data()
  suppressWarnings(
    x <- ackwards(d, k_max = 3, engine = "esem", cor = "polychoric")
  )

  expect_s3_class(x, "ackwards")
  validate_ackwards(x)
  expect_equal(x$cor, "polychoric")
  expect_equal(length(x$levels), 3L)

  # loadings_se populated
  for (ki in seq_len(x$k_max)) {
    lev <- x$levels[[as.character(ki)]]
    expect_false(is.null(lev$loadings_se),
      info = paste("loadings_se non-NULL level", ki)
    )
    expect_equal(x$levels[[as.character(ki)]]$scoring$basis, "polychoric",
      info = paste("basis polychoric level", ki)
    )
  }
})

# ── M31: cor = "polychoric" + incompatible estimator guard ───────────────────

test_that("cor = 'polychoric' with estimator = 'ML'/'MLR' errors loudly", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  expect_error(
    ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", estimator = "ML"),
    "incompatible"
  )
  expect_error(
    ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", estimator = "MLR"),
    "incompatible"
  )
})

test_that("cor = 'polychoric' with estimator = 'WLSMV'/'ULSMV' is unaffected by the guard", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  expect_no_error(
    suppressWarnings(
      ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", estimator = "WLSMV")
    )
  )
  expect_no_error(
    suppressWarnings(
      ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", estimator = "ULSMV")
    )
  )
})

test_that("estimator = 'WLSMV' with a continuous cor is allowed (not guarded)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  expect_no_error(
    suppressWarnings(
      ackwards(d, k_max = 2, engine = "esem", cor = "pearson", estimator = "WLSMV")
    )
  )
})

# ── Ordinal warning suppressed when cor = "polychoric" ───────────────────────

test_that("ordinal warning is NOT emitted when cor = 'polychoric'", {
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  expect_no_warning(
    suppressMessages(ackwards(d, k_max = 2, engine = "pca", cor = "polychoric")),
    message = "ordinal"
  )
})

# ── Heywood / improper-solution warning ───────────────────────────────────────

test_that("ESEM warns on improper solution (Heywood case) but still builds", {
  skip_if_not_installed("lavaan")
  # 8 variables with near-zero unique variance (eps=0.001), 2 true factors.
  # k=3 pushes one residual to the Heywood boundary (theta <= 0 for one variable).
  # lavaan clamps negative theta to 0 by default; our check fires at <= 0.
  # The skip_if below guards against future lavaan versions that may change this
  # bounding behaviour, which would make the fixture no longer trigger the warning.
  set.seed(42)
  n <- 200
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  eps <- 0.001
  d <- as.data.frame(cbind(
    x1 = sqrt(1 - eps) * f1 + sqrt(eps) * rnorm(n),
    x2 = sqrt(1 - eps) * f1 + sqrt(eps) * rnorm(n),
    x3 = sqrt(1 - eps) * f1 + sqrt(eps) * rnorm(n),
    x4 = sqrt(1 - eps) * f1 + sqrt(eps) * rnorm(n),
    x5 = sqrt(1 - eps) * f2 + sqrt(eps) * rnorm(n),
    x6 = sqrt(1 - eps) * f2 + sqrt(eps) * rnorm(n),
    x7 = sqrt(1 - eps) * f2 + sqrt(eps) * rnorm(n),
    x8 = sqrt(1 - eps) * f2 + sqrt(eps) * rnorm(n)
  ))
  # Pre-check: verify this version of lavaan still produces theta <= 0 here.
  lav_fit <- tryCatch(
    suppressWarnings({
      raw <- lavaan::efa(data = d, nfactors = 3L, rotation = "varimax")
      if (inherits(raw, "efaList")) raw[[1L]] else raw
    }),
    error = function(e) NULL
  )
  theta_ok <- !is.null(lav_fit) && tryCatch(
    {
      any(diag(lavaan::lavInspect(lav_fit, "theta")) <= 0)
    },
    error = function(e) FALSE
  )
  skip_if(!theta_ok, "lavaan no longer produces theta <= 0 for this fixture")

  expect_warning(
    x <- suppressMessages(ackwards(d, k_max = 3, engine = "esem")),
    "Heywood"
  )
  # Object still builds to the requested depth (not truncated)
  expect_s3_class(x, "ackwards")
  expect_equal(x$k_max, 3L)
})

# ── cor = "spearman" + engine = "esem" inconsistency warning ─────────────────

test_that("cor = 'spearman' with engine = 'esem' warns about inconsistent bases", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  expect_warning(
    suppressMessages(ackwards(d, k_max = 2, engine = "esem", cor = "spearman")),
    "inconsistent bases"
  )
})

# ── M26: cached sample statistics + parallel dispatch ─────────────────────────

test_that(".esem_lapply falls back to serial lapply when future.apply is absent", {
  # Force the is_installed() check to report future.apply missing.
  testthat::local_mocked_bindings(
    is_installed = function(...) FALSE,
    .package = "rlang"
  )
  out <- .esem_lapply(1:3, function(x) x^2)
  expect_identical(out, list(1, 4, 9))
})

test_that(".esem_lapply uses future.apply when available", {
  skip_if_not_installed("future.apply")
  out <- .esem_lapply(1:3, function(x) x^2)
  expect_identical(out, list(1, 4, 9))
})

test_that("ESEM results are identical across serial and parallel future plans", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_os("windows")
  skip_if(!future::supportsMulticore(), "multicore plan unavailable here")

  d <- .make_ordinal_data()

  future::plan(future::sequential)
  serial <- suppressMessages(suppressWarnings(
    ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", seed = 42)
  ))

  future::plan(future::multicore, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)
  par <- suppressMessages(suppressWarnings(
    ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", seed = 42)
  ))

  for (ki in seq_along(serial$levels)) {
    expect_equal(
      serial$levels[[ki]]$loadings,
      par$levels[[ki]]$loadings,
      info = paste("loadings level", ki)
    )
  }
  expect_equal(tidy(serial)$r, tidy(par)$r)
})

# ── glance() fit columns for ESEM ─────────────────────────────────────────────

test_that("glance() for ESEM has CFI/TLI/RMSEA/SRMR/BIC at deepest level (ML: BIC available)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # Default cor = "pearson" -> estimator = "ML", which has a proper
  # log-likelihood, so BIC is a real number here.
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem"))
  g <- generics::glance(x)
  expect_true(all(c("CFI", "TLI", "RMSEA", "SRMR", "BIC") %in% names(g)))
  expect_false(is.na(g$CFI))
  expect_false(is.na(g$TLI))
  expect_false(is.na(g$RMSEA))
  expect_false(is.na(g$SRMR))
  expect_false(is.na(g$BIC))
  expect_equal(g$deepest_converged, 3L)
})

test_that("ESEM under WLSMV: p_value uses the scaled test when df > 0; BIC genuinely NA", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("psych") # needed for detect_ordinal helper path
  d <- .make_ordinal_data()
  suppressWarnings(
    x <- ackwards(d, k_max = 3, engine = "esem", cor = "polychoric")
  )

  for (ki in seq_len(x$k_max)) {
    fv <- x$levels[[as.character(ki)]]$fit
    if (fv[["dof"]] > 0) {
      # lavaan's naive WLSMV p-value has no valid reference distribution
      # (lavaan's own summary() reports it as "Unknown"/NA); ackwards falls
      # back to the mean-and-variance-adjusted scaled test, which does.
      expect_false(is.na(fv[["p_value"]]), info = paste("p_value level", ki))
      expect_true(fv[["p_value"]] >= 0 && fv[["p_value"]] <= 1,
        info = paste("p_value range level", ki)
      )
    } else {
      # A saturated model (df = 0, e.g. k = 3 on 6 items here) has no
      # chi-square test to perform under any estimator -- NA is correct.
      expect_true(is.na(fv[["p_value"]]), info = paste("saturated level", ki))
    }
    # WLSMV has no proper log-likelihood -> BIC is genuinely unavailable.
    expect_true(is.na(fv[["BIC"]]), info = paste("BIC level", ki))
  }

  g <- generics::glance(x)
  expect_true(is.na(g$BIC))
})

test_that("ESEM under ULSMV: p_value uses the scaled test when df > 0; BIC genuinely NA", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  suppressWarnings(
    x <- ackwards(d, k_max = 3, engine = "esem", cor = "polychoric", estimator = "ULSMV")
  )
  for (ki in seq_len(x$k_max)) {
    fv <- x$levels[[as.character(ki)]]$fit
    if (fv[["dof"]] > 0) {
      expect_false(is.na(fv[["p_value"]]), info = paste("p_value level", ki))
      expect_true(fv[["p_value"]] >= 0 && fv[["p_value"]] <= 1)
    } else {
      expect_true(is.na(fv[["p_value"]]), info = paste("saturated level", ki))
    }
    # ULSMV, like WLSMV, is limited-information -> no BIC.
    expect_true(is.na(fv[["BIC"]]), info = paste("BIC level", ki))
  }
})

test_that("ESEM under MLR: fit row reports the scaled (not naive) test + indices; BIC real", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(x <- ackwards(d, k_max = 3, engine = "esem", estimator = "MLR"))

  # Reference: fit level 2 directly and confirm ackwards reports the *scaled*
  # variant, i.e. it matches chisq.scaled/cfi.scaled and (when they differ)
  # not the naive chisq/cfi.
  fit2 <- lavaan::efa(
    data = as.data.frame(d), nfactors = 2L, rotation = "varimax",
    estimator = "MLR", missing = "listwise"
  )[[1L]]
  fm <- lavaan::fitMeasures(
    fit2, c("chisq", "chisq.scaled", "cfi", "cfi.scaled", "bic")
  )
  fv2 <- x$levels[["2"]]$fit
  expect_equal(fv2[["chi"]], unname(fm[["chisq.scaled"]]), tolerance = 1e-6)
  expect_equal(fv2[["CFI"]], unname(fm[["cfi.scaled"]]), tolerance = 1e-6)
  # MLR retains a proper log-likelihood -> BIC is a real number.
  expect_false(is.na(fv2[["BIC"]]))
  expect_equal(fv2[["BIC"]], unname(fm[["bic"]]), tolerance = 1e-6)
  expect_false(is.na(generics::glance(x)$BIC))
})

test_that("ESEM with a continuous cor + WLSMV is allowed and yields a populated fit row", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  suppressWarnings(
    x <- ackwards(d, k_max = 3, engine = "esem", cor = "pearson", estimator = "WLSMV")
  )
  fv <- x$levels[["2"]]$fit
  expect_named(fv, c("chi", "dof", "p_value", "CFI", "TLI", "RMSEA", "SRMR", "BIC"))
  # A non-saturated level carries real fit indices.
  expect_false(is.na(fv[["CFI"]]))
  expect_false(is.na(fv[["RMSEA"]]))
  # Continuous WLSMV is still a limited-information estimator -> no BIC.
  expect_true(is.na(fv[["BIC"]]))
})
