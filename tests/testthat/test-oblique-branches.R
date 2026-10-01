# Oblique failure branches: the regression-weight fallbacks, a factor
# correlation that cannot be read, a non-finite factor correlation, and a
# failed rotation. Each is reached by planting the failure in the step that
# fails, and each test asserts which warning fired and what the level became.

# Warning messages with cli's styling and line wrapping removed.
.warnings_of <- function(expr) {
  gsub("\\s+", " ", cli::ansi_strip(testthat::capture_warnings(expr)))
}

# ── Regression-weight fallbacks ──────────────────────────────────────────────

test_that("efa: an oblique level falls back to psych's R^-1 L Phi weights", {
  skip_if_not_installed("GPArotation")
  local_mocked_bindings(.tenBerge_weights = function(R, L, Phi) stop("planted failure"))
  w <- .warnings_of(
    x <- suppressMessages(ackwards(sim16, k_max = 3, engine = "efa", rotation = "oblimin"))
  )
  expect_match(w, "tenBerge weights failed at k = 2: planted failure", fixed = TRUE, all = FALSE)
  expect_match(w, "Falling back to regression (Thurstone) weights.", fixed = TRUE, all = FALSE)
  for (k in c("2", "3")) {
    lev <- x$levels[[k]]
    expect_identical(lev$scoring$method, "regression")
    expect_gt(max(abs(lev$factor_cor[upper.tri(lev$factor_cor)])), 0.05)
    # The oblique regression rule, in the stored order and signs.
    expect_equal(unname(lev$scoring$weights),
      unname(solve(x$r) %*% lev$loadings %*% lev$factor_cor),
      tolerance = 1e-8
    )
  }
})

test_that("esem: an oblique level falls back to R^-1 L Phi weights", {
  skip_if_not_installed("lavaan")
  local_mocked_bindings(
    .tenBerge_weights = function(R, L, Phi) stop("planted failure"),
    .esem_lapply = function(X, FUN) lapply(X, FUN)
  )
  w <- .warnings_of(
    x <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", rotation = "oblimin", seed = 1))
  )
  expect_match(w, "tenBerge weights failed at k = 2; used regression (Thurstone) weights instead.",
    fixed = TRUE, all = FALSE
  )
  expect_length(x$levels, 3L)
  for (k in c("2", "3")) {
    lev <- x$levels[[k]]
    expect_identical(lev$scoring$method, "regression")
    expect_gt(max(abs(lev$factor_cor[upper.tri(lev$factor_cor)])), 0.05)
    expect_equal(unname(lev$scoring$weights),
      unname(solve(x$r) %*% lev$loadings %*% lev$factor_cor),
      tolerance = 1e-8
    )
  }
})

test_that("esem: a level with no usable weights truncates", {
  skip_if_not_installed("lavaan")
  real <- .regression_weights
  local_mocked_bindings(
    .tenBerge_weights = function(R, L, Phi) stop("planted failure"),
    .regression_weights = function(R, L, Phi) {
      if (ncol(L) >= 3L) stop("planted failure") else real(R, L, Phi)
    },
    .esem_lapply = function(X, FUN) lapply(X, FUN)
  )
  w <- .warnings_of(
    x <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", rotation = "oblimin", seed = 1))
  )
  expect_match(w, "ESEM failed at k = 3: could not compute scoring weights", fixed = TRUE, all = FALSE)
  expect_match(w, "Truncating hierarchy at level 2.", fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), c("1", "2"))
})

# ── lavaan's factor correlation (.esem_read_phi) ─────────────────────────────

test_that("esem: an unreadable factor correlation truncates an oblique fit only", {
  skip_if_not_installed("lavaan")
  # Serial dispatch for both fits: lavaan's random rotation starts draw from
  # the global stream here and from per-task streams under future.apply.
  local_mocked_bindings(.esem_lapply = function(X, FUN) lapply(X, FUN))
  ref <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", seed = 1))
  real <- .esem_read_phi
  local_mocked_bindings(
    .esem_read_phi = function(fit, factors_lav) {
      if (length(factors_lav) >= 3L) stop("planted failure") else real(fit, factors_lav)
    }
  )
  w <- .warnings_of(
    x <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", rotation = "oblimin", seed = 1))
  )
  expect_match(w, "ESEM failed at k = 3: could not extract the factor correlations",
    fixed = TRUE, all = FALSE
  )
  expect_identical(names(x$levels), c("1", "2"))

  # Under varimax the same failure takes the identity, and the fit matches an
  # unplanted one.
  w0 <- .warnings_of(x0 <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", seed = 1)))
  expect_false(any(grepl("factor correlations", w0, fixed = TRUE)))
  expect_identical(names(x0$levels), c("1", "2", "3"))
  expect_equal(unname(x0$levels[["3"]]$factor_cor), diag(3L))
  expect_equal(x0$edges$matrices, ref$edges$matrices, tolerance = 1e-10)
})

test_that(".esem_read_phi() matches lavaan's factors by name and rejects what it cannot match", {
  skip_if_not_installed("lavaan")
  x <- cached(ackwards(sim16,
    k_max = 3, engine = "esem", rotation = "oblimin",
    keep_fits = TRUE, seed = 1
  ))
  fit <- x$fits[["3"]]
  Phi <- lavaan::lavInspect(fit, "cor.lv")
  nm <- rownames(Phi)
  expect_gt(max(abs(Phi[upper.tri(Phi)])), 0.05)
  # Rows follow the names asked for, not their position in cor.lv.
  expect_equal(.esem_read_phi(fit, rev(nm)), unname(Phi[rev(nm), rev(nm)]))
  expect_error(.esem_read_phi(fit, c("g1", "g2", "g3")), "factor names", fixed = TRUE)

  local_mocked_bindings(lavInspect = function(object, what, ...) unname(Phi), .package = "lavaan")
  expect_error(.esem_read_phi(fit, nm), "factor names", fixed = TRUE)

  bad <- Phi
  bad[1L, 2L] <- NaN
  local_mocked_bindings(lavInspect = function(object, what, ...) bad, .package = "lavaan")
  expect_error(.esem_read_phi(fit, nm), "not finite", fixed = TRUE)
})

# ── A non-finite psych factor correlation ────────────────────────────────────

test_that("a non-finite factor correlation truncates PCA and EFA levels", {
  skip_if_not_installed("GPArotation")
  expect_false(.near_identity(matrix(NaN, 2L, 2L)))
  expect_error(
    .engine_phi(list(Phi = matrix(c(1, NaN, NaN, 1), 2L)), 2L),
    "the factor correlations are not finite.",
    fixed = TRUE
  )

  real_pca <- psych::pca
  real_fa <- psych::fa
  plant_nan <- function(f, nfactors) {
    if (nfactors == 3L) f$Phi[1L, 2L] <- NaN
    f
  }
  local_mocked_bindings(
    pca = function(r, nfactors = 1, ...) plant_nan(real_pca(r, nfactors = nfactors, ...), nfactors),
    fa = function(r, nfactors = 1, ...) plant_nan(real_fa(r, nfactors = nfactors, ...), nfactors),
    .package = "psych"
  )
  for (engine in c("pca", "efa")) {
    w <- .warnings_of(
      x <- suppressMessages(ackwards(sim16, k_max = 4, engine = engine, rotation = "oblimin"))
    )
    expect_match(w,
      paste(toupper(engine), "failed at k = 3: the factor correlations are not finite."),
      fixed = TRUE, all = FALSE
    )
    expect_identical(names(x$levels), c("1", "2"))
  }
})

# ── A failed rotation (start count decides) ──────────────────────────────────

test_that("pca: a failed oblique rotation truncates; other warnings pass through", {
  skip_if_not_installed("GPArotation")
  real_pca <- psych::pca
  plant <- function(msg) {
    local_mocked_bindings(
      pca = function(r, nfactors = 1, ...) {
        if (nfactors == 3L) warning(msg)
        real_pca(r, nfactors = nfactors, ...)
      },
      .package = "psych",
      .env = parent.frame()
    )
  }
  msgs <- c(
    "Convergence not obtained in GPFoblq. 1000 iterations used.",
    "The requested transformaton failed, Promax was used instead as an oblique transformation"
  )
  for (msg in msgs) {
    plant(msg)
    w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, rotation = "oblimin")))
    expect_match(w, paste0("The oblimin rotation failed at k = 3: ", msg), fixed = TRUE, all = FALSE)
    expect_identical(names(x$levels), c("1", "2"))
  }

  # Control: a warning that is not about the rotation reaches the user, and
  # the level is kept.
  plant("an unrelated psych warning")
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, rotation = "oblimin")))
  expect_true("an unrelated psych warning" %in% w)
  expect_false(any(grepl("rotation failed", w, fixed = TRUE)))
  expect_identical(names(x$levels), c("1", "2", "3", "4"))
})

test_that("efa: a rotation warning is shown and the level kept; a failed final step truncates", {
  skip_if_not_installed("GPArotation")
  real_fa <- psych::fa
  plant <- function(msg = NULL, drop_rot_mat = FALSE) {
    local_mocked_bindings(
      fa = function(r, nfactors = 1, ...) {
        if (nfactors == 3L && !is.null(msg)) warning(msg)
        f <- real_fa(r, nfactors = nfactors, ...)
        if (nfactors == 3L && drop_rot_mat) f$rot.mat <- NULL
        f
      },
      .package = "psych",
      .env = parent.frame()
    )
  }
  conv <- "Convergence not obtained in GPFoblq. 1000 iterations used."
  plant(conv)
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa", rotation = "oblimin")))
  expect_match(w, paste0("psych reported a rotation problem at k = 3: ", conv), fixed = TRUE, all = FALSE)
  expect_match(w, "The level is kept.", fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), c("1", "2", "3", "4"))

  # Control: under varimax the same warning stays muffled, as before.
  w0 <- .warnings_of(x0 <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa")))
  expect_false(any(grepl("rotation", w0, fixed = TRUE)))
  expect_identical(names(x0$levels), c("1", "2", "3", "4"))

  # A failed final step leaves no rotation matrix, with or without psych's
  # Promax notice.
  promax_note <- "The requested transformaton failed, Promax was used instead as an oblique transformation"
  plant(promax_note, drop_rot_mat = TRUE)
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa", rotation = "oblimin")))
  expect_match(w, paste0("The oblimin rotation failed at k = 3: ", promax_note), fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), c("1", "2"))

  plant(drop_rot_mat = TRUE)
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa", rotation = "oblimin")))
  expect_match(w, "The oblimin rotation failed at k = 3: psych returned no rotation matrix.",
    fixed = TRUE, all = FALSE
  )
  expect_identical(names(x$levels), c("1", "2"))
})

test_that("esem: lavaan's rotation warning is shown under an oblique rotation only", {
  skip_if_not_installed("lavaan")
  real_efa <- lavaan::efa
  conv <- "GP rotation algorithm did not converge after 10000 iterations"
  local_mocked_bindings(
    efa = function(..., nfactors) {
      if (nfactors == 3L) warning(conv)
      real_efa(..., nfactors = nfactors)
    },
    .package = "lavaan"
  )
  local_mocked_bindings(.esem_lapply = function(X, FUN) lapply(X, FUN))
  w <- .warnings_of(
    x <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", rotation = "oblimin", seed = 1))
  )
  expect_match(w, paste0("lavaan reported a rotation problem at k = 3: ", conv), fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), c("1", "2", "3"))

  # Control: under varimax lavaan's warnings stay muffled, as before.
  w0 <- .warnings_of(x0 <- suppressMessages(ackwards(sim16, k_max = 3, engine = "esem", seed = 1)))
  expect_false(any(grepl("rotation problem", w0, fixed = TRUE)))
  expect_identical(names(x0$levels), c("1", "2", "3"))
})
