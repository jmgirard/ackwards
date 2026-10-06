# Oblique failure branches: the regression-weight fallbacks, a factor
# correlation that cannot be read, a non-finite factor correlation, and a
# failed rotation. Each is reached by planting the failure in the step that
# fails. The fit tests assert which warning fired and what the level became,
# and the .esem_read_phi() and .esem_rotation_args() tests check the helpers
# directly. ESEM's counted rotation failures are tested in
# test-esem-rotation-starts.R.

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
  expect_match(w, "ESEM failed at k = 3: could not extract the factor correlations: planted failure",
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
    "Convergence not obtained in GPFoblq. 2000 iterations used.",
    # GPArotation's legacy algorithm writes the same warning in lowercase.
    "convergence not obtained in GPFoblq. 1000 iterations used.",
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

test_that("pca: an error inside an oblique fit truncates; under varimax it propagates", {
  skip_if_not_installed("GPArotation")
  # A real singular level: two redundant columns make the promax fit at
  # k = 9 fail inside psych::pca().
  d <- as.data.frame(sim16[, 1:8])
  d$dup1 <- d[[1]]
  d$dup2 <- d[[2]] + d[[3]]
  w <- .warnings_of(x <- suppressMessages(ackwards(d, k_max = 9, rotation = "promax")))
  expect_match(w, "PCA failed at k = 9:", fixed = TRUE, all = FALSE)
  expect_match(w, "Truncating hierarchy at level 8.", fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), as.character(1:8))

  # Varimax keeps its old behavior: an error inside psych::pca() propagates.
  real_pca <- psych::pca
  local_mocked_bindings(
    pca = function(r, nfactors = 1, ...) {
      if (nfactors == 3L) stop("planted psych error")
      real_pca(r, nfactors = nfactors, ...)
    },
    .package = "psych"
  )
  expect_error(suppressWarnings(ackwards(sim16, k_max = 4)), "planted psych error", fixed = TRUE)
})

test_that("efa: a rotation warning is shown and the level kept; a failed final step truncates", {
  skip_if_not_installed("GPArotation")
  real_fa <- psych::fa
  # The start count reads psych::fa's formals, so keep it on the real ones
  # while fa itself is replaced below (psych 2.6.5 rotates from 20 starts).
  real_starts <- .psych_fa_starts()
  skip_if(real_starts == 1L, "the installed psych::fa() rotates from one start")
  local_mocked_bindings(.psych_fa_starts = function(fa_formals = NULL) real_starts)
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
  conv <- "Convergence not obtained in GPFoblq. 2000 iterations used."
  plant(conv)
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa", rotation = "oblimin")))
  expect_match(w, paste0("psych reported a rotation problem at k = 3: ", conv), fixed = TRUE, all = FALSE)
  expect_match(w, "The level is kept.", fixed = TRUE, all = FALSE)
  expect_identical(names(x$levels), c("1", "2", "3", "4"))

  # Control: under varimax the same warning stays muffled, as before.
  w0 <- .warnings_of(x0 <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa")))
  expect_false(any(grepl("Convergence not obtained", w0, fixed = TRUE)))
  expect_false(any(grepl("rotation problem", w0, fixed = TRUE)))
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

test_that(".psych_fa_starts() reads psych::fa()'s n.rotations default", {
  expect_identical(.psych_fa_starts(list(n.rotations = 20)), 20L)
  expect_identical(.psych_fa_starts(list(n.rotations = 1)), 1L)
  # A psych without the argument rotates once.
  expect_identical(.psych_fa_starts(list(r = NULL)), 1L)
  expect_identical(.psych_fa_starts(), as.integer(formals(psych::fa)$n.rotations))
})

test_that("efa: with a one-start psych, a rotation warning truncates as it does for PCA", {
  skip_if_not_installed("GPArotation")
  local_mocked_bindings(.psych_fa_starts = function(fa_formals = NULL) 1L)
  real_fa <- psych::fa
  promax_note <- "The requested transformaton failed, Promax was used instead as an oblique transformation"
  local_mocked_bindings(
    fa = function(r, nfactors = 1, ...) {
      if (nfactors == 3L) warning(promax_note)
      real_fa(r, nfactors = nfactors, ...)
    },
    .package = "psych"
  )
  w <- .warnings_of(x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa", rotation = "oblimin")))
  # The rotation matrix is present, so only the start count sends it here.
  expect_match(w, paste0("The oblimin rotation failed at k = 3: ", promax_note), fixed = TRUE, all = FALSE)
  expect_false(any(grepl("rotation problem", w, fixed = TRUE)))
  expect_identical(names(x$levels), c("1", "2"))
})

test_that("efa: a net-negative one-factor solution is flipped in the loadings and the fallback weights", {
  skip_if_not_installed("GPArotation")
  real_fa <- psych::fa
  local_mocked_bindings(
    fa = function(r, nfactors = 1, ...) {
      f <- real_fa(r, nfactors = nfactors, ...)
      if (nfactors == 1L) {
        f$loadings <- -f$loadings
        f$weights <- -f$weights
      }
      f
    },
    .package = "psych"
  )
  local_mocked_bindings(.tenBerge_weights = function(R, L, Phi) stop("planted failure"))
  x <- suppressWarnings(suppressMessages(
    ackwards(sim16, k_max = 2, engine = "efa", rotation = "oblimin")
  ))
  lev <- x$levels[["1"]]
  expect_gt(sum(lev$loadings), 0)
  expect_identical(lev$scoring$method, "regression")
  # The regression weights were flipped with the loadings: R^-1 L at k = 1.
  expect_equal(unname(lev$scoring$weights), unname(solve(x$r) %*% lev$loadings), tolerance = 1e-8)
})

test_that(".esem_rotation_args() turns lavaan's rotation warnings on for every rotation", {
  # k = 1 is unrotated and passes the method name alone, in both lavaan forms.
  expect_identical(.esem_rotation_args("none", c("data", "rotation", "rotation_args")), list(rotation = "none"))
  expect_identical(.esem_rotation_args("none", c("data", "rotation", "rotation.args")), list(rotation = "none"))
  # lavaan >= 0.7: options inside `rotation`.
  expect_identical(
    .esem_rotation_args("varimax", c("data", "rotation", "rotation_args")),
    list(rotation = list("varimax", warn = TRUE))
  )
  expect_identical(
    .esem_rotation_args("oblimin", c("data", "rotation", "rotation_args")),
    list(rotation = list("oblimin", warn = TRUE))
  )
  # Earlier lavaan: options in `rotation.args`.
  expect_identical(
    .esem_rotation_args("varimax", c("data", "rotation", "rotation.args")),
    list(rotation = "varimax", rotation.args = list(warn = TRUE))
  )
  expect_identical(
    .esem_rotation_args("geomin", c("data", "rotation", "rotation.args")),
    list(rotation = "geomin", rotation.args = list(warn = TRUE))
  )
})
