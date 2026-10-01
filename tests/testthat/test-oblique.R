# Oblique rotation: the within-level factor correlation carried end to end,
# correlation-preserving weights, and the algebra-vs-scores oracle (IP2).

# Orthogonal ten Berge weights, written out here independently of
# .tenBerge_weights(): W = R^{-1} L (L' R^{-1} L)^{-1/2}.
.tenberge_orthogonal <- function(R, L) {
  A <- solve(R, L)
  e <- eigen(crossprod(L, A), symmetric = TRUE)
  A %*% e$vectors %*% diag(1 / sqrt(e$values), nrow = ncol(L)) %*% t(e$vectors)
}

# The +/-1 flips ackwards() applied to an engine's columns: stored = engine * s.
.applied_signs <- function(stored, engine) {
  sign(colSums(unname(stored) * unname(engine)))
}

# ── .tenBerge_weights() ──────────────────────────────────────────────────────

test_that(".tenBerge_weights() at Phi = I reproduces the orthogonal formula", {
  R <- cor(sim16)
  f <- psych::fa(R, nfactors = 3, rotate = "varimax", n.obs = nrow(sim16))
  L <- unclass(f$loadings)
  W <- .tenBerge_weights(R, L, diag(3L))
  expect_equal(unname(W), unname(.tenberge_orthogonal(R, L)), tolerance = 1e-12)
  # Live oracle on the same fit: psych's own tenBerge weights.
  W_psych <- unclass(psych::factor.scores(R, f, method = "tenBerge")$weights)
  expect_equal(unname(W), unname(W_psych), tolerance = 1e-8)
})

test_that(".tenBerge_weights() matches psych's oblique tenBerge weights", {
  skip_if_not_installed("GPArotation")
  R <- cor(sim16)
  for (rotation in c("oblimin", "promax")) {
    f <- psych::fa(R, nfactors = 3, rotate = rotation, n.obs = nrow(sim16))
    L <- unclass(f$loadings)
    Phi <- f$Phi
    expect_gt(max(abs(Phi[upper.tri(Phi)])), 0.05) # genuinely oblique
    W <- .tenBerge_weights(R, L, Phi)
    W_psych <- unclass(psych::factor.scores(R, f, method = "tenBerge")$weights)
    expect_equal(unname(W), unname(W_psych), tolerance = 1e-8, label = rotation)
    # Correlation preserving: the scores' covariance W'RW is Phi itself.
    expect_equal(unname(crossprod(W, R %*% W)), unname(Phi), tolerance = 1e-8)
    # The orthogonal formula on the same oblique pattern does not preserve it.
    W_orth <- .tenberge_orthogonal(R, L)
    expect_gt(max(abs(crossprod(W_orth, R %*% W_orth) - Phi)), 0.05)
  }
})

test_that(".tenBerge_weights() errors on a factor correlation that is not positive definite", {
  R <- cor(sim16)[1:4, 1:4]
  L <- matrix(c(0.8, 0.7, 0.1, 0.1, 0.1, 0.1, 0.8, 0.7), ncol = 2L)
  bad_phi <- matrix(c(1, 1.2, 1.2, 1), 2L)
  expect_error(
    .tenBerge_weights(R, L, bad_phi),
    "not positive definite"
  )
})

# ── PCA and EFA: the stored factor_cor is psych's Phi, in step with the loadings

for (engine in c("pca", "efa")) {
  test_that(paste0(engine, ": oblique factor_cor is psych's Phi, ordered and flipped with the loadings"), {
    skip_if_not_installed("GPArotation")
    x <- cached(ackwards(sim16,
      k_max = 4, engine = engine, rotation = "oblimin",
      keep_fits = TRUE
    ))
    for (ki in 2:4) {
      lev <- x$levels[[as.character(ki)]]
      fit <- x$fits[[as.character(ki)]]
      L_psych <- unclass(fit$loadings)
      s <- .applied_signs(lev$loadings, L_psych)
      # Same column order and signs as the stored loadings ...
      expect_equal(unname(lev$loadings), unname(sweep(L_psych, 2, s, "*")),
        tolerance = 1e-12
      )
      # ... and factor_cor is psych's Phi carried by those same flips.
      expect_equal(unname(lev$factor_cor), unname(fit$Phi * tcrossprod(s)),
        tolerance = 1e-12
      )
      expect_identical(dimnames(lev$factor_cor), list(lev$labels, lev$labels))
      off <- lev$factor_cor[upper.tri(lev$factor_cor)]
      expect_gt(max(abs(off)), 0.05) # not the identity
    }
  })

  test_that(paste0(engine, ": oblique scores correlate as factor_cor says"), {
    skip_if_not_installed("GPArotation")
    x <- cached(ackwards(sim16, k_max = 4, engine = engine, rotation = "oblimin"))
    R <- x$r
    for (ki in 2:4) {
      lev <- x$levels[[as.character(ki)]]
      W <- lev$scoring$weights
      expect_equal(unname(stats::cov2cor(crossprod(W, R %*% W))),
        unname(lev$factor_cor),
        tolerance = 1e-8
      )
    }
  })

  test_that(paste0(engine, ": oblique algebra and scores paths agree (IP2)"), {
    skip_if_not_installed("GPArotation")
    x <- cached(ackwards(sim16,
      k_max = 4, engine = engine, rotation = "oblimin",
      pairs = "all"
    ))
    E_scores <- compute_edges(
      levels = x$levels, R = x$r, edge_method = "scores",
      pairs = "all", data = sim16
    )$matrices
    for (key in names(x$edges$matrices)) {
      expect_equal(x$edges$matrices[[key]], E_scores[[key]],
        tolerance = 1e-8, label = paste(engine, "algebra vs scores", key)
      )
    }
  })
}

# ── ESEM: the stored factor_cor is lavaan's cor.lv, permuted and flipped ─────

# lavaan's standardized pattern loadings, columns in cor.lv's factor order.
.lavaan_pattern <- function(fit, items) {
  Phi <- lavaan::lavInspect(fit, "cor.lv")
  std <- lavaan::standardizedSolution(fit, type = "std.all")
  std <- std[std$op == "=~", , drop = FALSE]
  L <- matrix(0, length(items), nrow(Phi), dimnames = list(items, rownames(Phi)))
  for (f in rownames(Phi)) {
    rows <- std[std$lhs == f, , drop = FALSE]
    L[rows$rhs, f] <- rows$est.std
  }
  list(L = L, Phi = unname(Phi))
}

for (rotation in c("oblimin", "geomin")) {
  test_that(paste0("esem/", rotation, ": factor_cor is lavaan's cor.lv, ordered and flipped with the loadings"), {
    skip_if_not_installed("lavaan")
    x <- cached(ackwards(sim16,
      k_max = 4, engine = "esem", rotation = rotation,
      keep_fits = TRUE, seed = 1
    ))
    reordered <- FALSE
    for (ki in 2:4) {
      lev <- x$levels[[as.character(ki)]]
      lav <- .lavaan_pattern(x$fits[[as.character(ki)]], colnames(sim16))
      # Map each stored column to its lavaan column and sign.
      ord <- integer(ki)
      s <- numeric(ki)
      for (j in seq_len(ki)) {
        gap_pos <- colSums(abs(lev$loadings[, j] - lav$L))
        gap_neg <- colSums(abs(lev$loadings[, j] + lav$L))
        ord[j] <- which.min(pmin(gap_pos, gap_neg))
        s[j] <- if (gap_pos[ord[j]] <= gap_neg[ord[j]]) 1 else -1
      }
      expect_setequal(ord, seq_len(ki))
      reordered <- reordered || !identical(ord, seq_len(ki))
      expect_equal(unname(lev$loadings),
        unname(sweep(lav$L[, ord, drop = FALSE], 2, s, "*")),
        tolerance = 1e-10
      )
      expect_equal(unname(lev$factor_cor),
        lav$Phi[ord, ord] * tcrossprod(s),
        tolerance = 1e-10
      )
      off <- lev$factor_cor[upper.tri(lev$factor_cor)]
      expect_gt(max(abs(off)), 0.05) # not the identity
      # Correlation preserving: the stored weights' scores reproduce it.
      W <- lev$scoring$weights
      expect_identical(lev$scoring$method, "tenBerge")
      expect_equal(unname(stats::cov2cor(crossprod(W, x$r %*% W))),
        unname(lev$factor_cor),
        tolerance = 1e-8
      )
    }
    # Some level's stored order differs from lavaan's, so the checks above
    # saw the permutation carried and not only the identity.
    expect_true(reordered)
  })
}

test_that("oblique pattern and factor correlation reproduce the varimax common part", {
  skip_if_not_installed("GPArotation")
  skip_if_not_installed("lavaan")
  # A rotation leaves L Phi L' unchanged, so each oblique level's stored
  # loadings and factor_cor must give the varimax fit's L L'. This holds only
  # if factor_cor is in the loadings' column order and signs, and is checked
  # here against the varimax fit, not against psych's or lavaan's own Phi.
  specs <- list(
    c("pca", "oblimin"), c("pca", "promax"), c("efa", "oblimin"),
    c("efa", "promax"), c("esem", "oblimin"), c("esem", "geomin")
  )
  for (s in specs) {
    x1 <- cached(ackwards(sim16, k_max = 4, engine = s[[1]], rotation = s[[2]], seed = 1))
    x0 <- cached(ackwards(sim16, k_max = 4, engine = s[[1]], seed = 1))
    for (ki in c("2", "3", "4")) {
      l1 <- x1$levels[[ki]]
      l0 <- x0$levels[[ki]]
      expect_gt(max(abs(l1$factor_cor[upper.tri(l1$factor_cor)])), 0.05)
      expect_equal(
        unname(l1$loadings %*% l1$factor_cor %*% t(l1$loadings)),
        unname(tcrossprod(l0$loadings)),
        tolerance = 1e-6, label = paste(s[[1]], s[[2]], "level", ki)
      )
      # So the variance values sum to the varimax total.
      expect_equal(l1$variance[["cumulative"]], l0$variance[["cumulative"]],
        tolerance = 1e-6, label = paste(s[[1]], s[[2]], "level", ki)
      )
    }
  }
})

test_that("esem: oblique algebra and scores paths agree (IP2)", {
  skip_if_not_installed("lavaan")
  x <- cached(ackwards(sim16,
    k_max = 4, engine = "esem", rotation = "oblimin",
    pairs = "all", seed = 1
  ))
  E_scores <- compute_edges(
    levels = x$levels, R = x$r, edge_method = "scores",
    pairs = "all", data = sim16
  )$matrices
  for (key in names(x$edges$matrices)) {
    expect_equal(x$edges$matrices[[key]], E_scores[[key]],
      tolerance = 1e-8, label = paste("esem algebra vs scores", key)
    )
  }
})

# ── Variance explained under an oblique rotation ─────────────────────────────

test_that(".variance_explained() is colSums(L^2) / p at Phi = I", {
  L <- matrix(c(0.8, 0.7, 0.1, 0.2, 0.1, 0.2, 0.6, 0.7), ncol = 2L)
  v <- .variance_explained(L, 4L, c("a", "b"), diag(2L))
  expect_identical(names(v), c("a", "b", "cumulative"))
  expect_equal(unname(v[c("a", "b")]), colSums(L^2) / 4, tolerance = 1e-15)
  expect_equal(unname(v[["cumulative"]]), sum(L^2) / 4, tolerance = 1e-15)
})

for (engine in c("pca", "efa")) {
  test_that(paste0(engine, ": oblique variance equals psych's Vaccounted proportions"), {
    skip_if_not_installed("GPArotation")
    x <- cached(ackwards(sim16,
      k_max = 4, engine = engine, rotation = "oblimin",
      keep_fits = TRUE
    ))
    p <- ncol(sim16)
    for (ki in 2:4) {
      lev <- x$levels[[as.character(ki)]]
      fit <- x$fits[[as.character(ki)]]
      expect_equal(unname(lev$variance[lev$labels]),
        unname(fit$Vaccounted["Proportion Var", ]),
        tolerance = 1e-8
      )
      # Not vacuous: the squared-loading formula misses the oracle by far
      # more than its 1e-8 tolerance.
      expect_gt(max(abs(lev$variance[lev$labels] - colSums(lev$loadings^2) / p)), 1e-5)
    }
  })
}

test_that("esem: columns sort by the oblique variance key where it disagrees with squared loadings", {
  skip_if_not_installed("lavaan")
  x <- cached(ackwards(bfi25,
    k_max = 5, engine = "esem", rotation = "geomin",
    keep_fits = TRUE, seed = 1
  ))
  lev <- x$levels[["5"]]
  lav <- .lavaan_pattern(x$fits[["5"]], colnames(bfi25))
  key_oblique <- diag(lav$Phi %*% crossprod(lav$L))
  key_squares <- colSums(lav$L^2)
  ord <- order(key_oblique, decreasing = TRUE)
  # Precondition: on this level the two keys order the factors differently.
  expect_false(identical(ord, order(key_squares, decreasing = TRUE)))
  # The stored columns follow the oblique key ...
  s <- .applied_signs(lev$loadings, lav$L[, ord])
  expect_equal(unname(lev$loadings), unname(sweep(lav$L[, ord], 2, s, "*")),
    tolerance = 1e-10
  )
  # ... and the variance is that key over p, in descending order.
  expect_equal(unname(lev$variance[lev$labels]), unname(key_oblique[ord]) / ncol(bfi25),
    tolerance = 1e-10
  )
  expect_false(is.unsorted(rev(lev$variance[lev$labels])))
})

# ── Output surfaces under an oblique rotation (D-036) ────────────────────────

test_that("primary parents and signs use r under an oblique rotation", {
  skip_if_not_installed("GPArotation")
  x <- cached(ackwards(sim16, k_max = 4, rotation = "oblimin"))
  ed <- tidy(x, what = "edges")
  ed <- ed[ed$level_to == ed$level_from + 1L, , drop = FALSE]
  # beta is a separate quantity here: it differs from r on some edge.
  expect_gt(max(abs(ed$beta - ed$r)), 0.01)
  for (child in unique(ed$to)) {
    rows <- ed[ed$to == child, , drop = FALSE]
    primary <- rows[rows$is_primary, , drop = FALSE]
    expect_identical(nrow(primary), 1L)
    expect_identical(primary$from, rows$from[which.max(abs(rows$r))])
    expect_gt(primary$r, 0) # signs aligned on r
  }
})

test_that("print(), summary(), and autoplot() say that oblique edges are total correlations", {
  skip_if_not_installed("GPArotation")
  skip_if_not_installed("ggplot2")
  text_of <- function(obj) {
    cli::ansi_strip(paste(capture.output(print(obj), type = "message"), collapse = " "))
  }
  x <- cached(ackwards(sim16, k_max = 4, rotation = "oblimin"))
  for (txt in list(text_of(x), text_of(summary(x)))) {
    expect_match(txt, "Oblique rotation (oblimin): r is a total correlation", fixed = TRUE)
    expect_match(txt, "partialled coefficient is beta", fixed = TRUE)
  }
  p <- autoplot(x)
  expect_match(p$labels$caption, "total correlations (r) under the oblique oblimin", fixed = TRUE)

  # Control: a varimax fit carries none of these notes.
  x0 <- cached(ackwards(sim16, k_max = 4))
  expect_no_match(text_of(x0), "total correlation", fixed = TRUE)
  expect_no_match(text_of(summary(x0)), "total correlation", fixed = TRUE)
  expect_null(autoplot(x0)$labels$caption)
})

# ── Fit-time advisory and the prune() stance (D-036) ─────────────────────────

test_that("an oblique fit announces what its edges mean; a varimax fit does not", {
  skip_if_not_installed("GPArotation")
  msgs <- cli::ansi_strip(paste(
    testthat::capture_messages(
      suppressWarnings(ackwards(sim16, k_max = 3, rotation = "oblimin"))
    ),
    collapse = " "
  ))
  expect_match(msgs, "Oblique rotation (\"oblimin\")", fixed = TRUE)
  expect_match(msgs, "is a total correlation", fixed = TRUE)
  expect_match(msgs, "conventions were calibrated under varimax", fixed = TRUE)
  # The cost of reading lineage from r (the D-036 advisory).
  expect_match(gsub("\\s+", " ", msgs),
    "primary parent can then be a factor that only correlates with its real parent",
    fixed = TRUE
  )

  msgs0 <- cli::ansi_strip(paste(
    testthat::capture_messages(suppressWarnings(ackwards(sim16, k_max = 3))),
    collapse = " "
  ))
  expect_no_match(msgs0, "total correlation", fixed = TRUE)
})

test_that("prune() warns on an oblique object that its criterion assumes orthogonal levels", {
  skip_if_not_installed("GPArotation")
  x <- cached(ackwards(sim16, k_max = 4, rotation = "oblimin"))
  for (rules in list("redundant", "artifact", c("redundant", "artifact"))) {
    cnd <- expect_warning(
      xp <- prune(x, rules),
      "assumes? orthogonal levels",
      class = "rlang_warning"
    )
    # The warning names the rules that ran, and only those.
    msg <- gsub("\\s+", " ", cli::ansi_strip(conditionMessage(cnd)))
    for (r in c("redundant", "artifact")) {
      named <- grepl(paste0("\"", r, "\""), msg, fixed = TRUE)
      expect_identical(named, r %in% rules, label = paste(r, "named for", toString(rules)))
    }
    # The rules still ran, on r.
    expect_false(is.null(xp$prune))
  }
  # No rule runs, no warning: rules = "none" and manual-only flags.
  expect_no_warning(prune(x, "none"))
  expect_no_warning(prune(x, manual = "m4f1"))

  # Control: the same call on a varimax object stays silent on this point.
  x0 <- cached(ackwards(sim16, k_max = 4))
  expect_no_warning(prune(x0, "redundant"), message = "orthogonal levels")
})

test_that("EFA's oblique weights are labeled tenBerge and are not the orthogonal formula", {
  skip_if_not_installed("GPArotation")
  x <- cached(ackwards(sim16, k_max = 3, engine = "efa", rotation = "oblimin"))
  lev <- x$levels[["3"]]
  expect_identical(lev$scoring$method, "tenBerge")
  W_orth <- .tenberge_orthogonal(x$r, lev$loadings)
  expect_gt(max(abs(unname(lev$scoring$weights) - unname(W_orth))), 1e-3)
})
