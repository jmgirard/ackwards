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

test_that("EFA's oblique weights are labeled tenBerge and are not the orthogonal formula", {
  skip_if_not_installed("GPArotation")
  x <- cached(ackwards(sim16, k_max = 3, engine = "efa", rotation = "oblimin"))
  lev <- x$levels[["3"]]
  expect_identical(lev$scoring$method, "tenBerge")
  W_orth <- .tenberge_orthogonal(x$r, lev$loadings)
  expect_gt(max(abs(unname(lev$scoring$weights) - unname(W_orth))), 1e-3)
})
