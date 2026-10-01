# Phi-partialled edge reporting: `beta` on tidy(what = "edges") and `r2` on
# tidy(what = "variance"), both closed-form from the stored weights, R, and
# edge matrices (.partialled_edges).
#
# Oracles (cairn/ORACLES.md): the varimax identity beta == r and r2 ==
# colSums(r^2) on default fits of all three engines (invariant), and lm()
# standardized coefficients / R-squared on a hand-built correlated-composite
# case (live).

# ── Varimax identity on default fits ─────────────────────────────────────────

test_that("under varimax beta equals r and r2 equals the column sum of r^2 (three engines)", {
  skip_if_not_installed("lavaan")
  fits <- list(
    pca = cached(ackwards(sim16, k_max = 4)),
    efa = cached(ackwards(sim16, k_max = 4, engine = "efa")),
    esem = cached(ackwards(sim16, k_max = 4, engine = "esem", seed = 1))
  )
  for (nm in names(fits)) {
    x <- fits[[nm]]
    e <- tidy(x, what = "edges")
    expect_identical(
      names(e),
      c("from", "to", "level_from", "level_to", "r", "beta", "is_primary", "above_cut"),
      label = nm
    )
    expect_false(anyNA(e$beta), label = nm)
    expect_equal(e$beta, e$r, tolerance = 1e-10, label = paste(nm, "beta == r"))

    v <- tidy(x, what = "variance")
    expect_identical(names(v), c("level", "factor", "proportion", "cumulative", "r2"), label = nm)
    expect_true(all(is.na(v$r2[v$level == 1L])), label = paste(nm, "r2 NA at k = 1"))
    expect_false(anyNA(v$r2[v$level > 1L]), label = nm)
    # r2 for factor j at level k = sum over level k-1 of r^2 (adjacent pair).
    for (i in which(v$level > 1L)) {
      rows <- e[e$to == v$factor[i] & e$level_from == v$level[i] - 1L, ]
      expect_equal(v$r2[i], sum(rows$r^2), tolerance = 1e-10, label = paste(nm, v$factor[i]))
    }
  }
})

test_that("beta is present for skip-level pairs under pairs = 'all' and equals r under varimax", {
  x <- cached(ackwards(sim16, k_max = 4, pairs = "all"))
  e <- tidy(x, what = "edges")
  skip <- e[e$level_to - e$level_from > 1L, ]
  expect_gt(nrow(skip), 0L)
  expect_false(anyNA(skip$beta))
  expect_equal(skip$beta, skip$r, tolerance = 1e-10)
  # r2 still reads the adjacent level above only.
  v <- tidy(x, what = "variance")
  v_adj <- tidy(cached(ackwards(sim16, k_max = 4)), what = "variance")
  expect_equal(v$r2, v_adj$r2, tolerance = 1e-10)
})

# ── lm() oracle on a hand-built correlated-composite case ────────────────────

test_that("beta and r2 match lm() standardized coefficients and R-squared when composites correlate", {
  x <- cached(ackwards(sim16, k_max = 3))
  y <- x
  # Mix level 2's two weight columns so its score composites are correlated
  # (Phi_s is no longer the identity), then restore the 2:3 edge matrix to the
  # score correlations those mixed weights imply, so the stored E and the
  # stored weights describe the same composites.
  M <- matrix(c(1, 0.6, 0.3, 1), 2L)
  W2 <- y$levels[["2"]]$scoring$weights %*% M
  colnames(W2) <- colnames(y$levels[["2"]]$scoring$weights)
  y$levels[["2"]]$scoring$weights <- W2
  y$levels[["2"]]$scoring$score_var <- ackwards:::.score_var(W2, y$r)
  E <- ackwards:::compute_edges(
    levels = y$levels[c("2", "3")], R = y$r, edge_method = "algebra",
    pairs = "adjacent", build_tidy = FALSE
  )$matrices[["2:3"]]
  y$edges$matrices[["2:3"]] <- E
  key <- paste(y$edges$tidy$from, y$edges$tidy$to, sep = "\r")
  for (i in seq_len(nrow(E))) {
    for (j in seq_len(ncol(E))) {
      y$edges$tidy$r[key == paste(rownames(E)[i], colnames(E)[j], sep = "\r")] <- E[i, j]
    }
  }

  # Oracle: materialize the scores and regress each level-3 score on both
  # level-2 scores together. Pearson R and scale() share one standardization,
  # so the algebra and the regression describe the same composites exactly.
  Z <- scale(as.matrix(sim16))
  S2 <- Z %*% W2
  S3 <- Z %*% y$levels[["3"]]$scoring$weights
  expect_gt(abs(stats::cor(S2)[1L, 2L]), 0.3) # the composites really are correlated

  e <- tidy(y, what = "edges")
  v <- tidy(y, what = "variance")
  for (j in seq_len(ncol(S3))) {
    fit <- stats::lm(S3[, j] ~ S2)
    std_coef <- unname(stats::coef(fit)[-1L] * apply(S2, 2L, stats::sd) / stats::sd(S3[, j]))
    got <- e$beta[e$level_from == 2L & e$to == colnames(S3)[j]]
    expect_equal(got, std_coef, tolerance = 1e-8, label = colnames(S3)[j])
    expect_equal(
      v$r2[v$factor == colnames(S3)[j]], summary(fit)$r.squared,
      tolerance = 1e-8, label = paste("r2", colnames(S3)[j])
    )
    # And beta differs from r here, which is the point of the column.
    r_marg <- e$r[e$level_from == 2L & e$to == colnames(S3)[j]]
    expect_false(isTRUE(all.equal(got, r_marg, tolerance = 1e-3)))
  }
})

# ── Singular Phi_s ────────────────────────────────────────────────────────────

test_that("a non-invertible within-level score correlation gives NA beta/r2 and one warning naming the level", {
  x <- cached(ackwards(sim16, k_max = 3))
  y <- x
  # Two identical weight columns at level 2: Phi_s = [[1, 1], [1, 1]].
  W2 <- y$levels[["2"]]$scoring$weights
  W2[, 2L] <- W2[, 1L]
  y$levels[["2"]]$scoring$weights <- W2

  expect_warning(e <- tidy(y, what = "edges"), "k = 2", class = "rlang_warning")
  expect_true(all(is.na(e$beta[e$level_from == 2L])))
  expect_false(anyNA(e$beta[e$level_from == 1L]))
  expect_false(anyNA(e$r)) # the marginal edge is untouched

  expect_warning(v <- tidy(y, what = "variance"), "k = 2", class = "rlang_warning")
  expect_true(all(is.na(v$r2[v$level == 3L])))
  expect_false(anyNA(v$r2[v$level == 2L]))

  # Under pairs = "all" the singular level starts two stored pairs (2:3 and
  # 2:4) but tidy() warns once for the level, not once per pair.
  z <- cached(ackwards(sim16, k_max = 4, pairs = "all"))
  z$levels[["2"]]$scoring$weights <- W2
  n_tidy_warn <- 0L
  e_all <- withCallingHandlers(
    tidy(z, what = "edges"),
    warning = function(w) {
      n_tidy_warn <<- n_tidy_warn + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n_tidy_warn, 1L)
  expect_true(all(is.na(e_all$beta[e_all$level_from == 2L])))
  expect_false(anyNA(e_all$beta[e_all$level_from != 2L]))

  # The direct helper: exactly one warning, NA of the right shape.
  E <- y$edges$matrices[["2:3"]]
  res <- NULL
  n_warn <- 0L
  withCallingHandlers(
    res <- ackwards:::.partialled_edges(W2, y$r, E, level = 2L),
    warning = function(w) {
      n_warn <<- n_warn + 1L
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n_warn, 1L)
  expect_identical(dim(res$beta), dim(E))
  expect_true(all(is.na(res$beta)))
  expect_identical(names(res$r2), colnames(E))
  expect_true(all(is.na(res$r2)))
  expect_identical(res$status, "singular")
})

# ── Near-singular Phi_s ───────────────────────────────────────────────────────

# Smallest eigenvalue of the within-level score correlation implied by W and
# R, computed here independently of the package's helper.
phi_s_of <- function(W, R) stats::cov2cor(crossprod(W, R %*% W))
min_eig_of <- function(W, R) {
  min(eigen(phi_s_of(W, R), symmetric = TRUE, only.values = TRUE)$values)
}

# Replace level k's score weights and rebuild the score variances, every
# stored edge matrix (which beta and r2 read), and the tidy r column from
# them, so the planted weights and the stored edges describe the same
# composites.
# above_cut, is_primary, and factor_cor are left as fitted.
plant_weights <- function(x, k, W) {
  lev <- as.character(k)
  colnames(W) <- colnames(x$levels[[lev]]$scoring$weights)
  x$levels[[lev]]$scoring$weights <- W
  x$levels[[lev]]$scoring$score_var <- ackwards:::.score_var(W, x$r)
  mats <- ackwards:::compute_edges(
    levels = x$levels, R = x$r, edge_method = "algebra",
    pairs = x$meta$pairs, build_tidy = FALSE
  )$matrices
  x$edges$matrices <- mats
  key <- paste(x$edges$tidy$from, x$edges$tidy$to, sep = "\r")
  for (E in mats) {
    for (i in seq_len(nrow(E))) {
      for (j in seq_len(ncol(E))) {
        x$edges$tidy$r[key == paste(rownames(E)[i], colnames(E)[j], sep = "\r")] <- E[i, j]
      }
    }
  }
  x
}

# Level-2 weights whose second column nearly copies the first.
near_copy_l2 <- function(x, eps) {
  W <- x$levels[["2"]]$scoring$weights
  W[, 2L] <- W[, 1L] + eps * W[, 2L]
  W
}

# Every warning a call raises, muffled, as a list of conditions.
warnings_raised <- function(expr) {
  ws <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  ws
}

test_that("a near-singular level-2 Phi_s keeps beta and raises one warning naming the level, eigenvalue, and beta", {
  x <- cached(ackwards(sim16, k_max = 4))
  W2 <- near_copy_l2(x, 0.13)
  lam <- min_eig_of(W2, x$r)
  # Just under the 1e-2 cutoff, so a stricter eigenvalue cutoff (below 5e-3)
  # would stay silent. An rcond test at 1e-2 (rcond is about half the
  # eigenvalue for a 2 by 2) would also fire, and the control further down
  # catches it.
  expect_gte(lam, 5e-3)
  expect_lt(lam, 1e-2)
  y <- plant_weights(x, 2L, W2)

  ws <- warnings_raised(e <- tidy(y, what = "edges"))
  expect_length(ws, 1L)
  msg <- conditionMessage(ws[[1L]])
  expect_match(msg, "nearly singular", fixed = TRUE)
  expect_match(msg, "k = 2", fixed = TRUE)
  expect_match(msg, format(signif(lam, 2)), fixed = TRUE)
  expect_match(msg, "beta", fixed = TRUE)

  got <- e[e$level_from == 2L, ]
  expect_true(all(is.finite(got$beta)))
  E <- y$edges$matrices[["2:3"]]
  B <- solve(phi_s_of(W2, y$r), E)
  expect_equal(
    got$beta,
    B[cbind(match(got$from, rownames(B)), match(got$to, colnames(B)))],
    tolerance = 1e-8
  )
  # Edges from the other levels still get a finite beta.
  expect_true(all(is.finite(e$beta[e$level_from != 2L])))

  # The warning recurs on the next call: nothing is cached on the object.
  expect_length(warnings_raised(tidy(y, what = "edges")), 1L)
})

test_that("a near-singular level-3 Phi_s (3 by 3) keeps beta and raises one warning naming k = 3", {
  x <- cached(ackwards(sim16, k_max = 4))
  W3 <- x$levels[["3"]]$scoring$weights
  W3[, 3L] <- W3[, 1L] + W3[, 2L] + 0.001 * W3[, 3L]
  lam <- min_eig_of(W3, x$r)
  expect_gt(lam, 1e-8)
  expect_lt(lam, 1e-2)
  y <- plant_weights(x, 3L, W3)

  ws <- warnings_raised(e <- tidy(y, what = "edges"))
  expect_length(ws, 1L)
  msg <- conditionMessage(ws[[1L]])
  expect_match(msg, "nearly singular", fixed = TRUE)
  expect_match(msg, "k = 3", fixed = TRUE)
  expect_match(msg, format(signif(lam, 2)), fixed = TRUE)
  expect_match(msg, "beta", fixed = TRUE)

  got <- e[e$level_from == 3L, ]
  expect_true(all(is.finite(got$beta)))
  E <- y$edges$matrices[["3:4"]]
  B <- solve(phi_s_of(W3, y$r), E)
  expect_equal(
    got$beta,
    B[cbind(match(got$from, rownames(B)), match(got$to, colnames(B)))],
    tolerance = 1e-8
  )
})

test_that("the near-singular warning stays silent just above the cutoff and outside the edge table", {
  x <- cached(ackwards(sim16, k_max = 4))

  # Control just above the cutoff: no warning from the edge table.
  W_ctrl <- near_copy_l2(x, 0.15)
  lam_ctrl <- min_eig_of(W_ctrl, x$r)
  expect_gt(lam_ctrl, 1e-2)
  expect_lte(lam_ctrl, 1.5e-2)
  expect_length(warnings_raised(tidy(plant_weights(x, 2L, W_ctrl), what = "edges")), 0L)

  # On the near-singular plant, r2 and summary() stay quiet: only beta is
  # unstable.
  W2 <- near_copy_l2(x, 0.13)
  expect_lt(min_eig_of(W2, x$r), 1e-2)
  y <- plant_weights(x, 2L, W2)
  expect_length(warnings_raised(v <- tidy(y, what = "variance")), 0L)
  expect_true(all(is.finite(v$r2[v$level == 3L])))
  expect_length(warnings_raised(s <- summary(y)), 0L)
  expect_true(all(is.finite(s$variance$r2[s$variance$level == 3L])))
})

expect_no_tidy_warnings <- function(x) {
  expect_length(warnings_raised(tidy(x, what = "edges")), 0L)
  expect_length(warnings_raised(tidy(x, what = "variance")), 0L)
  expect_length(warnings_raised(summary(x)), 0L)
}

test_that("healthy varimax and PCA promax fits raise no near-singular warning", {
  expect_no_tidy_warnings(cached(ackwards(sim16, k_max = 4)))
  expect_no_tidy_warnings(cached(ackwards(bfi25, k_max = 8, rotation = "promax")))
})

test_that("a healthy EFA promax fit raises no near-singular warning", {
  skip_if_not_installed("GPArotation")
  expect_no_tidy_warnings(
    cached(ackwards(bfi25, k_max = 8, engine = "efa", rotation = "promax"))
  )
})

test_that("the singular path takes precedence over the near-singular one", {
  x <- cached(ackwards(sim16, k_max = 4))
  W2 <- x$levels[["2"]]$scoring$weights
  W2[, 2L] <- W2[, 1L] # Phi_s = [[1, 1], [1, 1]]: solve() rejects it
  y <- x
  y$levels[["2"]]$scoring$weights <- W2

  ws <- warnings_raised(e <- tidy(y, what = "edges"))
  expect_length(ws, 1L)
  expect_match(conditionMessage(ws[[1L]]), "cannot be inverted", fixed = TRUE)
  expect_no_match(conditionMessage(ws[[1L]]), "nearly singular", fixed = TRUE)
  expect_true(all(is.na(e$beta[e$level_from == 2L])))

  ws <- warnings_raised(v <- tidy(y, what = "variance"))
  expect_length(ws, 1L)
  expect_match(conditionMessage(ws[[1L]]), "cannot be inverted", fixed = TRUE)
  expect_true(all(is.na(v$r2[v$level == 3L])))
})

test_that("under pairs = 'all' a near-singular level warns once per call, not once per pair", {
  x <- cached(ackwards(sim16, k_max = 4, pairs = "all"))
  W2 <- near_copy_l2(x, 0.13)
  expect_lt(min_eig_of(W2, x$r), 1e-2)
  y <- plant_weights(x, 2L, W2)
  # Level 2 starts two stored pairs here.
  expect_true(all(c("2:3", "2:4") %in% names(y$edges$matrices)))

  for (call in 1:2) {
    ws <- warnings_raised(e <- tidy(y, what = "edges"))
    expect_length(ws, 1L)
    expect_match(conditionMessage(ws[[1L]]), "nearly singular", fixed = TRUE)
    expect_match(conditionMessage(ws[[1L]]), "k = 2", fixed = TRUE)
  }
  expect_true(all(is.finite(e$beta[e$level_from == 2L])))
})
