# Within-level factor correlations: engine routing, sign-alignment carry,
# tidy(what = "factor_cor"), and the summary() block.

# ── Routing: every engine builds factor_cor through .carry_factor_cor() ──────

test_that("every engine routes its factor_cor through .carry_factor_cor()", {
  skip_if_not_installed("lavaan")
  # Replace the helper with one that stamps its output. A level's factor_cor
  # carries the stamp only if the engine built it through the helper.
  # align_signs = FALSE keeps ackwards()'s own post-alignment carry out of the
  # picture, so the stamp can come from the engines alone.
  testthat::local_mocked_bindings(
    .carry_factor_cor = function(Phi, ord, signs) {
      out <- as.matrix(Phi)[ord, ord, drop = FALSE] * tcrossprod(as.numeric(signs))
      attr(out, "stamp") <- "carried"
      out
    }
  )
  for (engine in c("pca", "efa", "esem")) {
    x <- suppressWarnings(suppressMessages(
      ackwards(sim16, k_max = 3, engine = engine, align_signs = FALSE, seed = 1)
    ))
    for (ki in names(x$levels)) {
      expect_identical(
        attr(x$levels[[ki]]$factor_cor, "stamp"), "carried",
        label = paste(engine, "level", ki)
      )
    }
  }
})

# ── ackwards() flips factor_cor in step with align_signs ─────────────────────

test_that("ackwards() flips factor_cor rows and columns by the align_signs vectors", {
  # A default fit's factor_cor is the identity and its sign flips may all be
  # +1, so neither side of the carry is observable. Plant both: the engine
  # reports a non-identity Phi at k = 3, and the level-3 sign vector carries
  # exactly one -1 (the real alignment's vector with factor 2 negated, its
  # loadings negated to match, so lineage and edges stay consistent).
  Phi3 <- matrix(c(1, .35, -.2, .35, 1, .1, -.2, .1, 1), 3L)
  real_engine_phi <- ackwards:::.engine_phi
  real_align <- ackwards:::.align_signs
  testthat::local_mocked_bindings(
    .engine_phi = function(fit, k) if (k == 3L) Phi3 else real_engine_phi(fit, k),
    .align_signs = function(loadings_list, edges_list, lineage) {
      out <- real_align(loadings_list, edges_list, lineage)
      s <- rep(1, 3L)
      s[2L] <- -1
      out$loadings[[3L]] <- sweep(out$loadings[[3L]], 2L, s * out$signs[[3L]], "*")
      out$signs[[3L]] <- s
      out
    }
  )
  x <- suppressWarnings(ackwards(sim16, k_max = 3, engine = "pca"))
  signs3 <- c(1, -1, 1)
  expect_equal(unname(x$levels[["3"]]$factor_cor), Phi3 * tcrossprod(signs3))
  # The flipped factor's loadings really were negated relative to the
  # unaligned fit, so the sign vector used is the one the carry saw.
  x0 <- suppressWarnings(ackwards(sim16, k_max = 3, engine = "pca", align_signs = FALSE))
  ratio <- colSums(x$levels[["3"]]$loadings * x0$levels[["3"]]$loadings)
  expect_equal(unname(sign(ratio)), signs3)
  # Levels whose Phi is the identity are unchanged by the flip.
  expect_equal(unname(x$levels[["2"]]$factor_cor), diag(2L))
})

# ── tidy(what = "factor_cor") ────────────────────────────────────────────────

test_that("tidy(what = 'factor_cor') lists every unordered pair per level, all 0 under varimax", {
  x <- cached(ackwards(sim16, k_max = 4))
  fc <- tidy(x, what = "factor_cor")
  expect_identical(names(fc), c("level", "factor_a", "factor_b", "cor"))
  # choose(k, 2) rows per level, none for k = 1.
  expect_identical(as.vector(table(fc$level)), c(1L, 3L, 6L))
  expect_identical(fc$level[1L], 2L)
  # Upper-triangle order within a level: (1,2), (1,3), (2,3).
  l3 <- fc[fc$level == 3L, ]
  expect_identical(l3$factor_a, c("m3f1", "m3f1", "m3f2"))
  expect_identical(l3$factor_b, c("m3f2", "m3f3", "m3f3"))
  # Varimax: every within-level correlation is 0.
  expect_true(all(abs(fc$cor) < 1e-12))
})

test_that("tidy(what = 'factor_cor') returns zero rows with the same columns when only k = 1 exists", {
  x <- cached(ackwards(sim16, k_max = 3))
  x1 <- x
  x1$levels <- x1$levels["1"]
  fc <- tidy(x1, what = "factor_cor")
  expect_identical(names(fc), c("level", "factor_a", "factor_b", "cor"))
  expect_identical(nrow(fc), 0L)
  expect_type(fc$cor, "double")
})

test_that("tidy(what = 'factor_cor') returns a planted non-identity level-3 correlation", {
  x <- cached(ackwards(sim16, k_max = 3))
  y <- x
  Phi3 <- matrix(c(1, .4, -.25, .4, 1, .1, -.25, .1, 1), 3L)
  y$levels[["3"]]$factor_cor <- Phi3
  fc <- tidy(y, what = "factor_cor")
  l3 <- fc[fc$level == 3L, ]
  expect_equal(l3$cor, c(.4, -.25, .1))
  expect_identical(l3$factor_a, c("m3f1", "m3f1", "m3f2"))
  expect_identical(l3$factor_b, c("m3f2", "m3f3", "m3f3"))
  # Other levels are untouched.
  expect_true(all(abs(fc$cor[fc$level == 2L]) < 1e-12))
  # Label columns appear only when factor labels are set, keyed on the pair.
  y <- set_factor_labels(y, c(m3f1 = "Alpha"))
  fcl <- tidy(y, what = "factor_cor")
  expect_identical(names(fcl), c("level", "factor_a", "factor_b", "cor", "factor_a_label", "factor_b_label"))
  expect_identical(fcl$factor_a_label[fcl$level == 3L], c("Alpha", "Alpha", NA))
})
