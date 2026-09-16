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
