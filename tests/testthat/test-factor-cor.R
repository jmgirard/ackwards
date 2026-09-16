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
