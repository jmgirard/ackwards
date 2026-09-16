# Generate the default-output baseline fixture
#   tests/testthat/fixtures/baseline-m89.rds
#
# Purpose: prove that the factor-correlation carry and the partialled edge
# columns leave default-argument output unchanged. The script fits five
# default hierarchies and stores, for every level, the loadings, the score
# weights, the `variance` vector, and `factor_cor`, plus every stored edge
# matrix. The fixture was generated on master before any code change of the
# milestone that adds those features, so the test that reads it
# (tests/testthat/test-baseline-m89.R) compares the current code against the
# pre-change numbers. The comparison tolerance is 1e-12.
#
# The five fits (all default arguments except the ones named):
#   pca  on sim16, k_max = 4
#   efa  on sim16, k_max = 4
#   esem on sim16, k_max = 4, seed = 1
#   pca  on bfi25, k_max = 3, cor = "polychoric"
#   efa  on bfi25, k_max = 3, cor = "polychoric"
#
# The ESEM fit carries `seed = 1` because lavaan's rotation uses random
# starts: an unseeded ESEM fit differs run to run at about 1e-6 (observed
# 2026-09-16), so it cannot serve a 1e-12 oracle. A seeded fit reproduces to
# 0 across runs. The seed fixes the random starts only. It changes no method
# argument, so the fit is still the default ESEM hierarchy.
#
# Provenance: the `provenance` attribute names this generator, the data
# source (the package's own bundled datasets), the git commit the script ran
# on, and the date. Re-run only to re-baseline on purpose:
#   Rscript data-raw/baseline-m89.R
# `lavaan` must be installed for the ESEM fit.

pkgload::load_all(".", quiet = TRUE)

commit <- tryCatch(
  trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  error = function(e) NA_character_
)

specs <- list(
  pca_sim16 = list(data = "sim16", engine = "pca", k_max = 4L, cor = "pearson"),
  efa_sim16 = list(data = "sim16", engine = "efa", k_max = 4L, cor = "pearson"),
  esem_sim16 = list(data = "sim16", engine = "esem", k_max = 4L, cor = "pearson", seed = 1L),
  pca_bfi25_poly = list(data = "bfi25", engine = "pca", k_max = 3L, cor = "polychoric"),
  efa_bfi25_poly = list(data = "bfi25", engine = "efa", k_max = 3L, cor = "polychoric")
)

snapshot_fit <- function(spec) {
  dat <- get(spec$data)
  x <- suppressWarnings(suppressMessages(
    ackwards(dat, k_max = spec$k_max, engine = spec$engine, cor = spec$cor, seed = spec$seed)
  ))
  levels <- lapply(x$levels, function(lev) {
    list(
      loadings = lev$loadings,
      weights = lev$scoring$weights,
      variance = lev$variance,
      factor_cor = lev$factor_cor
    )
  })
  list(spec = spec, levels = levels, edges = x$edges$matrices)
}

baseline <- lapply(specs, snapshot_fit)

attr(baseline, "provenance") <- list(
  generator = "data-raw/baseline-m89.R",
  source = "ackwards() default fits on the bundled sim16 and bfi25 datasets",
  commit = commit,
  generated = as.character(Sys.Date()),
  pkg_version = as.character(utils::packageVersion("ackwards")),
  tolerance = 1e-12
)

saveRDS(baseline, "tests/testthat/fixtures/baseline-m89.rds", version = 2)
cat("Wrote tests/testthat/fixtures/baseline-m89.rds at commit", commit, "\n")
