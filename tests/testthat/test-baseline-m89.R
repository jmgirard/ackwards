# Default-output baseline (frozen oracle, cairn/ORACLES.md).
#
# The fixture holds the loadings, weights, variance, factor_cor, and every
# edge matrix of five default fits, generated on master before the
# factor-correlation carry and the partialled edge columns were added
# (data-raw/baseline-m89.R names the commit). Each fit is refit here fresh
# (never through cached(): a memoized second call asserts nothing) and every
# stored value must reproduce within 1e-12. The ESEM fit is seeded, because
# lavaan's rotation draws random starts (see the generator's header).

baseline <- readRDS(test_path("fixtures", "baseline-m89.rds"))

.refit_baseline <- function(spec) {
  dat <- get(spec$data, envir = asNamespace("ackwards"))
  suppressWarnings(suppressMessages(
    ackwards(dat, k_max = spec$k_max, engine = spec$engine, cor = spec$cor, seed = spec$seed)
  ))
}

for (name in names(baseline)) {
  test_that(paste0("default output is unchanged against the baseline: ", name), {
    entry <- baseline[[name]]
    if (entry$spec$engine == "esem") skip_if_not_installed("lavaan")
    x <- .refit_baseline(entry$spec)

    expect_identical(names(x$levels), names(entry$levels))
    for (ki in names(entry$levels)) {
      lev <- x$levels[[ki]]
      ref <- entry$levels[[ki]]
      expect_equal(lev$loadings, ref$loadings, tolerance = 1e-12, label = paste(name, ki, "loadings"))
      expect_equal(lev$scoring$weights, ref$weights, tolerance = 1e-12, label = paste(name, ki, "weights"))
      expect_equal(lev$variance, ref$variance, tolerance = 1e-12, label = paste(name, ki, "variance"))
      expect_equal(lev$factor_cor, ref$factor_cor, tolerance = 1e-12, label = paste(name, ki, "factor_cor"))
    }

    expect_identical(names(x$edges$matrices), names(entry$edges))
    for (key in names(entry$edges)) {
      expect_equal(x$edges$matrices[[key]], entry$edges[[key]], tolerance = 1e-12, label = paste(name, "edge", key))
    }
  })
}
