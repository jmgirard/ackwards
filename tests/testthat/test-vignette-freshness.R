# Guards the precomputed-vignette staleness invariant (M65): every generated
# vignettes/*.Rmd must carry a freshness stamp matching the md5 of the
# vignettes/*.Rmd.orig it was knitted from (via vignettes/precompute.R).
#
# The .Rmd.orig sources, precompute.R, and tools/check-vignette-freshness.R are
# all .Rbuildignore'd, so the built package that R CMD check inspects contains
# none of them -- this test therefore SKIPS there and runs only in the source
# checkout (devtools::test and the CI source-checkout step). The same guard is
# also run standalone in CI before the tarball is built, so the skip loses no
# coverage.

test_that("precomputed vignettes are fresh with their .Rmd.orig sources", {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  checker <- file.path(root, "tools", "check-vignette-freshness.R")
  vign_dir <- file.path(root, "vignettes")
  orig <- list.files(vign_dir, pattern = "\\.Rmd\\.orig$")

  skip_if_not(
    file.exists(checker) && length(orig) > 0L,
    "vignette sources absent (built package or non-source checkout)"
  )

  # Reuse the checker's function without triggering its script body (guarded by
  # sys.nframe() == 0L, false under sys.source).
  env <- new.env()
  sys.source(checker, envir = env)
  problems <- env$check_vignette_freshness(vign_dir)

  expect_equal(problems, character(0L))
})
