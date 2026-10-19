# The engines vignette's comparison table must not describe the ESEM WLSMV and
# ULSMV estimators as ordinal-only. `ackwards()` fits them on a continuous
# `cor` too (see `?ackwards`, `estimator`), and its help says so.
#
# Like test-docs-no-milestone-refs.R, this inspects source-tree files and
# skips when they are not reachable (an installed-package R CMD check).

.pkg_root_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) candidate else NULL
}

test_that("engines vignette does not call WLSMV/ULSMV ordinal-only", {
  paths <- Filter(Negate(is.null), list(
    .pkg_root_file("vignettes", "ackwards-engines.Rmd"),
    .pkg_root_file("vignettes", "ackwards-engines.Rmd.orig")
  ))
  skip_if(length(paths) == 0, "engines vignette not reachable from this test context")

  for (p in paths) {
    row <- grep("^\\| \\*\\*Estimators\\*\\* \\|", readLines(p), value = TRUE)
    expect_length(row, 1)
    cells <- trimws(strsplit(row, "|", fixed = TRUE)[[1]])
    esem_cell <- cells[length(cells)]
    wls_part <- sub("^.*?(WLSMV)", "\\1", esem_cell, perl = TRUE)
    expect_match(wls_part, "^WLSMV", info = basename(p))
    expect_match(wls_part, "continuous", info = basename(p))
  }
})
