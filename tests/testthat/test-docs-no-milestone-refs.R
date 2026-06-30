# User-facing docs (NEWS, README, vignettes) must not leak internal milestone
# numbers like "(M24)" -- they are meaningless outside this repo's dev process.
#
# This is a source-tree hygiene check: it inspects repo files that are not part
# of the installed package, so it runs under `devtools::test()` and the
# test-coverage CI job (both execute in the source tree) and skips gracefully
# when those files aren't reachable (e.g. an installed-package R CMD check). We
# deliberately do NOT ship copies into inst/ to force it under R CMD check --
# test-coverage CI already enforces the guard, so a reintroduced tag is caught.
# README.Rmd (the source) is scanned alongside README.md (the rendered/shipped
# file) so a leak is caught both at the source and in what users actually see.

.pkg_root_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) candidate else NULL
}

test_that("NEWS.md has no milestone-number tags", {
  path <- .pkg_root_file("NEWS.md")
  skip_if(is.null(path), "NEWS.md not reachable from this test context")

  hits <- grep("\\(M[0-9]+\\)", readLines(path), value = TRUE)
  expect_length(hits, 0)
})

test_that("README files have no milestone-number tags", {
  paths <- Filter(Negate(is.null), list(
    .pkg_root_file("README.md"),
    .pkg_root_file("README.Rmd")
  ))
  skip_if(length(paths) == 0, "README files not reachable from this test context")

  hits <- unlist(lapply(paths, function(p) {
    lines <- grep("\\(M[0-9]+\\)", readLines(p), value = TRUE)
    if (length(lines) > 0) paste0(basename(p), ": ", lines) else character(0)
  }))
  expect_length(hits, 0)
})

test_that("vignettes have no milestone-number tags", {
  dir <- .pkg_root_file("vignettes")
  skip_if(is.null(dir), "vignettes/ not reachable from this test context")

  rmd_files <- list.files(dir, pattern = "\\.Rmd$", full.names = TRUE)
  skip_if(length(rmd_files) == 0, "no .Rmd files found")

  hits <- unlist(lapply(rmd_files, function(f) {
    lines <- grep("\\(M[0-9]+\\)", readLines(f), value = TRUE)
    if (length(lines) > 0) paste0(basename(f), ": ", lines) else character(0)
  }))
  expect_length(hits, 0)
})
