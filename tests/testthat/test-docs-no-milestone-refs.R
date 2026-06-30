# User-facing docs (NEWS, README, vignettes) must not leak internal milestone
# numbers like "(M24)" -- they are meaningless outside this repo's dev process.
# Only runs against a source checkout (skips under an installed-package check
# where these files aren't present alongside tests/testthat).

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

test_that("README.md has no milestone-number tags", {
  path <- .pkg_root_file("README.md")
  skip_if(is.null(path), "README.md not reachable from this test context")

  hits <- grep("\\(M[0-9]+\\)", readLines(path), value = TRUE)
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
