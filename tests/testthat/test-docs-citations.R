# Citation hygiene (M30): the package citation (inst/CITATION) should name
# only Girard as the software author -- Goldberg is the method source, cited
# in DESCRIPTION and in ackwards()'s roxygen @references, not as the software
# citation. ackwards() implements the original Goldberg (2006) method, the
# Waller (2007) closed-form algebra, and the Forbes (2023) extension, so its
# \references block (man/ackwards.Rd, source-tree copy) must carry all three.
#
# This is a source-tree hygiene check in the same style as
# test-docs-no-milestone-refs.R: it inspects man/ackwards.Rd directly rather
# than round-tripping through tools::Rd2txt, so it runs wherever the source
# tree is reachable and skips gracefully otherwise (e.g. an installed-package
# R CMD check).

.pkg_root_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) candidate else NULL
}

test_that("ackwards() references block cites Goldberg, Waller, and Forbes", {
  path <- .pkg_root_file("man", "ackwards.Rd")
  skip_if(is.null(path), "man/ackwards.Rd not reachable from this test context")

  rd <- paste(readLines(path), collapse = "\n")
  references_start <- regexpr("\\\\references\\{", rd)
  skip_if(references_start < 0, "no \\references block found")
  references <- substring(rd, references_start)

  expect_match(references, "10.1016/j.jrp.2006.01.001", fixed = TRUE) # Goldberg 2006
  expect_match(references, "10.1016/j.jrp.2006.08.005", fixed = TRUE) # Waller 2007
  expect_match(references, "10.1037/met0000546", fixed = TRUE) # Forbes 2023
})
