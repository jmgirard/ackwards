# Citation hygiene (M30): the package citation (inst/CITATION) should name
# only Girard as the software author -- Goldberg is the method source, cited
# in DESCRIPTION and in ackwards()'s roxygen @references, not as the software
# citation. ackwards() implements the original Goldberg (2006) method, the
# Waller (2007) closed-form algebra, and the Forbes (2023) extension, so its
# \references block (man/ackwards.Rd, source-tree copy) must carry all three.
#
# These are source-tree hygiene checks in the same style as
# test-docs-no-milestone-refs.R: they inspect man/ackwards.Rd and inst/CITATION
# directly rather than round-tripping through tools::Rd2txt / citation(), so
# they run wherever the source tree is reachable and skip gracefully otherwise
# (e.g. an installed-package R CMD check).

.pkg_root_file <- function(...) {
  candidate <- testthat::test_path("..", "..", ...)
  if (file.exists(candidate)) candidate else NULL
}

test_that("ackwards() references block cites Goldberg, Waller, and Forbes", {
  path <- .pkg_root_file("man", "ackwards.Rd")
  skip_if(is.null(path), "man/ackwards.Rd not reachable from this test context")

  # Extract only the \references{...} block, bounded by its closing brace, so
  # a DOI that ever moved out of the block (e.g. into @details prose) cannot
  # satisfy these assertions. roxygen writes the closing brace on its own line
  # and \doi{...} always closes on the same line, so the first lone "}" after
  # \references{ is the block's terminator.
  lines <- readLines(path)
  start <- grep("^\\\\references\\{", lines)
  skip_if(length(start) == 0, "no \\references block found")
  after <- lines[(start + 1):length(lines)]
  close_rel <- match(TRUE, trimws(after) == "}")
  skip_if(is.na(close_rel), "\\references block not closed as expected")
  references <- paste(after[seq_len(close_rel - 1)], collapse = "\n")

  expect_match(references, "10.1016/j.jrp.2006.01.001", fixed = TRUE) # Goldberg 2006
  expect_match(references, "10.1016/j.jrp.2006.08.005", fixed = TRUE) # Waller 2007
  expect_match(references, "10.1037/met0000546", fixed = TRUE) # Forbes 2023
})

test_that("inst/CITATION source text names Girard, not Goldberg", {
  path <- .pkg_root_file("inst", "CITATION")
  skip_if(is.null(path), "inst/CITATION not reachable from this test context")

  citation_src <- paste(readLines(path), collapse = "\n")
  # Complements the runtime citation() check in test-utils.R: guards the source
  # file directly so a reintroduced Goldberg bibentry is caught even where the
  # package cannot be loaded. The method paper belongs in DESCRIPTION / roxygen
  # @references, not the software citation.
  expect_match(citation_src, "Girard", fixed = TRUE)
  expect_false(grepl("Goldberg", citation_src, fixed = TRUE))
})
