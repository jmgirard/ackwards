# Oracle-provenance guard (M57).
#
# Mechanizes CLAUDE.md Invariant 8 ("no unsourced or unreproducible reference
# value ships") at the fixture layer: every committed test fixture must carry a
# structured top-level `provenance` attr that names its `data-raw/` generator and
# its source. This is the automated "bar" a future un-sourced fixture cannot slip
# past. The registry that classifies every oracle is cairn/ORACLES.md.
#
# NB: data-raw/ is .Rbuildignore'd, so the generator scripts are absent at
# R CMD check time. The guard therefore checks that provenance NAMES a data-raw
# generator (the reproducible recipe), not that the file exists on disk.

.provenance_ok <- function(path) {
  obj <- readRDS(path)
  p <- attr(obj, "provenance")
  is.list(p) &&
    is.character(p$generator) && length(p$generator) == 1L && nzchar(p$generator) &&
    grepl("^data-raw/.*\\.R$", p$generator) &&
    is.character(p$source) && length(p$source) == 1L && nzchar(p$source)
}

test_that("every fixture carries a structured provenance attr (generator + source)", {
  fixtures <- list.files(test_path("fixtures"), pattern = "\\.rds$", full.names = TRUE)
  expect_gt(length(fixtures), 0L)
  for (f in fixtures) {
    expect_true(.provenance_ok(f), label = paste0("provenance for ", basename(f)))
  }
})

test_that("the provenance guard has teeth (rejects unsourced / partial fixtures)", {
  scratch <- tempfile(fileext = ".rds")
  on.exit(unlink(scratch), add = TRUE)

  # No provenance attr at all -> rejected.
  saveRDS(list(x = 1), scratch)
  expect_false(.provenance_ok(scratch))

  # Generator present but source missing -> still rejected.
  saveRDS(structure(list(x = 1), provenance = list(generator = "data-raw/x.R")), scratch)
  expect_false(.provenance_ok(scratch))

  # A generator that is not a data-raw/*.R recipe -> rejected.
  saveRDS(
    structure(list(x = 1), provenance = list(generator = "made it up", source = "trust me")),
    scratch
  )
  expect_false(.provenance_ok(scratch))
})
