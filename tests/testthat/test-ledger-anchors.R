# Guards the departures-ledger anchor invariant (M72): every D-entry, DESIGN
# principle, sibling reference note, and R source file cited in
# cairn/references/source-departures.md must still resolve.
#
# cairn/ and tools/ are .Rbuildignore'd, so the built package contains neither
# the ledger nor the checker -- this test SKIPS there and runs only in the
# source checkout (devtools::test, CI source step). CI also runs the checker
# standalone (see tools/check-ledger-anchors.R), so the skip loses no coverage.

test_that("source-departures.md anchors all resolve", {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  checker <- file.path(root, "tools", "check-ledger-anchors.R")
  ledger <- file.path(root, "cairn", "references", "source-departures.md")

  skip_if_not(
    file.exists(checker) && file.exists(ledger),
    "ledger or checker absent (built package or non-source checkout)"
  )

  # Reuse the checker's function without triggering its script body (guarded by
  # sys.nframe() == 0L, false under sys.source).
  env <- new.env()
  sys.source(checker, envir = env)
  problems <- env$check_ledger_anchors(root)

  expect_equal(problems, character(0L))
})
