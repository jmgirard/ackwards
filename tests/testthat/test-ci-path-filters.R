# Guards the CI paths-ignore filters (M82): the three check workflows skip the
# matrix on tracking-only commits, and that filter has two failure modes worth
# locking. It rots -- four entries pointed at files that had moved under cairn/
# months earlier -- and it can be widened too far: a blanket 'cairn/**' would
# let a commit that breaks a departures-ledger anchor skip the very check that
# catches it (tools/check-ledger-anchors.R reads cairn/DECISIONS.md,
# cairn/DESIGN.md and cairn/references/** from the source checkout).
#
# .github/ and tools/ are .Rbuildignore'd, so the built package contains neither
# the workflows nor the checker -- this test SKIPS there and runs only in the
# source checkout (devtools::test, CI source step). CI also runs the checker
# standalone (see tools/check-ci-path-filters.R), so the skip loses no coverage.

test_that("CI paths-ignore filters are current and preserve the ledger check", {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  checker <- file.path(root, "tools", "check-ci-path-filters.R")
  workflows <- file.path(root, ".github", "workflows")

  skip_if_not(
    file.exists(checker) && dir.exists(workflows),
    "workflows or checker absent (built package or non-source checkout)"
  )

  # Reuse the checker's function without triggering its script body (guarded by
  # sys.nframe() == 0L, false under sys.source).
  env <- new.env()
  sys.source(checker, envir = env)
  problems <- env$check_ci_path_filters(root)

  expect_equal(problems, character(0L))
})

test_that("the glob matcher follows GitHub's paths-ignore semantics", {
  checker <- file.path(
    normalizePath(test_path("..", ".."), mustWork = FALSE),
    "tools", "check-ci-path-filters.R"
  )
  skip_if_not(file.exists(checker), "checker absent (built package)")
  env <- new.env()
  sys.source(checker, envir = env)
  m <- env$ci_glob_matches

  # `**` spans directory separators; `*` does not.
  expect_true(m("cairn/**", "cairn/DECISIONS.md"))
  expect_true(m("cairn/**", "cairn/references/goldberg2006.md"))
  expect_false(m("cairn/*", "cairn/references/goldberg2006.md"))

  # Patterns match the whole path, not a prefix -- the check's clause 3 would
  # be vacuous if a narrow pattern matched a ledger-read path by prefix.
  expect_true(m("cairn/milestones/**", "cairn/milestones/archive/M77.md"))
  expect_false(m("cairn/milestones/**", "cairn/DECISIONS.md"))
  expect_false(m("cairn/ROADMAP.md", "cairn/ROADMAP.md.bak"))

  # A literal dot is a dot, not a regex wildcard.
  expect_false(m("cairn/ROADMAP.md", "cairn/ROADMAPXmd"))
})
