# Guards the CI paths-ignore filters (M82): the three check workflows skip the
# matrix on tracking-only commits, and that filter has two failure modes worth
# locking. It rots -- four entries pointed at files that had moved under cairn/
# months earlier -- and it can be widened too far: a blanket 'cairn/**' would
# let a commit that breaks a departures-ledger anchor skip the very check that
# catches it (tools/check-ledger-anchors.R reads cairn/DECISIONS.md,
# cairn/DESIGN.md and cairn/references/** from the source checkout).
#
# The golden test below (real repo state, expecting no problems) would pass
# against a checker that never reports anything, so the fixture tests carry the
# real weight: each mutates one thing and asserts the specific complaint. The
# comment-hidden fixture is a regression test for a fail-open bug the M82
# review caught -- the first parser ended a block at the first non-item line,
# so a 'cairn/**' written below a comment was never examined.
#
# .github/ and tools/ are .Rbuildignore'd, so the built package contains neither
# the workflows nor the checker -- these tests SKIP there and run only in the
# source checkout (devtools::test, CI source step). CI also runs the checker
# standalone (see tools/check-ci-path-filters.R), so the skip loses no coverage.

checker_env <- function() {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  checker <- file.path(root, "tools", "check-ci-path-filters.R")
  skip_if_not(file.exists(checker), "checker absent (built package)")
  env <- new.env()
  # Reuse the checker's functions without triggering its script body (guarded
  # by sys.nframe() == 0L, false under sys.source).
  sys.source(checker, envir = env)
  env
}

fixture <- function(case) test_path("fixtures", "ci-filters", case)

test_that("CI paths-ignore filters are current and preserve the ledger check", {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  skip_if_not(
    dir.exists(file.path(root, ".github", "workflows")),
    "workflows absent (built package or non-source checkout)"
  )
  expect_equal(checker_env()$check_ci_path_filters(root), character(0L))
})

test_that("a widened pattern is caught even when a comment precedes it", {
  # Regression: the blanket pattern sits below a comment line, which the
  # original parser treated as the end of the list.
  problems <- checker_env()$check_ci_path_filters(fixture("comment-hidden"))
  expect_true(any(grepl("match ledger-read path cairn/DECISIONS.md", problems)))
  expect_true(any(grepl("match ledger-read path cairn/references/", problems)))
})

test_that("a literal path that no longer exists is caught", {
  problems <- checker_env()$check_ci_path_filters(fixture("dead-literal"))
  expect_true(any(grepl("literal path\\(s\\) absent from repo: MILESTONES.md", problems)))
})

test_that("a dropped cairn tracking entry is caught", {
  problems <- checker_env()$check_ci_path_filters(fixture("dropped-entry"))
  expect_true(any(grepl("missing required cairn entries:.*cairn/milestones/\\*\\*", problems)))
})

test_that("a deleted block is caught", {
  # The fixture defines only R-CMD-check.yaml :: push; the other required
  # blocks are absent entirely, which must be reported rather than ignored.
  problems <- checker_env()$check_ci_path_filters(fixture("dead-literal"))
  expect_true(any(grepl("block\\(s\\) missing entirely", problems)))
})

test_that("shapes the parser cannot read are reported, never skipped", {
  env <- checker_env()

  flow <- env$check_ci_path_filters(fixture("flow-style"))
  expect_true(any(grepl("flow-style paths-ignore is not supported", flow)))

  dup <- env$check_ci_path_filters(fixture("dup-block"))
  expect_true(any(grepl("two paths-ignore blocks resolve to the same name", dup)))

  bad <- env$check_ci_path_filters(fixture("unparseable"))
  expect_true(any(grepl("unparseable line inside the block", bad)))
})

test_that("the parser recovers every entry, skipping comments and blanks", {
  parsed <- checker_env()$ci_parse_blocks(
    file.path(fixture("comment-hidden"), ".github", "workflows", "R-CMD-check.yaml")
  )
  pats <- parsed$blocks[["R-CMD-check.yaml :: push"]]
  # Nine real entries plus the blanket one written below the comment.
  expect_length(pats, 10L)
  expect_true("cairn/**" %in% pats)
  expect_true("cairn/reviews/**" %in% pats)
})

test_that("the glob matcher follows GitHub's paths-ignore semantics", {
  m <- checker_env()$ci_glob_matches

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
