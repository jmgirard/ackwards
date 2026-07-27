# Fork-safety guard for the test suite itself.
#
# R Core is explicit that forking and multi-threaded libraries do not mix.
# `?mclapply`: "It is _strongly discouraged_ to use these functions with
# multi-threaded libraries or packages (see 'mcfork' for more details). If in
# doubt, it is safer to use a non-FORK cluster". A forked child inherits only the
# calling thread, so a lock another thread held stays locked in the child
# forever -- and R's linear algebra runs through a BLAS/LAPACK that is threaded on
# many platforms.
#
# This package hit exactly that on CRAN at 0.1.1: `r-oldrel-macos-arm64` died
# with `*** caught segfault *** address 0x110, cause 'invalid permissions'`
# inside lavaan's `lav_inspect_vcov()`, in an `mcparallel()` child spawned by
# `future::plan(future::multicore)` -- while every other flavour passed the same
# test, `r-release-macos-arm64` and `r-oldrel-macos-x86_64` included.
# `future::supportsMulticore()` reports TRUE there: it answers "can this platform
# fork", never "is forking safe for this workload".
#
# The crash is platform-specific and intermittent, so no runnable test reproduces
# it here. This guard locks the rule instead of the symptom: any test block that
# establishes a forking plan must also skip macOS, the platform where it bit. The
# sweep reads the suite's own sources, which are installed alongside it and so
# are readable both under `devtools::test()` and under `R CMD check`.

# A forking plan, in either spelling the suite can write it: the namespaced
# `future::plan(future::multicore, ...)` and a bare `plan(multicore, ...)` after
# `library(future)`. Must NOT fire on `supportsMulticore()` (the capability probe
# sits beside every real fork site) nor on the non-forking `multisession`.
.forks_multicore <- function(lines) {
  any(grepl("future::multicore|plan\\([^)]*\\bmulticore\\b", lines))
}

# A macOS skip, in either spelling: skip_on_os("mac") and the combined
# skip_on_os(c("mac", "windows")) are equally valid guards.
.skips_mac <- function(lines) {
  any(grepl("skip_on_os\\([^)]*mac", lines))
}

# Split a test file into per-`test_that()` blocks. File-level matching would be
# false coverage here: test-esem.R already carries an unrelated
# `skip_on_os("windows")`, and a mac skip in some *other* block of the same file
# would satisfy a file-level check while the forking block stayed unguarded.
.fork_test_blocks <- function(path) {
  src <- readLines(path, warn = FALSE)
  starts <- grep("^test_that\\(", src)
  if (length(starts) == 0L) {
    return(list())
  }
  ends <- c(starts[-1L] - 1L, length(src))
  Map(function(a, b) src[a:b], starts, ends)
}

test_that("the fork/skip matchers see every spelling they must classify", {
  # Positive controls carried into the test: a sweep is only as strong as its
  # matcher, and the author of a matcher is exactly who cannot enumerate the
  # renderings it misses by inspection alone.
  expect_true(.forks_multicore("  future::plan(future::multicore, workers = 2)"))
  expect_true(.forks_multicore("  plan(multicore, workers = 2)"))
  expect_true(.forks_multicore("  plan(multicore)"))
  # Near-misses that must stay unflagged.
  expect_false(.forks_multicore('  skip_if(!future::supportsMulticore(), "no")'))
  expect_false(.forks_multicore("  future::plan(future::multisession, workers = 2)"))
  expect_false(.forks_multicore("  future::plan(future::sequential)"))

  expect_true(.skips_mac('  skip_on_os("mac")'))
  expect_true(.skips_mac('  skip_on_os(c("mac", "windows"))'))
  expect_false(.skips_mac('  skip_on_os("windows")'))
})

test_that("every test block that forks a future plan skips macOS", {
  files <- list.files(
    test_path("."),
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
  # This file spells out forking plans in its own matcher controls; exclude it so
  # the sweep measures the suite rather than itself. If it is ever renamed the
  # sweep self-hits and fails loudly, which is the safe direction.
  files <- files[basename(files) != "test-fork-safety.R"]
  expect_gt(length(files), 1L)

  forking <- character()
  unguarded <- character()

  for (f in files) {
    for (block in .fork_test_blocks(f)) {
      if (!.forks_multicore(block)) next
      label <- paste0(basename(f), ": ", sub("^test_that\\(", "", block[[1L]]))
      forking <- c(forking, label)
      if (!.skips_mac(block)) {
        unguarded <- c(unguarded, label)
      }
    }
  }

  # Positive control: a sweep that found nothing must not pass for free. These
  # are the blocks the rule exists to police.
  expect_gt(length(forking), 0L)
  expect_identical(unguarded, character())
})
