# Snapshot safety net for print.ackwards / print.summary_ackwards (M59).
# The M59 dedup routes both surfaces through shared helpers; these snapshots
# lock every user-visible string, count, and glyph so a behaviour-preserving
# refactor cannot silently reword output. Two output unifications are
# deliberate (cumulative-variance percent format; tick/cross glyph source) and
# are captured as reviewed, intentional snapshot updates.
#
# cli writes to stderr in non-interactive mode, so capture the message stream
# and strip ANSI (colour is also stripped by testthat's reproducible output,
# but ansi_strip keeps the snapshot glyph-only and locale-stable).

snap_print <- function(o) {
  cat(cli::ansi_strip(capture.output(print(o), type = "message")), sep = "\n")
}

test_that("print/summary snapshot: PCA, converged, no prune", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(sim16, k_max = 4))
  expect_snapshot(snap_print(x))
  expect_snapshot(snap_print(summary(x)))
})

test_that("print/summary snapshot: EFA fit-index glyph line", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(bfi25[, 1:10], k_max = 3, engine = "efa"))
  expect_snapshot(snap_print(summary(x)))
})

test_that("print/summary snapshot: pruning digest (redundant + artifact)", {
  skip_if_not_installed("psych")
  x <- cached(
    ackwards(bfi25[, 1:10], k_max = 4) |> prune(c("redundant", "artifact"))
  )
  expect_snapshot(snap_print(x))
  expect_snapshot(snap_print(summary(x)))
})

test_that("print snapshot: a non-converged level shows the red cross", {
  skip_if_not_installed("psych")
  # Convergence is data, not an error (Invariant 7): a non-converged level is
  # marked and printed, not dropped. Patch the flag (real non-convergence needs
  # a deep ESEM fit) to exercise print()'s red-cross branch.
  x <- cached(ackwards(sim16, k_max = 3))
  x$levels[["2"]]$converged <- FALSE
  expect_snapshot(snap_print(x))
})

test_that("print/summary snapshot: near-singular caution", {
  skip_if_not_installed("psych")
  x <- cached(ackwards(sim16, k_max = 3))
  x$meta$near_singular <- TRUE
  x$meta$min_eigenvalue <- 0.0012
  expect_snapshot(snap_print(x))
  expect_snapshot(snap_print(summary(x)))
})
