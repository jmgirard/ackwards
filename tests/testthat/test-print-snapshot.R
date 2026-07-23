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

# --- suggest_k() print table (M75) -----------------------------------------
# suggest_k()'s PA/MAP/VSS numbers are stochastic (Monte Carlo PA, engine
# fits), so a snapshot of a live run would churn. Build synthetic objects with
# fixed values — the same technique test-suggest_k.R's M42/m1 test uses — so the
# table layout is snapshotted deterministically.

.sk_synthetic <- function(requested, cd_available = FALSE, k_max = 5L) {
  structure(
    list(
      k_parallel_pc = 4L, k_parallel_fa = 3L,
      k_map = 4L, k_vss1 = 4L, k_vss2 = 4L,
      k_cd = if (cd_available) 2L else NA_integer_,
      cd_available = cd_available,
      cd_rmse = if (cd_available) {
        c(0.12, 0.08, 0.09, 0.10, 0.11)
      } else {
        NULL
      },
      criteria = data.frame(
        k = seq_len(k_max),
        ev_obs = c(5.2, 2.1, 1.4, 1.1, 0.6),
        ev_obs_fa = c(4.8, 1.7, 1.0, 0.7, 0.2),
        pa_pc_quant = c(1.3, 1.2, 1.15, 1.1, 1.05),
        pa_pc_suggested = seq_len(k_max) <= 4L,
        pa_fa_quant = c(0.9, 0.8, 0.75, 0.7, 0.65),
        pa_fa_suggested = seq_len(k_max) <= 3L,
        map = c(0.0668, 0.0495, 0.0402, 0.0230, 0.0320),
        vss1 = c(0.5505, 0.7439, 0.7685, 0.7884, 0.7881),
        vss2 = c(0.0000, 0.7834, 0.8641, 0.9002, 0.8882)
      ),
      criteria_requested = requested,
      k_max = k_max, n_obs = 1000L, n_vars = 16L,
      cor = "pearson", input_type = "data"
    ),
    class = "suggest_k"
  )
}

test_that("print.suggest_k() snapshot: full criteria table with CD", {
  expect_snapshot(
    snap_print(.sk_synthetic(
      c("pa_pc", "pa_fa", "map", "vss", "cd"),
      cd_available = TRUE
    ))
  )
})

test_that("print.suggest_k() snapshot: subset criteria (map + vss only)", {
  expect_snapshot(snap_print(.sk_synthetic(c("map", "vss"))))
})

test_that("print.suggest_k() renders a column-aligned table with a legend", {
  sk <- .sk_synthetic(c("pa_pc", "pa_fa", "map", "vss", "cd"), cd_available = TRUE)
  out <- cli::ansi_strip(capture.output(print(sk), type = "message"))

  # The table rows: a header line naming every requested criterion, then one
  # line per k. Locate them by the leading "k" header and the "k = "-free data
  # rows that carry the MAP values.
  header <- out[grepl("PA-PC", out) & grepl("VSS-1", out) & grepl("MAP", out)]
  expect_length(header, 1L)
  expect_match(header, "PA-PC")
  expect_match(header, "PA-FA")
  expect_match(header, "MAP")
  expect_match(header, "VSS-1")
  expect_match(header, "VSS-2")
  expect_match(header, "CD")

  # Data rows carry the fixed MAP values; every table line (header + data) must
  # share one width — the proof that columns are padded to a common grid.
  data_rows <- out[grepl("0\\.0668|0\\.0495|0\\.0402|0\\.0230|0\\.0320", out)]
  expect_length(data_rows, 5L)
  widths <- nchar(c(header, data_rows))
  expect_true(all(widths == widths[1]))

  # Numeric columns are right-aligned: the decimal point of the MAP value sits
  # at the same character index in every data row (star markers live outside the
  # numeric field, so they cannot shift it).
  map_dot <- vapply(
    data_rows,
    function(ln) regexpr("\\.0(668|495|402|230|320)", ln)[[1]],
    integer(1)
  )
  expect_true(all(map_dot == map_dot[1]))
  expect_true(all(map_dot > 0))

  # A legend explains the glyphs.
  expect_true(any(grepl("optimal", out) & grepl("retained", out)))
})
