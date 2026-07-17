# Verify every precomputed vignette is up to date with its source (M65).
#
# The heavy vignettes are authored as `vignettes/*.Rmd.orig` and knitted ahead
# of time into static `*.Rmd` by `vignettes/precompute.R`, which stamps each
# generated `.Rmd` with the md5 of the `.Rmd.orig` it came from. This checker
# re-derives that md5 and fails if a `.Rmd.orig` was edited without re-running
# precompute.R -- turning a stale-output slip from a CLAUDE.md prose warning into
# a mechanical failure (run in CI and at the DoD gate).
#
# Base R only (`tools::md5sum`), so it runs before any dependency install.
#
# USAGE (from the package root):
#   Rscript tools/check-vignette-freshness.R
#
# It is also sourced by tests/testthat/test-vignette-freshness.R, which reuses
# `check_vignette_freshness()` directly; sourcing defines the function without
# running the script body (guarded by `sys.nframe()`).

# Extract the machine-readable stamp line from a generated .Rmd, or NULL if the
# file carries no stamp (e.g. the live ackwards-interpret.Rmd).
.read_vignette_stamp <- function(rmd) {
  lines <- readLines(rmd, warn = FALSE)
  hit <- grep("precompute-stamp:", lines, value = TRUE, fixed = TRUE)
  if (length(hit) == 0L) {
    return(NULL)
  }
  line <- hit[[1L]]
  list(
    source = sub(".*source=(\\S+).*", "\\1", line),
    md5 = sub(".*md5=([0-9a-f]+).*", "\\1", line)
  )
}

# Returns a character vector of problems (empty = all fresh).
check_vignette_freshness <- function(vign_dir = "vignettes") {
  problems <- character(0L)

  # Forward: every .Rmd.orig must have a stamped, hash-matching .Rmd.
  orig_files <- list.files(vign_dir, pattern = "\\.Rmd\\.orig$", full.names = TRUE)
  for (orig in orig_files) {
    rmd <- sub("\\.Rmd\\.orig$", ".Rmd", orig)
    if (!file.exists(rmd)) {
      problems <- c(problems, sprintf(
        "%s has no generated %s (run Rscript vignettes/precompute.R).",
        basename(orig), basename(rmd)
      ))
      next
    }
    stamp <- .read_vignette_stamp(rmd)
    if (is.null(stamp)) {
      problems <- c(problems, sprintf(
        "%s carries no precompute-stamp (run Rscript vignettes/precompute.R).",
        basename(rmd)
      ))
      next
    }
    if (!identical(stamp$source, basename(orig))) {
      problems <- c(problems, sprintf(
        "%s stamp names source '%s' but sits beside '%s'.",
        basename(rmd), stamp$source, basename(orig)
      ))
      next
    }
    actual <- unname(tools::md5sum(orig))
    if (!identical(stamp$md5, actual)) {
      problems <- c(problems, sprintf(
        "%s is STALE: %s was edited since the last precompute (stamp md5 %s, actual %s). Re-run Rscript vignettes/precompute.R.",
        basename(rmd), basename(orig), stamp$md5, actual
      ))
    }
  }

  # Reverse: a stamped .Rmd whose named .Rmd.orig is gone is an orphan.
  rmd_files <- list.files(vign_dir, pattern = "\\.Rmd$", full.names = TRUE)
  for (rmd in rmd_files) {
    stamp <- .read_vignette_stamp(rmd)
    if (is.null(stamp)) {
      next
    }
    orig <- file.path(vign_dir, stamp$source)
    if (!file.exists(orig)) {
      problems <- c(problems, sprintf(
        "%s is stamped but its source %s is missing (orphan stamp).",
        basename(rmd), stamp$source
      ))
    }
  }

  problems
}

# Run only when invoked as a script (sys.nframe() == 0), not when sourced.
if (sys.nframe() == 0L) {
  if (!dir.exists("vignettes")) {
    stop("No vignettes/ directory found; run from the package root.", call. = FALSE)
  }
  problems <- check_vignette_freshness("vignettes")
  if (length(problems) > 0L) {
    message("Vignette freshness check FAILED:")
    for (p in problems) message("  - ", p)
    quit(status = 1L, save = "no")
  }
  message("Vignette freshness OK: every precomputed vignette matches its .Rmd.orig.")
}
