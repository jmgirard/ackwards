# Verify every anchor cited in cairn/references/source-departures.md resolves (M72).
#
# The departures ledger cites D-entries, DESIGN.md principles (IPn/GPn), sibling
# reference notes ([[citekey]]), and R source files. This checker asserts each
# of those anchors still resolves, turning silent citation-rot (a renumbered
# D-entry, a deleted note, a moved file) into a mechanical failure -- run in CI
# and at the DoD gate. It does NOT verify completeness: that a newly-added
# departure was recorded here is a process obligation (the ledger's Maintenance
# clause + IP9's D-entry requirement), not a mechanically-checkable property.
#
# cairn/ and tools/ are both .Rbuildignore'd, so the built package contains
# neither this script nor the ledger -- tests/testthat/test-ledger-anchors.R
# therefore SKIPS there and runs only in the source checkout; CI also runs this
# checker standalone, so the skip loses no coverage.
#
# Base R only. USAGE (from the package root):
#   Rscript tools/check-ledger-anchors.R
#
# Sourced by tests/testthat/test-ledger-anchors.R, which reuses
# check_ledger_anchors() directly; sourcing defines the function without running
# the script body (guarded by sys.nframe()).

check_ledger_anchors <- function(root = ".") {
  ledger <- file.path(root, "cairn", "references", "source-departures.md")
  if (!file.exists(ledger)) return(character(0)) # built package: nothing to check
  txt <- paste(readLines(ledger, warn = FALSE), collapse = "\n")
  problems <- character(0)

  hits <- function(pattern) {
    unique(unlist(regmatches(txt, gregexpr(pattern, txt, perl = TRUE))))
  }

  # D-0NN -> "### D-0NN" heading in DECISIONS.md
  dec <- file.path(root, "cairn", "DECISIONS.md")
  dec_lines <- if (file.exists(dec)) readLines(dec, warn = FALSE) else character(0)
  for (d in hits("D-[0-9]{3}")) {
    if (!any(grepl(paste0("^### ", d, "\\b"), dec_lines))) {
      problems <- c(problems, paste0("D-entry not found in DECISIONS.md: ", d))
    }
  }

  # IPn / GPn -> "- IPn:" / "- GPn:" bullet in DESIGN.md
  des <- file.path(root, "cairn", "DESIGN.md")
  des_lines <- if (file.exists(des)) readLines(des, warn = FALSE) else character(0)
  for (p in hits("\\b[IG]P[0-9]+\\b")) {
    if (!any(grepl(paste0("^- ", p, ":"), des_lines))) {
      problems <- c(problems, paste0("principle not found in DESIGN.md: ", p))
    }
  }

  # [[citekey]] -> cairn/references/<citekey>.md exists
  for (link in hits("\\[\\[[a-z0-9]+\\]\\]")) {
    key <- gsub("\\[|\\]", "", link)
    if (!file.exists(file.path(root, "cairn", "references", paste0(key, ".md")))) {
      problems <- c(problems, paste0("linked reference note missing: ", key, ".md"))
    }
  }

  # R/<file>.R -> file exists
  for (rf in hits("R/[A-Za-z0-9_]+\\.R")) {
    if (!file.exists(file.path(root, rf))) {
      problems <- c(problems, paste0("R source file missing: ", rf))
    }
  }

  problems
}

if (sys.nframe() == 0L) {
  problems <- check_ledger_anchors(".")
  if (length(problems)) {
    message("Ledger anchor check FAILED:")
    message(paste0("  - ", problems, collapse = "\n"))
    quit(status = 1L)
  }
  message("Ledger anchor check: all anchors in source-departures.md resolve.")
}
