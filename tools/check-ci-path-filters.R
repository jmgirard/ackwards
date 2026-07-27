#!/usr/bin/env Rscript
# Verify the CI `paths-ignore` filters (M82).
#
# Three clauses, over every paths-ignore block in .github/workflows/:
#
#   1. No listed literal (glob-free) path is absent from the repo. Four dead
#      entries -- DESIGN.md, MILESTONES.md, ROADMAP.md, nested-cv-guide.md --
#      survived here for months after those files moved under cairn/, which is
#      how the filter silently stopped covering the tracking files.
#   2. Every block ignores the same six cairn tracking paths. Copy-drift
#      between the identical blocks is what let clause 1 rot unnoticed.
#   3. No pattern matches cairn/DECISIONS.md, cairn/DESIGN.md, or anything
#      under cairn/references/. tools/check-ledger-anchors.R reads those from
#      the source checkout, so ignoring them would let a commit that breaks an
#      anchor skip the checker built to catch it.
#
# BASE R ONLY -- like check-vignette-freshness.R and check-ledger-anchors.R,
# this runs as a CI fail-fast step *before* setup-r-dependencies, so it cannot
# use the yaml package. The blocks it parses are a fixed, regular shape
# (`paths-ignore:` followed by `- 'path'` lines), which line parsing handles
# without a YAML reader.
#
# Follows the tools/check-ledger-anchors.R idiom: a function returning a
# character vector of problems, plus a main guard for direct Rscript use.

# Blocks that must exist. Without this, deleting a whole paths-ignore block
# would pass silently -- the clauses below only inspect blocks that are there.
# A new workflow that gains a filter belongs on this list.
REQUIRED_BLOCKS <- c(
  "R-CMD-check.yaml :: push", "R-CMD-check.yaml :: pull_request",
  "test-coverage.yaml :: push", "test-coverage.yaml :: pull_request",
  "macos-oldrel-check.yaml :: push"
)

REQUIRED_CAIRN <- c(
  "cairn/ROADMAP.md", "cairn/LESSONS.md", "cairn/PROFILE.md",
  "cairn/milestones/**", "cairn/legacy/**", "cairn/reviews/**"
)

# Paths read by tools/check-ledger-anchors.R from the source checkout.
LEDGER_READ <- c(
  "cairn/DECISIONS.md", "cairn/DESIGN.md",
  "cairn/references/source-departures.md", "cairn/references/goldberg2006.md"
)

# GitHub's paths-ignore globs: `**` spans `/`, `*` does not, and the pattern
# must match the whole path.
ci_glob_matches <- function(pattern, path) {
  rx <- gsub("\\*\\*", "\001", pattern, perl = TRUE)
  rx <- gsub("([.+^$(){}|\\[\\]\\\\])", "\\\\\\1", rx, perl = TRUE)
  rx <- gsub("\\*", "[^/]*", rx, perl = TRUE)
  rx <- gsub("\001", ".*", rx, perl = TRUE)
  grepl(paste0("^", rx, "$"), path, perl = TRUE)
}

# Returns a named list: "<file> :: <event>" -> character vector of patterns.
ci_parse_blocks <- function(file) {
  lines <- readLines(file, warn = FALSE)
  out <- list()
  event <- "?"
  i <- 1L
  while (i <= length(lines)) {
    ln <- lines[i]
    ev <- regmatches(ln, regexec("^  ([A-Za-z_]+):\\s*$", ln))[[1]]
    if (length(ev) == 2L) event <- ev[2]
    if (grepl("^\\s*paths-ignore:\\s*$", ln)) {
      pats <- character()
      j <- i + 1L
      while (j <= length(lines)) {
        m <- regmatches(
          lines[j],
          regexec("^\\s*-\\s*['\"]?([^'\"#]+?)['\"]?\\s*$", lines[j])
        )[[1]]
        if (length(m) != 2L) break
        pats <- c(pats, m[2])
        j <- j + 1L
      }
      out[[paste0(basename(file), " :: ", event)]] <- pats
      i <- j
      next
    }
    i <- i + 1L
  }
  out
}

check_ci_path_filters <- function(root = ".") {
  wf_dir <- file.path(root, ".github", "workflows")
  files <- list.files(wf_dir, pattern = "\\.ya?ml$", full.names = TRUE)
  if (length(files) == 0L) return("no workflow files found")

  blocks <- unlist(lapply(files, ci_parse_blocks), recursive = FALSE)
  if (length(blocks) == 0L) {
    return("no paths-ignore blocks found -- has the filter been deleted?")
  }

  problems <- character()

  gone <- setdiff(REQUIRED_BLOCKS, names(blocks))
  if (length(gone)) {
    problems <- c(problems, sprintf(
      "paths-ignore block(s) missing entirely: %s", paste(gone, collapse = ", ")
    ))
  }

  for (nm in names(blocks)) {
    pats <- blocks[[nm]]
    if (length(pats) == 0L) {
      problems <- c(problems, sprintf("%s: paths-ignore block is empty", nm))
      next
    }

    literals <- pats[!grepl("[*?\\[]", pats)]
    missing <- literals[!file.exists(file.path(root, literals))]
    if (length(missing)) {
      problems <- c(problems, sprintf(
        "%s: literal path(s) absent from repo: %s", nm,
        paste(missing, collapse = ", ")
      ))
    }

    absent <- setdiff(REQUIRED_CAIRN, pats)
    if (length(absent)) {
      problems <- c(problems, sprintf(
        "%s: missing required cairn entries: %s", nm,
        paste(absent, collapse = ", ")
      ))
    }

    for (g in LEDGER_READ) {
      hit <- pats[vapply(pats, ci_glob_matches, logical(1), path = g)]
      if (length(hit)) {
        problems <- c(problems, sprintf(
          "%s: pattern(s) %s match ledger-read path %s", nm,
          paste(hit, collapse = ", "), g
        ))
      }
    }
  }
  problems
}

if (sys.nframe() == 0L) {
  problems <- check_ci_path_filters(".")
  if (length(problems)) {
    message("CI path-filter check FAILED:")
    message(paste0("  - ", problems, collapse = "\n"))
    quit(status = 1L)
  }
  message("CI path-filter check: all paths-ignore blocks are current and ",
          "leave the ledger-read paths unfiltered.")
}
