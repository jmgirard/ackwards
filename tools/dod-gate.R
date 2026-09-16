#!/usr/bin/env Rscript
# Definition-of-done gate (M48). Runs the full CLAUDE.md gate sequence once,
# serially, in one process:
#   vignette freshness (M65; fail-fast, base R) ->
#   ledger anchors (M72) -> CI path filters (M82) -> prose check (M85) ->
#   devtools::check() (must be 0/0/0, vignettes included)
#   -> covr::package_coverage() (target 100%)
#   -> styler::style_pkg() -> lintr::lint_package()
#   -> pkgdown::check_pkgdown() (mirrors the pkgdown GHA; catches exported
#      topics missing from _pkgdown.yml, which R CMD check does not)
# Usage, from the package root:  Rscript tools/dod-gate.R
# Exits non-zero if any step fails, printing every failure it found.

# Parallel testthat (DESCRIPTION Config) defaults to 2 workers; use the
# machine. Applies to the check()'s test phase and the coverage run alike.
Sys.setenv(TESTTHAT_CPUS = max(1L, parallel::detectCores() - 1L))

failures <- character()
note <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")

# Vignette freshness (M65), fail-fast before the minutes-long check(). The
# test-vignette-freshness.R testthat wrapper SKIPS under check() (the tarball
# .Rbuildignore's every .Rmd.orig), so enforce the guard directly here against
# the source checkout. Base R only; sys.source blocks the script's own body.
fresh_env <- new.env()
sys.source("tools/check-vignette-freshness.R", envir = fresh_env)
fresh_problems <- fresh_env$check_vignette_freshness("vignettes")
if (length(fresh_problems) > 0) {
  for (p in fresh_problems) note("vignette-freshness: %s", p)
  failures <- c(failures, "vignette freshness (re-run Rscript vignettes/precompute.R)")
} else {
  note("vignette-freshness: clean")
}

# Departures-ledger anchor integrity (M72), fail-fast like the vignette check.
# The test-ledger-anchors.R wrapper SKIPS under check() (cairn/ + tools/ are
# .Rbuildignore'd), so enforce the guard directly here against the source
# checkout. Base R only; sys.source blocks the script's own body.
ledger_env <- new.env()
sys.source("tools/check-ledger-anchors.R", envir = ledger_env)
ledger_problems <- ledger_env$check_ledger_anchors(".")
if (length(ledger_problems) > 0) {
  for (p in ledger_problems) note("ledger-anchors: %s", p)
  failures <- c(failures, "ledger anchor integrity (see cairn/references/source-departures.md)")
} else {
  note("ledger-anchors: clean")
}

# CI paths-ignore filter integrity (M82), fail-fast like the two above. The
# filters skip the check matrix on tracking-only commits; four dead entries in
# them went unnoticed for months, and a blanket cairn/** would silently disable
# the ledger-anchor guard. Needs the source checkout (.github/ is
# .Rbuildignore'd). sys.source blocks the script's own body.
ci_env <- new.env()
sys.source("tools/check-ci-path-filters.R", envir = ci_env)
ci_problems <- ci_env$check_ci_path_filters(".")
if (length(ci_problems) > 0) {
  for (p in ci_problems) note("ci-path-filters: %s", p)
  failures <- c(failures, "CI paths-ignore filters (see tools/check-ci-path-filters.R)")
} else {
  note("ci-path-filters: clean")
}

# Plain-English prose check (M85), fail-fast like the three above. Sweeps the
# README, DESCRIPTION, the NEWS development section, and every roxygen line
# for dashes, semicolons, banned phrases, and sentences over 30 words. The
# vignette sources join this domain when their own rewrite lands. Base R
# only; sys.source blocks the script's own body.
prose_env <- new.env()
sys.source("tools/check-prose.R", envir = prose_env)
prose_reports <- prose_env$check_prose(
  c("README.Rmd", "DESCRIPTION", "NEWS.md", "R"),
  banned = prose_env$read_prose_list("prose-banned.txt", "tools"),
  abbrev = prose_env$read_prose_list("prose-abbrev.txt", "tools")
)
if (nrow(prose_reports) > 0) {
  for (i in seq_len(nrow(prose_reports))) {
    note(
      "prose: %s:%d [%s] %s", prose_reports$file[i], prose_reports$line[i],
      prose_reports$class[i], prose_reports$text[i]
    )
  }
  failures <- c(failures, sprintf("prose check: %d report(s) (see tools/check-prose.R)", nrow(prose_reports)))
} else {
  note("prose: clean")
}

t0 <- Sys.time()
chk <- devtools::check(error_on = "never", quiet = TRUE)
n_bad <- length(chk$errors) + length(chk$warnings) + length(chk$notes)
note(
  "check: %d errors | %d warnings | %d notes  [%.0fs]",
  length(chk$errors), length(chk$warnings), length(chk$notes),
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
)
if (n_bad > 0) {
  print(chk)
  failures <- c(failures, "devtools::check() not 0/0/0")
}

t0 <- Sys.time()
cov <- covr::package_coverage()
pct <- covr::percent_coverage(cov)
note("coverage: %.2f%%  [%.0fs]", pct, as.numeric(difftime(Sys.time(), t0, units = "secs")))
if (pct < 100) {
  print(cov)
  failures <- c(failures, sprintf("coverage %.2f%% < 100%%", pct))
}

styled <- styler::style_pkg()
if (any(styled$changed)) {
  note("styler: %d file(s) restyled -- review + commit them", sum(styled$changed))
  failures <- c(failures, "styler changed files (uncommitted)")
} else {
  note("styler: clean")
}

lints <- lintr::lint_package()
if (length(lints) > 0) {
  print(lints)
  failures <- c(failures, sprintf("%d lint(s)", length(lints)))
} else {
  note("lintr: clean")
}

if (rlang::is_installed("pkgdown")) {
  pkgdown_ok <- tryCatch(
    {
      pkgdown::check_pkgdown()
      TRUE
    },
    error = function(e) {
      note("pkgdown: %s", conditionMessage(e))
      FALSE
    }
  )
  if (!pkgdown_ok) failures <- c(failures, "pkgdown::check_pkgdown() failed")
  if (pkgdown_ok) note("pkgdown: reference index complete")
} else {
  note("pkgdown not installed -- eyeball NAMESPACE exports against _pkgdown.yml")
}

cat("\n")
if (length(failures) > 0) {
  note("GATE FAILED:\n- %s", paste(failures, collapse = "\n- "))
  quit(status = 1L)
}
note("GATE PASSED (prose clean, check 0/0/0, coverage 100%%, style/lint clean, pkgdown index complete)")
