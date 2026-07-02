#!/usr/bin/env Rscript
# Definition-of-done gate (M48). Runs the full CLAUDE.md gate sequence once,
# serially, in one process:
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
note("GATE PASSED (check 0/0/0, coverage 100%%, style/lint clean, pkgdown index complete)")
