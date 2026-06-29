# CRAN submission comments — ackwards 0.1.0

## R CMD check results

0 errors | 0 warnings on all platforms. 1 note on win-builder R-devel (see Notes).

| Platform | R version | Result |
|---|---|---|
| macOS 15.5 arm64 (local) | R 4.4.x | 0/0/0 |
| ubuntu-latest (R-hub) | R-* | 0/0/0 |
| macos-latest arm64 (R-hub) | R-* | 0/0/0 |
| windows-latest (R-hub) | R-* | 0/0/0 |
| win-builder R-devel | R-devel | 0/0/1 (see Notes) |

## Notes

### "Possibly misspelled words" in DESCRIPTION

win-builder flags: Ackwards (2:13), EFA (9:65), ESEM (9:73), ackwards (8:10).

These are all correctly spelled technical terms:
- `ackwards` / `Ackwards` — the package name, a deliberate play on "bass-ackwards"
- `EFA` — standard abbreviation for Exploratory Factor Analysis
- `ESEM` — standard abbreviation for Exploratory Structural Equation Modeling

No misspellings are present.

## Package scope

`ackwards` implements Goldberg's (2006) bass-ackwards hierarchical factor
analysis method and extensions (Waller 2007; Forbes 2023). It is the only CRAN
package providing an ESEM engine for this method, ordinal/polychoric support,
and the Forbes (2023) redundancy-pruning extension.

## Suggests dependencies

`psych` is in `Imports` as it is the engine substrate for the default PCA and
EFA paths; placing it in `Suggests` would require an install prompt for core
functionality. All other optional dependencies (`lavaan`, `ggplot2`, `EFAtools`,
`knitr`, `rmarkdown`, `testthat`, `covr`) are in `Suggests`, gated behind
`rlang::check_installed()` or `requireNamespace()` calls. The package installs
and loads without any of them; functionality degrades gracefully with informative
error messages when a required `Suggests` package is absent.

## Downstream dependencies

This is a new package (first submission). There are no downstream dependencies
on CRAN to check.

## Notes on `\donttest{}`

Two example blocks use `\donttest{}`:
- `suggest_k()` — runs `psych::fa.parallel()` with Monte Carlo simulation
  (stochastic, ~10–20 s per call; not appropriate for routine check timing).
- `autoplot.ackwards()` — fits two models and renders ~10 ggplot2 figures
  (comprehensive coverage but slow for default check timing).

All other examples run under `requireNamespace()` guards and execute fully
during `R CMD check`.
