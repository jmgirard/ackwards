# CRAN submission comments for ackwards 0.3.0

## Reason for this update

This is a feature release. It adds oblique rotation as an option (varimax
stays the default, and default output does not change), partialled edge
coefficients and within-level factor correlations in `tidy()`, fixes to the
ESEM engine, and revised documentation. GPArotation joins `Suggests` for the
oblique option. See NEWS.md.

The published 0.2.0 shows OK on all 13 CRAN check flavours (results page
read on 2026-10-19).

## R CMD check results

0 errors | 0 warnings | 0 notes on every platform below.

| Platform | R version | Errors / Warnings / Notes |
|---|---|---|
| macOS 27 arm64 (local, `--as-cran`) | R 4.6.1 | 0 / 0 / 0 |
| ubuntu-latest (GitHub Actions) | release / devel / oldrel-1 | PENDING |
| macos-latest (GitHub Actions) | release | PENDING |
| macos-latest arm64 (GitHub Actions) | oldrel-1 | PENDING |
| windows-latest (GitHub Actions) | release | PENDING |
| win-builder | R-devel | PENDING |
| R-hub `atlas` (ATLAS BLAS) | R-devel | PENDING |
| R-hub `nold` (no long double) | R-devel | PENDING |

The `atlas` and `nold` rows cover the two alternative-numerics flavours that
an earlier 0.2.0 tarball failed on.

## Package scope

`ackwards` implements Goldberg's (2006) bass-ackwards hierarchical structural
analysis method and its modern descendants (Waller 2007, Forbes 2023). It has
PCA, EFA, and ESEM engines and polychoric support for ordinal data. It also
provides an analysis toolkit: `suggest_k()` (depth range), `comparability()`
(split-half replicability gate), `prune()` (Forbes redundancy), `boot_edges()`
(bootstrap edge CIs), `predict()` (out-of-sample scoring), and `check_items()`
(pre-analysis item screening). It is the only CRAN package providing an ESEM
engine for this method and the Forbes (2023) redundancy-pruning extension.

## Suggests dependencies

`psych` is in `Imports` as it is the engine substrate for the default PCA and
EFA paths and for polychoric correlations. Moving it to `Suggests` means an
install prompt for core functionality. All other optional
dependencies (`covr`, `EFAtools`, `future`, `future.apply`, `ggplot2`,
`GPArotation`, `gt`, `knitr`, `lavaan`, `rmarkdown`, `testthat`) are in
`Suggests`, gated behind `rlang::is_installed()` / `requireNamespace()` calls.
The package installs and loads without any of them. When a required `Suggests`
package is absent, the function stops with an informative error message.

## Downstream dependencies

There are no reverse dependencies on CRAN to check.

## Notes on `\donttest{}`

The following example blocks use `\donttest{}` because they run stochastic
Monte-Carlo simulation, fit many models, or render several `ggplot2` figures,
which is slow for routine check timing:

- `suggest_k()` and `autoplot.suggest_k()`: parallel analysis / Comparison Data
  (Monte-Carlo, about 10 to 20 s per call).
- `comparability()` and `autoplot.comparability()`: repeated split-half refits.
- `boot_edges()`: nonparametric bootstrap over many replicates.
- `autoplot.ackwards()`: fits several models and renders about 15 figures.
- `forbes2023`: fits the 10-level hierarchy of the 155-variable applied example.

All other examples (including `check_items()`, `ackwards()`, `prune()`,
`top_items()`, `tidy()`/`glance()`/`augment()`/`predict()`) run in full during
`R CMD check`, guarded by `requireNamespace()` where they touch a `Suggests`
package.
