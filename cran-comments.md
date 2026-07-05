# CRAN submission comments — ackwards 0.1.0

## Resubmission

This resubmission addresses the reviewer's request to reduce the vignette
re-build time (previously "checking re-building of vignette outputs … [317s]").
Seven of the eight vignettes are now **precomputed**: they are authored as
`*.Rmd.orig` and knitted ahead of time into static `*.Rmd` with all results and
figures baked in (the rOpenSci precompute pattern), so `R CMD check` builds them
with pandoc only. The one remaining live vignette fits a single small model.
Local vignette re-build time dropped from ~64 s to ~6 s; the heavy per-vignette
model fits, parallel analysis, bootstrap, and split-half refits no longer run on
CRAN. No user-facing content changed — the rendered vignettes are identical.

Several slow `@examples` were also moved onto the bundled continuous `sim16`
dataset (from the ordinal `bfi25` polychoric path) where the dataset was
incidental to the example.

## R CMD check results

0 errors | 0 warnings on every platform below.

The only notes are:

* **New submission** — expected; the package is not yet on CRAN. Emitted by the
  incoming-feasibility check on every `--as-cran` run.
* **Possibly misspelled words** (win-builder R-devel only) — all correctly
  spelled; see Notes below.

| Platform | R version | Errors / Warnings |
|---|---|---|
| macOS 26 arm64 (local, `--as-cran`) | R 4.6.1 | 0 / 0 |
| ubuntu-latest (GitHub Actions) | release / devel / oldrel-1 | 0 / 0 |
| macos-latest (GitHub Actions) | release | 0 / 0 |
| windows-latest (GitHub Actions) | release | 0 / 0 |
| win-builder | R-devel | 0 / 0 |

## Notes

### "Possibly misspelled words" in DESCRIPTION

win-builder flags the following, all correctly spelled:

- `ackwards` / `Ackwards` — the package name, a deliberate play on "bass-ackwards".
- `PCA`, `EFA`, `ESEM` — standard abbreviations for Principal Component Analysis,
  Exploratory Factor Analysis, and Exploratory Structural Equation Modeling.
- `Goldberg`, `Waller`, `Forbes` — author surnames of the cited method papers.

No misspellings are present.

## Package scope

`ackwards` implements Goldberg's (2006) bass-ackwards hierarchical structural
analysis method and its modern descendants (Waller 2007; Forbes 2023). Alongside
the PCA / EFA / ESEM engines and polychoric support for ordinal data, it provides
a full analysis toolkit: `suggest_k()` (depth range), `comparability()`
(split-half replicability gate), `prune()` (Forbes redundancy), `boot_edges()`
(bootstrap edge CIs), `predict()` (out-of-sample scoring), and `check_items()`
(pre-analysis item screening). It is the only CRAN package providing an ESEM
engine for this method and the Forbes (2023) redundancy-pruning extension.

## Suggests dependencies

`psych` is in `Imports` as it is the engine substrate for the default PCA and
EFA paths and for polychoric correlations; placing it in `Suggests` would require
an install prompt for core functionality. All other optional dependencies
(`lavaan`, `ggplot2`, `EFAtools`, `future`, `future.apply`, `gt`, `knitr`,
`rmarkdown`, `testthat`, `covr`) are in `Suggests`, gated behind
`rlang::is_installed()` / `requireNamespace()` calls. The package installs and
loads without any of them; functionality degrades gracefully with informative
error messages when a required `Suggests` package is absent.

## Downstream dependencies

This is a new package (first submission). There are no reverse dependencies on
CRAN to check.

## Notes on `\donttest{}`

The following example blocks use `\donttest{}` because they run stochastic
Monte-Carlo simulation, fit many models, or render several `ggplot2` figures —
comprehensive coverage but slow for routine check timing:

- `suggest_k()` and `autoplot.suggest_k()` — parallel analysis / Comparison Data
  (Monte-Carlo, ~10–20 s per call).
- `comparability()` and `autoplot.comparability()` — repeated split-half refits.
- `boot_edges()` — nonparametric bootstrap over many replicates.
- `autoplot.ackwards()` — fits several models and renders ~15 figures.

All other examples (including `check_items()`, `ackwards()`, `prune()`,
`top_items()`, `tidy()`/`glance()`/`augment()`/`predict()`) run in full during
`R CMD check`, guarded by `requireNamespace()` where they touch a `Suggests`
package.
