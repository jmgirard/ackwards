# CRAN submission comments — ackwards 0.2.0

## Reason for this update

This release arrives sooner than the usual gap between updates because it
fixes a check failure on the currently published version. `ackwards` 0.1.1
shows an ERROR on `r-oldrel-macos-arm64`: two of the package's own tests set
a *forking* parallel plan (`future::multicore`) to verify that parallel and
serial fits agree exactly, and R segfaults inside the forked worker on that
flavour — forking does not mix with the multi-threaded numerical libraries R's
linear algebra runs through (`?mclapply` warns against the combination). Both
tests now skip on macOS and continue to run on every other platform.

Nothing user-facing was affected: `ackwards()` and `boot_edges()` never set a
parallel plan themselves, the default plan is sequential, and choosing one is
the caller's business. The failure was confined to the test suite.

The release also carries the user-visible additions made since 0.1.1 —
publication-figure controls for `autoplot()`, secondary correlation edges in
the pruned view, a near-redundant band in `prune()`, and cleaner layouts for
deep hierarchies (see NEWS.md). Hence the minor version bump rather than a
patch.

## R CMD check results

0 errors | 0 warnings on every platform below.

The only note is:

* **"Days since last update: 3"** — expected, and the reason is the section
  above: 0.1.1 was accepted on 2026-07-24 and then began erroring on
  `r-oldrel-macos-arm64`. This submission exists to clear that ERROR. I would
  otherwise have held the release for the usual interval, and I do not intend
  to submit again on this cadence.

| Platform | R version | Errors / Warnings |
|---|---|---|
| macOS 26 arm64 (local, `--as-cran`) | R 4.6.1 | 0 / 0 |
| ubuntu-latest (GitHub Actions) | release / devel / oldrel-1 | 0 / 0 |
| macos-latest (GitHub Actions) | release | 0 / 0 |
| macos-latest arm64 (GitHub Actions) | oldrel-1 | 0 / 0 |
| windows-latest (GitHub Actions) | release | 0 / 0 |
| win-builder | R-devel | 0 / 0 |

The `macos-latest arm64 / oldrel-1` row is new in this release: it reproduces
the CRAN flavour the 0.1.1 ERROR appeared on, which the previous check matrix
did not cover.

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

There are no reverse dependencies on CRAN to check.

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
