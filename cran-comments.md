# CRAN submission comments — ackwards 0.2.0

## Resubmission

This replaces the 0.2.0 tarball auto-rejected on 2026-07-27, whose noLD special
check reported 1 ERROR. That is the same failure recorded against the published
0.1.1 under both ATLAS and noLD (CRAN issue database, deadline 2026-08-21), and
it is fixed here.

The cause: the internal primary-parent assignment used
`apply(abs(E), 2, which.max)`. `which.max()` returns `integer(0)` for an all-NA
column, so `apply()` could not simplify and returned a list where an integer
vector was expected, and the caller aborted on "invalid subscript type 'list'".
An all-NA column means a factor whose scores carry no usable variance, which
arises when the requested number of factors exceeds what the data can identify.
Under most numerical libraries the affected fit is merely ill-conditioned;
under ATLAS and no-long-double it degenerates fully. Such a level is now
recorded, warned about, and skipped — matching how the package has always
treated a non-converged level — and the hierarchy keeps every level above it.

Verified on R-hub's `atlas` and `nold` containers; results in the table below.

## Reason for this update

This release arrives sooner than the usual gap between updates because it
fixes check failures on the currently published version — the ATLAS/noLD
issue described above, and the one below. `ackwards` 0.1.1
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

* **"Days since last update"** — expected, and the reason is the two sections
  above: 0.1.1 was accepted on 2026-07-24 and is now failing two of CRAN's
  check flavours. This submission exists to clear both. I would otherwise have
  held the release for the usual interval, and I do not intend to submit again
  on this cadence.

| Platform | R version | Errors / Warnings |
|---|---|---|
| macOS 26 arm64 (local, `--as-cran`) | R 4.6.1 | 0 / 0 |
| ubuntu-latest (GitHub Actions) | release / devel / oldrel-1 | 0 / 0 |
| macos-latest (GitHub Actions) | release | 0 / 0 |
| macos-latest arm64 (GitHub Actions) | oldrel-1 | 0 / 0 |
| windows-latest (GitHub Actions) | release | 0 / 0 |
| win-builder | R-devel | 0 / 0 |
| R-hub `atlas` (ATLAS BLAS) | R-devel | 0 / 0 |
| R-hub `nold` (no long double) | R-devel | 0 / 0 |

The `atlas` and `nold` rows are the two flavours the previous 0.2.0 tarball
failed on; both now report `Status: OK` with 0 test failures.

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
