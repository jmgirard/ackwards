# CRAN submission comments — ackwards 0.1.1

## Resubmission

This tarball replaces the 0.1.1 submission that was withdrawn from the
incoming queue on 2026-07-16: lavaan 0.7-2, published on 2026-07-16 (after
that submission), renamed an argument the ESEM engine relied on and made the
withdrawn tarball's test suite fail. This build detects the installed
lavaan's argument vocabulary and passes its full test suite under both
lavaan 0.6-21 and 0.7-2 (the ubuntu devel and win-builder results below ran
against lavaan 0.7-2).

The 0.1.1 version itself addresses both points from the CRAN review of the
0.1.0 submission (2026-07-12):

1. **"Please always explain all acronyms in the description text."**
   The Description field now spells out principal component analysis (PCA),
   exploratory factor analysis (EFA), and exploratory structural equation
   modeling (ESEM) at first use.

2. **"You write information messages to the console that cannot be easily
   suppressed" (R/label_template.R).** `label_template()` no longer writes
   to the console. It now returns its result visibly as a small classed
   object (`"ackwards_labels"`), and the copy-paste `c(...)` scaffold is
   rendered by that class's `print()` method — so it displays on a
   top-level call but is silent on assignment or when nested in another
   call. This was the only such call site: all other console output in the
   package already lives in `print()`/`summary()` methods, and
   advisory/progress messaging goes through suppressible `cli` conditions.

The version is bumped to 0.1.1, which also folds in two additions made
since the 0.1.0 submission: a bundled CC-BY 4.0 dataset (`forbes2023`,
attributed in `Authors@R` and `LICENSE.note`) and a corrected default for
the redundancy-pruning criterion (see NEWS.md).

## R CMD check results

0 errors | 0 warnings on every platform below.

The only notes are:

* **New submission** — expected; the package is not yet on CRAN. Emitted by
  the incoming-feasibility check on every `--as-cran` run.
* **Possibly misspelled words** (win-builder R-devel only) — `ackwards` /
  `Ackwards` (the package name, a deliberate play on "bass-ackwards") and
  `Goldberg` / `Waller` / `Forbes` (author surnames of the cited method
  papers). All correctly spelled. The acronyms flagged last round (PCA,
  EFA, ESEM) are now expanded in the Description text itself.

| Platform | R version | Errors / Warnings |
|---|---|---|
| macOS 26 arm64 (local, `--as-cran`) | R 4.6.1 | 0 / 0 |
| ubuntu-latest (GitHub Actions) | release / devel / oldrel-1 | 0 / 0 |
| macos-latest (GitHub Actions) | release | 0 / 0 |
| windows-latest (GitHub Actions) | release | 0 / 0 |
| win-builder | R-devel | 0 / 0 |

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
