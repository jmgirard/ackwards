# CRAN submission comments — ackwards 0.1.0

## R CMD check results

Local check (`devtools::check()`):

```
0 errors | 0 warnings | 0 notes
```

Tested on:
- macOS 15.5 (Darwin 25.5.0), R 4.4.x (local dev environment)

**Before submitting, also run:**

```r
# R-devel win-builder (the check CRAN actually runs):
devtools::check_win_devel()

# macOS release builder:
devtools::check_mac_release()

# R-hub v2 (via GitHub Actions — run rhub::rhub_setup() once first):
rhub::rhub_check()
```

## Package scope

`ackwards` implements Goldberg's (2006) bass-ackwards hierarchical factor
analysis method and extensions (Waller 2007; Forbes 2023). It is the only CRAN
package providing an ESEM engine for this method, ordinal/polychoric support,
and the Forbes (2023) redundancy-pruning extension.

## Suggests dependencies

All heavy dependencies (`psych`, `lavaan`, `GPArotation`, `ggplot2`, `EFAtools`)
are in `Suggests`, gated behind `rlang::check_installed()` calls. The package
installs and loads without any of them; functionality degrades gracefully with
informative error messages when a required Suggests package is absent.

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
