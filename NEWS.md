# ackwards 0.0.0.9000 (dev)

* Added ESEM engine (`method = "esem"`) using `lavaan::efa()` (requires lavaan
  >= 0.6-13) with WLSMV ordinal estimation, rotation-aware SEs in `loadings_se`,
  and per-level fit indices (CFI, TLI, RMSEA, SRMR). Self-computed tenBerge
  weights keep `compute_edges()` on the algebra path for ESEM as for EFA.
* Added `cor = "polychoric"` as a general basis option for all three engines.
  PCA/EFA use `psych::polychoric()` with NPD smoothing; ESEM uses lavaan's own
  polychoric matrix (extracted from WLSMV sampstat) to ensure consistency with
  the fitted model.
* Added `estimator` argument for the ESEM engine (`NULL` auto-selects `"WLSMV"`
  for `cor = "polychoric"`, `"ML"` otherwise; user-overridable).
* Added `loadings_se` field to the §4 level contract (p × k matrix; `NULL` for
  PCA and EFA engines that do not produce rotation-aware SEs).
* Updated ordinal-detection warning: no longer fires when `cor = "polychoric"`
  is already set; message updated to remove "future release" language.

* Added EFA engine (`method = "efa"`) using `psych::fa()` with tenBerge scoring weights,
  keeping the W'RW algebra path in `compute_edges()`. Convergence failures truncate the
  hierarchy at the deepest converged level (warn + valid object). Heywood cases warn but
  do not truncate. Orthogonal rotation only for now; `fm` controls extraction method
  (`"minres"` default, `"ml"`, `"pa"`).
* Added `ackwards()` implementing Goldberg's (2006) bass-ackwards method with a PCA engine.
* Added `compute_edges()` computing between-level factor-score correlations via W'RW algebra.
* Added `print.ackwards()`, `tidy.ackwards()`, and `glance.ackwards()` for output.
* Added `ba_layout()` computing a Sugiyama-style barycenter layout for bass-ackwards diagrams.
* Added `autoplot.ackwards()` and `plot.ackwards()` for ggplot2-based hierarchy diagrams
  (requires ggplot2; edges coloured by sign, scaled by |r|, solid/dashed by strength).
* Added `suggest_k()` with parallel analysis and MAP (Velicer) criteria to guide choice of `k`.
* Fixed: `scoring$basis` in both PCA and EFA engines now reflects the actual `cor=` argument
  passed by the user rather than always reporting `"pearson"`.
* Fixed: EFA tenBerge weight fallback now labels `scoring$method = "regression"` when
  `.tenBerge_weights()` fails and psych regression weights are used instead (Invariant 6).
* Added `tidy(x, what = "fit")` returning one row per fit index × level (`chi`, `dof`,
  `p_value`, `RMSEA`, `TLI`, `BIC` for EFA; eigenvalues for PCA).
* Fixed: `rotation = "cfQ"` now errors immediately rather than silently running oblimin.
* Fixed: `scores = TRUE` / `keep_fits = TRUE` now warn that storage is not yet implemented.
* Fixed: `print.ackwards()` cut threshold now reads from the stored object rather than a param.
* Fixed: `k = 1` now errors; at least two levels are required for a hierarchy.
* Fixed: `kappa` default corrected from `1 / (2 * p)` to `1 / p`, which is the
  varimax-equivalent Crawford-Ferguson value (Crawford & Ferguson, 1970; Browne, 2001;
  Kim & Eaton, 2015). The previous value (equamax-adjacent) was silently wrong.
