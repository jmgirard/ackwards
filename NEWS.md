# ackwards 0.1.0

First public release. `ackwards` implements Goldberg's (2006) bass-ackwards
method and its modern extensions: it fits factor solutions at every level from 1
to `k` and characterizes the hierarchy through the between-level factor-score
correlations that connect them. Initial features, roughly in order of importance:

* **`ackwards()`** — fit the hierarchy with a PCA, EFA, or ESEM engine;
  between-level edges from Waller's (2007) exact `W′RW` algebra.
* **`suggest_k()`** — bracket a plausible depth range from five retention
  criteria (parallel analysis, MAP, VSS, Comparison Data).
* **`comparability()`** — gate hierarchy depth on split-half replicability
  (Everett 1983; Goldberg 1990).
* **`factorability()`** — screen a dataset (or correlation matrix) before you
  fit: Kaiser-Meyer-Olkin sampling adequacy (overall and per item), Bartlett's
  test of sphericity, the `N:p` ratio, and the Ledermann bound on identifiable
  factors, reported as numbers-and-bands rather than pass/fail. `ackwards()`
  runs a light version internally and warns only at the consequential extreme
  (`k_max` above the Ledermann bound for EFA/ESEM, or poor sampling adequacy).
* **Forbes (2023) extension** — `pairs = "all"` skip-level edges, `prune()` for
  redundant/artifactual factors, and `boot_edges()` bootstrap edge CIs.
* **Ordinal data** — `cor = "polychoric"` (WLSMV for ESEM), `check_items()`
  pre-analysis screening, and near-singularity diagnostics.
* **Interpretation & scoring** — `top_items()` for reading factors;
  `label_template()`/`set_factor_labels()`/`factor_labels()` to attach persistent
  substantive names shown across `print()`, `summary()`, `autoplot()`, `tidy()`,
  and `top_items()`; `augment()`/`predict()` for factor scores, in and out of
  sample.
* **Output** — `autoplot()` hierarchy diagrams and `tidy()`/`glance()`/`summary()`.
* **Bundled data** — `bfi25` (ordinal Big Five) and `sim16` (continuous), plus
  eight vignettes.

Beyond base R, `psych` is the only hard dependency; `lavaan`, `ggplot2`, and
others are optional.
