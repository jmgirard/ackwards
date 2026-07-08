# Changelog

## ackwards (development version)

- [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  now reads Comparison Data (CD) results from **EFAtools \>= 0.8.0**,
  which restructured `CD()`’s return value (the per-iteration RMSE
  matrix moved from the top-level `RMSE_eigenvalues` field into
  `results[[1]]$rmse_eigenvalues`). Without this, the CD criterion and
  its
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  panel silently dropped out when a current EFAtools was installed.
  Older EFAtools versions still work.

## ackwards 0.1.0

First public release. `ackwards` implements Goldberg’s (2006)
bass-ackwards method and its modern extensions: it fits factor solutions
at every level from 1 to `k` and characterizes the hierarchy through the
between-level factor-score correlations that connect them. Initial
features, roughly in order of importance:

- **[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)**
  — fit the hierarchy with a PCA, EFA, or ESEM engine; between-level
  edges from Waller’s (2007) exact `W′RW` algebra.
- **[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)**
  — bracket a plausible depth range from five retention criteria
  (parallel analysis, MAP, VSS, Comparison Data).
- **[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)**
  — gate hierarchy depth on split-half replicability (Everett 1983;
  Goldberg 1990).
- **[`factorability()`](https://jmgirard.github.io/ackwards/reference/factorability.md)**
  — screen a dataset (or correlation matrix) before you fit:
  Kaiser-Meyer-Olkin sampling adequacy (overall and per item),
  Bartlett’s test of sphericity, the `N:p` ratio, and the Ledermann
  bound on identifiable factors, reported as numbers-and-bands rather
  than pass/fail.
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  runs a light version internally and warns only at the consequential
  extreme (`k_max` above the Ledermann bound for EFA/ESEM, or poor
  sampling adequacy).
- **Forbes (2023) extension** — `pairs = "all"` skip-level edges,
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
  for redundant/artifactual factors, and
  [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)
  bootstrap edge CIs.
- **Ordinal data** — `cor = "polychoric"` (WLSMV for ESEM),
  [`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md)
  pre-analysis screening, and near-singularity diagnostics.
- **Interpretation & scoring** —
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  for reading factors;
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)/[`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md)/[`factor_labels()`](https://jmgirard.github.io/ackwards/reference/factor_labels.md)
  to attach persistent substantive names shown across
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md);
  [`augment()`](https://generics.r-lib.org/reference/augment.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  for factor scores, in and out of sample.
- **Output** —
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  hierarchy diagrams and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)/[`summary()`](https://rdrr.io/r/base/summary.html).
- **Bundled data** — `bfi25` (ordinal Big Five) and `sim16`
  (continuous), plus eight vignettes.

Beyond base R, `psych` is the only hard dependency; `lavaan`, `ggplot2`,
and others are optional.
