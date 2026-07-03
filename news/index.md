# Changelog

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
  and
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  for reading and naming factors;
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
