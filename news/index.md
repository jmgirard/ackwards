# Changelog

## ackwards (development version)

- **Factorability diagnostics.** New
  [`factorability()`](https://jmgirard.github.io/ackwards/reference/factorability.md)
  screens a dataset (or a correlation matrix) before you fit: it reports
  the Kaiser-Meyer-Olkin measure of sampling adequacy (overall and per
  item), Bartlett’s test of sphericity, the sample size and `N:p` ratio,
  and the Ledermann bound on how many common factors `p` variables can
  identify. It reports the numbers and their conventional bands rather
  than a pass/fail flag — every cutoff here is a contested rule of
  thumb.
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  now runs a light version of this screen internally and warns (once)
  only at the consequential extreme: when `k_max` exceeds the Ledermann
  bound for an EFA/ESEM fit (components are exempt), or when sampling
  adequacy is poor (KMO below .5, or `N:p` below 5).

- **Persistent factor labels.** New
  [`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md)
  attaches substantive names to factors (e.g. `"Neuroticism"` for
  `"m5f1"`) on the object itself, and
  [`factor_labels()`](https://jmgirard.github.io/ackwards/reference/factor_labels.md)
  reads them back. Once set, the labels are shown by
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html), and
  [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  — no need to re-supply them each time. Labels are display only: factor
  IDs (`m{k}f{j}`) are never changed, an unlabeled factor falls back to
  its ID, and the labels ride along through
  [`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md),
  [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md),
  [`augment()`](https://generics.r-lib.org/reference/augment.html), and
  [`predict()`](https://rdrr.io/r/stats/predict.html). In
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) the label
  columns (`factor_label`, or `from_label`/`to_label` for edges) appear
  only when labels are set, so existing output is unchanged otherwise.
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  remains the convenient scaffold to fill in and pass to
  [`set_factor_labels()`](https://jmgirard.github.io/ackwards/reference/set_factor_labels.md).

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
