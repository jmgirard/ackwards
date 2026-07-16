# ackwards 0.1.1

CRAN resubmission of the first release, addressing reviewer feedback on the
0.1.0 submission; also picks up everything added since that submission.

* **Corrected the historical citations for `comparability()`'s split-half
  benchmarks.** The `.90` replication threshold traces to Everett (1983) and
  to its use in Goldberg's lexical research program by Saucier (1997) and
  Saucier, Georgiades, Tsaousis, and Goldberg (2005) — not to Goldberg
  (1990), which contains no split-half analyses and is no longer cited for
  this purpose. The `.95` reference line is now sourced to Lorenzo-Seva and
  ten Berge (2006). Affects the roxygen help page, `print()`/`autoplot()`
  footer text, the README, and the replicability-workflow vignette; no
  behavior changes.

* **`label_template()` now returns its scaffold visibly (behavior change).**
  The named character vector carries class `"ackwards_labels"`, and the
  editable `c(...)` literal is rendered by its `print()` method instead of
  being written to the console unconditionally. A top-level call looks the
  same as before; assigning the result (`labs <- label_template(x)`) or
  passing it inline (e.g. inside `autoplot()`) is now silent, per CRAN
  policy on console output.

* **DESCRIPTION** now spells out principal component analysis (PCA),
  exploratory factor analysis (EFA), and exploratory structural equation
  modeling (ESEM) per CRAN feedback.

* **New bundled dataset `forbes2023`.** The 155-variable "Assessing Mental Health"
  Spearman correlation matrix that forms Forbes's (2023) applied example is now
  exported, so `ackwards(forbes2023, k_max = 10)` reproduces her worked hierarchy
  directly. It joins `bfi25` and `sim16`. The matrix is redistributed under
  CC-BY 4.0 with attribution to M. K. Forbes (see `LICENSE.note`).

* **`prune("redundant")` gains `redundancy_criterion`, defaulting to `"direct"`
  (behavior change).** Redundancy chains are now traced by the **direct
  (skip-level)** correlation between a factor and each ancestor level — Forbes's
  (2023) actual `ChaseCorrPaths` rule — rather than the previous adjacent-hop
  walk. Because correlation is non-transitive, the two can differ in deep
  (many-level) hierarchies: on shallow ones (e.g. the bundled `sim16`) results
  are unchanged, but a factor can now be flagged redundant with an ancestor it
  correlates with directly even if an intermediate step is weak (and vice
  versa). Pass `redundancy_criterion = "adjacent"` for the old behavior. This
  makes `prune("redundant")` reproduce Forbes's published applied example
  exactly. `print()` and `summary()` name the active criterion. Note that under
  `"direct"`, a chain's `r_to_prev` column reports the adjacent-level
  correlation and can sit below `redundancy_r` (membership is set by the direct
  skip-level link; see `?prune`).

* **Validation.** The Forbes (2023) fidelity suite now also reproduces her
  155-variable "Assessing Mental Health" applied example (`k_max = 10`), not just
  the three simulation studies: between-level correlations match her reference
  implementation to 1.3e-14 across all 45 level-pairs, loading congruences agree
  within her two-decimal rounding, and her redundancy chase is reproduced for all
  54 components. The published matrix ships as a test fixture under CC-BY 4.0
  (see `LICENSE.note`).

* `suggest_k()` now reads Comparison Data (CD) results from **EFAtools >= 0.8.0**,
  which restructured `CD()`'s return value (the per-iteration RMSE matrix moved
  from the top-level `RMSE_eigenvalues` field into `results[[1]]$rmse_eigenvalues`).
  Without this, the CD criterion and its `autoplot()` panel silently dropped out
  when a current EFAtools was installed. Older EFAtools versions still work.

* **Console output consistency.** `summary()`'s per-level fit-index pass/fail
  mark is now the same terminal-adaptive `cli` glyph as `print()`'s convergence
  mark (a tick/cross that degrades to `v`/`x` in a non-UTF-8 console), replacing
  a hard-coded Unicode `✔`/`✘`. And `print()`'s cumulative-variance percentages
  now carry a fixed single decimal (e.g. `20.0%`, where a whole-number value
  previously printed as `20%`), matching `summary()`. Cosmetic only; the
  reported values are unchanged.

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
