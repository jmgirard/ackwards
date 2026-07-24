# ackwards 0.1.1

CRAN resubmission of the first release, addressing reviewer feedback on the
0.1.0 submission; also picks up everything added since that submission.

* **Cleaner deep-hierarchy diagrams (`k >= 10`).** `ba_layout()` now orders each
  level by a traversal of the primary-parent forest, laying every subtree out as
  a contiguous block. This drives primary-tree edge crossings to zero in deep
  hierarchies (the "bent levels" the previous single-pass ordering left behind —
  e.g. 3 crossings down to 0 on the 155-variable AMH example at `k_max = 10`),
  with no change to shallow layouts and the primary-child x-placement unchanged.
  Separately, `autoplot(show_r = TRUE)` now dodges overlapping edge-correlation
  labels apart so dense diagrams stay legible. The layout stays fully
  deterministic.

* **`autoplot(drop_pruned = TRUE)` can now draw secondary correlation edges.**
  A new `show_secondary = TRUE` adds the between-level correlations the pruned
  view otherwise hides — every kept cross-level factor pair with
  `|r| >= cut_show` that is not the single strongest-ancestor primary edge,
  including a factor's weaker second parents *and* direct skip-level
  correlations (a skip-level `|r|` is its own non-transitive fact). They render
  in a channel deliberately distinct from the primary edges — dimmed and thinner,
  drawn beneath them — while still inheriting the sign encoding (`sign_by`), so
  the sign color/linetype is never conflated with the secondary channel. Default
  `FALSE` leaves the pruned view unchanged. The visualization vignette
  illustrates the new argument.

* **`prune(x, "artifact")` now reports a near-redundant band.** A new
  `x$prune$near_redundant` table flags cross-level factor pairs that sit *just
  below* the redundancy thresholds — where `prune(x, "redundant")` drops full
  redundancy (`|r| >= redundancy_r`), the band surfaces the messier candidates a
  hair under it (e.g. `|r| = 0.89` / `phi = 0.94`), which Forbes (2023) treats as
  the main use of the artifact flags. A pair is flagged when its direct
  (skip-level) `|r|` **or** its Tucker `phi` falls within the new `near_margin`
  argument (default `0.1`) below the corresponding threshold, and the pair is not
  itself fully redundant. The band is report-only — no factor is dropped on its
  basis. Under EFA/ESEM, `redundancy_phi` now auto-resolves in artifact mode too
  (announced via cli), so the `phi` band is active; under PCA only the `|r|` band
  applies. The Forbes-extension vignette's artifact example now illustrates a
  genuinely near-redundant pair rather than one redundancy already drops.

* **Fixed compatibility with lavaan 0.7.** lavaan 0.7 renamed its
  sample-statistics slot argument (breaking the ESEM engine's multi-level
  reuse of anchor-level sample statistics — every level beyond k = 1 failed
  to build) and now requires an explicit `ordered = FALSE` to use WLSMV/ULSMV
  with continuous data. The ESEM engine detects the installed lavaan's
  argument vocabulary and works with both lavaan >= 0.7 and >= 0.6-13.

* **`print()` for `suggest_k()` now renders an aligned criteria table.** The
  per-criterion evidence prints as a column-aligned grid (one row per k, one
  column per requested criterion) with a header row and a glyph legend, instead
  of the previous per-k concatenated lines. Numeric columns are right-aligned so
  the optimal-k star never shifts the decimals, and the "retained" tick is drawn
  in text presentation so terminals that render it as a wide emoji keep the
  columns aligned. The returned object and all criteria values are unchanged;
  only the printed layout is new.

* **Clarified documentation (no behavior change).** The `?prune` help and the
  bass-ackward vignette now explain the default `redundancy_criterion = "direct"`
  as a *star anchored on the chain's deepest factor* — each ancestor correlates
  directly with that leaf — distinct from adjacent-hop chaining and from an
  all-pairs screen. The `?ackwards` rationale for varimax now states the reason
  correctly: the between-level edge algebra is exact for any linear scoring, so
  orthogonality is an interpretive choice (it keeps within-level factors
  uncorrelated, so between-level edges are not confounded by within-level factor
  intercorrelation), not a numerical necessity. The artifact-mode discussion now
  frames automated flags as *removing* investigator degrees of freedom, leaving
  only the substantive drop decision to the researcher.

* **New vignette: "Reproducing Forbes (2023): The AMH Applied Example"**
  (`vignette("ackwards-forbes2023")`). A full worked reproduction of the
  paper's 155-variable applied example on the bundled `forbes2023` dataset:
  the 10-level hierarchy, skip-level correlations, the redundancy chase
  (including where the default `direct` criterion and the `adjacent` opt-in
  disagree), and the pruned-factor diagram in the paper's publication style.

* Corrected the Forbes (2023) article title in the `forbes2023` help page
  ("bass-ackward method", per the published title).

* **Corrected the `n_obs` advisory for PCA on correlation-matrix input.** The
  message (and the `n_obs` help text) claimed supplying `n_obs` would enable
  chi-square/RMSEA/TLI, but the PCA engine's level fit is eigenvalue-based and
  never computes them; `n_obs` is recorded in the result metadata and feeds
  the N-based sampling-adequacy checks only. The message now says so. No
  behavior changes.

* Completed the Goldberg (2006) reference to its full published title ("Doing
  it all Bass-Ackwards: The development of hierarchical factor structures from
  the top down") in the `ackwards()` help page and four vignettes.

* **Corrected the historical citations for `comparability()`'s split-half
  benchmarks.** The `.90` replication threshold traces to Everett (1983) and
  to its use in Goldberg's lexical research program by Saucier (1997) and
  Saucier, Georgiades, Tsaousis, and Goldberg (2005) — not to Goldberg
  (1990), which contains no split-half analyses and is no longer cited for
  this purpose. The `.95` reference line is now sourced to Lorenzo-Seva and
  ten Berge (2006). Affects the roxygen help page, `print()`/`autoplot()`
  footer text, the README, and the replicability-workflow vignette; no
  behavior changes.

* **Sourced `suggest_k()`'s k-selection guidance.** The consensus-range stance
  now cites Lim and Jahng (2019) with Achim's (2021) counterpoint, and the
  "PA-PC tends to overextract" note cites Saucier (1997) alongside Forbes
  (2023), in both the help page and the suggest-k vignette; no behavior changes.

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
