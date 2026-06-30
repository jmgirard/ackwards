# DESIGN.md — Bass-Ackwards Hierarchical Structural Analysis (R package)

> Working design brief. Carry this into the repo as the seed for implementation.
> It records what we've settled, the contracts to implement, and recommended defaults.
> §14 logs the resolved decisions and the few small calls left for build time.

Package name: **`ackwards`** (distinctive, searchable, nods to the method's own name; avoids
colliding with `psych::bassAckward()`). Confirm with `available::available("ackwards")` before
committing. Main call: `ackwards::ackwards()`.

---

## 1. Purpose

Implement Goldberg's (2006) bass-ackwards method and its modern descendants: extract a
series of factor/component solutions from 1 to `k`, then characterize the hierarchy by the
correlations between factor/component scores at adjacent (and, for extensions, all) levels.
The package provides multiple extraction engines, principled tooling for choosing `k`,
sign/lineage-aware naming, tidy outputs, and attractive console + graphical reporting.

## 2. Scope & positioning

**Prior art.** `psych::bassAckward()` (Revelle) already implements the core loop for PCA and
EFA with rotation and `tenBerge` scores. We do **not** reimplement factor extraction numerics;
we wrap established engines and add what `psych` does not.

**Our contribution.**
- A first-class **ESEM engine** (lavaan), bringing the Mplus psychopathology/HiTOP workflow
  into R — per-level fit indices and SEs, ordinal/Likert support — which `psych::bassAckward`
  does not provide. (Kim & Eaton 2015; Forbush et al. 2024 are the reference workflows.)
- Proper **ordinal/Likert handling** (polychoric basis, ordinal estimators).
- The **Forbes (2023) extended method** (redundancy + artefact pruning, all-levels correlations).
- A clean, tidy, serializable **result object** with `print`/`summary`/`tidy`/`glance`/`autoplot`.
- A **lineage-aligned layered diagram** rather than a misleading tree.

**Honesty caveat to state in docs and `print`.** A bass-ackwards result is a *series of linked
solutions*, not a fitted hierarchical model. The "hierarchy" is descriptive/emergent (edges are
score correlations), not a higher-order SEM. `psych`'s own docs make this point; we should too.

## 3. Design philosophy & priorities

These are the project's standing priorities (owner-stated) and the recommended stance on each.

- **Manageable dependencies.** `psych` is in `Imports` (M21 — it is the engine substrate for the
  default PCA/EFA paths; requiring an install prompt for the default engine contradicted this
  principle). Everything else heavy or engine-specific remains in `Suggests` behind
  `rlang::check_installed()` guards. Installing the package must not pull in the SEM or plotting
  stacks. `GPArotation` was removed (M21 — varimax routes through `stats::varimax`). (See §12.)
- **Best practices.** testthat 3e, roxygen2, pkgdown, CI (`R CMD check`, lintr/styler), semantic
  versioning, snapshot tests against `psych` for the PCA path, a vignette reproducing a published
  structure.
- **Attractive outputs.** `{cli}` for all console output (`print`, `summary`, progress, warnings);
  `{ggraph}`/`{ggplot2}` for diagrams. Visuals and pretty printing are a feature, not an afterthought.
- **Efficiency / Rcpp.** *Recommendation: do not reach for Rcpp initially.* The heavy compute
  (eigen/ML/WLSMV, polychoric estimation) already lives in compiled code inside lavaan/psych/BLAS.
  The bass-ackwards-specific work is small matrix products (`W'RW`, where items `p` is modest and
  factors `k` is small) and small assignment/pruning steps — none are bottlenecks. The realistic
  hotspots are (a) **many bootstrap/permutation resamples** for stability/CIs, and (b) the **ESEM
  engine at large p**, where lavaan re-deriving the polychoric matrix + asymptotic weight matrix
  (NACOV) at every level dominates. Both are addressed with reuse + `{future}` parallelism rather
  than Rcpp: M26 computes the ESEM sample statistics once (data-derived, level-independent) and
  reuses them via `slotSampleStats=`, then dispatches the independent per-level fits through
  `future.apply` (opt-in via `future::plan()`; serial fallback). Only consider Rcpp if profiling
  proves a remaining C-level bottleneck. Measure before optimizing.
- **Safe, reproducible, and *loud* defaults.** (See §9.) The guiding rule: when the package makes
  a consequential automatic choice (correlation basis, `k`, rotation), it **announces it via cli**
  so the user who never reads an argument still sees what happened.

## 4. Engines

Three engines, all first-class, behind one user-facing function (`ackwards()`) that dispatches.
Rotation is a single fixed choice across all engines: **varimax** (orthogonal; the `T'=T^-1`
property enables the closed-form W'RW algebra). Oblique rotation is out of scope — see §9 and §14.1.

| Engine | Extraction | Typical use | Notes |
|---|---|---|---|
| `"pca"` | principal components | original Goldberg; Forbes extension | fastest, no Heywood, never fails to converge; Waller algebra exact |
| `"efa"` | exploratory factor analysis | factor-model rationale, continuous data | `psych::fa` / `stats::factanal` style |
| `"esem"` | EFA estimated in SEM framework (lavaan) | clinical/HiTOP workflow, ordinal data, downstream SEM linking | per-level fit + SEs; can fail to converge at deeper levels |

**Engine interface (contract).** Every engine, for each level `k = 1..K`, returns a standardized
per-level object:

```r
level <- list(
  k           = <int>,
  loadings    = <p x k matrix>,       # pattern matrix, item x factor
  loadings_se = <p x k matrix | NULL>,# rotation-aware SEs; NULL for PCA/EFA; populated by ESEM
  variance    = <named numeric>,      # per-factor + cumulative variance explained
  fit         = <named numeric>,      # scalar fit indices / eigenvalues (chi, df, CFI, TLI,
                                      # RMSEA, SRMR for ESEM; eigenvalues for PCA)
  converged   = <logical>,            # FALSE allowed; never fatal (see §13)
  factor_cor  = <k x k matrix>,       # within-level factor correlations (I for orthogonal)
  labels      = <character[k]>,       # default m{k}f{j}; user-overridable
  scoring     = <scoring descriptor>  # see §5
)
```

A level that fails to converge is **recorded, warned about, and skipped**, never thrown. The
hierarchy is built up to the deepest converged level.

## 5. Scoring descriptor + `compute_edges()` (the centerpiece)

The between-level edges are computed through **one shared code path** that uses exact algebra when
it can and falls back to materialized scores when it must. This unifies result semantics across
engines without forcing the expensive computation everywhere.

### 5.1 Why the algebra generalizes

For any scoring scheme where scores are a **linear** map of the standardized observed variables,
`S = Z W`, the between-level score correlations are exactly:

```
cor(S_a, S_b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}
           where  R   = input correlation matrix (pearson or polychoric)
                  D_x = diag(W_x' R W_x)   # score variances; NOT assumed 1
```

No scores need to be materialized — only `R` and the weight matrices. Waller (2007) derived this
for orthogonal components, but the identity holds for **any linear `W`**, which includes regression
(Thurstone, `W = R^{-1} Λ`), Bartlett, and tenBerge scores. So the standard PCA *and* EFA/ESEM
(regression-scored) paths are all algebra-eligible.

**Standardization is the trap.** Component/factor scores are not unit-variance in general (Bartlett
and oblique scores especially). Always divide by the real score SDs from `sqrt(diag(W'RW))`.
Assuming unit variance silently corrupts the correlations.

### 5.2 Scoring descriptor (returned per level by each engine)

```r
scoring <- list(
  linear    = TRUE,                # are scores a fixed linear map S = Z W?
  method    = "tenBerge",          # "components" | "regression" | "bartlett" | "tenBerge" | "EAP" | ...
  basis     = "pearson",           # "pearson" | "polychoric" | "spearman"
  weights   = W,                   # p x k score-coefficient matrix; NULL if !linear
  score_var = v                    # length-k vector = diag(W'RW); NULL if !linear
)
```

### 5.3 `compute_edges()` contract

```r
compute_edges(levels, R,
              edge_method = c("auto", "algebra", "scores"),
              pairs  = c("adjacent", "all"),   # "all" for Forbes extension
              data   = NULL,
              align  = TRUE,
              use    = "pairwise") {

  for (each pair (a, b) in chosen pairs) {

    algebra_ok <- levels[[a]]$scoring$linear &&
                  levels[[b]]$scoring$linear && !is.null(R)

    if (edge_method == "algebra" || (edge_method == "auto" && algebra_ok)) {
      Wa <- levels[[a]]$scoring$weights
      Wb <- levels[[b]]$scoring$weights
      C  <- crossprod(Wa, R %*% Wb)                       # W_a' R W_b
      sa <- sqrt(diag(crossprod(Wa, R %*% Wa)))           # score SDs
      sb <- sqrt(diag(crossprod(Wb, R %*% Wb)))
      E  <- sweep(sweep(C, 1, sa, "/"), 2, sb, "/")       # standardize
    } else {
      stopifnot(!is.null(data))                           # nonlinear/EAP or edge_method="scores"
      Sa <- score(levels[[a]], data); Sb <- score(levels[[b]], data)
      E  <- cor(Sa, Sb, use = use)
    }

    if (align) E <- align_signs(E, lineage)               # see §7
    edges[[pair]] <- E
  }
  # return: list of (k_a x k_b) matrices  +  tidy edge tibble (from, to, r, is_primary, above_cut)
}
```

**Two situations force the `scores` route:**
1. **Nonlinear scoring** — EAP/MAP/ML scores from a *categorical/ordinal* ESEM. (If the ordinal
   ESEM instead uses regression scores off polychoric loadings, it stays linear → algebra on the
   polychoric `R`.)
2. **User wants empirical, sample-realized correlations** (including FIML/missing-data realities)
   rather than the model-implied quantity. Under missingness these differ; document which one
   "the edge" represents (default: model-consistent / algebra where eligible).

### 5.4 Built-in cross-check (correctness oracle)

For a linear engine, algebra and materialized-scores must agree within sampling error. Keep the
`scores` route available even where `algebra` is the default, and add a test asserting
`max(abs(E_algebra - E_scores)) < tol` on a known dataset. This catches weight-convention,
standardization, sign, and column-ordering bugs cheaply.

## 6. Result object (S3) & storage rule

S3, not R6/S4 — matches `psych`/`lavaan`/broom conventions and serializes cleanly.

```r
structure(list(
  call, engine, rotation, cor, n_obs, k_max, seed, pkg_version,
  levels  = <list indexed by k; each is the §4 level object>,
  edges   = <list of adjacent (and optionally all-pairs) correlation matrices
             + tidy edge tibble>,
  lineage = <primary-parent matching: ordering + sign anchors; kept SEPARATE from IDs>,
  scores  = NULL,        # opt-in (large n; privacy); accessor recomputes on demand
  fits    = NULL,        # opt-in raw engine objects (keep_fits = TRUE)
  r       = <p x p input correlation matrix>,   # cheap; kept so scores/edges recomputable
  data    = NULL,        # opt-in raw data
  meta    = <suggest_k output, timestamps, convergence summary, chosen defaults,
             input_type = "data"|"cor_matrix">,
  prune   = NULL         # Forbes extension: node flags + chain table + phi table;
                         # populated by .apply_pruning() when prune != "none"
), class = "ackwards")
```

**Storage rule.** *Light tidy core always; heavy artifacts opt-in.*
- **Always:** loadings, variance, fit/convergence, weight matrices, within-level correlations,
  edges, lineage, the `p x p` input `R`, meta. This is enough for diagram + interpretation and is
  small and shareable.
- **Opt-in:** factor `scores` (O(n × Σk); also a privacy/sharing liability), raw `fits`
  (lavaan fits can be large), raw `data`. Defaults: `scores = FALSE`, `keep_fits = FALSE`,
  keep `R = TRUE`.

This directly answers the "giant list of fitted models vs. just loadings/scores/correlations"
question: neither — a structured light core, with the heavy bits nullable and recomputable.

## 7. Naming & sign conventions

- **IDs:** models numbered by factor count (`m1`…`mK`); factors within a model `m{k}f{j}`.
- **Deterministic within-level order = primary-parent (recursive), variance tiebreak.** Sort
  factors within a level by (primary parent's position in the level above, then descending variance).
  The single top factor and any ties fall back to variance. This places each child under its parent
  so splits read cleanly; the `f{j}` ID order follows this layout order so tables and plots agree.
  Limitation: ordering reduces but can't eliminate crossings (a child with a strong *secondary*
  parent elsewhere still throws a crossing edge); the layout's barycenter step (§11) refines
  x-position continuously beyond the discrete order.
- **Lineage lives in `edges`/`lineage`, never in the ID.** A child can have multiple strong
  parents (the point of the method), so don't encode parentage into `m{k}f{j}` — it breaks when the
  structure isn't a clean tree.
- **Matching:** assign each factor its primary parent via **greedy per-column argmax** on `|r|`
  (i.e. `apply(abs(E), 2, which.max)`). Multiple children sharing the same parent is normal and
  expected — adjacent levels always have n_b = n_a + 1 (level k has k factors vs k−1 above), so
  by pigeonhole at least two children share a parent. A bijection (LSAP/Hungarian) is therefore
  ill-posed and was removed. Tucker's φ serves as a supplementary congruence criterion for
  `redundancy_phi` filtering (Forbes extension, §14.19) but is not used for matching itself.
- **Sign alignment — anchor to the primary parent, NOT "all positive."** Sign is one DoF per
  factor; you cannot make every cross-level correlation positive. Rule: anchor `m1f1` to a defined
  orientation (positive sum of loadings, or a substantive marker item), then orient each factor so
  its correlation with its **primary parent** is positive, propagating top-down. Document that
  non-primary edges may legitimately be negative.

## 8. Suggesting k

`suggest_k()` returns **several complementary criteria and a consensus range**, never a single
number. The five criteria implemented (M12):

| Criterion | Source | Function | Notes |
|---|---|---|---|
| **PA-PC** | Horn (1965), PC basis | `psych::fa.parallel(fa = "both")` | Tends to overextract; use as upper bound |
| **PA-FA** | Horn (1965), FA basis | same call | More conservative; model-consistent for EFA/ESEM |
| **MAP** | Velicer (1976) | `psych::vss()` | Minimise average squared partial; usually conservative |
| **VSS-1/VSS-2** | Revelle & Rocklin (1979) | same call (already returned) | Maximise very-simple-structure fit at complexities 1 and 2 |
| **CD** | Ruscio & Roche (2012) | `EFAtools::CD()` (optional) | Resamples raw data; among top performers in simulation; skipped gracefully when `EFAtools` absent. M21: `cd_rmse = colMeans(RMSE_eigenvalues)` stored in object; dedicated 4th panel in `autoplot.suggest_k()` (2×2 grid when CD present; 3-panel single column otherwise). |

**Selective computation (M25).** A `criteria` argument (`rlang::arg_match(multiple = TRUE)`, default
all five) lets callers request any subset of `c("pa_pc", "pa_fa", "map", "vss", "cd")`. Skipping is a
real compute saving, not output filtering: `pa_pc`/`pa_fa` share one `fa.parallel(fa = "both")` call
(both or neither run) and `map`/`vss` share one `vss()` call, so e.g. `criteria = "map"` avoids the
PA simulation entirely; `vss` toggles VSS-1+VSS-2 as a unit. Non-requested `k_*` fields are `NA`; the
`criteria` data frame keeps a stable schema (NA-filled non-run columns); the object records
`criteria_requested`. `print()`/`autoplot()` render only requested criteria and the consensus range is
computed from requested criteria only.

Empirical Kaiser Criterion (EKC) and EGA (`{EGAnet}`) are **out of scope** — additional deps
without sufficient incremental benefit over the five criteria above. Note: `suggest_k()` accepts
`cor = "pearson"` (default) or `"spearman"` only; `"polychoric"` is not supported (PA and VSS do
not have a polychoric eigen-decomposition path). Users analysing ordinal data should run
`suggest_k()` on the Pearson matrix and then switch to `cor = "polychoric"` in `ackwards()`.
Report that `k` is a maximum *depth*; users often deliberately set it a level or two past the
consensus to watch factors fragment. Note the overextraction/non-replicability caution (Forbes 2023).
Add a `seed` argument for reproducibility of the **CD step only** — `psych::fa.parallel()` does not
respond reliably to `set.seed()`, so PA simulation results will vary across calls regardless of
`seed`. See `vignette("ackwards-suggest-k")` for a narrative/educational treatment of all five
criteria including pros/cons, bias direction, and engine-to-criterion pairing (M14).

## 9. Defaults (high-stakes — users will not override these)

Principle: **safe, robust, reproducible, and self-disclosing.** Every consequential auto-choice is
announced via cli and documented in roxygen with its rationale.

| Decision | Default | Rationale |
|---|---|---|
| `engine` | `"pca"` | original method; fastest; never fails to converge; algebra-exact. Docs steer to `efa`/`esem` when a measurement-model rationale exists. |
| `rotation` | **varimax — the only supported rotation; not a user argument** | **Only orthogonal rotation produces interpretable between-level factor score correlations** in the bass-ackwards method (Goldberg 2006; Kim & Eaton 2015): T'T = I so T' = T^-1, enabling closed-form W'RW algebra (Waller 2007) without materialising scores. CF(κ = 1/p) ≡ varimax (Crawford & Ferguson 1970; Browne 2001) — no reference paper varies κ. Oblique rotation is out of scope (resolved 2026-06, M13): it confounds the cross-level signal that is the method's whole point. `rotation` was removed as a user argument in M13 (varimax hardcoded internally). |
| `estimator` (ESEM only) | **`"WLSMV"`** for `cor = "polychoric"`; `"ML"` otherwise | WLSMV (mean-and-variance-adjusted WLS) is the standard limited-information ordinal estimator (matches Kim & Eaton 2015; Forbush et al. 2024); gives correct fit indices for categorical indicators without full-information ML cost. |
| `cor` (basis) | **`"pearson"`** (matches `psych`/`lavaan`); ordinal opt-in via `cor = "polychoric"` | No silent basis-switching (it can change the structure and break comparison to published work). Instead, **detect likely-ordinal columns and emit a suppressible cli warning** pointing to the polychoric option — loud *advice*, not silent action. |
| scores (method) | **`"tenBerge"`** on the active basis (pearson or polychoric) for factor engines; `"components"` for PCA; `"EAP"` opt-in only | tenBerge preserves factor correlations (the property bass-ackwards cares about) and stays linear → algebra-eligible. For ordinal ESEM, tenBerge-on-polychoric gives the clean model-implied edge; EAP's shrinkage attenuates cross-level correlations, so it's an opt-in (triggers the scores route + raw-data requirement), not the default. |
| `edge_method` | `"auto"` | algebra when linear, scores otherwise. |
| `pairs` | `"adjacent"` | classic Goldberg; `"all"` switched on with the Forbes extension. |
| `extension` (Forbes pruning) | **off** | pruning is an interpretive choice with thresholds; turning it on silently would change results. Opt-in with documented thresholds (|r| ≥ .9, congruence > .95). |
| `redundancy_phi` | **`NULL` (auto)** — PCA → no φ filter (W′RW algebra is exact; `|r|` alone is sufficient); EFA/ESEM → `0.95` (Lorenzo-Seva & ten Berge 2006; factor-score indeterminacy off-PCA makes `|r|`-only liberal; φ adds a congruence guard). Explicit number overrides on any engine. `NA` is the opt-out (no φ filter regardless of engine). Announces auto-resolve via cli (Invariant 6). **Added M25.** |
| sign `align_signs` | `TRUE` | unaligned signs make output unreadable. |
| `keep_fits` / `keep_scores` | `FALSE` / `FALSE` | memory + privacy. |
| `k_max` | required | force a deliberate choice; don't silently pick. |
| `seed` | `NULL` but captured | stochastic engines (rotation starts, ML) need reproducibility; encourage setting. |
| `missing` | **`"pairwise"`** | preserves existing behaviour (pairwise-complete correlations); warns when NAs present. `"listwise"` gives fully consistent N across fit and edges (reduces to complete cases pre-fit). `"fiml"` (ESEM ML/MLR only) uses Full Information ML and derives edge R from lavaan's FIML saturated model. FIML errors for PCA, EFA, and WLSMV/ULSMV. Added M16; **see §14 "ESEM ML/MLR pairwise" limitation entry** (the ML/MLR fit-vs-edges inconsistency under missingness is resolved for `"listwise"` and `"fiml"`; `"pairwise"` retains the existing minor inconsistency and now warns). Ignored (with a warning) when a correlation matrix is supplied — added M22. |
| `n_obs` | `NULL` (new M22) | required for `engine = "efa"` when a correlation matrix is supplied (psych needs N for chi-square/RMSEA); optional for `"pca"` (stored as `NA_integer_`). Ignored (with a warning) when raw data are supplied (N comes from `nrow(data)`). |

### Documentation standard (owner priority)

- Every default above is documented in roxygen with **why**, not just **what**.
- An `@details` section explains the algebra-vs-scores edge computation in plain language, and
  when a user would prefer `edge_method = "scores"`.
- Runnable `@examples` for each exported function; a vignette reproduces a published structure
  (e.g., a `psych::bfi` bass-ackwards) end to end.
- `@seealso` cross-links engines, `suggest_k()`, and the plotting/printing methods.

## 10. Console & tidy representations

- **`print.ackwards`** — compact cli "card": call, engine, rotation, cor basis, `n`, level
  range, per-level convergence, deepest usable level, count of edges above the show-cut. No matrix
  dumps. States the "series of linked solutions, not a fitted hierarchy" caveat once.
- **`summary.ackwards`** — per-level variance explained and fit indices (ESEM); a readable
  lineage list (`m1f1 → m2f1, m2f2 → …`); flagged redundant/artefact components when the extension
  is on.
- **`top_items(x, level, cut, n, sort)`** (M18) — per-factor salient-item listing, filtered to
  `|loading| >= cut` and sorted descending. Returns a `top_items` S3 object with a grouped cli print
  method. Loadings reflect primary-parent sign alignment (Inv. 4). The `$data` field is a subset of
  the `tidy(what = "loadings")` table.
- **broom-style tidiers:**
  - `tidy(x, what = "edges")` *(default)* — one row per directed edge `(from m_k f_i, to m_{k+1}
    f_j, r, is_primary, above_cut)`. This is the graph edge list and the primary plotting input.
  - `tidy(x, what = "loadings")` — long format `(level, factor, item, loading)`, sign-aligned;
    feeds loading heatmaps.
  - `tidy(x, what = "variance")` — per-factor variance by level.
  - `glance(x)` — one row of model-level meta (engine, k_max, n, deepest converged, cor basis).
  - `as_tibble(x)` → `tidy(x, "edges")`.

## 11. Visualization

Split **layout** (light, core) from **rendering** (optional, Suggests). This keeps the dependency
footprint sane and lets users plot the layout however they like.

- **`ba_layout(x)` (core, no heavy deps).** Returns tidy nodes `(id, level, x, y, label)` and edges.
  The structure is a **layered DAG**, not a tree (multiple parents allowed), so do **not** use a
  tree layout. Compute `y = -level`; compute `x` with a **two-pass barycenter algorithm**
  (Sugiyama-style crossing reduction):
  1. **Top-down pass** -- determines left-to-right *order* at each level. Each factor's rank is the
     |r|-weighted mean of its parents' ranks, grouping siblings that share a parent and minimising
     crossings.
  2. **Bottom-up pass** -- assigns actual x coordinates. The deepest level (level k) is evenly
     spaced; every upper-level factor is placed at the **simple mean x of its primary children**:
     a factor with one primary child lands directly above it; a factor with multiple primary children
     lands exactly halfway between them. Falls back to |r|-weighted mean of all children for factors
     with no primary children. Spreading enforces `min_sep` between adjacent nodes only when needed.
     The level-1 node is shifted to `x = 0` after both passes.
  This produces the "boxes aligned to their primary-child structure" look the owner wants.
- **`autoplot.ackwards()` / `plot()` (Suggests: ggraph/ggplot2/igraph/tidygraph).** Renders the
  layered diagram from `ba_layout()` + the edge tibble. Edge encodings: color by **sign**,
  width/alpha by **|r|**; draw **solid** above a "strong" cut and **dashed** between the show-cut
  and strong-cut (a real convention — strong between-level correlations solid, weaker primary
  associations dashed). Nodes labeled `m{k}f{j}`, optionally with substantive labels and top-loading
  items. Pyramid vs. tree-ish orientation as an option.
- **Forbes extension rendering (more complex — phase it).** Keep the same layout; **annotate**
  rather than re-layout: fade/strike pruned (redundant/artefact) nodes, and allow **skip-level**
  edges (longer curves) since the extension correlates *all* levels, not just adjacent. Default to
  showing only above-cut edges to control clutter. Treat as a later milestone; ship the clean
  adjacent-level Goldberg diagram first.
- **`label_template(x, style)` (M18, core, no heavy deps).** Generates the named-character-vector
  scaffold for `autoplot(node_labels=)`. Styles: `"id"` (round-trip identity), `"forbes"` (level-
  letter + within-level index), `"blank"` (empty strings). Returns the vector and also prints an
  editable `c(...)` literal for copy-paste. Factor IDs in canonical `ba_layout()` order.

## 12. Dependencies

| Tier | Packages | Purpose |
|---|---|---|
| Imports | `stats`, `utils`, `cli`, `rlang`, `generics`, `psych` | core, console output, guards, tidy/augment/glance generics; PCA/EFA engines + polychoric correlations (M21: moved from Suggests — engine substrate for default path) |
| Suggests — ESEM | `lavaan (>= 0.6-13)` | ESEM engine |
| Suggests — suggest_k | `EFAtools` (optional) | `EFAtools::CD()` for Comparison Data (skipped gracefully when absent); no `EGAnet`/`paran` dep |
| ~~Suggests — matching~~ | ~~`clue`~~ | ~~Hungarian assignment~~ — removed M5; greedy argmax (§7) requires no dep |
| ~~Suggests — rotations~~ | ~~`GPArotation`~~ | ~~CF-family rotations~~ — removed M21; varimax routes through `stats::varimax` and GPArotation never loaded on any supported path |
| Suggests — viz | `ggplot2`, `gt` | diagrams (uses `ggplot2` directly, not `ggraph`/`igraph`/`tidygraph`); `gt` for wide comparison tables in vignettes (M24) |
| Suggests — ESEM parallelism | `future.apply`, `future` (optional) | parallel backend for ESEM per-level fits (M26); gated by `rlang::is_installed()` with a serial `lapply` fallback; users opt in via `future::plan()` (`future` declared because the parallel test calls `plan()` directly) |
| Suggests — infra | `testthat (>= 3.0.0)`, `knitr`, `rmarkdown`, `covr` | testing, vignettes, coverage |

`methods` is **not** imported — no `methods::` usage in the package. `ggraph`, `igraph`,
`tidygraph`, `EGAnet`, `paran` are **not** in DESCRIPTION (earlier design considered them;
implementation chose leaner routes). `future`/`future.apply`/`parallel` are **not** in
Imports: M26 added `future.apply` + `future` to **Suggests** only (optional ESEM parallelism),
keeping the framework choice (sequential/multisession/multicore) entirely in the user's hands via
`future::plan()` — no `ncores` argument, no behaviour change by default. `EFAtools` (M12) and `future.apply` (M26) are
both gated behind `rlang::is_installed()` — never hard-required.

Gate every Suggests use with `rlang::check_installed()` and a helpful cli message naming the engine
that needs it. **No Rcpp dependency planned** (see §3).

## 13. Testing & infrastructure

- testthat 3e; snapshot tests of the **PCA path against `psych::bassAckward()`** on `psych::bfi`
  to prove numerical agreement on the original method.
- The §5.4 **algebra-vs-scores cross-check** as a standing test for every linear engine.
- Per-level **convergence handling** tested: a deliberately hard level must warn + skip, not error,
  and the result must still build to the deepest converged level.
- Sign-alignment and lineage-matching tested on a constructed case with a known multi-parent child.
- roxygen2 docs, pkgdown site, CI (`R CMD check`, lintr/styler), a reproduction vignette.

## 14. Decisions resolved & remaining

**Resolved (design round):**
1. Rotation → **varimax (orthogonal) across all engines, hardcoded internal constant since M13**;
   oblique rotation is **out of scope** (resolved 2026-06 — it confounds the cross-level signal that
   is the method's whole point). `rotation` removed as a user argument in M13.
2. Correlation basis → **`pearson` default, `polychoric` opt-in**, with an ordinal-detection cli
   **warning** (no silent switching).
3. Within-level order → **primary-parent (recursive), variance tiebreak**; `f{j}` IDs follow it.
4. Ordinal ESEM scoring → **tenBerge on polychoric (linear → algebra)** as default; **EAP opt-in**.
   Default Likert path does **not** need `data` at edge time.
5. Package name → **`ackwards`** (verify via `available::available("ackwards")`).

**Resolved during build (M1–M3):**
6. Ordinal-detection heuristic: **≤ 7 distinct integer values** triggers the warning.
7. CF `kappa`: **κ = 1/p** (varimax-equivalent; Crawford & Ferguson 1970). The `kappa` argument
   was accepted and stored but never wired to any engine — all engines hardcoded `cfT → "varimax"`.
   Removed entirely in M13: CF(κ=1/p) ≡ varimax; the literature never varies kappa; exposing it
   implied quartimax/equamax were reasonable alternatives (they are not for this method).
8. `cor = "spearman"` added alongside `"pearson"` as a non-polychoric rank-based option.
9. `fm` argument exposed for EFA engine (`"minres"` default, `"ml"`, `"pa"`).
10. `cut_show` argument exposed (default 0.3) to control which edges are flagged `above_cut`.

**Resolved for M4:**
11. ESEM engine → **`lavaan::efa()`** (EFA-in-SEM, not full ESEM with structural paths). Gives
    per-level fit indices + rotation-aware SEs + WLSMV ordinal estimation. Does **not** require
    full ESEM block syntax; the bass-ackwards levels are independent EFA solutions.
12. ESEM scoring → self-compute tenBerge weights from lavaan's estimated loadings `Λ` and latent
    correlation matrix `R` via the existing `.tenBerge_weights(R, Λ)`. `lavPredict()` is not used
    for the default path (it lacks tenBerge); it is the eventual hook for EAP opt-in.
13. EAP scoring → **out of scope (declined M28).** EAP's shrinkage attenuates cross-level
    correlations — the primary signal bass-ackwards measures — making it theoretically inferior to
    tenBerge for this method. An EAP request still returns `cli_abort("not yet implemented")`.
    The scores-route seam in `compute_edges()` is preserved, but implementing EAP is not planned.
14. `cor = "polychoric"` → **general basis** for all engines (not ESEM-only). PCA/EFA compute `R`
    via `psych::polychoric()` in `ackwards()` then feed it to the engine as usual. ESEM uses
    lavaan's own latent correlation matrix as `R` for edges (do not mix psych-polychoric `R` with
    lavaan-WLSMV loadings).
15. ESEM ordinal estimator → **WLSMV** default; ULSMV a documented option.
16. `loadings_se` → added to the §4 level contract (p×k matrix, NULL for PCA/EFA).

**Completed in M4 (in addition to items 11–16 above):**
- NPD polychoric matrices: warn + smooth via `psych::cor.smooth()` (implemented).
- Convergence truncation tested for the ESEM engine (k=4+ on p=6 triggers lavaan error → warn + truncate → object builds to k_eff=3).
- Algebra-vs-scores cross-check documented as Pearson/continuous paths only; polychoric ESEM edges (algebra uses lavaan polychoric R; scores route uses Pearson standardization) diverge by design and are excluded from the oracle.

**Resolved for M5 (Forbes extension):**
17. API → two orthogonal args: `pairs = c("adjacent","all")` and `prune = c("none","redundant","artefact")` (char vector; default `"none"`). `prune != "none"` auto-upgrades `pairs` to `"all"` with a loud `cli_inform()` (chains need all-levels edges to assess).
18. Prune action → **flag-only, never remove.** Adds `pruned`/`prune_reason` annotation columns to the edge/node tidy structures; the object retains all levels (preserves Invariant 5 and the algebra-vs-scores oracle). Pruning is interpretive relabeling, **not** re-estimation — say so in `print`/docs (extends the §2 honesty caveat).
19. Redundancy → faithful to Forbes (2023): score-correlation chains `|r| ≥ .9` (default, tunable), retention rule = keep the bottom node if the chain reaches level k (most-specific, best-defined), else keep the topmost node (broadest manifestation). Tucker's φ (`> .95`, Lorenzo-Seva & ten Berge 2006) computed on aligned loadings as an **optional conjunctive** criterion. φ formula `Σaᵢbᵢ / sqrt(Σaᵢ² · Σbᵢ²)`; base R, no new dependency.
20. **Additive enrichments over the paper** (default output still matches Forbes's examples): always *report* both `r` and `φ` for every redundancy candidate (report-first, flag-second — borderline cases like the paper's own `.89`/`.93` alcohol component stay visible); report endpoint `r` (direct, from all-levels edges) alongside the chain and flag where they disagree (correlation is non-transitive: a clean adjacent chain neither implies nor is implied by endpoint identity — the chain answers "perpetuates at every level," endpoint `r` answers "same construct").
21. Artefact → **never auto-flagged.** `prune = "artefact"` surfaces φ for inspection; removal is a documented researcher judgment (Forbes is explicit this introduces researcher DoF / confirmation bias; cf. Wicherts et al. 2016).

**Known limitations / deferred to future milestones:**
- `factor_cor` in the ESEM engine is not permuted by the variance-sort `ord` vector. Safe permanently: only orthogonal rotation is supported (`factor_cor = I`; permutation of I is I), and oblique rotation is out of scope (§9, §14.1). The guard comment in `engine_esem.R` documents what *would* be required if that decision were ever reversed.
- Algebra-vs-scores cross-check does not cover `cor = "polychoric"` paths (see above).
- `cor = "spearman"` + `engine = "esem"` is semantically inconsistent (lavaan fits Pearson ML on raw data while edges use Spearman R); a warning is now emitted (M10).
- ~~ESEM engine does not detect or warn on improper/Heywood solutions~~ — **resolved M10**: engine now inspects `lavaan::lavInspect(fit, "theta")` and warns when `any(diag(theta) <= 0)` (parity with EFA engine, Invariant 7).
- **ESEM ML/MLR with `missing = "pairwise"`**: lavaan uses listwise deletion for the model fit while edges are computed from a separately-computed pairwise `stats::cor()` — a minor inconsistency (fit statistics at complete-case N, edges at full pairwise N). Documented in `$meta$missing`; a per-call advisory warning fires whenever NAs are detected. Use `missing = "listwise"` or `"fiml"` to resolve. *(Added M16; §9 `missing` row cross-references this.)*
- **Forbes-extension improvements deferred past M5** (the published method has weaker spots worth strengthening later; M5 ships the faithful method + the §14.20 reporting enrichments):
  - ~~*Structural artefact signals.*~~ **Done M25 (Wave 2).** `prune = "artefact"` now populates `x$prune$structural` with per-factor `few_items` / `orphan` / `split_merge` signals (flag/report only; no auto-pruning). Args: `min_items = 3L`, `orphan_r = 0.5`. `split_merge` is `TRUE` when a factor's primary items came from multiple different primary factors at the preceding level. CLI and `print()`/`summary()` report the flagged count.
  - ~~*Factor-score-indeterminacy caveat for EFA/ESEM redundancy.*~~ **Done M25 (Wave 3).** `redundancy_phi = NULL` now auto-resolves: PCA → no φ filter; EFA/ESEM → `0.95`. `NA` is the explicit opt-out. Announces via cli when auto-applied (Invariant 6). See §9 defaults table.
  - *Selection bias in the "strongest" edge.* Plotting the max correlation across many all-levels pairs capitalizes on chance (85 → 1,320 correlations as levels grow). Add bootstrap CIs / SEs on edges (reuse the `loadings_se` infrastructure) so the strongest-edge claim is inferentially honest. **Deferred; warrants its own milestone** (perf-heavy: each resample re-runs full extraction at every level; needs `{future}` parallelisation + print/plot/vignette integration).

## 15. Milestones

The milestone history (M1–M26, completed) is **no longer duplicated here** — it lives in a single
detailed log, [`MILESTONES.md`](MILESTONES.md), with a one-line index in `CLAUDE.md`. This section
was originally a forward-looking roadmap; now that the roadmap is complete, keeping a second copy
of the log invited drift (see the git history around 2026-06-30). For what remains / is deferred,
see §14 (Decisions remaining) and `CLAUDE.md`'s "Out of scope" list.

### Key references
- Goldberg, L. R. (2006). Doing it all bass-ackwards. *J. Research in Personality*, 40(4), 347–358.
- Waller, N. (2007). A general method for computing hierarchical component structures... *JRP*, 41(4), 745–752.
- Kim, H., & Eaton, N. R. (2015). The hierarchical structure of common mental disorders. *J. Abnormal Psychology*, 124(4), 1064–1078.
- Forbush, K. T., et al. (2024). Integrating "Lumpers" vs "Splitters"... *Clinical Psychological Science*, 12(4), 625–643.
- Forbes, M. K. (2023). Improving hierarchical models of individual differences: An extension of Goldberg's bass-ackward method. *Psychological Methods*.
- Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation modeling. *Structural Equation Modeling*, 16(3), 397–438.
- Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal rotation. *Psychometrika*, 35(3), 321–332.
- Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. *Multivariate Behavioral Research*, 36(1), 111–150.
- Revelle, W. `psych::bassAckward()` documentation.
