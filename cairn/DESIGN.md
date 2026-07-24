# DESIGN.md — Bass-Ackwards Hierarchical Structural Analysis (R package)

> Working design brief. Carry this into the repo as the seed for implementation.
> It records what we've settled, the contracts to implement, and recommended defaults.
> Cross-cutting decisions are logged in [`DECISIONS.md`](DECISIONS.md); §14 points there.

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

**Audience (design interview 2026-07-17).** The primary audience is applied hierarchy
researchers — personality and psychopathology researchers running the published
Goldberg/Forbes/HiTOP workflows. Teaching and methodological extension are served but secondary;
when audiences conflict, the applied researcher running a published workflow wins.

## 2. Scope & positioning

**Prior art.** `psych::bassAckward()` (Revelle) already implements the core loop for PCA and
EFA with rotation and `tenBerge` scores. We do **not** reimplement factor extraction numerics;
we wrap established engines and add what `psych` does not.

**Our contribution.**
- A first-class **ESEM engine** (lavaan), bringing the Mplus psychopathology/HiTOP workflow
  into R — per-level fit indices and SEs, ordinal/Likert support — which `psych::bassAckward`
  does not provide. (Kim & Eaton 2015; Forbush et al. 2024 are the reference workflows.)
- Proper **ordinal/Likert handling** (polychoric basis, ordinal estimators).
- The **Forbes (2023) extended method** (redundancy + artifact pruning via the standalone
  `prune()` verb, all-levels correlations).
- A **replicability gate on hierarchy depth**: split-half factor comparability
  (Everett 1983; the split-half gate of Goldberg's own research program — Saucier 1997;
  Saucier et al. 2005, .90 threshold — dropped by the modern ESEM/HiTOP lineage; citation
  lineage verified against the PDFs 2026-07-16, see `cairn/references/`) as the standalone
  `comparability()` verb (M46; D-022), with the
  recommended end-to-end workflow documented in `vignette("ackwards-girard")`.
- A clean, tidy, serializable **result object** with `print`/`summary`/`tidy`/`glance`/`autoplot`.
- A **lineage-aligned layered diagram** rather than a misleading tree.

**Contract boundary (design interview 2026-07-17).** The package's job ends at interpretation
aids: fit, characterize, and help interpret the hierarchy (`tidy`/`autoplot`/`top_items`/
`comparability`). Publication-quality *table* exports (e.g. gt-based) are open future work;
generating manuscript *prose* is permanently out of scope. (Today gt appears only inside
vignettes; no table-export function is exported.)

**Ambition & maintenance posture (design interview 2026-07-17).** ackwards is maintained
indefinitely as the citable reference implementation of the extended bass-ackwards method
(Forbes 2023 footnote 3 already points readers here). Consequences: upstream-compatibility
vigilance is standing work, not incident response (the lavaan 0.7 episode, PR #66, is the
model), and the Forbes fidelity tests are a permanent contract. Upstream posture:
**track CRAN-current** psych/lavaan; the declared floors (`lavaan >= 0.6-13`, `R >= 4.1`) are
best-effort, droppable without ceremony when they become a burden — no old-version CI leg.

**Capability bar (design interview 2026-07-17).** A new capability earns a place by implementing
or directly serving a *published* method, with verifiable fidelity (the pattern set by the
Forbes M44/M53 oracle tests and `comparability()`'s Everett/Saucier lineage). Package-invented
conventions stay out (cf. D-020 dropping `n_obs = "effective"` for exactly this reason). Demand
can prioritize the queue, but publication lineage is the gate.

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
- **API stability: fluid until 1.0; 1.0 is paper-coupled** (design interview 2026-07-17). While
  0.x, exported behavior may break with a NEWS notice (the tracking-rules pre-1.0 waiver is
  standing). The 1.0 trigger is publication of the companion BRM methods paper: a published
  paper freezes the API its worked examples demonstrate. From 1.0, changes follow a deprecation
  cycle.

## Design principles (IP/GP)

Numbered per cairn's tracking rules: **IP** = inviolable (changing one requires a D-entry);
**GP** = guiding (tradeable with stated justification). Adopted at the 2026-07-17 design
interview (D-031), formalizing CLAUDE.md's former "Invariants" list. Numbers are never reused
or renumbered.

### Inviolable principles

- IP1: **One edge path.** All between-level correlations go through `compute_edges()`: exact
  `W'RW` algebra when scoring is linear; materialized scores only when nonlinear (EAP) or when
  the user asks. **Always** standardize by real score SDs `sqrt(diag(W'RW))` — never assume
  unit variance. (§5; D-004.)
- IP2: **Both routes, and they must agree.** The scores route stays available even where
  algebra is the default (`edge_method = "scores"`), and a standing test asserts
  algebra-vs-scores agreement within tolerance for every linear engine — the package's cheapest
  correctness oracle. (§5.4; D-004.)
- IP3: **Light core, heavy opt-in.** The object always carries loadings/variance/fit/weights/
  edges/lineage/`R`/meta; `scores`, raw `fits`, raw `data` are NULL by default and recomputable.
  Small-and-shareable by default is a privacy promise, not just a memory one. (§6; D-005.)
- IP4: **Sign alignment anchors to the primary parent**, never "all positive" (impossible —
  sign is one DoF per factor). Each child orients against its parent's *already-aligned* sign;
  secondary edges may legitimately be negative. (§7; D-010.)
- IP5: **Lineage lives in edges, never in IDs.** `m{k}f{j}` are stable labels; parentage is
  edge structure; no verb ever mutates an ID. (§7; D-009, D-029.)
- IP6: **Loud defaults.** Every consequential auto-choice announces itself via cli; the
  package advises loudly and never silently switches basis or any consequential setting.
  (§9; D-006, D-008.)
- IP7: **Convergence is data, not an error.** A non-converging level (or replicate, or split
  half) warns and is recorded/skipped; the result still builds to the deepest converged level.
  One bad level never aborts a run. (§4; D-003, D-023.)
- IP8: **Oracle-backed numerics.** Every numeric result is verified against ≥2 independent
  oracle *types* (published/closed-form, independent package, seeded simulation); no unsourced
  or unreproducible reference value ships; committed fixtures carry a structured `provenance`
  attr naming their `data-raw/` generator. Registry: `cairn/ORACLES.md`. Live independent-impl
  oracles are the stronger form — don't freeze them into fixtures without cause (M57).
- IP9: **Forbes (2023) reproducibility is a permanent *capability*, not a default lock-in.**
  The package must always be able to reproduce Forbes's published results exactly — test-backed
  (M44 sims + M53 AMH, 54/54 chase components) with the reproducing settings available and
  documented. Defaults, however, are free to adopt a better method when one arrives, with loud
  documentation (IP6) and a D-entry. (Adopted 2026-07-17; supersedes the earlier "default
  output must reproduce Forbes's examples exactly" contract wording.)

### Guiding principles

- GP1: **Published-method capability bar.** A new methodological capability earns its place by
  implementing or directly serving a *published* method with verifiable fidelity;
  package-invented conventions stay out (the D-020 precedent). Sound-engineering utilities are
  the tradeable exception, with stated justification.
- GP2: **Report-first, flag-second.** Computed flags (`above_cut`, prune flags) never appear
  without the underlying values beside them, and contested cutoffs never become returned
  verdicts — thresholds render as reference lines/annotations, not pass/fail columns.
  (D-014, D-017, D-028, D-030.)
- GP3: **Descriptive honesty.** A bass-ackwards result is presented as a series of linked
  solutions whose edges are score correlations — never as a fitted hierarchical model — and
  docs + `print()` keep saying so. (D-001.)
- GP4: **Wrap, don't reimplement.** Prefer wrapping established engines (psych, lavaan,
  stats) over reimplementing numerics; self-implement only where no wrapper exists (the D-016
  tenBerge-weights precedent shows the justified trade).
- GP5: **Lean install.** `psych` is the only heavy Imports; everything else heavy sits in
  Suggests behind `rlang` guards; no Rcpp. New Imports go through the dependency gate and a
  D-entry. (§3, §12; D-011.)
- GP6: **Reproducible by construction.** Stochastic results carry their seeds in the object;
  parallel and serial execution agree bit-for-bit (D-023); where an upstream defeats seeding
  (`psych::fa.parallel`), the exception is documented loudly, never hidden. (§8.)

### Lineage map (pre-interview numbering → principles)

CLAUDE.md's "Invariants" 1–8 — cited as `Invariant N` in ~22 `R/` and `tests/` files — map
**identically**: Invariant N → IPn for N = 1..8. No in-code repoint is required; existing
comments keep resolving. IP9 and all GPs are new numbers with no legacy citations.

## 4. Engines

Three engines, all first-class, behind one user-facing function (`ackwards()`) that dispatches.
Rotation is a single fixed choice across all engines: **varimax** (orthogonal; keeps within-level
factors uncorrelated so between-level edges are not confounded — the W'RW algebra itself is exact
for any linear `W`, see §5.1/§9). Oblique rotation is out of scope — see §9 and D-002.

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
             input_type = "data"|"cor_matrix",
             item_labels = <named chr of column "label" attributes, or NULL> (M36),
             item_means/item_sds = <fit-time item moments from the
             post-missing-handling data; NULL for cor_matrix input> (M45 --
             the frame of reference for out-of-sample scoring),
             min_eigenvalue = <smallest eigenvalue of R; NA where not computed>,
             near_singular = <TRUE when min_eigenvalue < 1e-4 -- a durable
             rank-deficiency signal re-surfaced by print()/summary() (M49)>>,
  prune   = NULL         # Forbes extension: node flags + chain table + phi table;
                         # populated by prune(x, ...) -- a standalone verb piped off
                         # ackwards(), not an ackwards() argument (M34; see s.14)
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
  `redundancy_phi` filtering (Forbes extension, D-017) but is not used for matching itself.
- **Sign alignment — anchor to the primary parent, NOT "all positive."** Sign is one DoF per
  factor; you cannot make every cross-level correlation positive. Rule: anchor `m1f1` to a defined
  orientation (positive sum of loadings, or a substantive marker item), then orient each factor so
  its correlation with its **primary parent** is positive, propagating top-down. Document that
  non-primary edges may legitimately be negative. *Implementation note (M35):* "propagating
  top-down" means each child's flip is chosen against the parent's **already-aligned** sign, not the
  parent's raw orientation — otherwise a flipped parent leaves its own primary edge displaying
  negative. So every primary-parent edge is non-negative; only secondary edges may be red.

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

`suggest_k()` answers "what depth range is *plausible*" from the eigenstructure. The companion
question — "which factors actually *replicate*" — is answered directly by the standalone
`comparability()` verb (M46, D-022): split-half factor comparability per level per factor, the
depth **floor** to `suggest_k()`'s ceiling. The two are deliberately separate functions:
`suggest_k()` is a pre-fit screen over eigenvalue criteria, `comparability()` is engine-specific
and fit-based.

## 9. Defaults (high-stakes — users will not override these)

Principle: **safe, robust, reproducible, and self-disclosing.** Every consequential auto-choice is
announced via cli and documented in roxygen with its rationale.

| Decision | Default | Rationale |
|---|---|---|
| `engine` | `"pca"` | original method; fastest; never fails to converge; algebra-exact. Docs steer to `efa`/`esem` when a measurement-model rationale exists. |
| `rotation` | **varimax — the only supported rotation; not a user argument** | **Only orthogonal rotation produces interpretable between-level factor score correlations** in the bass-ackwards method (Goldberg 2006; Kim & Eaton 2015): varimax makes the within-level factors orthogonal, so a between-level edge reflects only the cross-level relationship instead of being confounded by within-level factor intercorrelation. The varimax criterion itself originates with Kaiser (1958); CF(κ = 1/p) ≡ varimax (Crawford & Ferguson 1970; Browne 2001) — no reference paper varies κ. Oblique rotation is out of scope (resolved 2026-06, M13): it confounds the cross-level signal that is the method's whole point. `rotation` was removed as a user argument in M13 (varimax hardcoded internally). *(Rationale wording corrected M76 per RR01: the earlier "T'T = I … enabling closed-form W'RW algebra" phrasing conflated the interpretive Φ = I reason above with algebra-exactness — the W'RW identity is exact for **any fixed linear scoring**, oblique included, and Waller 2007 §3 gives the oblique closed form. Orthogonality is the interpretive choice, not a numerical prerequisite. See §5.1.)* |
| `estimator` (ESEM only) | **`"WLSMV"`** for `cor = "polychoric"`; `"ML"` otherwise | WLSMV (mean-and-variance-adjusted WLS) is the standard limited-information ordinal estimator (matches Kim & Eaton 2015; Forbush et al. 2024 use the ULSMV variant); gives correct fit indices for categorical indicators without full-information ML cost. |
| `cor` (basis) | **`"pearson"`** (matches `psych`/`lavaan`); ordinal opt-in via `cor = "polychoric"` | No silent basis-switching (it can change the structure and break comparison to published work). Instead, **detect likely-ordinal columns and emit a suppressible cli warning** pointing to the polychoric option — loud *advice*, not silent action. |
| scores (method) | **`"tenBerge"`** on the active basis (pearson or polychoric) for factor engines; `"components"` for PCA; `"EAP"` opt-in only | tenBerge preserves factor correlations (the property bass-ackwards cares about) and stays linear → algebra-eligible; the correlation-preservation-vs-determinacy trade-off this choice accepts is formalized by Grice (2001) and Beauducel, Hilger, & Kuhl (2024), and the factor-score hierarchy's criterion validity — with a categorical-indicator caveat — by Williams et al. (2025). For ordinal ESEM, tenBerge-on-polychoric gives the clean model-implied edge; EAP's shrinkage attenuates cross-level correlations, so it's an opt-in (triggers the scores route + raw-data requirement), not the default. |
| `edge_method` | `"auto"` | algebra when linear, scores otherwise. |
| `pairs` (`ackwards()`) | `"adjacent"` | classic Goldberg; `"all"` reveals skip-level correlations for inspection/plotting. Since M34, `prune()` recomputes its own all-pairs edges on demand regardless of this setting — pruning no longer requires (or auto-upgrades) `pairs = "all"` here. |
| `prune()` `rules` (standalone verb, M34 — not an `ackwards()` argument) | **`"none"`** | pruning is an interpretive choice with thresholds, kept as a separate, cheap, re-runnable step piped off an already-extracted object (`ackwards(...) |> prune(...)`) so new thresholds never require re-extraction. Turning it on silently would change results — opt-in with documented thresholds (|r| ≥ .9, congruence > .95). Canonical rule name is `"artifact"` (US spelling); `"artefact"` is accepted as an alias (nod to Commonwealth spelling and to Forbes). |
| `redundancy_phi` (`prune()` argument, M34) | **`NULL` (auto)** — PCA → no φ filter (component scores are **determinate** — exact linear functions of the data — so `|r|` is the true correlation between the components themselves and suffices alone); EFA/ESEM → `0.95` (Lorenzo-Seva & ten Berge 2006; factor-score indeterminacy (Grice 2001) makes `|r|`-only liberal; φ adds a congruence guard). Explicit number overrides on any engine. `NA` is the opt-out (no φ filter regardless of engine). Announces auto-resolve via cli (IP6). **Added M25**; moved from an `ackwards()` argument to a `prune()` argument in M34. *(Rationale wording corrected M43: the earlier "W′RW algebra is exact" phrasing conflated algebra-exactness — equally true of tenBerge EFA — with score determinacy, the actual reason.)* |
| sign `align_signs` | `TRUE` | unaligned signs make output unreadable. |
| `keep_fits` / `keep_scores` | `FALSE` / `FALSE` | memory + privacy. |
| `k_max` | required | force a deliberate choice; don't silently pick — bass-ackward stopping criteria are themselves an open question with only simulation-level guidance (Tong, Qu, & Zhang 2025). |
| `seed` | `NULL` but captured | stochastic engines (rotation starts, ML) need reproducibility; encourage setting. |
| `missing` | **`"pairwise"`** | preserves existing behaviour (pairwise-complete correlations); warns when NAs present. `"listwise"` gives fully consistent N across fit and edges (reduces to complete cases pre-fit). `"fiml"` uses Full Information ML: for ESEM (ML/MLR) it derives edge R from lavaan's FIML saturated model; for **PCA/EFA on the Pearson basis** (M38) it estimates R via `psych::corFiml()` and feeds it to the `W'RW` algebra (IP1 — one edge path; no new dep). FIML **errors** for WLSMV/ULSMV and for a **non-Pearson PCA/EFA basis** (corFiml is MVN-only). Added M16; PCA/EFA route added M38; **see the "Known limitations" section's ESEM ML/MLR pairwise entry** (the ML/MLR fit-vs-edges inconsistency under missingness is resolved for `"listwise"` and `"fiml"`; `"pairwise"` retains the existing minor inconsistency and now warns). Ignored (with a warning) when a correlation matrix is supplied — added M22. |
| `n_obs` | `NULL` (M22); **string `"total"`/`"complete"` on the raw-data FIML PCA/EFA path (M38)** | Correlation-matrix input: a positive integer, required for `engine = "efa"` (psych needs N for chi-square/RMSEA), optional for `"pca"` (stored as `NA_integer_`). Raw data: N normally comes from `nrow(data)` and a numeric `n_obs` is ignored (warning). Under `missing = "fiml"` (PCA/EFA), `n_obs` may be `"total"` (**default** — all rows contributing to the FIML likelihood, matching the FIML convention; Enders 2010) or `"complete"` (complete-case N, conservative lower bound). Point estimates are unaffected by the choice; only the EFA fit indices, which are *approximate* under this two-step (FIML matrix → normal-theory EFA) route regardless of N (Zhang & Savalei 2020). `"effective"` was considered and dropped — no canonical formula, so it would be a package-invented convention. A string `n_obs` is rejected off this path. |

Where any default above **departs** from Goldberg (2006) or Forbes (2023) — or
deliberately matches them — the divergence, its rationale, and its empirical/
mathematical support are catalogued in `cairn/references/source-departures.md`
(kept current per that page's Maintenance clause: a new/changed departure carries
a D-entry per IP9 and a ledger row in the same change).

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
  lineage list (`m1f1 → m2f1, m2f2 → …`); flagged redundant/artifact components when the object
  has been pruned via `prune()`.
- **`top_items(x, level, cut, n, sort, by, show_labels)`** (M18; `by`/`show_labels` M36) —
  salient-item listing, filtered to `|loading| >= cut` and sorted descending. Returns a `top_items`
  S3 object with a grouped cli print method. `by = "factor"` (default) groups items under each
  factor; `by = "item"` inverts to list, per item, the factors it loads on (cross-loading view);
  `n`/`sort` apply within the chosen unit. When the fit data carried per-column `"label"` attributes
  (captured into `meta$item_labels`), `show_labels = TRUE` (default) prints items as `code: label`
  with a per-item bare-code fallback. Loadings reflect primary-parent sign alignment (IP4). The
  `$data` field is a subset of the `tidy(what = "loadings")` table (plus a `label` column when
  labels are available).
- **broom-style tidiers:**
  - `tidy(x, what = "edges")` *(default)* — one row per directed edge `(from m_k f_i, to m_{k+1}
    f_j, r, is_primary, above_cut)`. This is the graph edge list and the primary plotting input.
  - `tidy(x, what = "loadings")` — long format `(level, factor, item, loading)`, sign-aligned;
    feeds loading heatmaps.
  - `tidy(x, what = "variance")` — per-factor variance by level: `(level, factor, proportion,
    cumulative)`, both on a 0-1 scale (M32; matches the engine's internal representation and
    broom/psych convention — percent formatting is a display concern for `print`/`summary`, not
    a `tidy()` concern).
  - `tidy(x, what = "fit")` — per-level fit statistic: `(level, statistic, value)` (M32: column
    renamed from `index`, which read like a row position rather than a fit-statistic name — it
    holds fit-index names for EFA/ESEM and eigenvalue positions for PCA). `format = "wide"` pivots
    to one row per level. No pass/fail `cutoffs=`/`meets` output (removed M32) — conventional Hu &
    Bentler (1999) thresholds are conventional and contested, so `tidy()` reports values only;
    `autoplot(what = "fit")` and `summary()` render them as reference lines/inline annotations.
  - `glance(x)` — one row of model-level meta (engine, k_max, n, deepest converged, cor basis).
  - `as_tibble(x)` → `tidy(x, "edges")`.
- **`augment(x, data, append, id_cols, scaling)` / `predict(x, newdata, scaling)` (M45)** —
  per-observation factor scores, including **out-of-sample scoring** (fit on a training split,
  score a test split without retraining). `scaling = "fit"` (default) standardizes the supplied
  data by the **fit-time item moments** stored in `meta` so all scores share the training metric;
  `"sample"` re-standardizes by the supplied data's own moments (the pre-M45 behaviour; the only
  option for cor-matrix objects, which carry no moments). `predict()` is a thin exported wrapper
  returning exactly `augment(x, data = newdata, append = FALSE, scaling = scaling)` — the
  idiomatic front door (`psych::predict.psych` / `lavPredict` precedent).

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
- **`autoplot.ackwards()` / `plot()` (Suggests: ggplot2).** Renders the
  layered diagram from `ba_layout()` + the edge tibble. Edge encodings are **user-assignable and
  always legended** (M35): `sign_by` chooses the channel for **sign** (color by default; also
  linetype, both, or none) and `magnitude_by` chooses the channel for **|r|** (linewidth by
  default, or none). No aesthetic is mapped without a matching legend. (The pre-M35 strong/weak
  `cut_strong` linetype split was retired — it double-encoded magnitude and was legend-less.)
  Nodes labeled `m{k}f{j}`, optionally with substantive labels and top-loading items. `direction`
  toggles vertical (default, level 1 at top) vs. horizontal (level 1 at left) orientation.
- **Forbes extension rendering (more complex — phase it).** Keep the same layout; **annotate**
  rather than re-layout: fade/strike pruned (redundant/artifact) nodes, and allow **skip-level**
  edges (longer curves) since the extension correlates *all* levels, not just adjacent. Default to
  showing only above-cut edges to control clutter. Treat as a later milestone; ship the clean
  adjacent-level Goldberg diagram first. **A fully-pruned level** (every factor flagged) has its
  left/bottom axis label rendered in *italic* (M40) in the normal render path, denoting its status
  alongside the grey node fill; partially-pruned levels keep a plain label, and under
  `drop_pruned = TRUE` a fully-pruned level's nodes are removed so no label is drawn.
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
- **Oracle discipline (IP8).** Every numeric result is pinned by ≥2 independent
  oracle *types*; no unsourced/unreproducible reference value ships. The full registry — every
  oracle classified as frozen / live / invariant / closed-form, with its asserting test and
  provenance — is `cairn/ORACLES.md`; committed fixtures carry a structured `provenance` attr
  (guarded by `test-oracle-provenance.R`).
- Per-level **convergence handling** tested: a deliberately hard level must warn + skip, not error,
  and the result must still build to the deepest converged level.
- Sign-alignment and lineage-matching tested on a constructed case with a known multi-parent child.
- roxygen2 docs, pkgdown site, CI (`R CMD check`, lintr/styler), a reproduction vignette.

## 14. Decisions resolved & remaining

The decision log formerly embedded here (46 numbered items, M1–M53) was extracted in M64
(2026-07-16). **Still-governing, cross-cutting decisions live in
[`DECISIONS.md`](DECISIONS.md)** (D-001 onward, each citing its source anchor); the
pre-extraction log is entombed verbatim as
[`legacy/DESIGN-s14-decision-log.md`](legacy/DESIGN-s14-decision-log.md), against which every
historical `§14.x` citation resolves. Live known limitations moved to the next section.

## Known limitations

- `factor_cor` in the ESEM engine is not permuted by the variance-sort `ord` vector. Safe
  permanently: only orthogonal rotation is supported (`factor_cor = I`; permutation of I is I), and
  oblique rotation is out of scope (§9, D-002). The guard comment in `engine_esem.R` documents what
  *would* be required if that decision were ever reversed.
- Algebra-vs-scores cross-check does not cover `cor = "polychoric"` paths (§5.4), nor the
  `missing = "fiml"` PCA/EFA path (D-020): there the algebra uses the `psych::corFiml()` matrix
  while the scores route standardizes the raw, NA-bearing data (pairwise Pearson SDs), so the two
  bases diverge under missingness by design — the same reason polychoric is excluded. The oracle
  tests therefore run only on complete-data linear engines; no FIML/polychoric object is fed to
  them, so there is no false-failure risk, but the cross-check does not *certify* those paths.
- `cor = "spearman"` + `engine = "esem"` is semantically inconsistent (lavaan fits Pearson ML on
  raw data while edges use Spearman R); a warning is emitted (M10).
- **ESEM convergence at depth on real ordinal data** is flakier than the calm warn-and-skip
  framing suggests (owner-reported, design interview 2026-07-17): the M49 polychoric robustness
  work (sparse cross-cells, `correct = 0.5`, `check_items()`) exists because real clinical data
  is messier than the bundled teaching sets. IP7's warn+skip is a load-bearing safety
  net, not a rarity.
- **ESEM ML/MLR with `missing = "pairwise"`**: lavaan uses listwise deletion for the model fit
  while edges are computed from a separately-computed pairwise `stats::cor()` — a minor
  inconsistency (fit statistics at complete-case N, edges at full pairwise N). Documented in
  `$meta$missing`; a per-call advisory warning fires whenever NAs are detected. Use
  `missing = "listwise"` or `"fiml"` to resolve. *(Added M16; §9 `missing` row cross-references
  this.)*

## 15. Milestones

The milestone history is **not duplicated here** — status lives in
[`ROADMAP.md`](ROADMAP.md), finished milestones in
[`milestones/archive/`](milestones/archive/) plus git, and the pre-migration
history (M1–M53) is entombed verbatim in
[`legacy/MILESTONES.md`](legacy/MILESTONES.md). This section was originally a
forward-looking roadmap; now that the roadmap is complete, keeping a second copy
of the log invited drift (see the git history around 2026-06-30). For what remains / is deferred,
see `DECISIONS.md`, the "Known limitations" section, and `CLAUDE.md`'s "Out of scope" list.
*(Corrected 2026-07-19: this section previously pointed at a repo-root
`MILESTONES.md` and a CLAUDE.md milestone index, both of which the cairn
migration removed.)*

### Key references
- Everett, J. E. (1983). Factor comparability as a means of determining the number of factors and their rotation. *Multivariate Behavioral Research*, 18(2), 197–218.
- Saucier, G. (1997). Effects of variable selection on the factor structure of person descriptors. *J. Personality and Social Psychology*, 73(6), 1296–1312.
- Saucier, G., Georgiades, S., Tsaousis, I., & Goldberg, L. R. (2005). The factor structure of Greek personality adjectives. *J. Personality and Social Psychology*, 88(5), 856–875.
- Goldberg, L. R. (2006). Doing it all bass-ackwards. *J. Research in Personality*, 40(4), 347–358.
- Waller, N. (2007). A general method for computing hierarchical component structures... *JRP*, 41(4), 745–752.
- Kim, H., & Eaton, N. R. (2015). The hierarchical structure of common mental disorders. *J. Abnormal Psychology*, 124(4), 1064–1078.
- Forbush, K. T., et al. (2024). Integrating "Lumpers" vs "Splitters"... *Clinical Psychological Science*, 12(4), 625–643.
- Forbes, M. K. (2023). Improving hierarchical models of individual differences: An extension of Goldberg's bass-ackward method. *Psychological Methods*.
- Asparouhov, T., & Muthén, B. (2009). Exploratory structural equation modeling. *Structural Equation Modeling*, 16(3), 397–438.
- Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal rotation. *Psychometrika*, 35(3), 321–332.
- Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. *Multivariate Behavioral Research*, 36(1), 111–150.
- Revelle, W. `psych::bassAckward()` documentation.
