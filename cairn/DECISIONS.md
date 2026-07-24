# Decisions

Append-only. Never renumber; supersede with a new entry. D-entries record
choices with rationale — never deferrals.

Pre-migration decision history lives in `cairn/legacy/` — since M64
(2026-07-16) the decision log formerly embedded at DESIGN.md §14 is entombed
verbatim as `cairn/legacy/DESIGN-s14-decision-log.md`, and every historical
`§14.x` anchor in this file (D-001–D-015 source lines) and in the milestone
archives resolves against that copy. This file re-records the still-governing,
cross-cutting decisions, each citing its source anchor (the legacy §14 log
and/or DESIGN §2–§12); new decisions continue from the next D-number.

### D-001 (date: see legacy §2): A result is a series of linked solutions, not a fitted hierarchy

**Context:** A bass-ackwards result is a set of linked per-level solutions whose edges are score correlations, not a fitted higher-order SEM.
**Decision:** Treat the "hierarchy" as descriptive/emergent and state this honesty caveat in docs and `print()`.
**Consequences:** `print()` states it once; `prune()` is interpretive relabeling, not re-estimation; frames every claim the package makes.
_Source: DESIGN §2._

### D-002 (2026-06): Varimax is the only rotation; oblique out of scope

**Context:** Only orthogonal rotation (`T' = T^-1`) yields interpretable between-level score correlations and the closed-form `W'RW` algebra; oblique confounds the cross-level signal that is the method's point. *(Wording corrected M76 per RR01: the "and the closed-form `W'RW` algebra" clause conflates — the `W'RW` identity is exact for any fixed linear scoring, oblique included; orthogonality's role is the interpretive `Φ = I` (uncorrelated within-level factors), not enabling the algebra. Decision unchanged.)*
**Decision:** Varimax hardcoded as an internal constant (not a user argument); oblique out of scope; the `kappa` arg removed since CF(κ=1/p) ≡ varimax.
**Consequences:** No `rotation` argument; algebra path valid for all engines; `factor_cor = I` everywhere.
_Source: DESIGN §14.1; M13._

### D-003 (date: see legacy §14.11): Three engines behind one dispatch; non-convergence never fatal

**Context:** Multiple extraction methods must share one user surface, and deep levels can fail to converge (especially ESEM/lavaan).
**Decision:** Engines `"pca"/"efa"/"esem"` dispatched by `ackwards()`, each returning the standardized per-level object (§4 contract); a non-converged level is recorded, warned, and skipped — never thrown (Invariant 7).
**Consequences:** Hierarchy builds to the deepest converged level; ESEM is `lavaan::efa()` (EFA-in-SEM with per-level fit + SEs), not full ESEM.
_Source: DESIGN §4; §14.11; M4._

### D-004 (date: see legacy §5): One shared edge path — algebra when linear, scores otherwise

**Context:** Edge correlations must be consistent across engines and scoring schemes without forcing expensive score materialization everywhere.
**Decision:** A single `compute_edges()` path uses exact `W'RW` algebra when scoring is linear (`edge_method = "auto"`) and falls back to materialized scores otherwise; always divide by real score SDs from `sqrt(diag(W'RW))`.
**Consequences:** One edge path (Invariant 1) reused by `prune()`/`comparability()`/`boot_edges()`; algebra-vs-scores cross-check serves as a standing correctness oracle for linear engines.
_Source: DESIGN §5; §9 (edge_method)._

### D-005 (date: see legacy §6): S3 object with a light core; heavy artifacts opt-in

**Context:** Answers the "giant list of fitted models vs. just loadings/scores" storage question.
**Decision:** S3 object always retains a light tidy core (loadings, variance, fit/convergence, weight matrices, within-level cor, edges, lineage, `p×p` R, meta); heavy artifacts (scores, raw fits, raw data) are opt-in, nullable, and recomputable.
**Consequences:** Small and shareable by default; keeping `R = TRUE` lets scores/edges recompute (Invariant 3); `scores`/`keep_fits` default `FALSE`.
_Source: DESIGN §6._

### D-006 (date: see legacy §14.2): Pearson default, polychoric opt-in, no silent basis switching

**Context:** Silently switching correlation basis can change the recovered structure and break comparison to published work.
**Decision:** `cor = "pearson"` default, `"polychoric"` (and `"spearman"`) opt-in; never switch silently — detect likely-ordinal columns and emit a suppressible cli advisory. Polychoric is a general basis for all engines.
**Consequences:** Ordinal users opt in explicitly; PCA/EFA compute polychoric R via `psych`, ESEM uses lavaan's own latent R (never mixed).
_Source: DESIGN §14.2; §14.14; M4._

### D-007 (date: see legacy §14.4): Default scoring is tenBerge / components; EAP out of scope

**Context:** The scoring map must preserve factor correlations (the measured signal) and stay linear for the algebra path.
**Decision:** Default scores = tenBerge on the active basis for factor engines, `"components"` for PCA; EAP is out of scope (declined M28 — its shrinkage attenuates cross-level correlations).
**Consequences:** Default paths stay algebra-eligible and need no raw data at edge time; an EAP request errors while the scores-route seam is preserved.
_Source: DESIGN §14.4; §14.13; M28._

### D-008 (date: see legacy §9): Safe, reproducible, self-disclosing defaults

**Context:** Users will not override high-stakes defaults and often never read an argument.
**Decision:** Defaults are safe/robust/reproducible, and every consequential auto-choice (basis, k, `redundancy_phi`, estimator) announces itself via cli and is documented with its rationale (Invariant 6).
**Consequences:** No silent consequential behavior; the package chooses loud advice over silent action across `ackwards()`/`suggest_k()`/`prune()`.
_Source: DESIGN §9; §3._

### D-009 (date: see legacy §7): Factor IDs, within-level order, and greedy-argmax matching

**Context:** A child can have several strong parents, so parentage is not a clean tree; adjacent levels always add exactly one factor.
**Decision:** IDs `m{k}f{j}` numbered by factor count; within-level order = primary-parent (recursive) with variance tiebreak; lineage lives in `edges`/`lineage`, never in the ID; primary parent assigned by greedy per-column argmax on `|r|` (a bijection/Hungarian match is ill-posed and was removed).
**Consequences:** No `clue` dependency; tables and plots share one order; multiple children per parent is expected and normal.
_Source: DESIGN §7; M5._

### D-010 (date: see legacy §7): Sign alignment anchors to the primary parent, not "all positive"

**Context:** Sign is one degree of freedom per factor; not every cross-level correlation can be made positive.
**Decision:** Anchor `m1f1` to a defined orientation, then orient each factor so its correlation with its already-aligned primary parent is positive, propagating top-down.
**Consequences:** Every primary-parent edge is non-negative; only secondary edges may be negative/red; all tidy loadings and edges reflect this alignment.
_Source: DESIGN §7; M35._

### D-011 (date: see legacy §12): Manageable dependencies — psych in Imports, everything else guarded Suggests

**Context:** Installing the package must not pull in the SEM or plotting stacks (owner priority).
**Decision:** `psych` in Imports as the default-path engine substrate (PCA/EFA + polychoric); all other heavy/engine-specific packages in Suggests behind `rlang::check_installed()` guards; no Rcpp. `GPArotation` and `clue` removed.
**Consequences:** Lean install; optional features (ESEM, CD, plotting, parallelism) skip gracefully when their Suggest is absent.
_Source: DESIGN §3; §12; M21._

### D-012 (date: see legacy §14.27): Pruning is a standalone, flag-only `prune()` verb

**Context:** Extraction is the expensive, deterministic step; pruning is cheap, interpretive, and threshold-dependent.
**Decision:** Pruning is a standalone pipeable S3 verb `prune()` (not an `ackwards()` argument), flag-only — it annotates nodes/edges and never removes a level (preserving Invariant 5 and the oracle); default `rules = "none"`, artifact never auto-flagged.
**Consequences:** Researchers can re-prune with new thresholds without re-extracting; object class is unchanged so all methods work; `pairs` no longer auto-upgrades.
_Source: DESIGN §14.27; §14.18; M34._

### D-013 (date: see legacy §8): `suggest_k()` returns several criteria and a consensus range, never one number

**Context:** A single "number of factors" is misleading; `k` is a maximum depth users often push past deliberately.
**Decision:** `suggest_k()` reports several complementary criteria (PA-PC, PA-FA, MAP, VSS, CD) plus a consensus range, never a single number; pearson/spearman only (no polychoric eigen path).
**Consequences:** Deliberately separate from the fit-based `comparability()` replicability floor; ordinal users screen on Pearson, then switch basis in `ackwards()`.
_Source: DESIGN §8; M12._

### D-014 (date: see legacy §14.24): Report fit values, never a returned pass/fail

**Context:** Hu & Bentler (1999) fit cutoffs are conventional and contested.
**Decision:** `tidy(what = "fit")` reports values only — no `meets`/pass-fail column; thresholds appear only as visual reference lines/annotations in `autoplot()`/`summary()`. The same stance governs `factorability()`'s conventional bands.
**Consequences:** The package never returns a computed pass/fail a user could mistake for a verdict.
_Source: DESIGN §14.24; M32._

### D-015 (date: see legacy §11): Split light layout from optional rendering; layered DAG, not a tree

**Context:** Diagram rendering needs ggplot2, but the layout must stay dependency-light, and the structure is a layered DAG (multiple parents), not a tree.
**Decision:** Split `ba_layout()` (core, no heavy deps, Sugiyama-style two-pass barycenter x-placement) from rendering (`autoplot()`/`plot()`, Suggests: ggplot2); never use a tree layout.
**Consequences:** Users can plot the layout however they like; edge encodings are user-assignable and always legended.
_Source: DESIGN §11; M35._

### D-016 (date: see legacy §14.12): ESEM tenBerge weights are self-computed, not lavPredict

**Context:** ESEM edge computation needs linear tenBerge weights, but `lavPredict()` offers no tenBerge method.
**Decision:** Compute tenBerge weights directly from lavaan's estimated loadings `Λ` and latent correlation `R` via the shared `.tenBerge_weights(R, Λ)`; `lavPredict()` is unused on the default path.
**Consequences:** The one-edge-path `W'RW` algebra (Invariant 1) applies unchanged to ESEM; scoring is consistent across all three engines.
_Source: legacy §14.12; M4._

### D-017 (M53; original rules see legacy §14.19–21): Redundancy chased by direct skip-level correlation; report-first

**Context:** Forbes (2023)'s `ChaseCorrPaths()` chases redundancy via the direct correlation to each ancestor level; an adjacent-hop walk diverges under non-transitivity (7/54 AMH components).
**Decision:** `prune("redundant")` default `redundancy_criterion = "direct"` (chase contiguously upward while `|r| ≥ redundancy_r`, default .9); `"adjacent"` retained opt-in; retention keeps the bottom node if the chain reaches level k, else the topmost; Tucker's φ (> .95) is an optional conjunctive filter; always report both `r` and `φ` plus endpoint `r` for every candidate (report-first, flag-second).
**Consequences:** Reproduces Forbes's published AMH chase exactly (54/54, test-backed); borderline cases stay visible; supersedes the pre-M53 adjacent-chain assumption.
_Source: legacy §14.19–21; M5, M53._

### D-018 (date: see legacy §14.22–26): M32 output-surface naming and scale

**Context:** Several tidy-output surfaces carried misleading names/scales (a `index` key that read as row position; percent-scaled variance; an unrecorded effective estimator).
**Decision:** `tidy(what = "fit")` long-format key column is `statistic`; `k_max` keeps its shared name in `ackwards()`/`suggest_k()` (same dial, different stage — disambiguated in roxygen); variance is proportion 0–1 (`proportion`/`cumulative`), percent formatting lives in the display layer; the effective ESEM estimator is recorded in `$meta$estimator` with a `summary()` scaled-variant footnote.
**Consequences:** Breaking pre-CRAN with no deprecation path; aligns with broom/psych convention.
_Source: legacy §14.22–26; M32._

### D-019 (date: see legacy §14.28–31): M34 pruning mechanics

**Context:** Extracting `prune()` as a verb (D-012) left mechanics to fix: edge basis, manual flags, rule naming, and the removed `ackwards()` args.
**Decision:** `prune()` recomputes all-pairs edges fresh via `compute_edges()` from stored `levels`/`r`, never writing back to `x$edges`; `manual =` flags user-named nodes (auto rule's `prune_reason` wins on overlap; unknown labels error); canonical rule name `"artifact"` with `"artefact"` alias (`"tucker"` rejected — the mode is more than φ); the five removed args are loudly rejected if passed via `...`.
**Consequences:** `pairs` is a pure display setting; pruning works regardless of fit-time `pairs`; no masked-argument footgun (Invariant 6).
_Source: legacy §14.28–31; M34._

### D-020 (date: see legacy §14.32–33): FIML as a first-class PCA/EFA route; n_obs "total" default

**Context:** M16's "FIML errors for PCA/EFA" was reversible once `psych::corFiml()` could feed the existing algebra; a user could already smuggle it in via the correlation-matrix seam.
**Decision:** `missing = "fiml"` + `engine = "pca"/"efa"` + `cor = "pearson"` routes R through `psych::corFiml()` into `W'RW` (still errors for non-Pearson bases and WLSMV/ULSMV); on this path `n_obs` accepts `"total"` (default, FIML convention — Enders 2010) or `"complete"`; `"effective"` dropped (no canonical formula — it would be a package-invented convention).
**Consequences:** One corFiml call per run (Invariant 1 clean); EFA fit indices are approximate under the two-step route regardless of N (Zhang & Savalei 2020) — announced via cli.
_Source: legacy §14.32–33; M38._

### D-021 (2026-07-01): Out-of-sample scoring standardizes by fit-time moments

**Context:** A test observation's score must not depend on which other observations share its split; train/test scores need one metric (cross-validation use case).
**Decision:** `augment()`/`predict.ackwards()` default `scaling = "fit"` (training means/SDs stored in `meta$item_means`/`item_sds`); `"sample"` is the explicit opt-in for re-standardizing a different population in its own metric (and the only option for cor-matrix objects); `predict()` exported for discoverability.
**Consequences:** Breaking pre-CRAN for subset/new-data scoring values; replicators find the standard S3 idiom.
_Source: legacy §14.34; M45._

### D-022 (2026-07-01): Split-half factor comparability as a standalone verb

**Context:** Goldberg's research program gates hierarchy depth on split-half factor comparability (Everett 1983); the package needed that replicability floor.
**Decision:** `comparability()` standalone verb: full-sample fit anchors labels; per split, each half's solution is matched to the anchor by greedy-with-removal max-|r| bijection; the coefficient is the correlation of matched half-solution scores on the pooled R via `compute_edges()` (Invariant 1), with Tucker's φ alongside; report-first (benchmarks are reference lines only); PCA/EFA + pearson/spearman only; `n_splits = 10` default.
**Consequences:** The replicability-gated (Girard) workflow is documented end-to-end; ESEM/polychoric variants stay demand-gated ROADMAP candidates.
_Source: legacy §14.35; M46._

### D-023 (2026-07-01): Bootstrap percentile CIs on every edge, as a standalone verb

**Context:** Point-estimate edges (and the Forbes strongest-edge practice) lacked any precision statement; analytic SEs don't exist for edges of varimax-rotated hierarchies.
**Decision:** `boot_edges()` standalone verb: nonparametric row resampling, refit per replicate with full-sample anchoring (matching + sign orientation), edges via `compute_edges()`; percentile CIs (respect the [-1,1] bound and skew near .9); all indices drawn upfront from `seed` so serial ≡ parallel bit-for-bit; failed replicates → NA, counted, never aborting (Invariant 7); PCA/EFA + pearson/spearman only.
**Consequences:** `tidy(what = "edges")` gains `se`/`lo`/`hi`/`n_boot_ok`; per-edge error bars only — explicitly not a familywise correction (documented, not silently "fixed").
_Source: legacy §14.36; M47._

### D-024 (2026-07-02): bfi25 ships IPIP item labels; suggest_k() warns on ordinal input

**Context:** `top_items()`'s label capture path (M36) had no zero-setup demonstration, and `suggest_k()` lacked the ordinal advisory its siblings emit (Invariant 6 symmetry).
**Decision:** bfi25 item columns carry public-domain IPIP marker stems as plain `label` attributes (dropped by base row-subsetting, so docs say fit directly rather than pre-filter); `suggest_k()` runs `detect_ordinal()` on raw data and warns once per session, wording pointed at the final `ackwards()` fit (it cannot switch basis itself; a cor matrix has no items to inspect).
**Consequences:** `top_items()` prints `code: label` out of the box; advisory parity across the three raw-data entry points.
_Source: legacy §14.37–38; M50._

### D-025 (2026-07-02): Dual EFA chi-squares declined

**Context:** A review enhancement proposed reporting psych's residual-based `chi_empirical` alongside the likelihood-ratio `chi` in the EFA fit row.
**Decision:** Declined as net-negative: it re-opens the exact chi-square/p-value mispairing M42/C1 fixed, adds an NA-heavy EFA-only column to a schema shared with ESEM, and has zero downstream consumer (readable off `keep_fits = TRUE`).
**Consequences:** A confusion surface for no capability stays out; same shape as the D-027 `categorical` decline.
_Source: legacy §14.40; M49 Phase A._

### D-026 (2026-07-02): Item-label ecosystem — one shipped, two declined, one deferred

**Context:** bfi25's item labels (D-024) attracted adjacent asks: a setter, a third dataset, and factor labels.
**Decision:** `label_items()` setter declined (duplicates `labelled::var_label()`, deepens the item-vs-factor "label" overload); a third teaching dataset declined (`sim16`/`bfi25` are a deliberate two-foil pair; M54 clarified `forbes2023` is a fidelity dataset outside this scope); a persistent factor-label pipeline deferred — later resolved as D-029 (M51). The item-vs-factor label vocabulary split is load-bearing for all future API.
**Consequences:** One documented labeling idiom (`labelled::var_label()` / `attr()`); no vocabulary drift.
_Source: legacy §14.41; M49 Phase A; M54._

### D-027 (2026-07-01): `categorical` convenience flag declined

**Context:** A proposed `categorical = TRUE` on `ackwards()` would flip `cor` and the ESEM estimator together.
**Decision:** Declined as redundant: `cor = "polychoric"` already auto-selects WLSMV (the `estimator = NULL` rule, DESIGN §9), so the flag would be a pure synonym adding a conflict surface (`categorical = TRUE` + `cor = "pearson"`?) for zero new capability. Discoverability handled in docs and the runtime ordinal advisory.
**Consequences:** No defaults change; the ordinal vignette and cli warning name the real option.
_Source: legacy §14 (M40 spin-off block); M40._

### D-028 (2026-07-02): Polychoric robustness — `correct` arg, `check_items()`, trust tiering

**Context:** Real ordinal data broke `psych::polychoric()` opaquely (sparse cross-cells under the default continuity correction; constant items silently deleted), and fit-time warnings are transient on shared objects.
**Decision:** `ackwards()` gains `correct = 0.5` forwarded to psych on the PCA/EFA polychoric path (failure message names the `correct = 0` remedy); new export `check_items()` + shared internal screen (`constant` errors, `near-constant` warns once, `sparse category` report-only under polychoric); a "When to trust the result" roxygen section tiers every diagnostic fatal/caution/informational; basis-agnostic `.near_singular_check()` records `meta$min_eigenvalue`/`near_singular` so `print()`/`summary()` re-surface the caution durably.
**Consequences:** Wrap-and-diagnose, not reimplementation; report-first; no estimate changes.
_Source: legacy §14.42–44; M49._

### D-029 (2026-07-02): Persistent factor labels — verb + getter, `label (id)` display

**Context:** D-026 deferred the factor-label pipeline to its own naming/storage/precedence design; factor labels must stay lexically distinct from item labels.
**Decision:** `set_factor_labels(x, labels)` (pipeable; merge/update semantics; `NULL` clears all, `NA`/`""` removes one; unknown IDs error) + `factor_labels()` getter; storage in `x$meta$factor_labels` (rides through every verb unmodified); every text surface shows `label (id)` (IDs stay the indexing handle); `autoplot()` nodes show the label only (a stored label is a persistent `node_labels` default; call-time overrides win); `tidy()` gains label columns only when labels are set (unlabeled objects byte-identical pre/post — Invariant 5: IDs never mutated); `augment()`/`predict()` column names stay `m{k}f{j}`.
**Consequences:** Zero churn for unlabeled objects; one pipeline idiom (no replacement-function mutation).
_Source: legacy §14.45; M51._

### D-030 (2026-07-02): `factorability()` export + advisory-only internal screen

**Context:** Neither `check_items()` (per-item) nor `suggest_k()` (depth) answered "is this matrix worth factoring, and is N adequate".
**Decision:** `factorability()` reports KMO (overall + per-item MSA), Bartlett, N/p/N:p, and the Ledermann bound on the chosen basis — values with conventional bands framed as rules of thumb, never a pass/fail column (extends D-014); `ackwards()` runs an advisory-only screen: a Ledermann under-identification warning (EFA/ESEM only — PCA exempt, no latent-model df constraint) against the requested `k_max`, plus one consolidated sampling-adequacy warning (KMO < .5, N:p < 5, or N < 100), never aborting (Invariants 6/7).
**Consequences:** No signature or default change; no new dependency (psych supplies KMO/Bartlett).
_Source: legacy §14.46; M52._

### D-031 (2026-07-17): Design-interview principle set adopted; Forbes fidelity reframed as capability

**Context:** The /design-interview run formalized CLAUDE.md's eight numbered "Invariants" plus interview-elicited candidates into DESIGN.md's numbered IP/GP principles (the ROADMAP candidate since 2026-07-11). CLAUDE.md had claimed "the default output must reproduce Forbes's examples exactly," which the owner judged an over-lock.
**Decision:** Invariants 1–8 adopted as IP1–IP8 with **identity numbering**, so the ~22 in-code `Invariant N` citations resolve without a repoint. IP9 reframes the Forbes contract as a permanent reproducibility *capability* — exact reproduction of Forbes (2023) stays test-backed and its settings available/documented, but defaults may adopt a better method with loud documentation (IP6) and a D-entry. Six GPs adopted: GP1 published-method capability bar, GP2 report-first/flag-second, GP3 descriptive honesty, GP4 wrap-don't-reimplement, GP5 lean install, GP6 reproducible-by-construction. Deliberately *not* numbered: D-002's varimax-only policy and D-029's byte-identical-unlabeled-objects guarantee stay D-entries. Phase-1 elicited posture (audience, contract boundary, ambition, capability bar, fluid-until-1.0 with the BRM paper as the 1.0 trigger, track-CRAN-current upstream) recorded as DESIGN §1–§3 prose.
**Consequences:** Milestone plans cite "Principles touched" against these numbers; changing an IP takes a D-entry; CLAUDE.md's Invariants section is now a pointer; the "default must reproduce Forbes exactly" wording is superseded.
_Source: design interview 2026-07-17; DESIGN "Design principles"._

### D-032 (2026-07-24; extends D-017): Forbes's `ChaseCorrPaths` is contiguous — the direct criterion's break-at-gap is faithful; gap-tolerant chase rejected

**Context:** M53 (D-017) adopted the direct/skip-level correlation as the redundancy chase, breaking *contiguously* at the first sub-threshold level. M78 asked whether Forbes's `ChaseCorrPaths` is instead *gap-tolerant* — skipping a sub-threshold level to reach a deeper ancestor the leaf still correlates with directly. Settled empirically from committed data (no `ackwards()` fit, no network): her chase output (`amh$corr_chase`, 54 endpoints in `tests/testthat/fixtures/forbes2023_amh.rds`) is reproduced 54/54 by a contiguous chase over her own `comp_corr`; the single AMH component where the two semantics diverge, `g2` (level 7), has a mid-chain gap (direct `|r|` to level 6 `< .9`, re-emerging `>= .9` at level 5) and she reports `g2--null` (no move) — contiguous, not the level-5 skip. 0/54 components match gap-tolerant.
**Decision:** `ChaseCorrPaths` is **contiguous**. The existing `redundancy_criterion = "direct"` behaviour (`.strong_links_direct` breaks at the first sub-threshold hop) is confirmed faithful and left unchanged; the gap-tolerant variant is rejected. A regression test (`test-prune.R`, "M78: direct chase stops at a mid-chain gap") locks the contiguous semantics with the `g2` case in miniature.
**Consequences:** No code or output change; D-017's default stands, the contiguity question now closed and test-guarded.
_Source: M78; OSF `pcwm8` AMH fixture (chase output); extends D-017._
