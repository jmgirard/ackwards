# Decisions

Append-only. Never renumber; supersede with a new entry. D-entries record
choices with rationale — never deferrals.

Pre-migration decision history lives in `cairn/DESIGN.md` §14 (embedded log,
kept verbatim) and `cairn/legacy/`. This file re-records only the
still-governing, cross-cutting decisions, each citing its DESIGN §14 anchor;
new decisions continue from the next D-number.

### D-001 (date: see legacy §2): A result is a series of linked solutions, not a fitted hierarchy

**Context:** A bass-ackwards result is a set of linked per-level solutions whose edges are score correlations, not a fitted higher-order SEM.
**Decision:** Treat the "hierarchy" as descriptive/emergent and state this honesty caveat in docs and `print()`.
**Consequences:** `print()` states it once; `prune()` is interpretive relabeling, not re-estimation; frames every claim the package makes.
_Source: DESIGN §2._

### D-002 (2026-06): Varimax is the only rotation; oblique out of scope

**Context:** Only orthogonal rotation (`T' = T^-1`) yields interpretable between-level score correlations and the closed-form `W'RW` algebra; oblique confounds the cross-level signal that is the method's point.
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
