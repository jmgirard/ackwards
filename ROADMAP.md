# ROADMAP.md — planned milestones (M41 review findings; follow-ups to be scoped)

Forward-looking counterpart to [`MILESTONES.md`](MILESTONES.md). `MILESTONES.md` is the source of
truth for **completed** milestones; this file captures the **intent and source notes** for the
*not-yet-built* milestones, so their context survives across planning sessions.

This file currently holds the **still-pending findings of the M41 independent review**
(statistical correctness, software design, vignette quality, and a defaults/decision audit;
conducted 2026-07-01 by Claude Fable 5). Each Critical/Major finding was verified with a runnable
reproduction before being recorded. **M42 (shipped 2026-07-01) fixed the code findings** — C1,
M1, m1, m2, m6, m7, m8, m10, e3; see its `MILESTONES.md` entry — and those items have been
removed below per this file's maintenance rule. What remains is the doc-fix set (proposed M43)
and the Forbes-fixture scoping question (proposed M44), each needing its own `/plan-milestone`
run.

---

## What was verified clean (no action needed)

The statistical core survives independent re-derivation. Verified numerically on
`na.omit(bfi25)` (n = 816) and `sim16` unless noted:

- **tenBerge weights** — `.tenBerge_weights()`'s `W = R⁻¹L(L'R⁻¹L)^{-1/2}` matches the ten Berge
  et al. (1999) orthogonal-case formula and agrees with `psych::factor.scores(method = "tenBerge")`
  to 1.6e-15; `W'RW = I` holds to 1e-8.
- **Edge algebra & standardization** — `compute_edges()`'s `W'RW` with `sqrt(diag(W'RW))`
  standardization re-derived independently; algebra vs. materialized-scores agreement 2.3e-15
  (PCA), and the existing suite carries the same oracle plus a `psych::bassAckward()` comparison
  for both PCA and EFA.
- **Sign alignment (M35 contract)** — every primary-parent edge ≥ 0 on PCA and EFA, `pairs = "all"`
  included; the skip-level recompute path (`ackwards.R` final `compute_edges()` from flipped
  weights) is algebraically identical to sweeping, so skip-level signs are correct.
- **Forbes machinery** — Tucker's φ formula exact; DFS chain enumeration and the retention rule
  (bottom iff chain reaches `k_max`, else top) verified on all five sim16 chains; the
  primary-parent-only restriction on chain links is **mathematically lossless** for thresholds
  > √.5: a child's squared correlations with its orthogonal parents sum to ≤ 1, so only the
  primary parent can reach |r| ≥ .9.
- **ESEM fit extraction** — engine fit rows match `lavaan::fitMeasures()` exactly: naive values
  under ML, scaled variants under WLSMV (chi/CFI verified to 1e-8), `BIC = NA` under WLSMV;
  polychoric edge `R` is bit-identical to lavaan's `sampstat$cov` (never psych's polychoric);
  `factor_cor = I` under varimax; per-level variance sorted descending.
- **`suggest_k` mappings** — MAP/VSS-1/VSS-2 recommendations match direct `psych::vss()` output;
  MAP correctly computed on components (`fm = "pc"`, per Velicer's definition).
- **Missing-data guard matrix (M16/M38)** — `.resolve_missing()` enforces exactly the documented
  combinations; PCA cumulative variance equals the top-k eigenvalue share; `detect_ordinal()`
  boundary behavior (8-level integer not flagged, binary flagged) as documented.

---

## Findings

### Critical

*(none pending — C1, the EFA chi/p-value pairing, was fixed in M42.)*

### Major

**M2. Engines vignette "Missing data" section documents pre-M38 behavior.**
`vignettes/ackwards-engines.Rmd` (Missing data + FIML subsection + "Which option to use?" table)
still states `"fiml"` is **ESEM-only** and "Errors for PCA/EFA", presents the manual
`psych::corFiml()` → cor-matrix route as the only PCA/EFA option, and says "A future release may
promote this pattern to a first-class `missing = 'fiml'` route" — **M38 shipped exactly that**.
The vignette now contradicts `?ackwards` and actual behavior (`missing = "fiml"` works for
PCA/EFA on the Pearson basis and announces itself via cli). M37 wrote this section with a forward
reference; M38 (code) never circled back to it, and M39's prose pass didn't cover the engines
vignette. **Fix:** rewrite the section around the first-class route; keep the corFiml seam as a
"what it does under the hood / non-Pearson caveat" note; update the recommendation table row.

**M3. suggest-k vignette misstates the Comparison Data mechanism.**
"CD resamples from the marginal item distributions **without preserving inter-item
correlations**, so on datasets with strong structure it can over-retain" (worked-recommendation
section). Ruscio & Roche (2012) comparison data are generated to reproduce the observed
**correlation structure under a known k-factor model** (plus the marginal distributions) — 
preserving the correlational structure is precisely the method's advance over PA. The earlier
five-criteria description in the same vignette is fine; this later sentence is wrong and the
"therefore it can over-retain" causal claim built on it is unsupported. **Fix:** correct the
mechanism sentence and re-derive (or drop) the over-retention aside.

**M4. Forbes vignette presents the artifact rule's by-construction zero as an empirical finding.**
"For the BFI, the artifact criterion flags `r n_artifact` factors — no factor at any level has a
loading pattern more similar to a factor from a non-adjacent level... This is a good result for a
well-validated instrument." But `prune(x, "artifact")` **never auto-flags** (DESIGN §14.21), so
`n_artifact` is 0 for every dataset; the same vignette's "Tuning the thresholds" section states
this correctly ("no factors are auto-flagged"), making the vignette internally contradictory. The
opening definition ("a factor whose loading pattern is more similar to a factor at a non-adjacent
level than to its own-level neighbors") also implies a specific automated comparison that is not
what the code computes (it reports φ for *all* cross-level pairs, plus structural signals, for
researcher judgment). **Fix:** rewrite the artifact section around report-and-judge semantics —
show `x$prune$phi` extremes instead of the vacuous zero-count table.

**M5. Forbes vignette still describes the retired `cut_strong` linetype split.**
"Solid lines indicate strong connections (|r| ≥ 0.5 by default); dashed lines indicate weaker
ones — the same `cut_strong` threshold used throughout" — `cut_strong` was retired in **M35**
(deprecated arg, warning, no effect), and the figures above that sentence are drawn with uniform
black lines (`edge_linewidth = 0.6`). Nothing in the current render produces a solid/dashed split.
**Fix:** delete/replace the sentence (and audit the paragraph for other pre-M35 remnants).

**M6. The Forbes fidelity contract is untested.**
CLAUDE.md states "the default output must reproduce Forbes's examples exactly" as the baseline
contract, but no test anywhere reproduces any Forbes (2023) example — `test-prune.R` asserts the
retention rule and chain structure on constructed/simulated cases only (good tests, wrong
question). The algorithm matches the paper's *description* (thresholds, conjunctive φ, retention
rule — verified in this review), but "reproduces her examples exactly" is an empirical claim no
one has checked. **Fix (scoping decision for the owner):** if Forbes's OSF materials/data are
obtainable, add a fixture test reproducing a published chain/pruning table; if not, soften the
CLAUDE.md contract wording to "faithful to the published algorithm" so the docs don't promise a
verification that doesn't exist.

### Minor

*(m1, m2, m6, m7, m8, and m10 were fixed in M42; the doc-facing minors below remain.)*

- **m3.** Forbes vignette chain example: "the intermediate nodes m3f2 and m4f2 are flagged" — for
  a chain reaching `k_max` the retention rule keeps only the bottom node, so the **top** node is
  flagged too. The illustrative sentence misstates the rule the code (correctly) implements.
- **m4.** Engines vignette "Choosing an engine" table recommends `"pca", fm = "pca"` for
  replicating Goldberg — `fm = "pca"` is not a legal value (`minres`/`ml`/`pa`) and errors;
  `fm` is EFA-only anyway. Should read plain `engine = "pca"`.
- **m5.** Engines vignette says "the BFI with > 2,000 participants" — the vignette fits
  `na.omit(bfi25)` (1,000 rows, 816 complete). Stale from a `psych::bfi` (n = 2,800) draft.
- **m9.** `data-raw/sim16.R` header comments use the pre-M34 API (`prune(artefact = TRUE)`,
  `prune(redundant = TRUE)`) and "phi >= .95" (the filter is strict `>`); `R/data.R` (`?sim16`)
  says "6-criteria consensus," counting VSS twice without saying so.
- **m11.** suggest-k vignette hardcodes stochastic PA outcomes in prose ("PA-FA exceeds PA-PC in
  this run (6 vs. 5)", "CD agrees with PA-FA (k = 6)"); `fa.parallel` ignores `seed`, so a
  vignette rebuild can contradict its own prose. Compute these inline (the sim16 section already
  does it right) or hedge harder.

### Enhancements / justification fixes (no behavior change)

- **e1.** §9's PCA rationale for `redundancy_phi = NULL` ("no φ filter — the W'RW algebra is
  exact") conflates two properties. The algebra is equally exact for tenBerge-scored EFA; the
  *actual* reason PCA needs no congruence guard is that **component scores are determinate**
  (no factor-score indeterminacy), so `|r|` between component scores is the true correlation
  between the components themselves. The EFA/ESEM half of the rationale (indeterminacy) is
  correct — the PCA half should say determinacy, not algebra-exactness. Fix wording in DESIGN §9,
  CLAUDE.md, and `prune()` roxygen.
- **e2.** Consider exposing both chi-squares for EFA (`chi` = likelihood — the M42/C1 fix,
  `chi_empirical` = psych's residual-based) — psych prints both for a reason (the empirical one
  is robust to non-normality/NPD matrices). Owner set this aside when approving M42 ("can be
  revisited later").
- **e4.** Bootstrap CIs on (skip-level) edges — the standing DESIGN §14 deferral; this review
  re-affirms it as the highest-value statistical addition (the selection-bias concern about
  "strongest edge" claims is real) and as its own perf-heavy milestone.

---

## Defaults & decision audit (Phase 4 verdicts)

Every DESIGN §9 defaults-table row and §14 numbered decision was audited for both the choice and
its stated justification. Verdicts: **sound** unless listed below.

**§9 rows — all sound**, with two annotations: `redundancy_phi` auto-rule is
**sound-but-misjustified** on the PCA side (see e1; the 0.95 value itself is a defensible borrow
of Lorenzo-Seva & ten Berge's equivalence threshold, applied conjunctively — conservative by
construction); `missing = "pairwise"` is sound *given* the documented ESEM ML/MLR inconsistency
and its per-call warning (the M16 disclosure regime is the right call).

**§14 decisions 1–33 — all sound.** Spot-checked in detail: item 6 (≤ 7 distinct integer values)
is an honest, documented heuristic with the expected false-positive (integer-coded counts) and
false-negative (8+ category Likert) edges — acceptable for a warning-only signal; item 7/§14.1
(CF(κ=1/p) ≡ varimax, kappa removal) verified against Crawford & Ferguson (1970)/Browne (2001);
item 19 (|r| ≥ .9, retention rule, φ > .95 conjunctive) verified faithful to Forbes both in code
and output; items 27–31 (prune verb) sound as a design but shipped a `drop_pruned` regression (fixed in
M42) — the *decision* was right, the migration missed one consumer; items 32–33 (M38 FIML +
`n_obs` strings) sound and well-cited (Enders 2010; Zhang & Savalei 2020), but shipped the M2
doc gap.

**Declined decisions — all declines hold.** EAP (shrinkage attenuates the cross-level signal —
statistically correct reasoning); oblique rotation (T′ = T⁻¹ is load-bearing for the algebra;
oblique would confound the method's core signal); `categorical` flag (pure synonym, correctly
refused); EKC/EGA (dependency cost, defensible); Hungarian matching removal (bijection is
ill-posed under the pigeonhole argument — and this review adds the stronger Σr² ≤ 1 argument that
greedy argmax and "the ≥ .9 link" necessarily coincide).

**Arbitrary-constant inventory.** Documented with rationale: `cut_show = 0.3`, `redundancy_r =
0.9`, `redundancy_phi = 0.95`, `min_items = 3`, `orphan_r = 0.5`, ordinal ≤ 7, `k_max` default
`min(p−1, 8)`, `n_iter = 20`, PA quantile 0.95, Hu-Bentler reference lines, CD α = 0.30 (inherited
from EFAtools, named in the vignette). Undocumented but acceptable (cosmetic/numeric-hygiene):
symmetry/diagonal tolerances `1e-8` (utils.R), the eigenvalue floor `.Machine$double.eps` in
`.tenBerge_weights()`, linewidth legend range `c(0.4, 1.8)`, edge-label nudge `0.15`, level-label
offset `0.8`, arrowhead `0.15 cm`. None warrants promotion to an argument; the two numeric-hygiene
constants deserve a one-line comment at most.

**Rationale drift.** Three stale in-code comments found (M1's layout.R comment, m6's summary.R
comment, m9's data-raw comments); CLAUDE.md / DESIGN.md / MILESTONES.md / roxygen otherwise tell
the same story for every audited decision.

---

## Proposed follow-up milestones (each needs its own /plan-milestone)

- **M43 — review fixes, docs:** M2 (engines missing-data rewrite), M3 (CD mechanism), M4
  (artifact section rewrite), M5 (cut_strong remnant), m3, m4, m5, m9, m11, e1 (§9/roxygen
  justification wording). Doc-only; no export/signature change. Note: the engines vignette's
  rendered EFA fit tables will also pick up the M42 chi values on rebuild.
- **M44 — Forbes fixture (scoping):** M6 — owner decides: obtain Forbes (2023) materials and add
  an exact-reproduction test, or amend the CLAUDE.md contract wording.
- **(Unscheduled)** e2 (dual chi-squares for EFA) and e4 (bootstrap-CI milestone, deferred in
  DESIGN §14).

## Provenance

The M31–M40 arc was decomposed from a page-by-page pkgdown-site review the owner did on
**2026-06-30** (origin transcript `9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl`). The M41 findings
above were produced by the model-led review milestone (branch `m41-fable-review`, 2026-07-01);
its scope and acceptance criteria are logged in `MILESTONES.md` (M41 entry).

**How to maintain this file:** when a future milestone is scoped, capture its intent and source
notes here (raw review questions + banked decisions — inputs to planning, not a finished spec);
when it ships, delete its section and rely on its `MILESTONES.md` entry. Keep this file scoped to
*pending* work only.
