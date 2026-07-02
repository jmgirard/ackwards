# ROADMAP.md — planned milestones (M41 review findings; follow-ups to be scoped)

Forward-looking counterpart to [`MILESTONES.md`](MILESTONES.md). `MILESTONES.md` is the source of
truth for **completed** milestones; this file captures the **intent and source notes** for the
*not-yet-built* milestones, so their context survives across planning sessions.

This file currently holds the **still-pending findings of the M41 independent review**
(statistical correctness, software design, vignette quality, and a defaults/decision audit;
conducted 2026-07-01 by Claude Fable 5). **M42 fixed the code findings** (C1, M1, m1, m2, m6,
m7, m8, m10, e3) and **M43 fixed the doc findings** (M2–M5, m3, m4, m5, m9, m11, e1) — see their
`MILESTONES.md` entries; fixed items are removed below per this file's maintenance rule. What
remains pending is the Forbes-fixture scoping question (proposed M44) and two unscheduled
enhancements. The "verified clean" and defaults-audit sections at the bottom are **reference
results** of the completed review (the M41 `MILESTONES.md` entry points here for their detail),
not pending work.

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

*(M2–M5, the vignette findings, were fixed in M43; M6 below remains.)*

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

*(all fixed: m1, m2, m6, m7, m8, m10 in M42; m3, m4, m5, m9, m11 in M43.)*

### Enhancements / justification fixes (no behavior change)

*(e1 fixed in M43; e3 in M42.)*

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
