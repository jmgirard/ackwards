# ROADMAP.md — planned milestones (M41 review findings; follow-ups to be scoped)

Forward-looking counterpart to
[`MILESTONES.md`](https://jmgirard.github.io/ackwards/MILESTONES.md).
`MILESTONES.md` is the source of truth for **completed** milestones;
this file captures the **intent and source notes** for the
*not-yet-built* milestones, so their context survives across planning
sessions.

**The M41 review remediation is complete.** M42 fixed the code findings
(C1, M1, m1, m2, m6, m7, m8, m10, e3), M43 fixed the doc findings
(M2–M5, m3, m4, m5, m9, m11, e1), and M44 resolved the last finding (M6)
by shipping a Forbes-fidelity fixture test — see their `MILESTONES.md`
entries. **There are currently no pending milestones**; the unscheduled
ideas below would each need a scoping discussion before a
`/plan-milestone` run. The “verified clean” and defaults-audit sections
at the bottom are **reference results** of the completed review (the M41
`MILESTONES.md` entry points here for their detail), not pending work.

------------------------------------------------------------------------

## What was verified clean (no action needed)

The statistical core survives independent re-derivation. Verified
numerically on `na.omit(bfi25)` (n = 816) and `sim16` unless noted:

- **tenBerge weights** — `.tenBerge_weights()`’s
  `W = R⁻¹L(L'R⁻¹L)^{-1/2}` matches the ten Berge et al. (1999)
  orthogonal-case formula and agrees with
  `psych::factor.scores(method = "tenBerge")` to 1.6e-15; `W'RW = I`
  holds to 1e-8.
- **Edge algebra & standardization** —
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)’s
  `W'RW` with `sqrt(diag(W'RW))` standardization re-derived
  independently; algebra vs. materialized-scores agreement 2.3e-15
  (PCA), and the existing suite carries the same oracle plus a
  [`psych::bassAckward()`](https://rdrr.io/pkg/psych/man/bassAckward.html)
  comparison for both PCA and EFA.
- **Sign alignment (M35 contract)** — every primary-parent edge ≥ 0 on
  PCA and EFA, `pairs = "all"` included; the skip-level recompute path
  (`ackwards.R` final
  [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  from flipped weights) is algebraically identical to sweeping, so
  skip-level signs are correct.
- **Forbes machinery** — Tucker’s φ formula exact; DFS chain enumeration
  and the retention rule (bottom iff chain reaches `k_max`, else top)
  verified on all five sim16 chains; the primary-parent-only restriction
  on chain links is **mathematically lossless** for thresholds \> √.5: a
  child’s squared correlations with its orthogonal parents sum to ≤ 1,
  so only the primary parent can reach \|r\| ≥ .9.
- **ESEM fit extraction** — engine fit rows match
  [`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html)
  exactly: naive values under ML, scaled variants under WLSMV (chi/CFI
  verified to 1e-8), `BIC = NA` under WLSMV; polychoric edge `R` is
  bit-identical to lavaan’s `sampstat$cov` (never psych’s polychoric);
  `factor_cor = I` under varimax; per-level variance sorted descending.
- **`suggest_k` mappings** — MAP/VSS-1/VSS-2 recommendations match
  direct [`psych::vss()`](https://rdrr.io/pkg/psych/man/VSS.html)
  output; MAP correctly computed on components (`fm = "pc"`, per
  Velicer’s definition).
- **Missing-data guard matrix (M16/M38)** — `.resolve_missing()`
  enforces exactly the documented combinations; PCA cumulative variance
  equals the top-k eigenvalue share; `detect_ordinal()` boundary
  behavior (8-level integer not flagged, binary flagged) as documented.

------------------------------------------------------------------------

## Findings

**All findings are resolved.** Critical: C1 fixed in M42. Major: M2–M5
fixed in M43; M6 (the untested Forbes-fidelity contract) fixed in M44
via a fixture test reproducing the paper’s three simulation studies
against Forbes’s own reference implementation. Minor: m1, m2, m6, m7,
m8, m10 in M42; m3, m4, m5, m9, m11 in M43. Enhancements: e3 in M42; e1
in M43; e4 (bootstrap edge CIs) shipped as M47; e2 (dual EFA
chi-squares) **declined in M49** (DESIGN §14 item 40).

------------------------------------------------------------------------

## Defaults & decision audit (Phase 4 verdicts)

Every DESIGN §9 defaults-table row and §14 numbered decision was audited
for both the choice and its stated justification. Verdicts: **sound**
unless listed below.

**§9 rows — all sound**, with two annotations: `redundancy_phi`
auto-rule is **sound-but-misjustified** on the PCA side (see e1; the
0.95 value itself is a defensible borrow of Lorenzo-Seva & ten Berge’s
equivalence threshold, applied conjunctively — conservative by
construction); `missing = "pairwise"` is sound *given* the documented
ESEM ML/MLR inconsistency and its per-call warning (the M16 disclosure
regime is the right call).

**§14 decisions 1–33 — all sound.** Spot-checked in detail: item 6 (≤ 7
distinct integer values) is an honest, documented heuristic with the
expected false-positive (integer-coded counts) and false-negative (8+
category Likert) edges — acceptable for a warning-only signal; item
7/§14.1 (CF(κ=1/p) ≡ varimax, kappa removal) verified against Crawford &
Ferguson (1970)/Browne (2001); item 19 (\|r\| ≥ .9, retention rule, φ \>
.95 conjunctive) verified faithful to Forbes both in code and output;
items 27–31 (prune verb) sound as a design but shipped a `drop_pruned`
regression (fixed in M42) — the *decision* was right, the migration
missed one consumer; items 32–33 (M38 FIML + `n_obs` strings) sound and
well-cited (Enders 2010; Zhang & Savalei 2020), but shipped the M2 doc
gap.

**Declined decisions — all declines hold.** EAP (shrinkage attenuates
the cross-level signal — statistically correct reasoning); oblique
rotation (T′ = T⁻¹ is load-bearing for the algebra; oblique would
confound the method’s core signal); `categorical` flag (pure synonym,
correctly refused); EKC/EGA (dependency cost, defensible); Hungarian
matching removal (bijection is ill-posed under the pigeonhole argument —
and this review adds the stronger Σr² ≤ 1 argument that greedy argmax
and “the ≥ .9 link” necessarily coincide).

**Arbitrary-constant inventory.** Documented with rationale:
`cut_show = 0.3`, `redundancy_r = 0.9`, `redundancy_phi = 0.95`,
`min_items = 3`, `orphan_r = 0.5`, ordinal ≤ 7, `k_max` default
`min(p−1, 8)`, `n_iter = 20`, PA quantile 0.95, Hu-Bentler reference
lines, CD α = 0.30 (inherited from EFAtools, named in the vignette).
Undocumented but acceptable (cosmetic/numeric-hygiene):
symmetry/diagonal tolerances `1e-8` (utils.R), the eigenvalue floor
`.Machine$double.eps` in `.tenBerge_weights()`, linewidth legend range
`c(0.4, 1.8)`, edge-label nudge `0.15`, level-label offset `0.8`,
arrowhead `0.15 cm`. None warrants promotion to an argument; the two
numeric-hygiene constants deserve a one-line comment at most.

**Rationale drift.** Three stale in-code comments found (M1’s layout.R
comment, m6’s summary.R comment, m9’s data-raw comments); CLAUDE.md /
DESIGN.md / MILESTONES.md / roxygen otherwise tell the same story for
every audited decision.

------------------------------------------------------------------------

## Unscheduled ideas (each needs a scoping discussion before /plan-milestone)

- **AMH applied-example fidelity extension.** The M44 feasibility study
  verified `ackwards` against Forbes’s reference implementation on her
  full 155-variable AMH applied example (max deviation 3.9e-14 across
  all 45 level-pairs at `k_max = 10`), but the AMH Spearman matrix
  (`corSpearman_AMH.csv`, OSF `pcwm8`) carries **no license**, so it was
  not committed. Options when the owner contacts Forbes: (a) she adds a
  license/permission → commit a fixture like the simulation one; (b)
  keep it uncommitted and add a `skip_if_offline()`/`skip_on_cran()`
  test that downloads the matrix from OSF at test time (no
  redistribution; expected values from the M44 feasibility run). The
  download URLs and comparison script survive in the M44 milestone
  entry.
- **[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
  and
  [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)
  engine/basis extensions** (deferred from M46/M47, DESIGN
  §14.35/§14.36) — **feasible but demand-gated; keep off the schedule
  until asked.** Feasibility verdict (M49 Phase A):
  - **ESEM
    [`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)**
    is the plausible one: `2 * n_splits` (~20) lavaan hierarchies per
    call, minutes with the M26 cached-sample-stats + `future.apply`
    machinery; the matching / anchoring code is already engine-agnostic.
    Main work is per-half convergence handling (Invariant 7).
  - **ESEM
    [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)**
    is painful: `n_boot` (~1000) × (k_max − 1) WLSMV fits — tens of
    minutes to hours even parallel; users would have to cut `n_boot`,
    undermining the percentile CIs, and resample non-convergence /
    Heywood cases inflate the dropped-replicate count.
  - **Polychoric basis** (both verbs) is the least realistic — not just
    slow: a bootstrap resample can drop an entire response category,
    changing the threshold structure mid-replicate (a drop-or-merge
    policy question — a genuine statistical wrinkle, not just
    performance). Mirror
    [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
    Pearson/Spearman-only scope; screen on Pearson, fit the final model
    polychoric.

  Net: only
  ESEM-[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
  is worth reconsidering if demand appears; the EFA-on-Pearson screen
  covers the common case for each verb.
- **Factor-label pipeline** (headline 0.2.0 candidate; owner-approved
  deferral 2026-07-02) — persistent *factor* labels (distinct from M50’s
  *item* / variable labels): a `set_factor_labels()`-style setter
  storing user names on the `ackwards` object, honored by
  [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  (node labels), `print`/`summary` (cli), and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html);
  [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  remains the scaffold that seeds them. **Purely additive** (no
  signature/default change), so nothing is lost by shipping 0.1.0
  without it; the naming/storage/precedence design gets proper DESIGN
  §14 treatment when scoped. Keep the item-vs-factor “label” vocabulary
  split (DESIGN §14 item 39) front of mind — the two concepts must not
  blur.

## Provenance

The M31–M40 arc was decomposed from a page-by-page pkgdown-site review
the owner did on **2026-06-30** (origin transcript
`9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl`). The M41 findings above
were produced by the model-led review milestone (branch
`m41-fable-review`, 2026-07-01); its scope and acceptance criteria are
logged in `MILESTONES.md` (M41 entry).

**How to maintain this file:** when a future milestone is scoped,
capture its intent and source notes here (raw review questions + banked
decisions — inputs to planning, not a finished spec); when it ships,
delete its section and rely on its `MILESTONES.md` entry. Keep this file
scoped to *pending* work only.
