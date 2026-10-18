# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-10-18 (triage pass. Compressed four Forbes rows, retitled "D-032's premise contradicted by its source" → "D-032 premise contradicted by its source", "Promote Tucker's φ…" → "Tucker's φ as the rotation-consistency diagnostic", "D-017's retention rule vs her published Figure 6B" → "D-017 retention vs Figure 6B", "Unbroken chase to level a" kept. Retitled six over-long titles: "Oblique `boot_edges()` intervals for beta and r2", "Oblique `comparability()`", "Indefinite R near-singular warning", "Manuscript varimax and oblimin example", "Exported factor-matching helper", "Partial ESEM rotation failure on the object". Kept the `.esem_rotation_args()` row with no further extension. No drops, no decision entry. Validate green.)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M099 | Chase-gap reason corrected and a manual bridge documented | planned | — | high | milestones/M099-chase-gap-supersede.md |
| M098 | Counted lavaan rotation failures on every ESEM rotation | done | — | normal | milestones/archive/M098-esem-rotation-start-count.md |
| M097 | Cross-check coverage and the matrix behind each edge | done | — | normal | milestones/archive/M097-cross-check-coverage.md |
| M84 | Cross-branch-only secondary edges in the pruned view | blocked | — | normal | milestones/M84-cross-branch-secondary-edges.md |
| M96 | Varimax and oblimin side by side on bfi25 in the engines vignette | done | — | normal | milestones/archive/M96-rotation-comparison.md |
<!-- M01–M80 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 3 most recent terminal rows. -->

## Candidates

- Further CRAN macOS coverage beyond M82: a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its push-to-master job) and a macOS x86_64 CI row. The 0.1.1 failure was arm64-only. Promote either on a CRAN macOS failure M82's job could not have caught. — added 2026-07-27, trimmed 2026-09-07, compressed 2026-10-01
- Untested axis from M83: disabling testthat parallelism outright (`Config/testthat/parallel: false` / `TESTTHAT_PARALLEL=false`), not only lowering the worker count. At `TESTTHAT_CPUS=1` the crashes still read `testthat subprocess exited`, so no sweep removed the subprocess. Costs M48's speedup and changes DESCRIPTION, so it needs release re-verification. Promote if a CRAN Windows flavour fails with the -1073741819 signature. — added 2026-07-27
- [low] ESEM `comparability()`: split-half per level per factor is feasible but runs 2·n_splits lavaan hierarchies per call and needs per-half convergence handling (D-022 / M46). Demand-gated: promote when a user asks. — added 2026-07-11, merged 2026-07-16, split 2026-10-01
- [low] WLSMV/polychoric `boot_edges()`: costs n_boot × (k_max−1) fits, and a resample can drop a response category (D-023 / M47). Demand-gated: promote when a user asks. — added 2026-07-11, merged 2026-07-16, split 2026-10-01
- [low] Oblique `boot_edges()` intervals for beta and r2: Oblique `boot_edges()` (RR02 Q5 item 12, held out of M90's Scope): it refits with the object's rotation but bootstraps only `r`, so `beta` and `r2` get no intervals. Promote when a user needs intervals on the partialled coefficients. — added 2026-09-30, split 2026-10-01
- [low] Oblique `comparability()`: Oblique `comparability()` (RR02 Q5 item 12, held out of M90's Scope): it has no `rotation` argument and always fits varimax. Promote when a user needs split-half replicability of an oblique hierarchy. — added 2026-09-30, split 2026-10-01
- [low] `.esem_rotation_args()` detects lavaan's list form through the deprecated `rotation_args` formal, and its pre-0.7 branch is only unit-tested (M90 review). M098's count of one warning per failed start was measured only on lavaan 0.7.2, not on the pre-0.7 `rotation.args` path (M098 review). Promote when lavaan drops `rotation_args`. — added 2026-09-30, split 2026-10-01, extended 2026-10-06
- [low] Indefinite R near-singular warning: Indefinite R and the near-singular warning (M91 review, F4). With pairwise missing data or a non-PD user matrix, and EFA or ESEM falling back to regression weights, Φ_s can have a negative eigenvalue. `tidy()` then calls it "nearly singular" with a negative value, and `r2` can leave [0, 1] with no warning. Promote on a real fit that shows either. — added 2026-09-30
- [low] User-facing sample-realized edges: an option to build edges from materialized scores, so that under missing data an edge shows the sample-realized correlation, not the model-implied one. DESIGN §5.3 described it, but no exported function offers it. Demand-gated: promote when a user asks for sample-realized edges. — added 2026-10-01 — M93 plan
- [low] Manuscript varimax and oblimin example: Manuscript worked example of varimax and oblimin on one dataset. The manuscript covers rotation in prose only. M96 puts the bfi25 comparison in the engines vignette. If a manuscript revision or a reviewer asks for an example, promote it. — added 2026-10-01 — M96 plan
- [low] Exported factor-matching helper: Exported helper that matches factors across two fits by loading congruence. Two fits can number the same factor differently (oblimin and varimax swap two level-5 IDs on bfi25), and M96's vignette matches them with `psych::factor.congruence()`. If a user asks for it, promote it. — added 2026-10-01 — M96 plan
- [low] Partial ESEM rotation failure on the object: Record a partial ESEM rotation failure on the object. A kept level where some lavaan starts failed shows only a warning (M098 review), so `print()`, `glance()`, and a run under `suppressWarnings()` lose the count. Promote when a user asks to read the count after the fit. — added 2026-10-06 — M098 review
### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23

### Forbes correspondence (2026-07-30)

Batch from her reply to the M76–M81 write-up. The cross-branch secondary-edge request → M84 (blocked on her). Her publication-figure and near-redundant-band feedback needed no action. The rotation question went to RB02/RR02 and produced D-034; the rows below are the residue, all gated on the same design session as [H] above.

- [high] D-032 premise contradicted by its source: D-032 rejected gap-tolerant chains partly on reading Forbes's contiguous `ChaseCorrPaths` as deliberate. She confirms the finding and denies the inference: the contiguity is a coding limitation, and she handles a dead intervening level by hand, still on the redundancy criterion. M53's 54/54 reproduction and M78's `g2` test stand. **Gated:** a supersede takes the call. — added 2026-07-30 — started M099
- Tucker's φ as the rotation-consistency diagnostic: RR02 rec 5 (consider). RR02 Q3 finds that reading the edges both as structure and as a rotation-robustness check is a role conflict, and gives the robustness role to Tucker's φ, which `prune()` already reports. Docs-only if adopted. **Gated:** design-session agenda item. — added 2026-07-30
- D-017 retention vs Figure 6B: on the AMH chain `E1-F1-G1-H1-I1-J1`, D-017's retention keeps only `J1`, and her Figure 6B keeps `E1` and `J1` (measured at M84 on her data, her node set supplied by hand). Her 6A colouring does not map 1:1 onto 6B's node set either. The chase still reproduces 54/54 (M53). The divergence is in the retention step, which her code does not specify. **Gated:** ask her retention rule at the design session. — added 2026-07-30
- Unbroken chase to level a: her `ChaseCorrPaths()` returns "X--null" for a level-3+ component whose chase is unbroken to level a (`which.min()` on a vector with no `FALSE`), and `prune("redundant")` chases it to a1. Found by M90's oblique fixture, `references/source-departures.md` row E6 (corrected M90). **Gated:** confirm her intent at the design session, then record the departure's D-entry (IP9) or change `prune()`. — added 2026-09-30
