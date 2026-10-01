# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-10-01 (M94 done and archived. Added the edge-wording sweep row from its review, extended the DESIGN §5 row, pruned the M91 row. Validate green.)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M95 | Edge wording sweep across shipped docs and DESIGN §5 | planned | — | normal | milestones/M95-edge-wording-sweep.md |
| M84 | Cross-branch-only secondary edges in the pruned view | blocked | — | normal | milestones/M84-cross-branch-secondary-edges.md |
| M94 | Manuscript accuracy pass for four older claims | done | M93 | normal | milestones/archive/M94-manuscript-accuracy-pass.md |
| M93 | Reconcile DESIGN's edge_method text with the code | done | — | normal | milestones/archive/M93-design-edge-method-reconcile.md |
| M92 | Manuscript wording for the oblique rotation option | done | — | normal | milestones/archive/M92-manuscript-oblique-wording.md |
<!-- M01–M80 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 3 most recent terminal rows. -->

## Candidates

- Further CRAN macOS coverage beyond M82: a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its push-to-master job) and a macOS x86_64 CI row. The 0.1.1 failure was arm64-only. Promote either on a CRAN macOS failure M82's job could not have caught. — added 2026-07-27, trimmed 2026-09-07, compressed 2026-10-01
- Untested axis from M83: disabling testthat parallelism outright (`Config/testthat/parallel: false` / `TESTTHAT_PARALLEL=false`), not only lowering the worker count. At `TESTTHAT_CPUS=1` the crashes still read `testthat subprocess exited`, so no sweep removed the subprocess. Costs M48's speedup and changes DESCRIPTION, so it needs release re-verification. Promote if a CRAN Windows flavour fails with the -1073741819 signature. — added 2026-07-27
- [low] ESEM `comparability()`: split-half per level per factor is feasible but runs 2·n_splits lavaan hierarchies per call and needs per-half convergence handling (D-022 / M46). Demand-gated: promote when a user asks. — added 2026-07-11, merged 2026-07-16, split 2026-10-01
- [low] WLSMV/polychoric `boot_edges()`: costs n_boot × (k_max−1) fits, and a resample can drop a response category (D-023 / M47). Demand-gated: promote when a user asks. — added 2026-07-11, merged 2026-07-16, split 2026-10-01
- [low] Oblique `boot_edges()` (RR02 Q5 item 12, held out of M90's Scope): it refits with the object's rotation but bootstraps only `r`, so `beta` and `r2` get no intervals. Promote when a user needs intervals on the partialled coefficients. — added 2026-09-30, split 2026-10-01
- [low] Oblique `comparability()` (RR02 Q5 item 12, held out of M90's Scope): it has no `rotation` argument and always fits varimax. Promote when a user needs split-half replicability of an oblique hierarchy. — added 2026-09-30, split 2026-10-01
- [low] ESEM varimax hides lavaan's rotation non-convergence, because lavaan's `rotation.args$warn` is off by default. This predates M90 (M90 review). Promote on a report of a silently non-converged varimax ESEM level. — added 2026-09-30, split 2026-10-01
- [low] `.esem_rotation_args()` detects lavaan's list form through the deprecated `rotation_args` formal, and its pre-0.7 branch is only unit-tested (M90 review). Promote when lavaan drops `rotation_args`. — added 2026-09-30, split 2026-10-01
- [low] Indefinite R and the near-singular warning (M91 review, F4). With pairwise missing data or a non-PD user matrix, and EFA or ESEM falling back to regression weights, Φ_s can have a negative eigenvalue. `tidy()` then calls it "nearly singular" with a negative value, and `r2` can leave [0, 1] with no warning. Promote on a real fit that shows either. — added 2026-09-30
- [low] User-facing sample-realized edges: an option to build edges from materialized scores, so that under missing data an edge shows the sample-realized correlation, not the model-implied one. DESIGN §5.3 described it, but no exported function offers it. Demand-gated: promote when a user asks for sample-realized edges. — added 2026-10-01 — M93 plan
- [low] The algebra-vs-scores cross-check misses more paths than DESIGN's first Known limitations entry names (M93 review, F3). Under `cor = "spearman"` the algebra uses Spearman R, but the scores branch of `compute_edges()` takes a Pearson `cor()` of standardized raw data. ESEM FIML and pairwise-missing Pearson data also differ in basis. Add them to the entry, or test them. Promote on a wrong edge from one of those paths. — added 2026-10-01 — M93 review
### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23

### Forbes correspondence (2026-07-30)

Batch from her reply to the M76–M81 write-up. The cross-branch secondary-edge request → M84 (blocked on her). Her publication-figure and near-redundant-band feedback needed no action. The rotation question went to RB02/RR02 and produced D-034; the rows below are the residue, all gated on the same design session as [H] above.

- [high] **D-032's premise contradicted by its source.** D-032 (2026-07-24) rejected gap-tolerant redundancy chains partly on the inference that Forbes's contiguous `ChaseCorrPaths` was deliberate. She confirms the empirical finding and denies the inference — the contiguity is a coding limitation, and she handles a dead intervening level by hand at the artefact stage, still on the redundancy criterion. M53's 54/54 reproduction and M78's `g2` regression test stand; only the reading of intent falls. **Gated:** a supersede takes the call. — added 2026-07-30
- **Promote Tucker's φ to the documented rotation-consistency diagnostic (RR02 rec 5, consider).** RR02 Q3 finds Forbes's "between-level correlations tell us whether the rotations are consistent/robust" reading is in role-conflict with reading the same edges as structure, and that the clean resolution is giving the robustness duty to Tucker's φ — which `prune()` already computes and reports. Docs-only if adopted. **Gated:** design-session agenda item. — added 2026-07-30
- **D-017's retention rule vs her published Figure 6B.** For the AMH chain `E1-F1-G1-H1-I1-J1`, D-017's retention (keep the bottom when the chain reaches `k_max`, else the topmost) retains only `J1`; Figure 6B retains `E1` **and** `J1`. Measured at M84 against her own data with her node set supplied by hand; her 6A green/grey colouring does not map 1:1 onto 6B's node set either. The chase itself still reproduces 54/54 (M53) — the divergence is in the *retention* step, which is the package's construction on top of her chase and not something her code specifies. **Gated:** what her retention rule actually is, is a question for her — take it to the design session. — added 2026-07-30
- **Unbroken chase to level a.** Her `ChaseCorrPaths()` returns "X--null" for a level-3+ component whose chase is unbroken to level a. The cause is `which.min()` on a vector with no `FALSE`, which counts zero links. `prune("redundant")` chases such a component to a1. Found by M90's oblique fixture (sim 3, minres, promax) and recorded in `references/source-departures.md` row E6 (corrected M90: was row M2). **Gated:** confirm her intent at the design session, then record the departure's D-entry (IP9) or change `prune()`. — added 2026-09-30
