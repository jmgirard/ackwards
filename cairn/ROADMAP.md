# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-09-30 (M91 done and archived; two candidate rows added from its review; M88's terminal row pruned; two LESSONS lines added)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M92 | Manuscript wording for the oblique rotation option | review | — | normal | milestones/M92-manuscript-oblique-wording.md |
| M84 | Cross-branch-only secondary edges in the pruned view | blocked | — | normal | milestones/M84-cross-branch-secondary-edges.md |
| M91 | Near-singular guard for the within-level score correlation | done | — | normal | milestones/archive/M91-near-singular-score-cor-guard.md |
| M89 | Real within-level factor correlations and Φ-partialled edge reporting | done | — | normal | milestones/archive/M89-factor-cor-partialled-edges.md |
| M90 | Oblique rotation as a documented non-default option | done | M89 | normal | milestones/archive/M90-oblique-rotation-option.md |
<!-- M01–M80 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 3 most recent terminal rows. -->

## Candidates

- Further CRAN **macOS**-flavour coverage beyond M82: a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its own push-to-master job) and a macOS **x86_64** row (`macos-13`) — the 0.1.1 failure was arm64-specific and `r-oldrel-macos-x86_64` passed. Promote either on any CRAN macOS failure M82's job could not have caught. The row's alternative-numerics half shipped 2026-09-07 (R-hub `atlas` + `nold` mandated in the release-walk, no longer "as applicable"). — added 2026-07-27, trimmed 2026-09-07
- Untested axis from M83: **disabling testthat parallelism outright** (`Config/testthat/parallel: false` / `TESTTHAT_PARALLEL=false`), as distinct from the worker counts M83 measured. Evidence it matters: at `TESTTHAT_CPUS=1` the crashes still read `testthat subprocess exited`, so the parallel machinery still spawned a subprocess — the crashing component was never removed in any sweep. Plausibly the actual fix, and the only candidate that eliminates rather than reduces. Costs M48's speedup (27s parallel vs 81s serial locally) and, unlike M83's mitigation, changes **tarball content** (DESCRIPTION), so it needs release re-verification. Promote if a CRAN Windows flavour fails with the -1073741819 signature. — added 2026-07-27
- DESIGN still describes `edge_method = "scores"` as a user choice (IP2 text, the §9 `edge_method` row, and the line near 121), but `ackwards()` has no such argument and always builds `r` by algebra. The Known-limitations entry corrected at M91 now contradicts them. Reconcile, with a D-entry if IP2's wording changes. Found by the M91 review (F5). — added 2026-09-30
- Manuscript accuracy pass, four older claims from the M92 review. "Exact" needs a qualifier for a polychoric R, where no observed scores reproduce the edges. The methods sentence on a fallback to scores for nonlinear scoring describes a path no shipped engine takes (see the `edge_method` row above). Line 154 credits Waller with the general weight-matrix form, which `references/waller2007.md` calls ours. The AI-use disclosure (front matter and closing section) says June to July 2026. — added 2026-09-30
- [low] Upstream report for the M83 Windows crash, if it isolates to a dependency's compiled code (`EFAtools` or `mnormt`; `psych`/`lavaan`/`GPArotation` are pure R): file with a minimal reproducer. Held out of M83 at the 2026-07-27 plan gate because a maintainer's timeline must not gate our CI health. Promote when M83's bisection names a component. — added 2026-07-27
- [low] ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16
- [low] Oblique `boot_edges()` and `comparability()` semantics (RR02 Q5 item 12, held out of M90's Scope). `boot_edges()` refits an oblique object with its rotation but bootstraps only `r`, so `beta` and `r2` get no intervals. `comparability()` has no `rotation` argument and always fits varimax. Promote when a user needs intervals on the partialled coefficients or split-half replicability of an oblique hierarchy. — added 2026-09-30
- [low] lavaan rotation warnings, two gaps from the M90 review. ESEM varimax still hides lavaan's rotation non-convergence, because lavaan's `rotation.args$warn` is off by default (this predates M90). `.esem_rotation_args()` detects lavaan's list form through the deprecated `rotation_args` formal, and its pre-0.7 branch is only unit-tested. Promote on a report of a silently non-converged varimax ESEM level, or when lavaan drops `rotation_args`. — added 2026-09-30
- [low] Indefinite R and the near-singular warning (M91 review, F4). With pairwise missing data or a non-PD user matrix, and EFA or ESEM falling back to regression weights, Φ_s can have a negative eigenvalue. `tidy()` then calls it "nearly singular" with a negative value, and `r2` can leave [0, 1] with no warning. Promote on a real fit that shows either. — added 2026-09-30

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23

### Forbes correspondence (2026-07-30)

Batch from her reply to the M76–M81 write-up. The cross-branch secondary-edge request → M84 (blocked on her). Her publication-figure and near-redundant-band feedback needed no action. The rotation question went to RB02/RR02 and produced D-034; the rows below are the residue, all gated on the same design session as [H] above.

- [high] **D-032's premise contradicted by its source.** D-032 (2026-07-24) rejected gap-tolerant redundancy chains partly on the inference that Forbes's contiguous `ChaseCorrPaths` was deliberate. She confirms the empirical finding and denies the inference — the contiguity is a coding limitation, and she handles a dead intervening level by hand at the artefact stage, still on the redundancy criterion. M53's 54/54 reproduction and M78's `g2` regression test stand; only the reading of intent falls. **Gated:** a supersede takes the call. — added 2026-07-30
- **Promote Tucker's φ to the documented rotation-consistency diagnostic (RR02 rec 5, consider).** RR02 Q3 finds Forbes's "between-level correlations tell us whether the rotations are consistent/robust" reading is in role-conflict with reading the same edges as structure, and that the clean resolution is giving the robustness duty to Tucker's φ — which `prune()` already computes and reports. Docs-only if adopted. **Gated:** design-session agenda item. — added 2026-07-30
- **D-017's retention rule vs her published Figure 6B.** For the AMH chain `E1-F1-G1-H1-I1-J1`, D-017's retention (keep the bottom when the chain reaches `k_max`, else the topmost) retains only `J1`; Figure 6B retains `E1` **and** `J1`. Measured at M84 against her own data with her node set supplied by hand; her 6A green/grey colouring does not map 1:1 onto 6B's node set either. The chase itself still reproduces 54/54 (M53) — the divergence is in the *retention* step, which is the package's construction on top of her chase and not something her code specifies. **Gated:** what her retention rule actually is, is a question for her — take it to the design session. — added 2026-07-30
- **Unbroken chase to level a.** Her `ChaseCorrPaths()` returns "X--null" for a level-3+ component whose chase is unbroken to level a. The cause is `which.min()` on a vector with no `FALSE`, which counts zero links. `prune("redundant")` chases such a component to a1. Found by M90's oblique fixture (sim 3, minres, promax) and recorded in `references/source-departures.md` row E6 (corrected M90: was row M2). **Gated:** confirm her intent at the design session, then record the departure's D-entry (IP9) or change `prune()`. — added 2026-09-30
