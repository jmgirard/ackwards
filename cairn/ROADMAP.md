# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-27 (0.2.0 release tail closed out — tag, GitHub release, dev-version bump, `CRAN-SUBMISSION` marker removed — and its candidate row retired; the 91 dangling-id-token advisories are all pre-migration M≤53 ids entombed in `cairn/legacy/`)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M80 | Deep-hierarchy layout quality at k>=10 — crossing reduction + edge-label dodging | done | — | normal | milestones/archive/M80-deep-hierarchy-layout-quality.md |
| M81 | Publication-figure polish — item lists, per-node box sizes, manual factor ordering | done | M80 | normal | milestones/archive/M81-publication-figure-polish.md |
| M82 | macOS oldrel (arm64) CI visibility + skip-filter repair | done | — | normal | milestones/archive/M82-macos-oldrel-ci-visibility.md |
| M83 | Windows parallel-testthat crash — measure, diagnose, mitigate | done | — | high | milestones/archive/M83-windows-parallel-testthat-crash.md |
| M84 | Cross-branch-only secondary edges in the pruned view | blocked | — | normal | milestones/M84-cross-branch-secondary-edges.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent done rows. -->

## Candidates

- Further CRAN-flavour coverage beyond M82 (grouped): a standing R-hub macOS run in `PROFILE.md`'s release-walk (M82 declined it as subsumed by its own push-to-master job) and a macOS **x86_64** row (`macos-13`) — the 0.1.1 failure was arm64-specific and `r-oldrel-macos-x86_64` passed. Promote either on any CRAN macOS failure M82's job could not have caught. **Promoted-by-evidence 2026-07-27: the alternative-numerics gap fired** — 0.2.0 was auto-rejected on a noLD ERROR no green macOS/Windows/Ubuntu matrix could have caught, because the defect needed degenerate rather than merely ill-conditioned arithmetic. Add R-hub `atlas` + `nold` to `PROFILE.md`'s release-walk as a pre-submission step (the `rhub.yaml` workflow already exists; `rhub::rhub_check(platforms = c("atlas", "nold"))` is the invocation). — added 2026-07-27, updated 2026-07-27
- Untested axis from M83: **disabling testthat parallelism outright** (`Config/testthat/parallel: false` / `TESTTHAT_PARALLEL=false`), as distinct from the worker counts M83 measured. Evidence it matters: at `TESTTHAT_CPUS=1` the crashes still read `testthat subprocess exited`, so the parallel machinery still spawned a subprocess — the crashing component was never removed in any sweep. Plausibly the actual fix, and the only candidate that eliminates rather than reduces. Costs M48's speedup (27s parallel vs 81s serial locally) and, unlike M83's mitigation, changes **tarball content** (DESCRIPTION), so it needs release re-verification. Promote if a CRAN Windows flavour fails with the -1073741819 signature. — added 2026-07-27
- Upstream report for the M83 Windows crash, if it isolates to a dependency's compiled code (`EFAtools` or `mnormt`; `psych`/`lavaan`/`GPArotation` are pure R): file with a minimal reproducer. Held out of M83 at the 2026-07-27 plan gate because a maintainer's timeline must not gate our CI health. Promote when M83's bisection names a component. — added 2026-07-27
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. **A, B → M76; D → M77; C → M78; E → M79; G → M80; F → M81 (all done).** H remains below.

- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23

### Forbes correspondence (2026-07-30)

Batch from her reply to the M76–M81 write-up. The cross-branch secondary-edge request → M84 (blocked on her). Her publication-figure and near-redundant-band feedback needed no action. The rotation question went to RB02/RR02 and produced D-034; the rows below are the residue, all gated on the same design session as [H] above.

- **Oblique rotation support (D-034; implementation gated).** D-002 superseded 2026-07-30 after RB02/RR02: varimax stays the sole default and no current output changes, but oblique is no longer out of scope — it becomes a documented non-default option. **Gated:** implementation waits on the Forbes design session settling the marginal-vs-partialled lineage semantics — which quantity drives primary-parent matching, sign anchoring, and the diagram (RR02 Q5 item 7, Q6). **RR02's BC1–BC8 bind the future milestone** and are ingested verbatim into its acceptance criteria when planned, after the fresh-context [O] criteria audit the ingest protocol requires (deferred with them — no milestone existed at ingest time). RR02 rules out plumbing-only enablement: `factor_cor`, oblique tenBerge scoring, `.variance_explained()`, and the varimax-calibrated `|r| >= .9` chase would go silently wrong rather than error. — added 2026-07-30
- **Φ-partialled edge decomposition as reporting (RR02 rec 4, consider).** Ship `B = Φ_s^-1 E` and per-node `R² = E' Φ_s^-1 E` as output even without oblique rotation: under varimax it is an identity check (`B ≡ E`), and shipping it first makes the eventual oblique semantics legible rather than novel. Closed-form from `W`, `R`, and `factor_cor` — all already stored, so no new object state (RR02 Q2). Related finding: `factor_cor` is currently produced by every engine and consumed by nothing (RR02 Beyond-the-brief 4). — added 2026-07-30
- **Promote Tucker's φ to the documented rotation-consistency diagnostic (RR02 rec 5, consider).** RR02 Q3 finds Forbes's "between-level correlations tell us whether the rotations are consistent/robust" reading is in role-conflict with reading the same edges as structure, and that the clean resolution is giving the robustness duty to Tucker's φ — which `prune()` already computes and reports. Docs-only if adopted. **Gated:** design-session agenda item. — added 2026-07-30
- **D-017's retention rule vs her published Figure 6B.** For the AMH chain `E1-F1-G1-H1-I1-J1`, D-017's retention (keep the bottom when the chain reaches `k_max`, else the topmost) retains only `J1`; Figure 6B retains `E1` **and** `J1`. Measured at M84 against her own data with her node set supplied by hand; her 6A green/grey colouring does not map 1:1 onto 6B's node set either. The chase itself still reproduces 54/54 (M53) — the divergence is in the *retention* step, which is the package's construction on top of her chase and not something her code specifies. **Gated:** what her retention rule actually is, is a question for her — take it to the design session. — added 2026-07-30
- **D-032's premise contradicted by its source.** D-032 (2026-07-24) rejected gap-tolerant redundancy chains partly on the inference that Forbes's contiguous `ChaseCorrPaths` was deliberate. She confirms the empirical finding and denies the inference — the contiguity is a coding limitation, and she handles a dead intervening level by hand at the artefact stage, still on the redundancy criterion. M53's 54/54 reproduction and M78's `g2` regression test stand; only the reading of intent falls. **Gated:** a supersede takes the call. — added 2026-07-30
