# Roadmap

_The only authority on milestone status. Grouped by status, not ID._
_Last hygiene check: 2026-07-23 (M75 merged + archived; M70 row pruned under terminal-row retention; Forbes website-review feedback batch A–H added as candidates)_

Pre-migration history: see `cairn/legacy/` (MILESTONES.md, ROADMAP.md, skills)
and git log. Milestone IDs run through M53; new work continues from M54.

## Milestones

| ID | Title | Status | Depends on | Priority | File/Archive |
|---|---|---|---|---|---|
| M71 | Author + verify the 5 application source notes as citation precedents + Forbes drift-watch | done | M70 | normal | milestones/archive/M71-reference-notes-applications.md |
| M72 | Correct the stale Forbes contract line + author the Goldberg/Forbes departures ledger | done | — | normal | milestones/archive/M72-source-departures-ledger.md |
| M73 | Draft the manuscript Introduction with the verified framing + application sources | done | — | normal | milestones/archive/M73-manuscript-introduction.md |
| M74 | Draft the manuscript Discussion + enrich seeded sections with method-backer citations | done | M73 | normal | milestones/archive/M74-manuscript-discussion-citations.md |
| M75 | Beautify suggest_k() print output as an aligned criteria table | done | — | normal | milestones/archive/M75-suggest-k-print-table.md |
<!-- M01–M70 done/dropped (entombed in cairn/legacy/MILESTONES.md + milestones/archive/); terminal-row retention keeps the 5 most recent. -->

## Candidates

- Owner-only post-M55 release tail: 0.1.1 **resubmitted 2026-07-17** (tarball from master `0a5da58`, after the 2026-07-16 withdrawal; CRAN-SUBMISSION committed) — **owner: confirm the CRAN submission email if pending; on acceptance tag `v0.1.1` + update README "on CRAN" phrasing; if CRAN bounces again, plan the next resubmission milestone** — added 2026-07-12, updated 2026-07-17
- ESEM engine/basis extensions (grouped, demand-gated — keep off schedule until asked): `comparability()` split-half per level per factor (feasible; 2·n_splits lavaan hierarchies per call, per-half convergence handling — D-022 / M46) and `boot_edges()` WLSMV/polychoric bootstrap edges (expensive, n_boot × (k_max−1) fits; resample can drop a response category — D-023 / M47) — added 2026-07-11, merged 2026-07-16

### Forbes website-review feedback (2026-07-23)

Batch from Forbes's hands-on review of the package website/vignettes. Grouped; register/plan individually as picked up.

- **[A] doc — clarify redundancy-criterion + oblique-rotation prose.** Rewrite the `redundancy_criterion` explanation (Forbes found the "direct/non-transitive" paragraph hard to follow) with a worked micro-example, and make explicit *why* oblique rotation confounds between-level signal: varimax gives `T'=T⁻¹` so the `W'RW` edge algebra is exact and all cross-level correlation is read cleanly between levels; oblique's within-level `Φ≠I` leaks into every between-level edge, confounding the within-vs-between question that is the method's point (DESIGN §5.1, §9, D-002). Also clarify the default chain semantics: `direct` anchors on the deepest leaf and requires each ancestor to correlate |r|≥.9 *directly with the leaf* (star pattern), NOT every-pair-with-every-pair and NOT adjacent-hop chaining (`prune.R:618-627`, `.strong_links_direct`). Ready. — added 2026-07-23
- **[B] doc — reframe "investigator degrees of freedom" wording.** Forbes (the source author) says the current artifact-mode framing is backwards: researcher *decisions* create DoF; standardized automated *flags* remove them (`ackwards-forbes.Rmd.orig:307-313`). Soften to: package declines to *auto-drop* because the drop *decision* is substantive, not because flagging adds DoF. Note in the M72 departures ledger. Ready. — added 2026-07-23
- **[C] behavior+design — gap-tolerant (skip-level) direct redundancy chains.** Forbes: a chain should skip a whole intermediate artifact level (e.g. b2→d2→e2→f2 skipping c). Current `.strong_links_direct` (`prune.R:152-197`) uses a `for(j in L-1:1)` loop with `break`, so a sub-threshold level *terminates* the chain even when the leaf still correlates ≥.9 with a grandparent directly. Change `break`→skip-and-continue to be truly direct/skip-tolerant. **Gated:** confirm her `ChaseCorrPaths` semantics (contiguous or not) and re-verify the M44/M53 fidelity suite (AMH 54/54) before changing; needs a D-entry (redundancy-semantics change). — added 2026-07-23
- **[D] feature+doc — "near-redundant" band flag in artifact mode.** Forbes uses φ/r flags mainly for the *near-redundant* band (just below thresholds, e.g. r=.89/φ=.93 — a messy re-rotation), since full redundancy is already dropped in the redundancy stage. Add a near-redundant band signal (within a margin of `redundancy_r`/`redundancy_phi`) to `prune("artifact")`, and fix the artifact vignette example (φ=.99 pair would already be flagged by `prune("redundant")` — poor illustration). Ready. — added 2026-07-23
- **[E] feature — secondary-correlation edges in the pruned/publication view.** After pruning, plot secondary edges |r|≥.3 that are NOT already implied by the primary hierarchy path (e.g. f1→d2 when the primary path is f1→d1→b1) as dashed lines. `drop_pruned=TRUE` currently draws primary-parent edges only (`layout.R:157-173`, `autoplot.R:398-400`); note "dashed" already encodes *sign*, so this needs a distinct visual channel. Ready. — added 2026-07-23
- **[F] feature — publication-figure polish.** Three asks: (a) list the top items under the lowest-level factors; (b) adjustable box sizes to fit manual node labels; (c) manual factor ordering per level to untangle/arrange the plot (override `ba_layout`'s ordinal pass). Ready; may split. — added 2026-07-23
- **[G] bug/robustness — deep-hierarchy layout at k≥10 (PARKED).** Forbes: "levels got a bit bent" at ten levels. `ba_layout` puts each level on a clean horizontal row (`y=-level`), so this is likely edge-tangle/crowding or aspect-ratio squish, not literal non-horizontal levels. **Gated:** need her k=10 screenshot to tell layout bug from density before scoping. — added 2026-07-23
- **[H] collaboration — replicability-gated hierarchies (PARKED).** Forbes offered to co-develop this. Overlaps existing `comparability()` (split-half per level) + `boot_edges()`. **Gated:** design-interview territory with Forbes in the room — do not spec unilaterally; schedule a design session before planning. — added 2026-07-23
