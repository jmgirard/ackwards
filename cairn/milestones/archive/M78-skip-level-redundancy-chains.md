# M78: Gap-tolerant (skip-level) direct redundancy chains — gated on ChaseCorrPaths semantics

**Status:** done (2026-07-24, PR #83 https://github.com/jmgirard/ackwards/pull/83)

**Goal:** Determine whether Forbes's `ChaseCorrPaths` is contiguous or gap-tolerant, and change the direct redundancy criterion only if gap-tolerant.

**Outcome:** Determined empirically, from the committed AMH fixture alone (her `comp_corr` + `corr_chase` endpoints — no `ackwards()` fit, no network), that `ChaseCorrPaths` is **contiguous**: a contiguous direct-chase reproduces her 54 endpoints 54/54, and the single AMH component where the two semantics diverge (`g2`, a level-7 mid-chain gap — direct `|r|` to level 6 = 0.8981 < .9, re-emerging 0.9084 at level 5) reports `g2--null` (no move); 0/54 match gap-tolerant. The existing `.strong_links_direct` break-at-first-sub-threshold-hop (`R/prune.R`) is therefore already faithful — no code change. Added a planted-gap regression test in `test-prune.R` ("M78: direct chase stops at a mid-chain gap") locking the contiguous no-skip semantics, with a positive control that flips NULL→top-reaching chain when only the gap edge is filled. No behavior change, no NEWS entry.

**Decisions:** D-032 (extends D-017) — `ChaseCorrPaths` is contiguous; the gap-tolerant variant is rejected; the current `redundancy_criterion = "direct"` default is confirmed faithful and test-guarded.

**Review:** Three fresh-context lenses (diff-bug/Opus, blame-history/Sonnet, prior-review/Sonnet) — zero findings; one benign non-actionable aside (cosmetic `E_1_2` fill in the positive control). Scorer no-op. dod-gate 0/0/0, coverage 100%, `cairn_validate` clean.
