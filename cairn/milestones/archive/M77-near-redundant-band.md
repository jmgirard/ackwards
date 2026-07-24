# M77: Near-redundant band flag in artifact mode + fix the artifact vignette example

**Status:** done (2026-07-24, PR #82 https://github.com/jmgirard/ackwards/pull/82)

**Goal:** Add a report-only "near-redundant" band signal to `prune("artifact")` flagging factor pairs just below the redundancy thresholds, and fix the artifact-mode vignette so its illustrated pair is genuinely near-redundant.

**Outcome:** `prune(x, "artifact")` now populates `x$prune$near_redundant` (new internal `.near_redundant_pairs()`): cross-level pairs that are not fully redundant yet have direct `|r|` **or** Tucker `phi` within the new `near_margin` argument (default `0.1`) below `redundancy_r`/`redundancy_phi`. Columns `from,to,level_from,level_to,r,phi,near_r,near_phi`. `redundancy_phi` auto-resolve (0.95, EFA/ESEM) now also fires in artifact mode, cli-announced (IP6). Report-only — never mutates the kept set (GP2). The Forbes-extension vignette's artifact example now features the genuine near-miss m1f1↔m2f1 (r=.89/φ=.94), replacing the top-φ pair m3f2↔m5f2 (φ=.987) that `prune("redundant")` already drops.

**Decisions:** M77-D1 (milestone-local): band = OR of the just-below `|r|`/φ windows; "not already flagged" judged at the pair level (not the chain/retention node flags, since correlation is non-transitive); `near_margin` fixed at 0.1 (not engine-dependent); PCA (`redundancy_phi` NULL) → `|r|`-band only.

**Review:** Three-lens fan-out (diff-bug/Opus, blame-history/Sonnet, prior-review/Sonnet) — zero surviving findings; scorer/triage no-oped. Surfaced out-of-scope (IP3, not actioned): `print()`/`summary()` don't yet display the band count. check() 0/0/0, coverage 100%, vignette freshness OK, CI green.
