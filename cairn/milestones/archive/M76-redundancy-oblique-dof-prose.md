# M76: Clarify redundancy-criterion + oblique-rotation prose and reframe artifact-mode DoF wording

**Status:** done (2026-07-23, PR #81 https://github.com/jmgirard/ackwards/pull/81)

**Goal:** Rewrite the `redundancy_criterion` explanation with a worked micro-example, correct the oblique-rotation rationale (the `W'RW` algebra is exact for any linear scoring — orthogonality is interpretive, not numerical), and reframe the artifact-mode DoF wording per Forbes's own correction.

**Outcome:** Docs-only, no runtime change. (1) `?prune` `@details` + bass-ackward vignette now describe `redundancy_criterion = "direct"` as a *star anchored on the chain's deepest leaf* (each ancestor correlates directly with the leaf), not adjacent-hop and not all-pairs. (2) The oblique-rotation rationale corrected at six sites — intro vignette, `?ackwards` roxygen, DESIGN §4 + §9 rotation row, D-002 Context: the `W'RW` identity is exact for any fixed linear scoring (oblique included), so varimax is the interpretive Φ=I choice, not a numerical necessity; D-002 (varimax only) unchanged. (3) Artifact-mode DoF wording reframed (standardized flags *remove* investigator DoF; the drop *decision* is substantive) + `match` row M3 in source-departures ledger.

**Decisions:** Two milestone-local (both under RR01): the `W'RW` identity is exact for any fixed linear scoring, orthogonality interpretive not numerical; the exact-vs-approximate axis is component-vs-estimated-score determinacy (rotation-independent), never the orthogonal/oblique axis. No cross-cutting D-entry (D-002 unchanged, only its rationale corrected via M43-style parenthetical).

**Review:** Escalated the math claim to Fable (RB01→RR01, archived); ingested BC1–BC7 verbatim as ACs. Three-lens fan-out: prior-review clean; blame-history + diff-bug each found one missed conflation copy (F1 second forbes-vignette DoF paragraph, score 93; F2 fifth site DESIGN §4, score 92) — both fixed. No sub-80 findings.
