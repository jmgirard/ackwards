# M72: Correct the stale Forbes contract line + author the Goldberg/Forbes departures ledger

**Status:** done (2026-07-23, PR #75 https://github.com/jmgirard/ackwards/pull/75)

**Goal:** Correct the pre-D-031 contract wording in forbes2023.md and author a consolidated cairn/references/ synthesis note cataloguing every ackwards departure from Goldberg (2006) / Forbes (2023) with rationale + empirical support, auditable going forward.

**Outcome:** forbes2023.md's "Why this is the primary source" corrected from the superseded "default output must reproduce Forbes's examples exactly" to the IP9/D-031 capability framing. New `cairn/references/source-departures.md`: 5 default-level departures (E1 required `k_max` vs auto-stop [tong2025 58%/71%, p.14]; E2 W′RW algebra vs score-then-correlate [Waller + IP2]; E3 tenBerge vs components [D-007]; E4 primary-parent sign alignment vs unaligned `comp.corr` [1.3e-14]; E5 exact Tucker φ vs rounded congruence [0.005]) + 2 source-matches (`cut_show=.30`; `redundancy_criterion="direct"`), each with rationale + support; additive extensions noted GP1-governed. Maintenance hook keeping it current: Maintenance clause (new departure ⇒ D-entry per IP9 ⇒ ledger row, same change) + DESIGN §9 pointer + `tools/check-ledger-anchors.R` (anchor-integrity, wired into dod-gate + CI, testthat wrapper skips in built pkg). tong2025.md gained the 58%/71% figure at p.14.

**Decisions:** none promoted (works under IP9/GP1; no principle text changed; no NEWS — cairn/ internal + a buildignored guard).

**Review:** 7/7 ACs fresh-evidenced; DoD gate 0/0/0 + 100% cov + `ledger-anchors: clean`; cairn_validate 0. Two passes: pass-1 3-lens found 2 citation-hygiene fixes (E3 "see M1", M1 dangling anchor) + a tong2025 p.14 traceability fix (rendering corrected the reviewer's stated p.11); pass-2 [O] delta-review caught the checker's CI/dod-gate wiring gap (docstrings claimed enforcement not actually wired) — fixed. F3 (E3 categorization, 35) logged/not-actioned.
