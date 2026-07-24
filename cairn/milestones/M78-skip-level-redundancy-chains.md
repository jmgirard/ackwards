# M78: Gap-tolerant (skip-level) direct redundancy chains — gated on ChaseCorrPaths semantics

- **Status:** planned
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP8, IP9
- **Branch/PR:** —

## Goal

Determine whether Forbes's `ChaseCorrPaths` is contiguous or gap-tolerant, and — only if gap-tolerant — change `.strong_links_direct` to skip sub-threshold levels while keeping the AMH fidelity oracle exact (candidate C).

## Scope

**In:**
- **Investigation first (gate):** determine `ChaseCorrPaths`'s contiguity — does the chase *break* at the first sub-threshold level (as current code does), or *skip and continue* while the leaf still correlates `|r| >= redundancy_r` directly with a deeper ancestor? Her chase **output** already ships as a fixture (`amh$corr_chase`, 54 endpoints) and the current contiguous code reproduces it 54/54 — so settle contiguity **empirically first**: find any AMH component whose chain has a mid-chain gap where the two semantics diverge, and read her committed endpoint. Only if AMH is silent (no distinguishing gap) fetch her source from OSF `pcwm8` (only the AMH matrix + her chase output ship, not her code).
- **Conditional change:** if gap-tolerant, replace the `break` at `prune.R:164/167/173` in `.strong_links_direct` (`prune.R:152-197`) with skip-and-continue; keep consecutive emitted links one level apart so they slot into the shared chain machinery unchanged.
- Re-verify the M44/M53 fidelity suite (`test-forbes-fidelity.R`): AMH 54/54 components matched to its stated tolerance (1.3e-14). Add a regression test with a planted mid-chain gap where contiguous and skip-tolerant diverge.
- Append a D-entry extending/superseding **D-017** (redundancy chased by direct skip-level correlation) recording the semantics finding and the decision.

**Out:** near-redundant band → M77; prose → M76.

## Acceptance criteria

- [ ] AC1: `ChaseCorrPaths`'s contiguity is determined and recorded in the D-entry. Settle it empirically first — whether any AMH component's chain has a mid-chain gap where contiguous vs gap-tolerant chase diverge, read against her committed `corr_chase` endpoint (`fixtures/forbes2023_amh.rds`); cite the fixture. Fall back to the OSF `pcwm8` source (file + line) only if AMH is silent. (RB tripwire: ip-touching)
- [ ] AC2: If the finding is gap-tolerant, `.strong_links_direct` skips sub-threshold levels and continues while the leaf's direct `|r|` to a considered deeper ancestor `>= redundancy_r`; if contiguous, code is unchanged and the D-entry records the current behavior as already faithful.
- [ ] AC3: `test-forbes-fidelity.R` still passes — AMH 54/54 components matched to the fixture's stated tolerance (1.3e-14).
- [ ] AC4: A regression test exercises a planted hierarchy with a mid-chain gap where contiguous vs skip-tolerant chains diverge, asserting the chosen semantics.
- [ ] AC5: A D-entry extending/superseding D-017 is appended to `cairn/DECISIONS.md`.
- [ ] AC6: `devtools::test()` clean; `devtools::check()` clean.

## Coverage

- AC1 → T1
- AC2 → T2, T3
- AC3 → T4
- AC4 → T2
- AC5 → T5
- AC6 → T6

## Tasks

- [ ] T1: Determine `ChaseCorrPaths` contiguity — empirically first: for each AMH component walk our all-levels edges, detect any mid-chain gap where contiguous and gap-tolerant chase diverge, and compare her committed `corr_chase` endpoint to both (her endpoint reveals the semantics). Only if no AMH component distinguishes them, fetch `ChaseCorrPaths` from OSF `pcwm8`. Write up the finding with its citation (fixture, or source file+line). (RB tripwire: ip-touching)
- [ ] T2: (test-first, conditional on gap-tolerant) regression test with a planted mid-chain gap where the two semantics diverge.
- [ ] T3: If warranted, change `break`→skip-and-continue in `.strong_links_direct` (`prune.R:152-197`), preserving one-level-apart emitted links.
- [ ] T4: Re-run `test-forbes-fidelity.R`; confirm AMH 54/54 to tolerance.
- [ ] T5: Append the D-entry extending D-017.
- [ ] T6: NEWS.md entry if behavior changed; `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-24: /milestone-plan re-run — reframed T1/AC1 investigation to empirical-first (her ChaseCorrPaths output already ships as `amh$corr_chase`, reproduced 54/54 by the current contiguous code; settle contiguity from committed data, fetch OSF source only if AMH is silent). Owner-confirmed; docs-only.

## Decisions

## Review
