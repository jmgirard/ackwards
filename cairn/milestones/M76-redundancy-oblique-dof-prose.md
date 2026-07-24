# M76: Clarify redundancy-criterion + oblique-rotation prose and reframe artifact-mode DoF wording

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** m76-redundancy-oblique-dof-prose

## Goal

Rewrite the `redundancy_criterion` and oblique-rotation explanations with a worked micro-example, and reframe the artifact-mode "investigator degrees of freedom" wording per Forbes's own correction (candidates A + B).

## Scope

**In:**
- Rewrite the `redundancy_criterion = "direct"` explanation (the `@details` block at `prune.R:618-627` and the matching redundancy prose in `vignettes/ackwards-forbes.Rmd.orig`) with a worked micro-example that states the **star/anchor-on-leaf** semantics — `direct` anchors on the deepest leaf and requires each ancestor to correlate `|r| >= redundancy_r` *directly with that leaf* — and explicitly denies both the every-pair-with-every-pair reading and adjacent-hop chaining.
- Add a paragraph making explicit **why oblique rotation confounds between-level signal**: varimax gives `T' = T^-1`, so the `W'RW` edge algebra is exact and all cross-level correlation is read cleanly between levels; oblique's within-level `Phi != I` leaks into every between-level edge, confounding the within-vs-between question the method exists to answer (DESIGN §5.1, §9, D-002).
- Reframe the artifact-mode DoF paragraph (`ackwards-forbes.Rmd.orig:307-313`): researcher *decisions* create degrees of freedom and standardized *flags* remove them, so the package declines to *auto-drop* because the drop *decision* is substantive — not because flagging manufactures DoF.
- Record the corrected framing as a row in the M72 departures ledger (`cairn/references/source-departures.md`).
- Re-precompute the touched vignette(s).

**Out:** the near-redundant band *feature* → M77; the skip-level chain *behavior change* → M78. Both are behavior, not prose.

## Acceptance criteria

- [ ] AC1: The `redundancy_criterion` `@details` in `prune.R` and the matching forbes-vignette prose carry a worked micro-example stating the star/anchor-on-leaf rule (each ancestor `|r| >= redundancy_r` direct-to-leaf) and explicitly denying both all-pairs and adjacent-hop readings.
- [ ] AC2: The oblique-rotation paragraph names the `T'=T^-1` (varimax) vs `Phi != I` (oblique) mechanism and cites DESIGN §5.1/§9/D-002.
- [ ] AC3: The artifact-mode DoF prose no longer states that flagging manufactures investigator DoF; it states the auto-drop *decision* is substantive.
- [ ] AC4: `cairn/references/source-departures.md` carries a ledger row recording the corrected DoF framing.
- [ ] AC5: Touched `*.Rmd.orig` re-precomputed; vignette-freshness stamp check passes; `devtools::test()` clean; `devtools::check()` clean (0/0, NOTEs justified).

## Coverage

- AC1 → T1, T2
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T5, T6

## Tasks

- [ ] T1: Rewrite the `redundancy_criterion` `@details` (`prune.R:618-627`) with the worked micro-example; `devtools::document()`.
- [ ] T2: Mirror the rewritten explanation into the redundancy section of `ackwards-forbes.Rmd.orig`; add the oblique-rotation "why" paragraph (DESIGN §5.1/§9/D-002).
- [ ] T3: Reframe the investigator-DoF paragraph (`ackwards-forbes.Rmd.orig:307-313`).
- [ ] T4: Add the ledger row to `cairn/references/source-departures.md`.
- [ ] T5: Re-run `Rscript vignettes/precompute.R`; `git checkout --` untouched `vignettes/assets`; confirm freshness stamps (`tools/check-vignette-freshness.R`).
- [ ] T6: NEWS.md entry (doc clarifications); `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: in-progress; branch m76-redundancy-oblique-dof-prose cut from master.
- 2026-07-23: T1/T2 investigation found the oblique "why" already lives in ackwards-intro.Rmd.orig:123-132 + R/ackwards.R:11-17, both overstating orthogonality's role (claiming it makes W'RW "exact", which DESIGN §5.1 says holds for any linear scoring). Proposed correcting both in place (a substantive scope amendment). At the amendment gate the owner chose to escalate the math claim (algebra exact under oblique; orthogonality interpretive not numerical, per D-002) to Fable via /milestone-brief before rewriting. Paused pending RR.

## Decisions

## Review
