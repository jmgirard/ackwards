# M76: Clarify redundancy-criterion + oblique-rotation prose and reframe artifact-mode DoF wording

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** RR01
- **Principles touched:** IP1
- **Branch/PR:** m76-redundancy-oblique-dof-prose

## Goal

Rewrite the `redundancy_criterion` explanation with a worked micro-example, correct the oblique-rotation rationale (the `W'RW` algebra is exact for any linear scoring — orthogonality is interpretive, not numerical; RR01), and reframe the artifact-mode "investigator degrees of freedom" wording per Forbes's own correction (candidates A + B).

## Scope

**In:**
- Rewrite the `redundancy_criterion = "direct"` explanation (`prune.R:618-627` `@details` + the matching redundancy prose in `vignettes/ackwards-forbes.Rmd.orig`) with a worked micro-example stating the **star/anchor-on-leaf** semantics — each ancestor correlates `|r| >= redundancy_r` *directly with the deepest leaf*, not every-pair and not adjacent-hop.
- **Correct the oblique-rotation rationale at all four sites RR01 identified** (the `W'RW` identity is exact for any fixed linear scoring, oblique included; varimax is chosen for the interpretive Φ = I reason, not a numerical one): the intro-vignette aside (`ackwards-intro.Rmd.orig:123-132`), the roxygen `rotation` bullet (`R/ackwards.R:11-17`), the DESIGN §9 `rotation` row, and D-002's Context line (`cairn/DECISIONS.md`, append-only M43-style parenthetical). Apply RR01's proposed replacements + BC1–BC7.
- Reframe the artifact-mode DoF paragraph (`ackwards-forbes.Rmd.orig:307-313`): researcher *decisions* create DoF and standardized *flags* remove them, so the package declines to *auto-drop* because the drop *decision* is substantive.
- Record the corrected DoF framing as a `match` row in the M72 departures ledger (`cairn/references/source-departures.md`).
- Re-precompute the touched vignette(s).

**Out:** the near-redundant band *feature* → M77; the skip-level chain *behavior change* → M78. Changing D-002 itself (supporting oblique) is explicitly rejected — RR01 corrects only the stated rationale.

## Acceptance criteria

- [ ] AC1: The `redundancy_criterion` `@details` in `prune.R` and the matching forbes-vignette prose carry a worked micro-example stating the star/anchor-on-leaf rule (each ancestor `|r| >= redundancy_r` direct-to-leaf) and explicitly denying both all-pairs and adjacent-hop readings.
- [ ] AC2: The artifact-mode DoF prose no longer states that flagging manufactures investigator DoF; it states the auto-drop *decision* is substantive.
- [ ] AC3: `cairn/references/source-departures.md` carries a ledger row recording the corrected DoF framing.
- [ ] AC4: Touched `*.Rmd.orig` re-precomputed; vignette-freshness stamp check passes; `devtools::document()` no diff; `devtools::test()` clean; `devtools::check()` clean (0/0, NOTEs justified).
- [ ] AC5 (BC1): The rewritten vignette aside contains no claim that orthogonality (or `T' = T^{-1}`) is necessary for, enables, or is what makes the `W'RW` closed form exact; it states affirmatively that the identity is exact for any (fixed) linear scoring, oblique included.
- [ ] AC6 (BC2): The rewritten roxygen `rotation` bullet contains no "enables the closed-form" (or equivalent necessity) claim about orthogonality, and states the interpretive Φ = I rationale (within-level factors uncorrelated, so between-level edges are not confounded by within-level factor intercorrelation).
- [ ] AC7 (BC3): Both rewritten texts retain the substantive point that an oblique rotation's within-level factor correlation would contaminate/confound the between-level edges — the interpretive rationale must not be weakened while removing the numerical one.
- [ ] AC8 (BC4): Neither rewritten text suggests oblique rotation is supported, supportable, or desirable (D-002 unchanged).
- [ ] AC9 (BC5): Uncorrelatedness claims in the rewritten texts are attributed to the factors (Φ = I), not asserted unconditionally of estimated factor scores.
- [ ] AC10 (BC6): The rewrite does not remove or contradict (i) DESIGN §5.1's "holds for any linear W" and "standardization is the trap" notes, or (ii) the `forbes2023.md` correspondence note that `W'RW = I` for varimax PCA is what makes Forbes's unstandardized `comp.corr` equal the standardized edges.
- [ ] AC11 (BC7): The DESIGN §9 `rotation` row no longer asserts that `T' = T^{-1}` enables the closed-form algebra; the corrected row carries an M43-style parenthetical noting the wording correction and this review (RR01) as its source.

## Coverage

- AC1 → T1, T2
- AC2 → T7
- AC3 → T8
- AC4 → T9, T10
- AC5 → T3
- AC6 → T4
- AC7 → T3, T4
- AC8 → T3, T4
- AC9 → T3, T4
- AC10 → T3, T4, T5
- AC11 → T5

## Tasks

- [x] T1: Rewrite the `redundancy_criterion` `@details` (`prune.R:618-627`) with the worked micro-example; `devtools::document()`.
- [x] T2: Mirror the rewritten redundancy explanation into the redundancy section of `ackwards-forbes.Rmd.orig`.
- [x] T3: Rewrite the oblique "why" aside in `ackwards-intro.Rmd.orig:123-132` per RR01's proposed intro replacement (BC1/BC3/BC4/BC5); preserve the surrounding Goldberg/CF-VARIMAX/confound sentences.
- [x] T4: Rewrite the roxygen `rotation` bullet (`R/ackwards.R:11-17`) per RR01's proposed roxygen replacement (BC2/BC3/BC4/BC5); keep the literature-match + "only supported rotation" sentences; `devtools::document()`.
- [x] T5: Correct the DESIGN §9 `rotation` row with an M43-style parenthetical citing RR01 (BC7); confirm §5.1 + `forbes2023.md` correspondence notes untouched (BC6).
- [x] T6: Append a correction parenthetical to D-002's Context in `cairn/DECISIONS.md` (append-only, per RR01 rec 4).
- [x] T7: Reframe the investigator-DoF paragraph (`ackwards-forbes.Rmd.orig:307-313`).
- [x] T8: Add the `match` ledger row to `cairn/references/source-departures.md`.
- [x] T9: Re-run `Rscript vignettes/precompute.R`; `git checkout --` untouched `vignettes/assets`; confirm freshness stamps (`tools/check-vignette-freshness.R`).
- [ ] T10: NEWS.md entry (doc clarifications); `Rscript tools/dod-gate.R`.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: in-progress; branch m76-redundancy-oblique-dof-prose cut from master.
- 2026-07-23: T1/T2 investigation found the oblique "why" already lives in ackwards-intro.Rmd.orig:123-132 + R/ackwards.R:11-17, both overstating orthogonality's role (claiming it makes W'RW "exact", which DESIGN §5.1 says holds for any linear scoring). Proposed correcting both in place (a substantive scope amendment). At the amendment gate the owner chose to escalate the math claim (algebra exact under oblique; orthogonality interpretive not numerical, per D-002) to Fable via /milestone-brief before rewriting. Paused pending RR.
- 2026-07-23: blocked on RB01 (oblique-rotation algebra claim).
- 2026-07-23: RR01 ingested — concern confirmed (identity exact for any fixed linear scoring; Waller 2007 §3 gives the oblique closed form; orthogonality buys D=I convenience only). Scope broadened to correct 4 sites (intro vignette, roxygen, DESIGN §9 rotation row, D-002 Context); BC1–BC7 ingested verbatim as AC5–AC11 (Driving RR: RR01, no deviations). Beyond-brief triage: DESIGN §9 = apply (BC7); D-002 Context parenthetical upgraded consider→apply (append-only, completes the 4-site correction, T6). Status back to in-progress.
- 2026-07-23: T1/T2 done — star/anchor-on-leaf worked micro-example added to prune.R @details + forbes vignette redundancy prose (document()/precompute deferred to T9/T10).
- 2026-07-23: T3–T6 done — oblique rationale corrected at all four sites (intro-vignette aside, roxygen `rotation` bullet, DESIGN §9 rotation row w/ M43-style parenthetical, D-002 Context parenthetical). BC6 preserved: §5.1 + forbes2023.md correspondence untouched. document()/precompute deferred to T9.
- 2026-07-23: T7/T8 done — DoF paragraph reframed (flagging removes DoF; the auto-drop decision is what's substantive); `match` row M3 added to source-departures.md ledger + Disposition line updated.
- 2026-07-23: T9 done — document() regenerated man/ackwards.Rd + man/prune.Rd (no NAMESPACE change); precompute regenerated vignettes. Reverted the 4 untouched vignettes + 3 asset PNGs (run noise); restored intro.Rmd from master and re-applied stamp + prose only, dropping suggest_k timing churn AND the ✔→✔︎ tick-glyph drift that #80/2a0122d fixed (M75 lesson). Freshness check passes.

## Decisions

- 2026-07-23 (RR01): The `W'RW` between-level correlation identity is algebraically exact for **any fixed linear scoring map** `S = ZW` (regression/Bartlett/tenBerge, orthogonal *or* oblique) — no property of `W` enters the derivation. Orthogonality (varimax) is an **interpretive** choice: it makes within-level factors uncorrelated (Φ = I) so between-level edges are not confounded by within-level factor intercorrelation; it is *not* what makes the algebra exact. The only numerical thing orthogonality buys is `W'RW = I` for varimax PCA (a standardization convenience `compute_edges()` never relies on, since it divides by real score SDs unconditionally) — which is load-bearing only for the Forbes unstandardized-`comp.corr` fidelity correspondence. D-002 (varimax only) is **unchanged**; only its stated rationale is corrected. Applies to intro vignette, roxygen `rotation` bullet, DESIGN §9 rotation row, and D-002 Context.
- 2026-07-23 (RR01 Q1/Q5): The genuine exact-vs-approximate distinction is the component-vs-estimated-factor-score (determinacy; Grice 2001) axis, rotation-independent and engine-level — the same conflation M43 already corrected in the `redundancy_phi` row. The rewrite must NOT move that caveat onto the orthogonal/oblique axis inside these rotation passages (RR01 rec 6, reject).

## Review
