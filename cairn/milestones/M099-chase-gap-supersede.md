# M099: Chase-gap reason corrected and a manual bridge documented

- **Status:** review
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP9, GP2
- **Resolves:** —
- **Surface tier:** user-facing — it changes the `?prune` help page and NEWS.md
- **Branch/PR:** `m099-chase-gap-supersede`

## Goal

Replace D-032 with a decision whose reason matches what Forbes says about her chase, and show
users how to flag a factor across a chase gap by hand.

## Scope

**In:** one D-entry that supersedes D-032. The `redundancy_criterion = "direct"` chase stays
contiguous, with no output change. The `?prune` help explains the gap stop and the `manual =`
route. One AMH test covers both. The Forbes source note, the departures ledger, and the M78
lesson record her 2026-07-30 statement. One NEWS entry.

**Out:** a gap report in `prune("artifact")` and an opt-in bridged chase. The plan question
set declined both, and the new D-entry records them as alternatives. User docs do not name
Forbes's private email (question set). The `ackwards-forbes2023` vignette names her code
`ChaseCorrPaths()` correctly, so it stays unchanged. The retention divergence from Figure 6B
stays with the candidate row "D-017 retention vs Figure 6B". The `"X--null"` edge case stays
with the candidate row "Unbroken chase to level a".

## Acceptance criteria

- [ ] AC1: `cairn/DECISIONS.md` gains one entry whose heading says it supersedes D-032. The
      entry keeps `redundancy_criterion = "direct"` and its contiguous break unchanged. It
      states that the contiguity matches her reference code `ChaseCorrPaths`, not a rule of
      her method. It gives two reasons. First, the contiguous chase reproduces her AMH chase
      output 54/54 (D-017 from M53) and her `g2--null` gap case (D-032 from M78). Second, she bridges a gap by hand when she reviews
      artifacts. The nearest package step is `prune(manual = )`, which records reason
      `"manual"` and reports no chain.
- [ ] AC2: In the `?prune` help (roxygen in `R/prune.R`, regenerated into `man/prune.Rd`), the
      `"direct"` entry of `redundancy_criterion` says three things. First, the chase stops at
      the first ancestor level that fails the criterion. The criterion is `|r| >=
      redundancy_r`, plus phi above `redundancy_phi` when `redundancy_phi` is non-NULL,
      its automatic value off PCA included. A level above the stop
      that meets the criterion again does not restart the chase. Second, the contiguity
      follows Forbes's reference code `ChaseCorrPaths`. Third, `manual =` flags such a factor
      by the user's judgment, with reason `"manual"`.
- [ ] AC3: A test in `tests/testthat/test-prune.R` fits
      `ackwards(forbes2023, k_max = 10, pairs = "all")` and checks `m7f2`. Its largest direct
      `|r|` to a level-6 factor is below .9, and its largest direct `|r|` to a level-5 factor
      is .9 or more. After `prune(x, "redundant")`, `m7f2` has `pruned = FALSE` in
      `x$prune$nodes` and does not appear in the `id` column of `x$prune$chains`. After
      `prune(x, "redundant", manual = "m7f2")`, it has `pruned = TRUE` and
      `prune_reason = "manual"`.
- [ ] AC4: `cairn/references/forbes2023.md` records her 2026-07-30 statement on the chase. The
      contiguity is a coding limitation. She bridges a dead intervening level by hand when
      she reviews artifacts, on the redundancy criterion. The M2 row of
      `cairn/references/source-departures.md` cites that record. The M78 line of
      `cairn/LESSONS.md` is corrected in place and marked `corrected M099`.
- [ ] AC5: The development section of NEWS.md gains one documentation entry that describes
      the AC2 help text.

## Coverage

- AC1 → T3
- AC2 → T2
- AC3 → T1
- AC4 → T4
- AC5 → T5

## Tasks

- [x] T1: Write the AC3 test first, after the M78 test in `test-prune.R` (~line 422). Reuse
      the `cached()` AMH fit. It passes on current behavior. To show that it can fail, plant
      a gap skip in `.strong_links_direct` (`R/prune.R` ~line 232): a `next` in place of the
      `break` at the `|r|` check, plus a direct link from the found ancestor to the leaf
      (stats read from the `j:L` matrix). A bare `next` errors at `E_adj[Pnode, prev]`
      (~line 263). Confirm the "not in the `id` column" expectation itself fails,
      not an error, then revert.
- [x] T2: Edit the roxygen of `prune()`. The `rules` entry (~line 662) says "faithful to
      Forbes's reference code". The `"direct"` entry (~lines 688-693) gets the AC2 text in
      place of "Forbes's (2023) published `ChaseCorrPaths` rule". Make the same change to
      the internal comments at ~line 14 and ~line 221. Run `devtools::document()`.
- [x] T3: Write the superseding D-entry (next free id) to the AC1 wording. Name the declined
      alternatives (gap report, opt-in bridged chase). Name the class of evidence that would
      reopen them: Forbes asks for the bridge in code, or a user needs the hand bridge
      recorded as redundancy.
- [x] T4: Add the 2026-07-30 chase statement to the correspondence part of
      `cairn/references/forbes2023.md`, and cite it in the M2 row of `source-departures.md`.
      The M2 tag stays `match`, because the automated chase matches her code and the hand
      step is the user's judgment. Correct the M78 line of `cairn/LESSONS.md` in place: the
      committed output settled what the code does, not what its author intends
      (`corrected M099`).
- [x] T5: Add the NEWS entry. Then run `Rscript tools/dod-gate.R` and get exit 0.

## Work log

- 2026-10-18: created by /milestone-plan.
- 2026-10-18: tier check: milestone — superseding D-032 is a design decision.
- 2026-10-18: collision check: candidate row "D-032 premise contradicted by its source" absorbed. D-032 superseded (T3), D-017 kept. The rows "D-017 retention vs Figure 6B" and "Unbroken chase to level a" are distinct and stay. M78 (done) shipped the contiguity test that T1 extends.
- 2026-10-18: inbox sweep: 1 open issue (#111, upstream re-oracle watch), no overlap. No open PRs.
- 2026-10-18: lessons harvest: the M78 lesson (committed output settles semantics) bears on scope. T4 narrows it to code behavior, not intent.
- 2026-10-18: measured on the AMH fit: `m7f2` largest direct |r| .8981 to level 6, .9084 to level 5. It is in no chain, and `manual = "m7f2"` flags it with reason "manual".
- 2026-10-18: criteria audit (full mode, fresh Opus reader) returned 13 findings. Fixed: AC1 overclaimed `manual =` as serving her step, now "nearest step" with its limits. AC1 citations corrected to D-017/M53 and D-032/M78. AC2 gained the phi stop and the reference-code attribution. AC3 gained `pairs = "all"`, "largest", and defined flag and chain terms. AC4 "artefact stage" glossed. The AC5 gate run moved to T5. Answered: the M2 row stays `match`, see T4. The vignette needs no edit, see Out.
- 2026-10-18: question set: package behavior — keep the contiguous chase and add docs. User-doc attribution — package only, no citation of her email.
- 2026-10-18: plan gate chose docs plus `manual =` over a gap report in `prune("artifact")` and over an opt-in bridged chase, because her code has no automated bridge and the user chose the smallest change. Falsified by Forbes asking for the bridge in code, or by a user who needs a hand-bridged factor recorded as redundant with its chain values.
- 2026-10-18: the question set changed AC2 (no citation of her email), so the same reader re-audited the final criteria (full mode) and returned 4 findings, all fixed. T1's bare `next` plant would error, not fail, so it now adds the direct link. AC1's first reason narrowed to the AMH 54/54 output and the `g2` case, because row E6 departs. AC2 "when set" now reads non-NULL, automatic value included. AC4 now covers the LESSONS correction.
- 2026-10-18: plan gate chose an AMH test over a planted mock because AC3 needs a real gap, and `m7f2` is one in shipped data. Falsified by an AMH refit where `m7f2` no longer has the gap.
- 2026-10-18: implement started on branch `m099-chase-gap-supersede`, cut from master at `f136a00`.
- 2026-10-18: T1 done. The AMH `m7f2` test passes on current code. A planted gap skip (`next` plus a cross-gap link in `.strong_links_direct`) made 4 expectations fail as assertions with no error, and the plant was reverted (`git diff R/prune.R` empty). The test also checks the strong-link table.
- 2026-10-18: substantive amendment: AC3 said "absent from both label columns of `x$prune$chains`", but that table has one label column, `id`. AC3 now reads "has `pruned = FALSE` in `x$prune$nodes` and does not appear in the `id` column of `x$prune$chains`". The deliverable is unchanged, so no stop.
- 2026-10-18: re-audit: AC3 (full) — nothing. The fresh reader confirmed `id` is the only label column and that the clause still catches a gap-tolerant chase.
- 2026-10-18: T2 done. The `"direct"` help entry states the stop rule (phi included), the reference-code attribution, and the `manual =` route with reason "manual" and no chain row. The `rules` entry and the comments at the top of `R/prune.R` and above `.strong_links_direct` drop "faithful to Forbes", and the top comment cites D-040. The AMH test gained a check that the manual flag adds no chain row, so that help claim is test-backed. `document()` regenerated `man/prune.Rd`. Prune tests: 63 run, 0 failed.
- 2026-10-18: T3 done. D-040 supersedes D-032 to the AC1 wording. It cites the 54/54 figure that D-017 records and adds no new measurement. It names the two declined alternatives and their reopening evidence.
- 2026-10-18: T4 done. `forbes2023.md` gained a bullet under "Correspondence conventions" with her 2026-07-30 chase statement. The `source-departures.md` M2 row cites it and stays `match`. The M78 lesson is corrected in place (`corrected M099`). The simple-english hook flags older text in both reference pages. That text is outside this milestone and was left as it is.
- 2026-10-18: T5 in progress. The NEWS entry is in the development section. `tools/dod-gate.R` is running, and the claim audit reader was started over the added lines outside `cairn/`.
- 2026-10-18: first gate run passed (check 0/0/0, coverage 100%, style, lint, pkgdown). The claim audit then found the `manual` sentence overstated: a factor already flagged keeps "redundant", and a gap-stopped factor can have chain rows from the levels it passed. The help and NEWS now say that a factor no rule flagged gets "manual" and that `manual` adds no chain row. The test now checks that the chains table is unchanged by `manual`. The overlap case is held by the existing test "manual pruning combines with an auto rule". The gate is running again.
- 2026-10-18: claim audit: 22 claims read, 1 corrected — R/prune.R, man/prune.Rd, NEWS.md. The same reader re-read the corrected claim once and found it correct.
- 2026-10-18: T5 done. Second gate run on the corrected text passed: `tools/dod-gate.R` exit 0 (prose clean, check 0/0/0, coverage 100%, style and lint clean, pkgdown index complete). Status set to review.

## Decisions

## Review
