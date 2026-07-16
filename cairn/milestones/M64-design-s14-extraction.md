<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M64: DESIGN §14 → DECISIONS.md extraction (hybrid + entomb)

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** m64-design-s14-extraction · https://github.com/jmgirard/ackwards/pull/65   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create -->

Complete the M20-deferred migration step: retire DESIGN.md's embedded §14 decision log by
entombing it verbatim in `cairn/legacy/`, appending D-entries for the still-governing decisions
D-001–D-015 don't cover, and repointing every live inline `§14.x` reference.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:**
- Verbatim entombment of DESIGN.md §14 (lines 461–876 at plan time) as
  `cairn/legacy/DESIGN-s14-decision-log.md` with a provenance header (extracted M64, source commit).
- New D-entries (D-016+) for still-governing, cross-cutting §14 decisions with no existing D-entry
  (≈12–16: M53 redundancy criterion; M32 tidy/variance/meta decisions; M34 manual-prune/spelling/
  removed-args details; M38 FIML route + `n_obs`; M45 `scaling = "fit"` + `predict()`; M46
  `comparability()`; M47 `boot_edges()`; M49 `correct`/`check_items()`/trust-tiering; M50 bfi25
  labels + `suggest_k()` ordinal warning; M51 factor labels; M52 `factorability()`; the genuine
  declines — `categorical` flag M40, dual chi-square M49A, `label_items()`/third-dataset M49A).
- DECISIONS.md header gains an anchor-resolution note: `§14.x` citations resolve against the
  legacy copy. D-001–D-015 bodies untouched (append-only).
- DESIGN.md §14 reduced to a short pointer stub (section number kept — §9/§12/§15 numbering must
  not shift); the "Known limitations / deferred" block becomes its own live DESIGN section,
  trimmed to live items only (struck-through done items stay findable in the legacy copy).
- Repoint live inline `§14`/`§14.x` references: DESIGN.md internal (13 lines), CLAUDE.md (2),
  ROADMAP.md ESEM-extensions candidate row (§14.35/§14.36), `cairn/references/`
  (everett1983 ×2, asparouhov2009 ×1, lorenzoseva2006 ×3, tenberge1999 ×3),
  `tests/testthat/test-vignette-m24.R:200` comment.

**Out:**
- `cairn/legacy/` and `cairn/milestones/archive/` citations — entombed history, stay verbatim.
- Rewriting D-001–D-015 source lines — append-only; the header note covers resolution.
- IP/GP principle formalization → stays a ROADMAP candidate (`/design-interview`, M20 migration
  gate sibling); this milestone only lightens DESIGN.md for it.
- Full §14 → 46 individual D-entries (literal conversion) — declined at the plan gate 2026-07-16
  in favor of hybrid + entomb; the verbatim text is conserved in legacy, not duplicated.
- Any R behavior change — the sole code-file edit is one test comment.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [x] AC1: `cairn/legacy/DESIGN-s14-decision-log.md` exists, carries a provenance header, and its
      body is byte-identical to the §14 block of the pre-M64 DESIGN.md
      (`git show <base>:cairn/DESIGN.md` diff over the extracted range).
- [x] AC2: A triage ledger (in the PR description) dispositions every §14 numbered item (1–46),
      the build-time items, and each limitations bullet as: covered by D-00x / new D-entry
      D-0xx / stays as live limitation / history-only (entombed). No item is silently absent;
      no new D-entry duplicates D-001–D-015.
- [x] AC3: Every new D-entry has Context/Decision/Consequences plus a source citation (legacy
      anchor + originating milestone); DECISIONS.md header carries the anchor-resolution note;
      D-001–D-015 bodies show zero diff.
- [x] AC4: DESIGN.md §14 is a ≤10-line pointer stub; a live "Known limitations" section exists
      listing only unresolved limitations; no other DESIGN section is renumbered
      (grep: `## 15. Milestones` heading unchanged).
- [x] AC5: `grep -rn "§14" DESIGN.md CLAUDE.md ROADMAP.md cairn/references/ R/ tests/ vignettes/`
      (live surface) returns only the §14 stub's own pointer line(s) — zero stale `§14.x` targets;
      `cairn/legacy/`, `cairn/milestones/archive/`, and DECISIONS.md are excluded as
      history-bearing by design.
- [x] AC6: `cairn_validate.py` all checks pass, and the test suite is green
      (`TESTTHAT_CPUS=8 devtools::test()` — one test file is touched, comment-only).

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T2
- AC2 → T1
- AC3 → T3
- AC4 → T4
- AC5 → T5, T6
- AC6 → T6

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1: Triage pass — enumerate every §14 item (numbered 1–46, build-time 6–10 list, limitations
      bullets) against D-001–D-015; write the disposition ledger (destined for the PR description;
      not this file — line cap).
- [x] T2: Entomb — create `cairn/legacy/DESIGN-s14-decision-log.md`: provenance header + verbatim
      §14 body copied from the pre-branch DESIGN.md.
- [x] T3: Append D-016+ entries per the T1 ledger (Context/Decision/Consequences, legacy-anchor +
      milestone citations, real dates where known else "date: see legacy" per existing style);
      add the header resolution note.
- [x] T4: Rewrite DESIGN.md — §14 → pointer stub (to DECISIONS.md + the legacy file); new live
      "Known limitations" section (trimmed per T1); keep all section numbers stable.
- [x] T5: Repoint the live reference surface enumerated in Scope-In (DESIGN internal, CLAUDE.md,
      ROADMAP ESEM row, four references pages, the test comment) to D-numbers / the new section /
      legacy anchors as appropriate.
- [x] T6: Verify — scoped `§14` grep (AC5), `cairn_validate.py`, `TESTTHAT_CPUS=8` suite.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-16: created by /milestone-plan; promotes the 2026-07-11 ROADMAP candidate (M20 migration
  gate); gate resolved: hybrid + entomb, history left verbatim + resolution note, limitations stay
  in DESIGN trimmed.
- 2026-07-16: T1 done — ledger drafted (scratchpad → PR body): 15 new D-entries D-016–D-030
  (one more than the plan's upper estimate: item 12, ESEM self-computed tenBerge weights, is cited
  by two references pages and covered by no D-entry); 4 live limitations; rest covered/arch/history.
- 2026-07-16: T2 done — §14 body (416 lines) entombed verbatim at
  cairn/legacy/DESIGN-s14-decision-log.md; byte-identity vs base 0ce1095 diff-verified.
- 2026-07-16: T3 done — D-016–D-030 appended (15 entries, Context/Decision/Consequences + legacy
  anchors); DECISIONS.md preamble now carries the anchor-resolution note; D-001–D-015 bodies
  untouched.
- 2026-07-16: T4+T5 done — DESIGN.md 897→511 lines: §14 = 9-line stub, unnumbered live "Known
  limitations" section (4 entries), §15 heading unchanged; all live §14.x refs repointed
  (DESIGN ×7, CLAUDE.md ×2, ROADMAP ESEM row ×2, 4 references pages ×9, test comment ×1).
- 2026-07-16: T6 done — scoped grep clean (only stub pointer lines + the M64 row title);
  cairn_validate all-pass (dangling-token advisories 87→83); suite FAIL 0/PASS 2303 (one
  parallel-run-only WARN in test-pca.R ordinal advisory, absent in isolation, unrelated to the
  comment-only diff). Status → review.

## Decisions
<!-- owner: implement / review · append-only; milestone-local; promote
     cross-cutting ones to cairn/DECISIONS.md -->

## Review
<!-- owner: review · exclusive; evidence per criterion, consistency-gate
     results, review findings + triage. EXEMPT from the 150-line cap (M55). -->

Evidence gathered fresh 2026-07-16 on branch m64-design-s14-extraction (PR #65, base 0ce1095):

- **AC1** ✓ — `diff <(tail -n +8 cairn/legacy/DESIGN-s14-decision-log.md) <(git show 0ce1095:cairn/DESIGN.md | sed -n '461,876p')` → zero diff ("BYTE-IDENTICAL"); provenance header present (line 1 names M64, date, source commit).
- **AC2** ✓ — PR #65 body carries the full triage ledger (all 46 numbered items + build-time list + limitations bullets dispositioned; `gh pr view 65` confirmed); dedup constraint stated in ledger footer.
- **AC3** ✓ — D-016–D-030 all present, each with Context/Decision/Consequences + `_Source:` (counts 30/30/30 across the 30 total entries); `git diff master..HEAD -- cairn/DECISIONS.md` removes only the 5 preamble lines (resolution note replaces them); D-001–D-015 bodies zero diff.
- **AC4** ✓ — §14 stub = 9 lines; live `## Known limitations` section with 4 entries; `## 15. Milestones` at DESIGN.md:491 unchanged.
- **AC5** ✓ — scoped grep returns only DESIGN.md:5 + :468 (the stub's own pointers) and the M64 ROADMAP row title (the milestone's own name, not a stale anchor); zero stale `§14.x` targets.
- **AC6** ✓ — `cairn_validate.py` exit 0 all-pass (advisories 87→83); suite `TESTTHAT_CPUS=8` FAIL 0 / WARN 1 / PASS 2303 — the WARN is a parallel-run-only session-state artifact in test-pca.R (ordinal `.frequency="once"` advisory), absent when the file runs in isolation, unrelated to the comment-only test diff.

Consistency gate (r-package profile): `devtools::document()` → no diff (NAMESPACE/man clean); `pkgdown::check_pkgdown()` → no problems; full `devtools::check()` → 0 err / 0 warn / 0 note; NEWS.md justified-skip (no user-visible change — docs/tracking only); no new top-level files needing .Rbuildignore (legacy file is under cairn/, already ignored).

Independent review (3 fresh-context lenses + scorer, 2026-07-16):
- [O] diff-bug: 1 finding — D-016 Context's added parenthetical "(only EBM/regression-style scores)" misenumerated `lavPredict()`'s methods (omits Bartlett, misapplies EBM; invented during extraction — legacy §14.12 says only "it lacks tenBerge"). Scored 85 → **fixed**: parenthetical dropped, entry now matches the source claim exactly. All other axes clean (fidelity of all 15 entries, all 14 repoints, limitations triage, internal consistency).
- [S] blame-history: no findings — M61 citation fix preserved; entombment captured the freshest DESIGN.md state; repoints verified against D-entry content.
- [S] prior-PR-comments: no prior-PR evidence (solo repo, zero GitHub review comments) — clean no-op.
- Sub-80 findings: none (the single finding scored 85 and was actioned).
