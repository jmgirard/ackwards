# M64: DESIGN §14 → DECISIONS.md extraction (hybrid + entomb) — done 2026-07-16

**Goal:** complete the M20-deferred migration step — retire DESIGN.md's embedded §14 decision
log (416 lines, 46 numbered items).

**Outcome:** DESIGN.md 897→511 lines, now pure architecture. §14 body entombed byte-identically
(diff-verified vs base 0ce1095) as `cairn/legacy/DESIGN-s14-decision-log.md`, against which all
historical `§14.x` anchors (D-001–D-015 source lines, milestone archives) resolve — per the new
DECISIONS.md preamble note. Fifteen new D-entries (D-016–D-030) record the still-governing
decisions D-001–D-015 didn't cover (ESEM tenBerge mechanics; M53 redundancy criterion; M32
naming/scale; M34 prune mechanics; M38 FIML; M45 scoring; M46/M47 verbs; M50 labels; M49
robustness; M51 factor labels; M52 factorability; declines: dual chi-square, label ecosystem,
categorical flag). §14 is a 9-line pointer stub; a live unnumbered "Known limitations" section
carries the 4 unresolved limitations; all live `§14.x` refs repointed (DESIGN ×7, CLAUDE.md ×2,
ROADMAP ×2, four references pages ×9, one test comment). Triage ledger in the PR body.

**Key decisions (plan gate 2026-07-16):** hybrid + entomb over literal 46-entry conversion;
history files left verbatim + resolution note over repoint-everything; limitations kept in
DESIGN, trimmed to live items.

**Review:** all 6 ACs fresh-verified; check() 0/0/0, suite FAIL 0/PASS 2303, document() no-diff,
pkgdown + cairn_validate clean. One finding (scored 85, fixed): D-016 had invented an inaccurate
`lavPredict()` method enumeration absent from the source — dropped. Other two lenses clean.

**PR:** https://github.com/jmgirard/ackwards/pull/65 (squash 708e916).
