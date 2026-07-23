# M73: Draft the manuscript Introduction with the verified framing + application sources

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** `m73-manuscript-introduction` · https://github.com/jmgirard/ackwards/pull/77

## Goal

Replace the `[AUTHOR TO DRAFT]` Introduction stub in `manuscript/manuscript.qmd` with a complete, citation-backed first draft that frames the paper using the now-verified HiTOP and published-application source notes.

## Scope

**In:** Draft the `# Introduction` prose covering the stub's four suggested beats
(the misframed "how many factors?" problem; the hierarchical view of
personality *and* psychopathology structure; the gap — no general, tested,
documented tool despite decades of bespoke scripts; a one-paragraph roadmap).
Cite `kotov2017` for the HiTOP framing and a **curated, representative** set of
published bass-ackwards applications spanning the personality domain
(`markon2005`, `wright2014a`, `partsch2022`) and the psychopathology domain
(`kim2015`, `forbush2018`/`forbush2024`, `cowan2024`, `michelini2019`,
`carmichael2025`, `forbes2025`/`forbes2025a`) — exemplars per domain, not the
whole shelf — alongside the method lineage already cited (`goldberg2006`,
`waller2007`, `forbes2023`). Add each new key to `manuscript/references.bib`
with a Crossref-verified entry.

**Out:** Discussion drafting + Method/Package citation enrichment → M74 (depends
on M73). Any package R code change → not here (manuscript is `.Rbuildignore`d).
Authoring new `cairn/references/` notes → not needed; every cited source already
has a committed, extraction-verified note (M69–M72).

## Acceptance criteria

- [x] The `# Introduction` `[AUTHOR TO DRAFT]` blockquote stub is gone, replaced
      by drafted prose that addresses each of the four suggested beats (verified
      against a beat→paragraph checklist recorded in the Review section).
- [x] Every `@citekey` newly used in the Introduction traces to a committed,
      extraction-verified `cairn/references/` note (source note or a member of
      `applications.md`/`background.md`); the citekey→note map is recorded as
      evidence.
- [x] The Introduction cites `kotov2017` for HiTOP framing and ≥1 personality-domain
      and ≥1 psychopathology-domain application exemplar; no domain claim is made
      without a citation tracing to a verified note.
- [x] Every new `@citekey` resolves in `manuscript/references.bib`, and its
      bibliographic entry (author/year/title/journal/DOI) is Crossref-checked
      (per the M63 lesson: re-verify before propagating).
- [x] `manuscript/manuscript.qmd` renders to **both** PDF and docx with no
      unresolved-citation warnings and no new LaTeX errors (render log recorded).
- [x] `git diff` for the milestone touches only `manuscript/` and `cairn/` — no
      package `R/`, `tests/`, `NAMESPACE`, or `DESCRIPTION` change.

## Coverage

- AC1 → T3
- AC2 → T1, T3
- AC3 → T1, T3
- AC4 → T2
- AC5 → T4
- AC6 → T4

## Tasks

- [x] T1: Triage the shelf's application + framing notes into the curated Intro
      citation set — list each chosen citekey, the beat it supports, and confirm
      it has a committed, extraction-verified note (`INDEX.md` + the note's
      provenance block). Record the citekey→note map.
- [x] T2: Add the selected keys to `manuscript/references.bib` with
      Crossref-verified entries (author, year, title, journal, volume/issue,
      pages, DOI).
- [x] T3: Draft the Introduction prose replacing the stub — the four beats,
      integrating the citations; keep the scholarly argument the author's to
      refine (a complete draft, not the final word).
- [x] T4: Render PDF + docx (`quarto render`); confirm no unresolved citations /
      new LaTeX errors; confirm the diff scope is manuscript + tracking only.

## Work log

- 2026-07-23: created by /milestone-plan (split from the "Draft author-owned Intro + Discussion" candidate; M56 lineage). Intro half; Discussion + citation enrichment is M74.
- 2026-07-23: T1+T2 — curated 10-source Intro set (kotov2017; markon2005/wright2014a/partsch2022; kim2015/forbush2018/cowan2024/michelini2019/carmichael2025/forbes2025), each with a committed verified note; all 10 DOIs Crossref-verified (full author lists, version-of-record issue years) and added to references.bib. Dropped forbush2024/forbes2025a as redundant exemplars.
- 2026-07-23: T3+T4 — drafted 4-paragraph Introduction replacing the stub (all 4 beats, 10 new + 3 lineage citations). `quarto render` clean: PDF 16pp + docx, 0 unresolved-citation warnings, 0 LaTeX errors, all 10 new keys resolve in-text and in the reference list. Diff scope = manuscript/ + cairn/ only (no package files). Updated manuscript/README.md status (Intro now a complete draft; Discussion still stub → M74).

## Decisions

## Review

**PR:** https://github.com/jmgirard/ackwards/pull/77 (draft) · branch `m73-manuscript-introduction`, 3 commits ahead of `master` (unmoved since cut).

**Consistency gate.** `cairn_validate` exit 0 (82 pre-existing advisories, none introduced). `devtools::document()` no generated-file drift. No `IPn/GPn` changed → `cairn_impact` skipped. Remaining `consistency-gate` slot checks (README.Rmd knit, pkgdown, NEWS, `.Rbuildignore`, full `devtools::check()`) are package-artifact checks **not implicated**: the branch diff touches zero package files (`git diff --name-only master..HEAD` = `manuscript/{manuscript.qmd,references.bib,README.md}` + `cairn/`), and `manuscript/` is `.Rbuildignore`d — so the package is byte-identical to `master`, already green.

**AC evidence (fresh):**
- AC1 — `grep "AUTHOR TO DRAFT"` = 1 hit (line 352, the Discussion stub, untouched); Intro stub gone. Beat→paragraph map: ¶1 how-many-factors reframe (@goldberg2006); ¶2 hierarchical view across personality (@markon2005, @wright2014a, @partsch2022) + psychopathology (@kotov2017; @kim2015, @forbush2018, @cowan2024, @michelini2019, @carmichael2025, @forbes2025); ¶3 bass-ackwards as shared engine + the tooling gap (@goldberg2006, @waller2007, @forbes2023); ¶4 paper roadmap. All four beats present.
- AC2 — the 13 Intro `@keys` (10 new + 3 lineage) each carry an `INDEX.md` line and a committed extraction-verified note (kotov2017→background.md; markon2005/wright2014a/kim2015/forbush2018/cowan2024→applications.md; partsch2022/michelini2019/carmichael2025/forbes2025→own notes; all provenance-verified 2026-07-16/19/23).
- AC3 — kotov2017 (HiTOP) ✓; personality exemplars markon2005/wright2014a/partsch2022 ✓; psychopathology exemplars kim2015/forbush2018/cowan2024/michelini2019/carmichael2025/forbes2025 ✓.
- AC4 — all 10 new keys have a `references.bib` entry; all 10 DOIs Crossref-verified (title/journal/vol/issue/pages/full author list; version-of-record issue years for the 3 online-first-year cases: partsch2022, cowan2024, forbes2025).
- AC5 — fresh `quarto render`: exit 0, PDF (194 KB, 16pp) + docx (1.1 MB); 0 "not found"/LaTeX-error lines; `pdftotext` shows 0 `??` markers; all 10 new surnames appear both in-text and in the reference list.
- AC6 — `git diff --name-only master..HEAD` = manuscript/ + cairn/ only; no `R/`, `tests/`, `man/`, `NAMESPACE`, `DESCRIPTION`.
