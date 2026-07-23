# M73: Draft the manuscript Introduction with the verified framing + application sources

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** —

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

- [ ] The `# Introduction` `[AUTHOR TO DRAFT]` blockquote stub is gone, replaced
      by drafted prose that addresses each of the four suggested beats (verified
      against a beat→paragraph checklist recorded in the Review section).
- [ ] Every `@citekey` newly used in the Introduction traces to a committed,
      extraction-verified `cairn/references/` note (source note or a member of
      `applications.md`/`background.md`); the citekey→note map is recorded as
      evidence.
- [ ] The Introduction cites `kotov2017` for HiTOP framing and ≥1 personality-domain
      and ≥1 psychopathology-domain application exemplar; no domain claim is made
      without a citation tracing to a verified note.
- [ ] Every new `@citekey` resolves in `manuscript/references.bib`, and its
      bibliographic entry (author/year/title/journal/DOI) is Crossref-checked
      (per the M63 lesson: re-verify before propagating).
- [ ] `manuscript/manuscript.qmd` renders to **both** PDF and docx with no
      unresolved-citation warnings and no new LaTeX errors (render log recorded).
- [ ] `git diff` for the milestone touches only `manuscript/` and `cairn/` — no
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
- [ ] T3: Draft the Introduction prose replacing the stub — the four beats,
      integrating the citations; keep the scholarly argument the author's to
      refine (a complete draft, not the final word).
- [ ] T4: Render PDF + docx (`quarto render`); confirm no unresolved citations /
      new LaTeX errors; confirm the diff scope is manuscript + tracking only.

## Work log

- 2026-07-23: created by /milestone-plan (split from the "Draft author-owned Intro + Discussion" candidate; M56 lineage). Intro half; Discussion + citation enrichment is M74.
- 2026-07-23: T1+T2 — curated 10-source Intro set (kotov2017; markon2005/wright2014a/partsch2022; kim2015/forbush2018/cowan2024/michelini2019/carmichael2025/forbes2025), each with a committed verified note; all 10 DOIs Crossref-verified (full author lists, version-of-record issue years) and added to references.bib. Dropped forbush2024/forbes2025a as redundant exemplars.

## Decisions

## Review
