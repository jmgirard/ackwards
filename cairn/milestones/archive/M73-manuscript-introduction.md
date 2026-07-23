# M73: Draft the manuscript Introduction with the verified framing + application sources

**Status:** done (2026-07-23, PR #77 https://github.com/jmgirard/ackwards/pull/77)

**Goal:** Replace the `[AUTHOR TO DRAFT]` Introduction stub in `manuscript/manuscript.qmd` with a complete, citation-backed first draft that frames the paper using the now-verified HiTOP and published-application source notes.

**Outcome:** The `# Introduction` stub became four drafted paragraphs — the how-many-factors reframe, the hierarchical view across personality (markon2005, wright2014a, partsch2022) and psychopathology (kotov2017/HiTOP; kim2015, forbush2018, cowan2024, michelini2019, carmichael2025, forbes2025), bass-ackwards as the shared engine + the tooling gap, and a paper roadmap. Added 10 Crossref-verified entries to `manuscript/references.bib` (full author lists; version-of-record issue years for partsch2022/cowan2024/forbes2025), each tracing to a committed extraction-verified `cairn/references/` note. Updated `manuscript/README.md` status. `quarto render` clean (PDF 16pp + docx, 0 unresolved citations / LaTeX errors). No package code changed (manuscript is `.Rbuildignore`d). Discussion + method-backer citation enrichment remain M74 (depends on this).

**Decisions:** none (curated 10-source exemplar set, dropping forbush2024/forbes2025a as redundant, is a plan-scoped judgment, not a cross-cutting decision).

**Review:** 3 fresh-context lenses + scorer. Blame-history + prior-review: no findings (M56's Waller-title/bfi25-n/54-vs-55 findings none regressed). Diff-bug [O] verified all citations + 10 bib entries clean; one actioned finding F1 (score 85): the Intro gap paragraph near-verbatim duplicated the Package Statement of Need — fixed by rewriting Intro ¶3 to motivate (not re-specify) the gap. F2 (score 64, "three decades" overstatement) logged, not actioned; the F1 rewrite incidentally removed the introduced instance.
