# M74: Draft the manuscript Discussion + enrich the seeded sections with method-backer citations

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M73
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** m74-manuscript-discussion-citations

## Goal

Replace the `[AUTHOR TO DRAFT]` Discussion stub with a complete, citation-backed
first draft, and backfill the seeded Method/Package sections with the verified
method-backer citations their currently-uncited claims are owed.

## Scope

**In:** (1) Draft the `# Discussion` prose covering the stub's four beats (what
the package enables over bespoke scripts; interpretive cautions the package
enforces; scope/limitations — including the descriptive-not-confirmatory and
"sequential, not truly hierarchical" contrast against bottom-up
Schmid–Leiman/higher-order models; future directions + availability). Cite
`schmid1957` and `yung1999` at the scope-boundary/honest-caveat beat (the
citations `background.md` flags for exactly this section), plus application
context and `beauducel2024`/`williams2025` for the factor-score-validity caution
where apt. (2) Add method-backer citations to the seeded Method/Package sections
for their uncited claims — `kaiser1958` (varimax origin), `asparouhov2009`
(ESEM), `grice2001` (factor-score indeterminacy → tenBerge), and the
`suggest_k()` criteria (`revelle1979`, `ruscio2012a`) — **without altering any
technical claim** in those sections. Add each new key to `references.bib`
(Crossref-verified).

**Out:** Introduction drafting → M73 (this milestone depends on it). Rewriting
the seeded sections' technical substance → not here; this is citation backfill
only. Package R code change → not here.

## Acceptance criteria

- [ ] The `# Discussion` `[AUTHOR TO DRAFT]` blockquote stub is gone, replaced by
      drafted prose addressing each of the four suggested beats (verified against
      a beat→paragraph checklist recorded in the Review section).
- [ ] The Discussion cites `schmid1957` and `yung1999` at the scope-boundary beat;
      every `@citekey` newly used in the Discussion traces to a committed,
      extraction-verified `cairn/references/` note (citekey→note map recorded).
- [ ] The seeded Method and Package sections gain method-backer citations for
      their previously-uncited claims (≥ the varimax, ESEM, tenBerge/factor-score,
      and `suggest_k`-criteria claims named In-scope), each tracing to a verified
      note; a diff review confirms no technical claim in those sections was
      changed — only citations added.
- [ ] Every new `@citekey` resolves in `manuscript/references.bib` with a
      Crossref-checked entry.
- [ ] `manuscript/manuscript.qmd` renders to **both** PDF and docx with no
      unresolved-citation warnings and no new LaTeX errors (render log recorded).
- [ ] `git diff` for the milestone touches only `manuscript/` and `cairn/` — no
      package `R/`, `tests/`, `NAMESPACE`, or `DESCRIPTION` change.

## Coverage

- AC1 → T1
- AC2 → T1, T4
- AC3 → T2
- AC4 → T3, T4
- AC5 → T5
- AC6 → T5

## Tasks

- [x] T1: Draft the Discussion prose replacing the stub — the four beats,
      integrating `schmid1957`/`yung1999` at the scope-boundary caveat and
      application/factor-score context where apt (a complete draft to refine).
- [x] T2: Backfill method-backer citations into the seeded Method/Package
      sections at each named uncited claim (`kaiser1958`, `asparouhov2009`,
      `grice2001`, `revelle1979`, `ruscio2012a`), citations only — no technical
      wording change.
- [ ] T3: Add all new keys to `manuscript/references.bib` with Crossref-verified
      entries.
- [ ] T4: Record the citekey→note map for every newly cited source (Discussion +
      seeded-section backfill), confirming each has a committed, verified note.
- [ ] T5: Render PDF + docx; confirm no unresolved citations / new LaTeX errors;
      confirm the diff scope is manuscript + tracking only.

## Work log

- 2026-07-23: created by /milestone-plan (split from the "Draft author-owned Intro + Discussion" candidate; M56 lineage). Discussion + citation-enrichment half; depends on M73 (shared manuscript file, sequenced to avoid conflicts).
- 2026-07-23: question gate — author chose full Discussion citations (add beauducel2024/williams2025 factor-score caution + kotov2017/michelini2019 context) and concrete future-directions prose (ESEM/polychoric comparability + boot_edges).
- 2026-07-23: T1 done — replaced the `[AUTHOR TO DRAFT]` Discussion stub with a five-paragraph draft covering the four beats; cites hu1999, beauducel2024, williams2025, yung1999, schmid1957, kotov2017, michelini2019.
- 2026-07-23: T2 done — four citations-only backfill edits (kaiser1958 at varimax; asparouhov2009 at ESEM engine; grice2001 at factor-score materialization; revelle1979+ruscio2012a at the suggest_k retention-criteria group). `git diff` confirms no technical wording changed.

## Decisions

## Review
