# M71: Author + verify the 5 application source notes as citation precedents + Forbes drift-watch

- **Status:** review
- **Priority:** normal
- **Depends on:** M70
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** m71-reference-notes-applications

## Goal

Author a verified `cairn/references/` source note for each of the 5 application PDFs (carmichael2025, michelini2019, partsch2022, forbes2025 [DSM-5 reconstruction], forbes2025a [youth]) as manuscript/vignette citation precedents, with the two Forbes notes doubling as extended-bass-ackward drift-watch against the forbes2023 fidelity contract.

## Scope

**In:**
- 5 standalone `references/<citekey>.md` notes authored from `templates/source-note.md`, each verified against its rendered shelf PDF (standing facts anchored, provenance block, INDEX line).
- Each note's Role states the domain and whether the source uses bass-ackwards specifically or an adjacent hierarchical method.
- Confirm partsch2022 actually uses bass-ackwards against the now-on-shelf PDF, with a page anchor — the sweep's claim rested on abstract/secondary sources because the PDF was access-blocked then.
- The two Forbes notes (forbes2025, forbes2025a) each record, as a dated observation, whether the paper's extended-BA / hierarchical-PCA procedure diverges from the forbes2023 fidelity contract (stopping criterion, redundancy pruning), with anchors.

**Out:**
- Any DESIGN/roxygen wiring — applications back manuscript/vignette citations, not code defaults; docs-only, `cairn/`-only diff.
- If a Forbes note surfaces a genuine fidelity-contract divergence, it becomes a new `candidate` ROADMAP row, not silent scope creep here.
- The 5 backer notes → M70.
- BRM manuscript prose → the existing "Intro + Discussion" candidate row.

## Acceptance criteria

- [ ] AC1: 5 notes exist at `cairn/references/{carmichael2025,michelini2019,partsch2022,forbes2025,forbes2025a}.md`; each provenance `Extraction:` line begins its own line (M60 lesson) and reads `verified <YYYY-MM-DD> … — observed <YYYY-MM-DD>`; every standing fact carries a page/table anchor; each verbatim quote transcribed from the rendered page, not flattened `pdftotext` (M67/M69 lesson).
- [ ] AC2: `INDEX.md` gains one filename-first line (`- <name>.md — …`) per note under the Applications section, each naming its domain.
- [ ] AC3: partsch2022's bass-ackwards use is confirmed against the shelf PDF with a page anchor; each note states bass-ackwards vs adjacent method.
- [ ] AC4: forbes2025 and forbes2025a each carry a dated observation on whether their extended-BA/hPCA procedure diverges from the forbes2023 fidelity contract, with anchors.
- [ ] AC5: docs-only `cairn/`-only diff (no `R/`, no `.Rd`); `cairn_validate` exits 0 (references check: every page has its INDEX line).

## Coverage

- AC1 → T1, T2
- AC2 → T3
- AC3 → T1
- AC4 → T2
- AC5 → T4

## Tasks

- [x] T1: Author + verify the 3 new-domain notes (carmichael2025 — TBI; michelini2019 — ABCD child/adolescent; partsch2022 — VIA character strengths) against rendered shelf pages, each stating domain + method; for partsch2022 confirm bass-ackwards use with a page anchor.
- [x] T2: Author + verify the 2 Forbes continuation notes (forbes2025 — DSM-5 reconstruction; forbes2025a — youth), each recording the extended-BA/hPCA procedure and a dated observation on forbes2023 fidelity-contract divergence.
- [x] T3: Add the 5 filename-first `INDEX.md` Applications lines.
- [x] T4: Run `cairn_validate` (exit 0) and confirm the diff is `cairn/`-only; fix any fallout.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: /milestone-implement — set in-progress; branch `m71-reference-notes-applications` cut from master; cairn_validate baseline exits 0.
- 2026-07-23: T1 — authored carmichael2025.md (TBI; Forbes' extended BA, redundancy r≥.90 & φ>.95), michelini2019.md (ABCD child/adult; Goldberg BA on EFA/geomin, edge ≥.65), partsch2022.md (VIA strengths; Goldberg BA via Waller-2007 code, oblique Promax). partsch bass-ackwards confirmed p.832. All verbatim quotes read from rendered pages + cross-checked against pdftotext.
- 2026-07-23: minor amendment (INDEX placement) — the 5 standalone application notes go under a new `## Applications (one note per source)` INDEX subsection paralleling the existing "Method sources"/"Default-rationale backers" standalone subsections, not the "Applications (collapsed in applications.md)" block (those are collapsed, these are standalone). Executed in T3.
- 2026-07-23: T2 — authored forbes2025.md (DSM-5; hPCA/extended-BA, oblique oblimin, redundancy r>.9 & φ>.95 = matches contract) and forbes2025a.md (youth; extended-BA, deliberate divergence — dropped φ conjunction, r≥.9 only, p.287). AC4 drift-watch section in each. Divergence disposition recorded as MD-1 (no candidate row: already parameterized by `redundancy_phi`). All quotes read from rendered pages + cross-checked against pdftotext.
- 2026-07-23: T3 — added the 5 filename-first INDEX.md lines under a new `## Applications (one note per source)` subsection, each naming its domain (TBI, child/adolescent, personality, DSM-5, youth).
- 2026-07-23: T4 — cairn_validate exits 0 (references index<->disk PASS; the 83 advisories are pre-existing legacy-M-number + kaiser1958 M70 staleness, none from M71 pages); branch diff confirmed cairn/-only (no R/, .Rd, man/, NAMESPACE, data/). Status → review.

## Decisions

- MD-1 (2026-07-23): The forbes2025a redundancy divergence gets no ROADMAP candidate row. **Context:** the plan's Out scope said a genuine forbes2023 fidelity-contract divergence surfaced by a Forbes note "becomes a new `candidate` ROADMAP row." forbes2025a (p. 287) states a deliberate divergence — it dropped forbes2023's Tucker-φ conjunction, using `r ≥ .9` only. **Decision:** no candidate row is created — the looser cutoff is already expressible as `prune("redundant")` with the φ filter off (`redundancy_phi = NA` / the PCA auto-default), and the stricter forbes2023 conjunction is `redundancy_phi = 0.95` (DESIGN §9; D-017); there is no capability gap, so a row would be a no-op. The divergence is recorded as the AC4 drift-watch observation in `forbes2025a.md`. **Consequences:** AC4 satisfied by the observation, not a row; surfaced to the owner at the completion gate for override.

## Review
