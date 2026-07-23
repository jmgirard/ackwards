# M70: Author + verify the 5 default-rationale backer notes and wire them into DESIGN §9 + roxygen

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** —

## Goal

Author a verified `cairn/references/` source note for each of the 5 default-rationale backer PDFs (grice2001, beauducel2024, tong2025, williams2025, kaiser1958) and wire each citation into the DESIGN §9 default it supports and the relevant roxygen `@references`.

## Scope

**In:**
- Rename shelf `williams2025a.pdf → williams2025.pdf`; use the `williams2025` citekey throughout.
- 5 new `references/<citekey>.md` notes authored from `templates/source-note.md`, each verified against its rendered shelf PDF (standing facts anchored, provenance block, INDEX line).
- Wire each citation into DESIGN §9's rationale for the default it backs — kaiser1958 → the `rotation` row; grice2001/beauducel2024/williams2025 → the `scores (method)` + `redundancy_phi` rows (tenBerge / factor-score indeterminacy rationale); tong2025 → the `k_max` row — evidence citations only, no default value changes.
- Add each as a roxygen `@references` entry on the relevant function (`R/ackwards.R` for the rotation/scoring backers; `R/suggest_k.R` for tong2025); regenerate docs.
- Full DoD gate (this milestone touches exported `.Rd` docs, so it is **not** documentary-only).

**Out:**
- The 5 application notes (carmichael2025, michelini2019, partsch2022, forbes2025, forbes2025a) → M71.
- Any change to a default *value* or to D-002/D-007 — this adds supporting evidence, it does not reopen the decisions.
- BRM manuscript prose → the existing "Intro + Discussion" candidate row.

## Acceptance criteria

- [ ] AC1: `williams2025a.pdf` renamed to `williams2025.pdf` on the shelf; no INDEX line, note heading/body, or roxygen entry uses the `williams2025a` citekey (grep hits only this milestone file's own rename description).
- [ ] AC2: 5 notes exist at `cairn/references/{grice2001,beauducel2024,tong2025,williams2025,kaiser1958}.md`; each provenance `Extraction:` line begins its own line (M60 lesson) and reads `verified <YYYY-MM-DD> … — observed <YYYY-MM-DD>`; every extracted standing fact carries a page/table anchor; each verbatim quote transcribed from the rendered page, not flattened `pdftotext` (M67/M69 lesson).
- [ ] AC3: `INDEX.md` gains one filename-first line (`- <name>.md — …`) per new note under the appropriate section.
- [ ] AC4: each backer citation appears in the DESIGN §9 rationale for the default it supports (kaiser1958→rotation; grice2001+beauducel2024+williams2025→scores/tenBerge & redundancy_phi; tong2025→k_max).
- [ ] AC5: each backer appears as a roxygen `@references` entry on the relevant function (`R/ackwards.R`: kaiser1958, grice2001, beauducel2024, williams2025; `R/suggest_k.R`: tong2025); `devtools::document()` re-run, `.Rd` updated, NAMESPACE unchanged.
- [ ] AC6: `Rscript tools/dod-gate.R` passes (check 0 err/0 warn/0 note, coverage maintained, style/lint/pkgdown clean) and `cairn_validate` exits 0.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T5
- AC6 → T6

## Tasks

- [x] T1: Rename `cairn/references/sources/williams2025a.pdf` → `williams2025.pdf`; grep the repo for `williams2025a` and fix any stray reference.
- [x] T2: Author + verify the 5 backer notes against rendered shelf pages — grice2001 (factor-score indeterminacy; regression/Bartlett vs correlation-preserving), beauducel2024 (determinacy vs inter-factor-correlation trade-off), tong2025 (bass-ackward stopping-criterion simulation; the two BA criteria), williams2025 (latent-variable vs factor-score criterion validity; note the provided R function), kaiser1958 (the varimax criterion). Each note names in its Role which DESIGN §9 default it backs.
- [x] T3: Add the 5 filename-first `INDEX.md` lines under the right section.
- [x] T4: Wire each citation into the matching DESIGN §9 rationale row (rotation / scores / redundancy_phi / k_max) — citations only, no value change.
- [x] T5: Add the `@references` entries (`R/ackwards.R`, `R/suggest_k.R`) in the existing citation style; run `devtools::document()`.
- [ ] T6: Run `Rscript tools/dod-gate.R` + `cairn_validate`; fix any fallout.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: set in-progress; branch m70-reference-notes-backers cut from master.
- 2026-07-23: T1 — renamed shelf williams2025a.pdf → williams2025.pdf (gitignored). Minor AC1 wording refinement: scoped the grep to citekey usages (INDEX/note/roxygen), since the milestone file itself names the rename. (Reverted an over-eager AC1 tick — AC boxes are review-owned.)
- 2026-07-23: T5 — added @references (ackwards.R: kaiser1958, grice2001, beauducel2024, williams2025; suggest_k.R: tong2025); devtools::document() regenerated ackwards.Rd + suggest_k.Rd; NAMESPACE unchanged.
- 2026-07-23: T4 — wired the 5 citations into DESIGN §9: rotation row (kaiser1958), scores row (grice2001, beauducel2024, williams2025), redundancy_phi row (grice2001), k_max row (tong2025). Citations only; no default value changed.
- 2026-07-23: T3 — added a new INDEX.md section "Default-rationale backers" with the 5 filename-first lines.
- 2026-07-23: T2 — authored 5 backer notes (kaiser1958, grice2001, beauducel2024, tong2025, williams2025), each verified against rendered shelf pages. DOIs Crossref-confirmed for williams2025 (21677026231225414) + tong2025 (s40647-024-00423-2); beauducel third author confirmed "Kuhl" (not "Kühl"). tong2025 shelf copy is the preprint — pagination basis recorded as preprint pages. kaiser1958 prints no DOI (registered BF02289233).

## Decisions

## Review
