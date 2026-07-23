# M74: Draft the manuscript Discussion + enrich the seeded sections with method-backer citations

- **Status:** review
- **Priority:** normal
- **Depends on:** M73
- **Driving RR:** —
- **Principles touched:** —
- **Branch/PR:** m74-manuscript-discussion-citations · PR #78 (https://github.com/jmgirard/ackwards/pull/78)

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

- [x] The `# Discussion` `[AUTHOR TO DRAFT]` blockquote stub is gone, replaced by
      drafted prose addressing each of the four suggested beats (verified against
      a beat→paragraph checklist recorded in the Review section).
- [x] The Discussion cites `schmid1957` and `yung1999` at the scope-boundary beat;
      every `@citekey` newly used in the Discussion traces to a committed,
      extraction-verified `cairn/references/` note (citekey→note map recorded).
- [x] The seeded Method and Package sections gain method-backer citations for
      their previously-uncited claims (≥ the varimax, ESEM, tenBerge/factor-score,
      and `suggest_k`-criteria claims named In-scope), each tracing to a verified
      note; a diff review confirms no technical claim in those sections was
      changed — only citations added.
- [x] Every new `@citekey` resolves in `manuscript/references.bib` with a
      Crossref-checked entry.
- [x] `manuscript/manuscript.qmd` renders to **both** PDF and docx with no
      unresolved-citation warnings and no new LaTeX errors (render log recorded).
- [x] `git diff` for the milestone touches only `manuscript/` and `cairn/` — no
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
- [x] T3: Add all new keys to `manuscript/references.bib` with Crossref-verified
      entries.
- [x] T4: Record the citekey→note map for every newly cited source (Discussion +
      seeded-section backfill), confirming each has a committed, verified note.
- [x] T5: Render PDF + docx; confirm no unresolved citations / new LaTeX errors;
      confirm the diff scope is manuscript + tracking only.

## Work log

- 2026-07-23: created by /milestone-plan (split from the "Draft author-owned Intro + Discussion" candidate; M56 lineage). Discussion + citation-enrichment half; depends on M73 (shared manuscript file, sequenced to avoid conflicts).
- 2026-07-23: question gate — author chose full Discussion citations (add beauducel2024/williams2025 factor-score caution + kotov2017/michelini2019 context) and concrete future-directions prose (ESEM/polychoric comparability + boot_edges).
- 2026-07-23: T1 done — replaced the `[AUTHOR TO DRAFT]` Discussion stub with a five-paragraph draft covering the four beats; cites hu1999, beauducel2024, williams2025, yung1999, schmid1957, kotov2017, michelini2019.
- 2026-07-23: T2 done — four citations-only backfill edits (kaiser1958 at varimax; asparouhov2009 at ESEM engine; grice2001 at factor-score materialization; revelle1979+ruscio2012a at the suggest_k retention-criteria group). `git diff` confirms no technical wording changed.
- 2026-07-23: T3 done — added 7 new bib keys (kaiser1958, asparouhov2009, grice2001, revelle1979, ruscio2012a, beauducel2024, williams2025). Each DOI verified live against the Crossref API (title/authors/journal/volume/issue/pages/year); revelle1979's DOI 10.1207/s15327906mbr1404_2 recovered from Crossref (source prints none); full given names pulled from Crossref for beauducel2024 + williams2025.
- 2026-07-23: T4 done — recorded the citekey→note map (12 keys) in the Decisions section; every added citekey traces to an extraction-verified reference note.
- 2026-07-23: T5 done — `quarto render manuscript.qmd` produced both manuscript.pdf and manuscript.docx (Quarto 1.9.38, ackwards 0.1.1); render log clean (no citeproc "not found", no LaTeX errors); docx text confirms all 7 new references resolved and 0 literal `@citekey` leaks; `git diff master...HEAD` touches only `manuscript/` + `cairn/` (no R/tests/NAMESPACE/DESCRIPTION); render outputs are gitignored. All tasks complete → status review.

## Decisions

**Citekey→note map (T4 citation-provenance evidence).** Every `@citekey` this
milestone adds — Discussion (D) + seeded-section backfill (B) — traces to a
committed, extraction-verified `cairn/references/` note, and every new-to-bib key
resolves in `manuscript/references.bib`. `bib?` = new bib entry added by M74 (T3);
others were already in the bib.

| citekey | where | note (verified) | bib? |
|---|---|---|---|
| schmid1957 | D (scope beat) | background.md — verified 2026-07-19 (M68) | pre-existing |
| yung1999 | D (scope beat) | background.md — verified 2026-07-19 (M68) | pre-existing |
| beauducel2024 | D (factor-score caution) | beauducel2024.md — verified 2026-07-23 (M70) | new (M74) |
| williams2025 | D (factor-score caution) | williams2025.md — verified 2026-07-23 (M70) | new (M74) |
| hu1999 | D (fit-values caution) | hu1999.md — verified 2026-07-19 (M69) | pre-existing |
| kotov2017 | D (taxonomy context) | background.md — verified 2026-07-19 (M68) | pre-existing |
| michelini2019 | D (taxonomy context) | michelini2019.md — verified 2026-07-23 (M71) | pre-existing |
| kaiser1958 | B (varimax) | kaiser1958.md — verified 2026-07-23 (M70) | new (M74) |
| asparouhov2009 | B (ESEM engine) | asparouhov2009.md — verified 2026-07-19 (M67) | new (M74) |
| grice2001 | B (factor scores) | grice2001.md — verified 2026-07-23 (M70) | new (M74) |
| revelle1979 | B (suggest_k VSS) | revelle1979.md — verified 2026-07-19 (M69) | new (M74) |
| ruscio2012a | B (suggest_k CD) | ruscio2012a.md — verified 2026-07-19 (M69) | new (M74) |

## Review

**Reviewed 2026-07-23** · PR #78 · branch `m74-manuscript-discussion-citations` (base `master` @ `1629d26`, no divergence).

### Consistency gate
- `cairn_validate`: **PASS** (exit 0; "all checks passed"). 82 advisory warnings, all pre-existing historical milestone-ID cross-references (M46/M53/…) in DESIGN/references/archive — not introduced by M74, advisory not gate.
- Coverage completeness: **PASS** (bundled in `cairn_validate`; AC1–AC6 all map to existing tasks per the Coverage section).
- `cairn_impact`: **N/A** — no IP/GP principle changed (`Principles touched: —`).
- R-toolchain gate (`r-package` `consistency-gate` slot): **vacuously satisfied — built package untouched.** `git diff master...HEAD` changes only `manuscript/` + `cairn/`, both `.Rbuildignore`d (`^manuscript$` L23, `^cairn$` L2); no `R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `README`, `_pkgdown`, or `NEWS` change. Hence `document()` is no-diff (no roxygen touched) and `devtools::check()` is byte-identical to M73's green state — re-running it on an unchanged package is the needless-suite-run CLAUDE.md warns against, so justified-skipped on provable non-impact. `pkgdown::check_pkgdown()` ran anyway: "✔ No problems found." NEWS entry not required (the manuscript is `.Rbuildignore`d, not a shipped package surface).

### Acceptance criteria (fresh evidence)
- **AC1 — stub replaced, four beats covered.** `grep -c "AUTHOR TO DRAFT"` = 0; the Discussion is five paragraphs. Beat→paragraph checklist: ¶1 = what the package enables over bespoke scripts (engine breadth, single verified edge path, machine-precision reproduction, tested/documented); ¶2 = interpretive cautions (descriptive-not-confirmatory, fit-values-not-verdicts, `prune()` flags-not-deletions, precision beyond the single strongest edge); ¶3 = factor-score-validity caution; ¶4 = scope/limitations (orthogonal-only, sequential-not-hierarchical vs. Schmid–Leiman/higher-order, EAP + higher-order SEM out of scope, deep-level convergence); ¶5 = future directions (concrete: ESEM/polychoric comparability + `boot_edges`) + availability (CRAN, GitHub, docs, MIT). All four beats present. **PASS.**
- **AC2 — scope-beat cites + Discussion citekey traceability.** `@schmid1957` (×1) and `@yung1999` (×2) appear in the scope-boundary paragraph (¶4). All Discussion citekeys — schmid1957, yung1999, beauducel2024, hu1999, kotov2017, michelini2019, williams2025 — are in the citekey→note map (Decisions section) and each traces to an extraction-verified `cairn/references/` note. **PASS.**
- **AC3 — seeded-section citations only, no wording change.** `git diff master...HEAD -- manuscript.qmd` outside the Discussion shows exactly four hunks, each an inserted `[@key]` on an otherwise byte-identical line: `grice2001` (factor-score materialization), `kaiser1958` (varimax), `asparouhov2009` (ESEM engine), `revelle1979`+`ruscio2012a` (suggest_k criteria). No technical claim altered. Each cited note extraction-verified. **PASS.**
- **AC4 — new keys resolve, Crossref-checked.** All 7 new keys (`kaiser1958`, `asparouhov2009`, `grice2001`, `revelle1979`, `ruscio2012a`, `beauducel2024`, `williams2025`) defined in `references.bib`; each DOI verified live against the Crossref API during T3 (title/authors/journal/volume/issue/pages/year matched); `revelle1979` DOI recovered from Crossref. **PASS.**
- **AC5 — renders both formats, clean.** `quarto render manuscript.qmd` (Quarto 1.9.38, ackwards 0.1.1) produced `manuscript.pdf` + `manuscript.docx`; render-log scan found no citeproc "not found", no LaTeX error/undefined/missing; docx text confirms all 7 new references resolve with 0 literal `@citekey` leaks. **PASS.**
- **AC6 — diff scope.** `git diff --name-only master...HEAD` = `cairn/ROADMAP.md`, milestone file, `manuscript/manuscript.qmd`, `manuscript/references.bib` — no `R/`, `tests/`, `NAMESPACE`, or `DESCRIPTION`. **PASS.**

### Independent review (three lenses + scorer)
