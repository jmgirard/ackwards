<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M56: Behavior Research Methods manuscript (reproducible Quarto scaffold + worked example)

- **Status:** in-progress   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M54 (forbes2023 dataset) already done -->
- **Branch/PR:** `m56-brm-manuscript`   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create; a wrong goal returns to plan, never edited in place -->

Stand up a reproducible, in-repo Quarto/apaquarto manuscript for a *Behavior Research Methods* tool
paper on `ackwards` — rendering to PDF + Word, computing its worked examples live from the package,
with the verifiable prose sections seeded and the scholarly-argument sections left as author stubs.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:** a `manuscript/` directory (`.Rbuildignore`d) holding an APA-7 Quarto document (apaquarto
extension) that renders to **both PDF and docx**; correct BRM/APA front-matter (title, author +
ORCID + affiliation); a `references.bib` resolving all citations (primary method sources + software);
a **simple bfi25 walkthrough** and the **Forbes-AMH centerpiece**, both computing live from
`ackwards`/`forbes2023` (no hardcoded numbers), the AMH value fidelity-anchored to Forbes's published
reference; drafted **Abstract + Statement of Need + Method** prose (seeded from DESIGN/vignettes);
**Intro + Discussion** as structured stubs; a `manuscript/README.md` documenting the render command,
toolchain, and package version. `devtools::check()` stays clean throughout.

**Out:**
- Polished scholarly prose for Intro/Discussion → author-owned follow-on (candidate row, not a milestone).
- Actual journal submission (Editorial Manager upload, cover letter, author agreements) → owner-only, out.
- Wiring `forbes2023` into a *package vignette* → stays the separate candidate ROADMAP row (package docs ≠ journal manuscript).
- BRM house-style tweaks beyond APA-7 (journal-specific template polish) → deferred until submission time.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [ ] AC1. `manuscript/` exists and is excluded from the R build: `^manuscript$` is in `.Rbuildignore`, and `devtools::check()` is 0 err / 0 warn / 0 note (manuscript never enters the tarball).
- [ ] AC2. The apaquarto (APA-7) document renders reproducibly to **both PDF and docx** via a single documented `quarto render` from a clean state (toolchain: Quarto ≥ 1.4, tinytex — both present locally).
- [ ] AC3. Front-matter is correct (title; author "Jeffrey M. Girard" with ORCID 0000-0002-7359-3746 + affiliation; APA-7/BRM options), and the `.bib` bibliography resolves with **zero unresolved citations**, including the primary sources — Goldberg (2006) doi:10.1016/j.jrp.2006.01.001, Forbes (2023) doi:10.1037/met0000546, Waller (2007) doi:10.1016/j.jrp.2006.08.005 — plus software cites (R, psych, lavaan, ackwards).
- [ ] AC4. A bfi25 walkthrough section computes **live from `ackwards`** (no hardcoded numbers), producing ≥ 1 figure and ≥ 1 table from package output.
- [ ] AC5. The Forbes-AMH centerpiece computes **live from the bundled `forbes2023` dataset**, and ≥ 1 reported value is checked in-document against Forbes's published reference — anchored to the existing fidelity oracle (`tests/testthat/test-forbes-fidelity.R`, OSF `pcwm8`; Forbes 2023 doi:10.1037/met0000546).
- [ ] AC6. Abstract, Statement of Need, and Method sections are drafted (seeded from DESIGN.md + existing vignettes); Intro and Discussion are present as structured stubs (headings + bullet guidance), clearly marked as author-owned.
- [ ] AC7. `manuscript/README.md` documents the render command, required toolchain, and records the package version / `sessionInfo` used, so the manuscript reproduces standalone.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1, T8
- AC2 → T2, T8
- AC3 → T3
- AC4 → T4
- AC5 → T5
- AC6 → T6
- AC7 → T7

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1. Create `manuscript/`; add `^manuscript$` to `.Rbuildignore`; confirm the dir is excluded from `R CMD build`.
- [x] T2. Vendor the apaquarto extension into `manuscript/_extensions`; scaffold `manuscript.qmd` with APA-7 format targeting `apaquarto-pdf` + `apaquarto-docx`; confirm a trivial render of both formats.
- [x] T3. Author YAML front-matter (title, author/ORCID/affiliation, APA-7/BRM options); create `references.bib` with the primary method sources + software citations (`citation()` / `citation("psych")` / `citation("lavaan")`); wire `bibliography:` and confirm zero unresolved keys.
- [x] T4. Write the bfi25 walkthrough section with live `ackwards` chunks (≥ 1 figure + ≥ 1 table), reusing computations from `vignettes/ackwards-intro.Rmd` / `ackwards-girard.Rmd`.
- [x] T5. Write the Forbes-AMH centerpiece section computing from `forbes2023` (reusing `vignettes/ackwards-forbes.Rmd` computations), with an in-document fidelity check of ≥ 1 value against Forbes's reference via the existing oracle; ingest a `forbes2023` summary into `cairn/references/` if a citeable anchor is needed.
- [x] T6. Draft Abstract + Statement of Need + Method prose (seed from DESIGN.md §2/§3 + vignettes); add Intro + Discussion structured stubs with guidance bullets, marked author-owned.
- [x] T7. Write `manuscript/README.md`: render command, toolchain (Quarto ≥ 1.4, tinytex), and pinned package version / `sessionInfo`.
- [ ] T8. Full clean render to PDF + docx; run `Rscript tools/dod-gate.R` (or at minimum `devtools::check()`) to confirm 0/0/0 with `manuscript/` excluded; commit rendered outputs (or document their reproduction).

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan. Venue BRM; Quarto+apaquarto (user pref); in-repo `manuscript/`, `.Rbuildignore`d; deliverable = scaffold + wired example (bfi25 walkthrough + AMH centerpiece) + seeded prose; PDF+docx; reuse vignette computations, vignette candidate row stays separate.
- 2026-07-12: gate — sole author (Girard, Univ. of Kansas); gitignore built outputs; title "ackwards: An R Package for Bass-Ackwards Hierarchical Structural Analysis". apaquarto 5.0.18 vendored.
- 2026-07-12: T1–T7 done. Manuscript renders to PDF (13pp) + docx; walkthrough (bfi25) Fig 1 + Table 1 and AMH centerpiece (forbes2023, k=10) Fig 2 compute live; all citations resolve; AMH fidelity anchored to test-forbes-fidelity.R (~1e-14). Callout blocks replaced with blockquote stubs to drop the fontawesome5 LaTeX dep; fontawesome5/newtx/etc. installed manually via CTAN tarballs (local TinyTeX was a sync ahead of all mirrors, so tlmgr auto-install refused). T8 (check + final render) pending.

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->
