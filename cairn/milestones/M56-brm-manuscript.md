<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M56: Behavior Research Methods manuscript (reproducible Quarto scaffold + worked example)

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M54 (forbes2023 dataset) already done -->
- **Branch/PR:** `m56-brm-manuscript` · [PR #56](https://github.com/jmgirard/ackwards/pull/56)   <!-- owner: implement (branch) / review (PR URL) · create -->

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

- [x] AC1. `manuscript/` exists and is excluded from the R build: `^manuscript$` is in `.Rbuildignore`, and `devtools::check()` is 0 err / 0 warn / 0 note (manuscript never enters the tarball).
- [x] AC2. The apaquarto (APA-7) document renders reproducibly to **both PDF and docx** via a single documented `quarto render` from a clean state (toolchain: Quarto ≥ 1.4, tinytex — both present locally).
- [x] AC3. Front-matter is correct (title; author "Jeffrey M. Girard" with ORCID 0000-0002-7359-3746 + affiliation; APA-7/BRM options), and the `.bib` bibliography resolves with **zero unresolved citations**, including the primary sources — Goldberg (2006) doi:10.1016/j.jrp.2006.01.001, Forbes (2023) doi:10.1037/met0000546, Waller (2007) doi:10.1016/j.jrp.2006.08.005 — plus software cites (R, psych, lavaan, ackwards).
- [x] AC4. A bfi25 walkthrough section computes **live from `ackwards`** (no hardcoded numbers), producing ≥ 1 figure and ≥ 1 table from package output.
- [x] AC5. The Forbes-AMH centerpiece computes **live from the bundled `forbes2023` dataset**, and ≥ 1 reported value is checked in-document against Forbes's published reference — anchored to the existing fidelity oracle (`tests/testthat/test-forbes-fidelity.R`, OSF `pcwm8`; Forbes 2023 doi:10.1037/met0000546).
- [x] AC6. Abstract, Statement of Need, and Method sections are drafted (seeded from DESIGN.md + existing vignettes); Intro and Discussion are present as structured stubs (headings + bullet guidance), clearly marked as author-owned.
- [x] AC7. `manuscript/README.md` documents the render command, required toolchain, and records the package version / `sessionInfo` used, so the manuscript reproduces standalone.

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
- [x] T8. Full clean render to PDF + docx; run `Rscript tools/dod-gate.R` (or at minimum `devtools::check()`) to confirm 0 err/0 warn/0 note with `manuscript/` excluded; commit rendered outputs (or document their reproduction).

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan. Venue BRM; Quarto+apaquarto (user pref); in-repo `manuscript/`, `.Rbuildignore`d; deliverable = scaffold + wired example (bfi25 walkthrough + AMH centerpiece) + seeded prose; PDF+docx; reuse vignette computations, vignette candidate row stays separate.
- 2026-07-12: gate — sole author (Girard, Univ. of Kansas); gitignore built outputs; title "ackwards: An R Package for Bass-Ackwards Hierarchical Structural Analysis". apaquarto 5.0.18 vendored.
- 2026-07-12: T1–T7 done. Manuscript renders to PDF (13pp) + docx; walkthrough (bfi25) Fig 1 + Table 1 and AMH centerpiece (forbes2023, k=10) Fig 2 compute live; all citations resolve; AMH fidelity anchored to test-forbes-fidelity.R (~1e-14). Callout blocks replaced with blockquote stubs to drop the fontawesome5 LaTeX dep; fontawesome5/newtx/etc. installed manually via CTAN tarballs (local TinyTeX was a sync ahead of all mirrors, so tlmgr auto-install refused). T8 (check + final render) pending.
- 2026-07-12: T8 done — `devtools::check()` 0 err/0 warn/0 note (manuscript excluded, no top-level-files NOTE); final clean render to PDF + docx, zero unresolved citations. Status → review. (No R code changed, so style/lint/coverage/pkgdown are unaffected.)

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->

**AC evidence (fresh, 2026-07-12).**
- AC1: `^manuscript$` at `.Rbuildignore:23`; fresh `devtools::check()` = 0 err/0 warn/0 note, Status OK (no top-level-files NOTE → manuscript excluded).
- AC2: `quarto render` from a `git clean -fdX` state exits 0, producing `manuscript.pdf` (13 pp) + `manuscript.docx`.
- AC3: PDF p.1 shows title, author, ORCID 0000-0002-7359-3746, "Department of Psychology, University of Kansas"; `bibliography: references.bib` wired; all 3 primary DOIs in `.bib`; text extraction finds zero `??`/`[?]`; full reference list renders.
- AC4: bfi25 walkthrough computes live (`dim` 875×25); Figure 1 (`autoplot`) + Table 1 (edges via `kable`) render (PDF pp.11–12).
- AC5: AMH centerpiece computes live (`ackwards(forbes2023, k_max=10, pairs="all")`, 155×155; prune 55 total/37 flagged); Figure 2 renders (PDF p.13); fidelity paragraph anchored to `test-forbes-fidelity.R` (~1e-14, OSF pcwm8).
- AC6: Abstract + Method + package sections drafted; Introduction + Discussion carry `[AUTHOR TO DRAFT]` stubs (2 occurrences).
- AC7: `manuscript/README.md` documents render command, toolchain (Quarto ≥1.4, TinyTeX), pinned versions (R 4.6.1, ackwards 0.1.1, Quarto 1.9.38, apaquarto 5.0.18).

**Consistency gate.** `cairn_validate.py` exit 0 (after fixing a `0/0/0` non-ISO token in T8). `devtools::document()` no diff. No R/ code, NAMESPACE, exports, or `_pkgdown.yml` changes → pkgdown/README/NEWS gates N/A (docs-only; manuscript is not user-facing package content). Coverage map complete: every AC → ≥1 existing task.

**Review-side fix.** Untracked stray endfloat build artifacts (`manuscript.fff`, `manuscript.ttt`) the initial gitignore missed; broadened `manuscript/` ignore to all LaTeX aux.

**Independent review (2 lenses + scorer).** [O] diff-bug + [S] blame-history, then [S] scorer.
- Blame-history: no findings (pure additions; no `.Rbuildignore`/`.gitignore` regression; no decision contradicted).
- Finding 1 (score 95, actioned+fixed): `references.bib` Waller (2007) title was wrong ("…by Bass-Ackward factor analysis"); Crossref + the package's own roxygen (`R/ackwards.R:148`) give "…by Goldberg's Bass-Ackwards method". Corrected + brace-protected; re-rendered.
- Finding 2 (score 78, below threshold but fixed as trivial): prose said "1,000 respondents" then the chunk printed 875 after `na.omit`; added a complete-case sentence.
- Finding 3 (score 55, logged, not actioned): "54 components" sits near a printed `total_factors = 55`; both correct (55 total factors vs 54 chase-eligible after excluding root m1f1) — clarity nit, left as-is.
- Out-of-scope discovery (follow-up, not M56): the same wrong Waller title exists pre-existing in `vignettes/ackwards-intro.Rmd.orig` / `ackwards-forbes.Rmd.orig`, inconsistent with roxygen — spawned as a separate task.
