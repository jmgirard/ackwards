<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M63: User-facing text accuracy pass (Goldberg 2006 citations + PCA n_obs message)

- **Status:** planned   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** —   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create -->

Correct the two user-facing text defects the M62 review surfaced: bring every
Goldberg (2006) full reference to the Crossref-verified published title, and
stop the PCA matrix-input `n_obs` message (and its roxygen twin) promising fit
statistics PCA never computes.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:**
- Goldberg (2006) full-reference lines ("Doing it all …") in roxygen
  (`R/ackwards.R:145`) and four vignette sources
  (`ackwards-engines.Rmd.orig:583`, `ackwards-forbes.Rmd.orig:527` — subtitle
  missing; `ackwards-intro.Rmd.orig:395`, `ackwards-girard.Rmd.orig:336` —
  lowercase "bass-ackwards:" vs published casing), unified to the exact form in
  `cairn/references/goldberg2006.md` (Crossref-verified,
  doi:10.1016/j.jrp.2006.01.001).
- The PCA matrix-input message `R/ackwards.R:429-435` and the `@param n_obs`
  roxygen `R/ackwards.R:45-46`: remove the false claim that `n_obs`
  enables/disables chi-square/RMSEA/TLI for PCA (`pca_levels()` never receives
  `n_obs`; its fit slot is eigenvalues only — `R/engine_pca.R:10-71`). Keep the
  message, corrected: `n_obs` stored as `NA`, recorded as metadata only
  (gate decision 2026-07-16).
- Regeneration in sync: `devtools::document()` → `man/ackwards.Rd`;
  `Rscript vignettes/precompute.R` → the four touched `.Rmd` + assets.

**Out:**
- In-text cites (README.md, NEWS.md, vignette prose) — no subtitle expected in
  author-year cites; nothing to fix.
- BRM manuscript references — author-owned prose candidate (stays on ROADMAP).
- Uncommitted `vignettes/*.html` build leftovers — not tracked, nothing to do.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [ ] AC1 — Every full-reference line matching `Doing it all` under `R/`,
      `man/`, and `vignettes/*.Rmd{,.orig}` carries the exact published title
      "Doing it all Bass-Ackwards: The development of hierarchical factor
      structures from the top down" (source: goldberg2006.md, Crossref-verified);
      grep evidence enumerating every match. (Anchor is the defective-context
      phrase "Doing it all", which cannot false-flag Waller 2007 — LESSONS
      2026-07-16/M62.)
- [ ] AC2 — Neither the cli message nor the `@param n_obs` roxygen claims
      `n_obs` enables/disables chi-square/RMSEA/TLI for PCA; a test asserts the
      corrected message fires (and states metadata-only storage) on PCA
      matrix input without `n_obs`.
- [ ] AC3 — Regenerated artifacts committed in sync: `man/ackwards.Rd` matches
      roxygen; regenerated vignette diff scoped to the four touched vignettes +
      their assets (untouched vignettes restored per LESSONS 2026-07-16/M61).
- [ ] AC4 — NEWS.md entry for the user-visible message correction.
- [ ] AC5 — `Rscript tools/dod-gate.R` passes (check 0 err/0 warn/0 note,
      coverage, style, lint, pkgdown).

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1, T3, T4
- AC2 → T1, T2
- AC3 → T1, T4
- AC4 → T5
- AC5 → T6

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [ ] T1 — Fix `R/ackwards.R`: full reference at :145 to the published title;
      `@param n_obs` at :45-46 (and the matrix-input Details bullet at :196-197
      if it repeats the claim) reworded to metadata-only for PCA; run
      `devtools::document()` and commit `man/ackwards.Rd` with it.
- [ ] T2 — Correct the cli message at `R/ackwards.R:429-435` (drop the
      enable-fit-stats promise; state metadata-only storage); add a test
      asserting the corrected message on PCA matrix-input without `n_obs`
      (none exists today — new coverage).
- [ ] T3 — Fix the four `.Rmd.orig` reference lines (engines:583, forbes:527
      add subtitle; intro:395, girard:336 casing) to the goldberg2006.md form.
- [ ] T4 — `Rscript vignettes/precompute.R`; `git checkout --` the untouched
      vignettes/assets (forbes2023, ordinal, suggest-k, visualization); commit
      regenerated `.Rmd` + assets.
- [ ] T5 — NEWS.md entry (message correction; citation completion can share
      the bullet).
- [ ] T6 — `Rscript tools/dod-gate.R` clean at the gate.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-16: created by /milestone-plan; absorbs the two M62-review candidate
  rows (citation hygiene + matrix-input message); bundling + keep-corrected-message
  decided at the plan gate.

## Decisions
<!-- owner: implement / review · append-only -->

## Review
<!-- owner: review · exclusive; EXEMPT from the 150-line cap (M55) -->
