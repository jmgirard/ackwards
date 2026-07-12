<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M54: Export `amh_cor` as a bundled dataset

- **Status:** planned   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** —   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal

Ship Forbes's (2023) 155-variable AMH Spearman correlation matrix as an exported,
documented, CC-BY-attributed user dataset `amh_cor`, extending the M53 test-only fixture.

## Scope

**In:** A new bundled dataset `data/amh_cor.rda` (the 155×155 AMH matrix, row/col
names preserved); its roxygen `@source`/`@references` doc block in `R/data.R`; a
`data-raw/amh_cor.R` generator that downloads the CC-BY matrix from OSF (guid
`s9bjz`, md5-pinned) and writes **both** `data/amh_cor.rda` and the slimmed
fidelity fixture from one download; single-sourcing the fidelity test's matrix
onto the exported `amh_cor` (fixture slimmed to expected-values-only); Forbes
added to DESCRIPTION as data-scoped `cph`; `LICENSE.note` extended to cover
`data/amh_cor.rda`; `_pkgdown.yml` Data section + NEWS entry; dataset validity
tests.

**Out:**
- Wiring `amh_cor` into a vignette (the Forbes vignette runs on the seed-regenerated
  simulations) → ROADMAP candidate.
- Any change to the fidelity oracle values or `prune()`/edge algorithms — M53's
  numbers stand unchanged; this milestone only relocates where the matrix is read from.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [ ] `data/amh_cor.rda` exists and loads lazily; `amh_cor` is a 155×155 numeric
      matrix, symmetric (`isSymmetric(unname(amh_cor))`), unit diagonal, with
      identical non-null row and column names. (evidence: a test asserting all four.)
- [ ] `amh_cor` is exported and documented: `?amh_cor` renders, `@source` names the
      OSF project + CC-BY 4.0, and `amh_cor` appears in the `_pkgdown.yml` Data
      section (`pkgdown::check_pkgdown()` clean).
- [ ] The Forbes fidelity test reads its matrix from the exported `amh_cor` (not a
      matrix embedded in the fixture) and still reproduces Forbes to ≤ 1e-12 across
      all 45 level-pairs at `k_max = 10` — matching the M53 result — per
      `Forbes (2023) doi:10.1037/met0000546` reference values in `forbes2023_amh.rds`.
- [ ] `data-raw/amh_cor.R` regenerates `data/amh_cor.rda` and the slimmed fixture
      from a single md5-pinned OSF download; the shipped `amh_cor` and the fixture's
      expected values trace to the same source matrix. (evidence: script runs clean;
      md5 guard present.)
- [ ] Licensing is complete and correct: DESCRIPTION lists Forbes as `cph` scoped by
      `comment` to the bundled matrix; `LICENSE.note` covers `data/amh_cor.rda`;
      `devtools::check()` clean (0 errors/warnings/notes), including no license NOTE.

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T2, T6
- AC2 → T3, T5
- AC3 → T4, T6
- AC4 → T1, T6
- AC5 → T5, T6

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [ ] T1: Write `data-raw/amh_cor.R` — adapt master's `data-raw/forbes2023_amh.R`
      to, from one md5-pinned OSF download, (a) `usethis::use_data(amh_cor, compress = "xz")`
      and (b) write the slimmed `forbes2023_amh.rds` (expected `comp_corr`/`cong`/
      `corr_chase`/`k_max` + provenance attr, **no** `R` matrix). Run it to produce
      both artifacts.
- [ ] T2: Regenerate `data/amh_cor.rda` via T1; confirm dims/symmetry/dimnames by hand
      before committing (never hand-edit `.rda`).
- [ ] T3: Add the `amh_cor` roxygen block to `R/data.R` (adapt the branch's block:
      `@format`, `@details` with `n_obs = 3175` note, `@source` CC-BY + OSF,
      `@references` both Forbes papers, `\donttest{}` example). `devtools::document()`.
- [ ] T4: Repoint `test-forbes-fidelity.R` AMH block to read the matrix from `amh_cor`;
      drop the `$amh$R` read; keep all expected-value assertions. Add/keep the
      direct-vs-adjacent prune pins unchanged.
- [ ] T5: DESCRIPTION — add Forbes `cph` (comment-scoped to the matrix); extend
      `LICENSE.note` to name `data/amh_cor.rda` under CC-BY; add `amh_cor` to
      `_pkgdown.yml` Data section; NEWS.md dataset entry.
- [ ] T6: Add `test-data-amh.R` (dims/symmetry/diag/dimnames validity); run
      `Rscript tools/dod-gate.R` (check → coverage → style → lint → pkgdown).

<!-- (RB tripwire: irreversible-api) The exported object name `amh_cor` is a
     permanent public-API contract; name resolved at the plan gate (2026-07-12),
     no Fable escalation warranted for a conventional data export. -->

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan. Extends M53 (fixture → exported dataset).
  Gate decisions: single-source the matrix; Forbes as data-scoped `cph` + LICENSE.note;
  name `amh_cor`. Source material exists on local `amh-fidelity` branch (adapt, don't
  cherry-pick — master's LICENSE.note/data-raw are newer).

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->
