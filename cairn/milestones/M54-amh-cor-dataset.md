<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M54: Export `forbes2023` as a bundled dataset

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Branch/PR:** m54-amh-cor-dataset · https://github.com/jmgirard/ackwards/pull/54   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal

Ship Forbes's (2023) 155-variable AMH Spearman correlation matrix as an exported,
documented, CC-BY-attributed user dataset `forbes2023`, extending the M53 test-only fixture.

## Scope

**In:** New bundled dataset `data/forbes2023.rda` (155×155 AMH matrix, names preserved) + roxygen
doc block; a `data-raw/forbes2023.R` generator (md5-pinned OSF download) writing **both** the
dataset and the slimmed fidelity fixture; single-sourcing the fidelity test's matrix onto exported
`forbes2023`; Forbes as data-scoped `cph`; `LICENSE.note`, `_pkgdown.yml` Data, NEWS; validity tests.

**Out:** Wiring `forbes2023` into a vignette → ROADMAP candidate. Any change to the fidelity oracle
or `prune()`/edge algorithms — M53's numbers stand; this only relocates where the matrix is read from.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [x] `data/forbes2023.rda` exists and loads lazily; `forbes2023` is a 155×155 numeric
      matrix, symmetric (`isSymmetric(unname(forbes2023))`), unit diagonal, with
      identical non-null row and column names. (evidence: a test asserting all four.)
- [x] `forbes2023` is exported and documented: `?forbes2023` renders, `@source` names the
      OSF project + CC-BY 4.0, and `forbes2023` appears in the `_pkgdown.yml` Data
      section (`pkgdown::check_pkgdown()` clean).
- [x] The Forbes fidelity test reads its matrix from the exported `forbes2023` (not a
      matrix embedded in the fixture) and still reproduces Forbes to ≤ 1e-12 across
      all 45 level-pairs at `k_max = 10` — matching the M53 result — per
      `Forbes (2023) doi:10.1037/met0000546` reference values in `forbes2023_amh.rds`.
- [x] `data-raw/forbes2023.R` regenerates `data/forbes2023.rda` and the slimmed fixture
      from a single md5-pinned OSF download; the shipped `forbes2023` and the fixture's
      expected values trace to the same source matrix. (evidence: script runs clean;
      md5 guard present.)
- [x] Licensing is complete and correct: DESCRIPTION lists Forbes as `cph` scoped by
      `comment` to the bundled matrix; `LICENSE.note` covers `data/forbes2023.rda`;
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

- [x] T1: Write `data-raw/forbes2023.R` — adapt master's `data-raw/forbes2023_amh.R`
      to, from one md5-pinned OSF download, (a) `usethis::use_data(forbes2023, compress = "xz")`
      and (b) write the slimmed `forbes2023_amh.rds` (expected `comp_corr`/`cong`/
      `corr_chase`/`k_max` + provenance attr, **no** `R` matrix). Run it to produce
      both artifacts.
- [x] T2: Regenerate `data/forbes2023.rda` via T1; confirm dims/symmetry/dimnames by hand
      before committing (never hand-edit `.rda`).
- [x] T3: Add the `forbes2023` roxygen block to `R/data.R` (adapt the branch's block:
      `@format`, `@details` with `n_obs = 3175` note, `@source` CC-BY + OSF,
      `@references` both Forbes papers, `\donttest{}` example). `devtools::document()`.
- [x] T4: Repoint `test-forbes-fidelity.R` AMH block to read the matrix from `forbes2023`;
      drop the `$amh$R` read; keep all expected-value assertions. Add/keep the
      direct-vs-adjacent prune pins unchanged.
- [x] T5: DESCRIPTION — add Forbes `cph` (comment-scoped to the matrix); extend
      `LICENSE.note` to name `data/forbes2023.rda` under CC-BY; add `forbes2023` to
      `_pkgdown.yml` Data section; NEWS.md dataset entry.
- [x] T6: Add `forbes2023` validity test to `test-data.R` (dims/symmetry/diag/dimnames); run
      `Rscript tools/dod-gate.R` (check → coverage → style → lint → pkgdown).

<!-- (RB tripwire: irreversible-api) The exported object name `forbes2023` is a
     permanent public-API contract; name resolved at the plan gate (2026-07-12),
     no Fable escalation warranted for a conventional data export. -->

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan. Extends M53 (fixture → exported dataset).
  Gate decisions: single-source the matrix; Forbes as data-scoped `cph` + LICENSE.note;
  name `forbes2023`. Source material exists on local `amh-fidelity` branch (adapt, don't
  cherry-pick — master's LICENSE.note/data-raw are newer).
- 2026-07-12: T1/T2 — consolidated `data-raw/amh_cor.R` writes both `data/amh_cor.rda`
  and the slimmed fixture from one md5-pinned OSF download (`c1dd9eca…`, verified
  byte-identical to the M53-validated matrix); removed superseded `data-raw/forbes2023_amh.R`;
  fixture 117 KB → 14 KB (matrix dropped).
- 2026-07-12: T3 — `forbes2023` roxygen block added to `R/data.R` (documented; `man/amh_cor.Rd`
  generated; NAMESPACE unchanged — bundled data is user-visible via LazyData, not `export()`).
- 2026-07-12: T4 — fidelity test AMH block reads matrix from exported `forbes2023` (not `$amh$R`);
  full `forbes-fidelity` suite green, still reproduces Forbes ≤1e-12 across all 45 pairs.
- 2026-07-12: T5 — Forbes added to DESCRIPTION as data-scoped `cph` (Authors@R parses, verified);
  `LICENSE.note` now covers `data/amh_cor.rda`; `forbes2023` in `_pkgdown.yml` Data (`check_pkgdown()`
  clean); NEWS dataset bullet added.
- 2026-07-12: T6 — `forbes2023` validity test added to `test-data.R`; `Rscript tools/dod-gate.R`
  green (check clean, 0 err/0 warn/0 note incl. no license NOTE, coverage 100%, style/lint clean, pkgdown complete).
  All tasks done → status `review`.
- 2026-07-12: review-gate "adjust first" — user renamed the dataset `amh_cor` → `forbes2023`
  (author-year style). Swept every source/doc/test/DESCRIPTION/LICENSE.note/pkgdown/NEWS ref;
  `data-raw/amh_cor.R` → `data-raw/forbes2023.R`; regenerated `data/forbes2023.rda`; re-ran the
  full gate — green. (Earlier work-log/task lines keep the original `amh_cor` name as history;
  the fixture `forbes2023_amh.rds` and the file slug are unchanged.)

## Decisions

- 2026-07-12 (T5): Forbes recorded in `Authors@R` as `role = "cph"` **scoped by `comment`** to
  `data/amh_cor.rda` only — she holds copyright in the bundled CC-BY matrix, not the MIT-licensed
  package code. Standard R idiom for vendoring third-party CC-BY data; keeps CRAN attribution
  machine-readable. Not promoted to a D-entry (milestone-local, no cross-cutting effect).
- 2026-07-12 (post-implement, user override): dataset renamed `amh_cor` → `forbes2023` (author-year),
  superseding the plan-gate name and the T5 entry's `data/amh_cor.rda` path (now `data/forbes2023.rda`).
- 2026-07-12 (review F2, owner-approved): `forbes2023` classified as a fidelity/reproduction
  dataset outside DESIGN §14.41(b)'s declined-third-*teaching*-dataset scope. Added a §14.41(b)
  clarification and reworded the roxygen so it no longer reads as a third teaching foil.
<!-- owner: implement / review · append-only; milestone-local -->

## Review
<!-- owner: review · exclusive -->

**Reviewed 2026-07-12 · PR #54 · branch up to date with master.** AC evidence (fresh):
- **AC1** — 155×155 numeric matrix, symmetric, unit diagonal, rownames==colnames (non-null/unique);
  `test-data.R` asserts all four; `data|forbes-fidelity` suites 16/16 pass (0 fail/warn/skip).
- **AC2** — `?forbes2023` renders, `@source` names OSF pcwm8 + CC-BY 4.0; in `_pkgdown.yml` Data;
  `check_pkgdown()` clean; `--run-donttest` example ran OK under check.
- **AC3** — fidelity test reads matrix from exported `forbes2023` (fixture no longer carries `$amh$R`);
  reproduces Forbes ≤1e-12 across all 45 pairs at k_max=10; oracle unchanged (Forbes's own ref impl).
- **AC4** — `data-raw/forbes2023.R` ran clean, wrote both artifacts; md5 guard present; fresh download
  byte-identical to the M53-validated matrix.
- **AC5** — `devtools::check(vignettes = TRUE)`: Status OK, zero errors/warnings/notes (no license
  NOTE); `Authors@R` parses with Forbes `cph` (comment-scoped); `LICENSE.note` covers the dataset.

Consistency gate: `cairn_validate.py` exit 0; `document()` no diff; `check_pkgdown()` clean; NEWS has
a user-facing entry; no new top-level files; no DESIGN principle changed. CI (PR #54): full matrix green.

Independent review — [O] diff-bug + [S] blame-history + [S] scorer. Core design cleared: blame
reviewer verified byte-for-byte that the slimmed fixture's oracle values are `identical()` to M53's
and the shipped matrix is `identical()` to the previously-embedded one (single-sourcing neither
weakened the oracle nor broke test self-containment). Findings:
- **F1 (85, fixed):** `ROADMAP.md` candidate bullet still said `amh_cor` — rename miss. Fixed.
- **F3 (85, fixed):** `LICENSE.note` "Both are the AMH correlation matrix" → "Both derive from …"
  (one file holds derived expected values, not the matrix).
- **F2 (50, logged):** §14.41(b) declines a third *teaching* dataset; `forbes2023` is a
  *fidelity/reproduction* dataset (arguably out of scope) but roxygen "complements sim16/bfi25"
  echoes the pattern — surfaced to owner at the merge gate.
