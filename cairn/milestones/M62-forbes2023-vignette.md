<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M62: Worked Forbes (2023) AMH example vignette

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —
- **Branch/PR:** m62-forbes2023-vignette · https://github.com/jmgirard/ackwards/pull/63

## Goal
<!-- owner: plan · create -->

Ship a precomputed vignette that reproduces the Forbes (2023) AMH applied example end-to-end on
the bundled `forbes2023` dataset.

## Scope
<!-- owner: plan · create/amend-via-gate -->

**In:**
- New standalone vignette `vignettes/ackwards-forbes2023.Rmd.orig` → precomputed
  `ackwards-forbes2023.Rmd` + `vignettes/assets/` figures. Full applied example:
  `ackwards(forbes2023, k_max = 10, pairs = "all")`, hierarchy overview, skip-level edge
  highlights, the redundancy chase via `prune("redundant")` (default `direct` criterion; note the
  `adjacent` opt-in divergence, 37 vs 36 flagged), the pruned-factor diagram in Forbes's
  publication style, and a fidelity note pointing at `test-forbes-fidelity.R`.
- Cross-links: `ackwards-forbes.Rmd.orig` ↔ new vignette; `@seealso` in `forbes2023` roxygen
  (`R/data.R`).
- Riding along (plan-gate decision 2026-07-16): fix the `R/data.R:161` citation typo
  ("bass-ackward**s** method" → "bass-ackward method", matching the other six cites and Crossref).
- `_pkgdown.yml` articles index + `pkgdown::check_pkgdown()`; NEWS.md entry.

**Out:**
- Extended-tutorial extras on AMH (`suggest_k()` on the matrix, `n_obs = 3175` fit statistics,
  threshold-sensitivity sweep) → declined at the plan gate 2026-07-16 (scope = reproduction, not a
  second tutorial); revisit conversationally if wanted later.
- BRM manuscript Intro/Discussion prose → existing ROADMAP candidate row (M56 Out).
- Rewriting `ackwards-forbes.Rmd` onto AMH → declined at the plan gate (bfi25 concept vignette
  stays intact); only a cross-link is added.
- Any new exported behavior or plotting feature — if the diagram needs one, that is an amendment
  gate, not silent scope growth.

**Test scope note:** docs-only milestone — the vignette's headline numbers are already locked by
`test-forbes-fidelity.R` (fixture `forbes2023_amh.rds`, provenance: Forbes's reference
implementation, OSF `pcwm8`); no new tests expected unless a code change proves necessary
(→ amendment + test).

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [x] AC1: `vignettes/ackwards-forbes2023.Rmd.orig` exists; the precomputed
      `ackwards-forbes2023.Rmd` + its `vignettes/assets/ackwards-forbes2023-*` figures are
      committed; `R CMD check` builds the vignette cleanly.
- [x] AC2: The rendered `.Rmd` shows the fixture-pinned reproduction numbers — 10-level hierarchy
      over 155 variables, 45 level-pairs under `pairs = "all"`, **37 of 55** factors flagged
      redundant under the default `direct` criterion (and 36 under `adjacent`, if shown) —
      grep-verifiable against `tests/testthat/test-forbes-fidelity.R:254-266`. Source: Forbes
      (2023), \doi{10.1037/met0000546}, via the `forbes2023_amh.rds` fixture (provenance:
      `data-raw/forbes2023.R`).
- [x] AC3: The pruned-factor diagram (`autoplot(…, drop_pruned = TRUE)`, publication style per
      `ackwards-forbes.Rmd.orig:246-281`) renders in the vignette with its figure asset committed.
- [x] AC4: Bidirectional cross-links present (`ackwards-forbes` ↔ `ackwards-forbes2023`);
      `forbes2023` roxygen gains a `@seealso`/vignette pointer; `grep -rn "extension of
      Goldberg's bass-ackwards" R/` returns no matches (the Forbes-citation typo fixed; the
      Waller 2007 title legitimately uses "Bass-Ackwards" and stays); `devtools::document()` run.
- [x] AC5: `_pkgdown.yml` articles index lists the new vignette; `pkgdown::check_pkgdown()` clean.
- [x] AC6: NEWS.md entry for the new vignette.
- [x] AC7: `Rscript tools/dod-gate.R` clean, and the PR diff contains no churn-only regenerated
      vignettes/assets (M61 lesson: `git checkout --` untouched precompute outputs).

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1, T3
- AC2 → T1, T3
- AC3 → T1, T3
- AC4 → T2
- AC5 → T4
- AC6 → T4
- AC7 → T3, T5

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1: Draft `vignettes/ackwards-forbes2023.Rmd.orig` — the full applied example: dataset
      intro (AMH, N = 3,175, CC-BY provenance), fit `ackwards(forbes2023, k_max = 10,
      pairs = "all")`, hierarchy overview, strongest skip-level edges, `prune("redundant")`
      with the direct-vs-adjacent note (37 vs 36), pruned-factor diagram (labeled + compressed
      variants), fidelity note, references (verify titles against Crossref — M56 lesson).
- [x] T2: Cross-links + typo: pointer from `ackwards-forbes.Rmd.orig` to the new vignette (and
      back), `@seealso`/vignette pointer in the `forbes2023` roxygen (`R/data.R`), fix
      `R/data.R:161` Forbes-citation "bass-ackwards"; `devtools::document()`.
- [x] T3: Run `Rscript vignettes/precompute.R`; verify the pinned numbers appear in the generated
      `.Rmd`; `git checkout --` all untouched regenerated vignettes/assets so the diff stays
      scoped (M61 lesson).
- [x] T4: Add the article to `_pkgdown.yml` ("The workflow" group, after `ackwards-forbes`); run
      `pkgdown::check_pkgdown()`; add the NEWS.md entry.
- [x] T5: Final gate: `Rscript tools/dod-gate.R`.

## Work log
<!-- owner: any skill · append-only -->

- 2026-07-16: created by /milestone-plan. Promotes the 2026-07-12 candidate row "Wire `forbes2023`
  into a vignette" (lineage: M54 Out). Plan-gate decisions: new standalone vignette; full applied
  example; ACs pinned to fixture values; data.R typo fixed in-milestone.
- 2026-07-16: /milestone-implement start — Forbes (2023) title verified via Crossref API
  ("bass-ackward method", singular; confirms the T2 fix). Branch m62-forbes2023-vignette.
- 2026-07-16: T1 done — vignette drafted against live numbers (fit 0.6s; 1320 edges/45 pairs;
  37/55 direct vs 36 adjacent; d4 chain m4f4→…→m10f4), test-knit clean, both diagrams verified
  visually. Crossref check found repo-wide Goldberg (2006) cites omit the published subtitle
  ("…The development of hierarchical factor structures from the top down") — full title used in
  the new vignette; repo-wide fix left out of scope (would regenerate all vignettes), surface at
  recap as a candidate.
- 2026-07-16: T2 done. Minor AC4 wording amendment: the planned repo-wide grep for
  "bass-ackwards method" would flag two *correct* uses (Waller 2007's published title, M56-verified;
  a prose mention in comparability.R) — recipe narrowed to the Forbes-citation pattern
  ("extension of Goldberg's bass-ackwards"), which was the actual defect. Intent unchanged.
- 2026-07-16: T3 done — precompute run; 5 churn-only .Rmd + 3 PNGs reverted (M61 lesson);
  ackwards-forbes.Rmd diff is exactly the 4-line cross-link; pinned numbers verified in the
  generated forbes2023 vignette (45 pairs, 37/55 vs 36, 13 chain-retained + 5 non-chain).
- 2026-07-16: T4 done — pkgdown articles index updated (check_pkgdown clean); NEWS entry added
  under a new "(development version)" heading (0.1.1 is submitted, its section is frozen).
- 2026-07-16: T5 done — DoD gate PASSED (check 0/0/0, coverage 100%, style/lint clean, pkgdown
  index complete). All tasks complete; status → review.

## Decisions
<!-- owner: implement / review · append-only -->

## Review
<!-- owner: review · exclusive -->

### Acceptance-criteria evidence (2026-07-16, all by command on the branch)

- AC1: `ls` — `.Rmd.orig` (10,350 B), generated `.Rmd` (24,032 B), both PNGs committed;
  `devtools::check()` 0/0/0 via `tools/dod-gate.R` this session (vignette built cleanly). ✔
- AC2: greps on the generated `.Rmd` — "45 level pairs" ×1, "37 are [flagged]" ×1, "flags 36
  components" ×1; values match `test-forbes-fidelity.R:254` (37L direct) and `:266` (36L
  adjacent); 45 pairs = choose(10,2) per the fixture loop. ✔
- AC3: both `assets/ackwards-forbes2023-pruned-diagram*` PNGs referenced at `.Rmd:740`/`:755`;
  rendered diagrams visually verified during implement (publication style, level gaps visible). ✔
- AC4: cross-link greps — `ackwards-forbes2023` in `ackwards-forbes.Rmd.orig` + generated `.Rmd`
  (1 each); back-link `vignette("ackwards-forbes")` in new `.orig` (1); `\seealso` in
  `man/forbes2023.Rd` (1); `grep "extension of Goldberg's bass-ackwards" R/` → 0. `document()`
  re-run at review: no diff. ✔
- AC5: `_pkgdown.yml:33` lists `ackwards-forbes2023`; `pkgdown::check_pkgdown()` fresh at review:
  "No problems found". ✔
- AC6: NEWS.md:4 — entry under new "(development version)" heading. ✔
- AC7: `Rscript tools/dod-gate.R` PASSED this session (check 0/0/0, coverage 100%, style/lint
  clean, pkgdown index complete); `git diff --name-only master..HEAD` = 12 files, all intentional
  (no churn-only regenerated vignettes/assets — 5 .Rmd + 3 PNGs reverted at T3). ✔

### Consistency gate (2026-07-16)

- `cairn_validate.py`: all checks passed (87 pre-existing advisory dangling-ID warnings, legacy
  M1–M53 refs).
- No principle change (`Principles touched: —`) → `cairn_impact` skipped.
- Profile consistency-gate slot: `document()` no diff ✔; generated files regenerated not
  hand-edited (no-diff check) ✔; README untouched by the diff (0 files) ✔;
  `pkgdown::check_pkgdown()` ✔; NEWS entry present ✔; no new top-level files (new `.orig`
  covered by `.Rbuildignore:21` pattern; check 0 NOTEs) ✔; full `check()` clean via gate ✔.

### Independent review (2026-07-16, three fresh-context lenses + scorer protocol)

- **[O] diff-bug**: no findings. Independently re-knit the `.orig` (committed `.Rmd`
  byte-identical modulo gt's random div id; both PNGs byte-identical); verified every vignette
  claim against test-forbes-fidelity.R (37/36/45/1320/13-retained ids/d4 chain/7-of-54) and all
  chunk API usage against R/ sources; citations match master-verified references.
- **[S] blame-history**: no findings. data.R plural was M54 copy-drift (cedcd37/541eb5f), not a
  deliberate choice — the fix aligns with the repo-wide verified form; NEWS dev-heading is the
  correct convention with 0.1.1 pending at CRAN.
- **[S] prior-PR-comments**: no prior-PR evidence (PRs #34–#55 on these files carry only codecov
  bot comments) — clean no-op, zero findings.
- Scorer: skipped — zero surviving findings to score.
- Logged sub-threshold/self-excluded observations (surfaced, not actioned): (1) pre-existing
  `R/ackwards.R:429-435` message tells matrix-input users `n_obs` "enables" chi-square/RMSEA/TLI,
  but the PCA engine never computes those — pre-existing, not introduced here; candidate-worthy.
  (2) repo-wide Goldberg (2006) cites omit the published subtitle (T1 work-log) — candidate-worthy.
