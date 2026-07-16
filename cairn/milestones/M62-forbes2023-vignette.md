<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M62: Worked Forbes (2023) AMH example vignette

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —
- **Branch/PR:** m62-forbes2023-vignette

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

- [ ] AC1: `vignettes/ackwards-forbes2023.Rmd.orig` exists; the precomputed
      `ackwards-forbes2023.Rmd` + its `vignettes/assets/ackwards-forbes2023-*` figures are
      committed; `R CMD check` builds the vignette cleanly.
- [ ] AC2: The rendered `.Rmd` shows the fixture-pinned reproduction numbers — 10-level hierarchy
      over 155 variables, 45 level-pairs under `pairs = "all"`, **37 of 55** factors flagged
      redundant under the default `direct` criterion (and 36 under `adjacent`, if shown) —
      grep-verifiable against `tests/testthat/test-forbes-fidelity.R:254-266`. Source: Forbes
      (2023), \doi{10.1037/met0000546}, via the `forbes2023_amh.rds` fixture (provenance:
      `data-raw/forbes2023.R`).
- [ ] AC3: The pruned-factor diagram (`autoplot(…, drop_pruned = TRUE)`, publication style per
      `ackwards-forbes.Rmd.orig:246-281`) renders in the vignette with its figure asset committed.
- [ ] AC4: Bidirectional cross-links present (`ackwards-forbes` ↔ `ackwards-forbes2023`);
      `forbes2023` roxygen gains a `@seealso`/vignette pointer; `grep -rn "bass-ackwards method"
      R/` returns no matches (typo fixed); `devtools::document()` run.
- [ ] AC5: `_pkgdown.yml` articles index lists the new vignette; `pkgdown::check_pkgdown()` clean.
- [ ] AC6: NEWS.md entry for the new vignette.
- [ ] AC7: `Rscript tools/dod-gate.R` clean, and the PR diff contains no churn-only regenerated
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

- [ ] T1: Draft `vignettes/ackwards-forbes2023.Rmd.orig` — the full applied example: dataset
      intro (AMH, N = 3,175, CC-BY provenance), fit `ackwards(forbes2023, k_max = 10,
      pairs = "all")`, hierarchy overview, strongest skip-level edges, `prune("redundant")`
      with the direct-vs-adjacent note (37 vs 36), pruned-factor diagram (labeled + compressed
      variants), fidelity note, references (verify titles against Crossref — M56 lesson).
- [ ] T2: Cross-links + typo: pointer from `ackwards-forbes.Rmd.orig` to the new vignette (and
      back), `@seealso`/vignette pointer in the `forbes2023` roxygen (`R/data.R`), fix
      `R/data.R:161` "bass-ackwards method"; `devtools::document()`.
- [ ] T3: Run `Rscript vignettes/precompute.R`; verify the pinned numbers appear in the generated
      `.Rmd`; `git checkout --` all untouched regenerated vignettes/assets so the diff stays
      scoped (M61 lesson).
- [ ] T4: Add the article to `_pkgdown.yml` ("The workflow" group, after `ackwards-forbes`); run
      `pkgdown::check_pkgdown()`; add the NEWS.md entry.
- [ ] T5: Final gate: `Rscript tools/dod-gate.R`.

## Work log
<!-- owner: any skill · append-only -->

- 2026-07-16: created by /milestone-plan. Promotes the 2026-07-12 candidate row "Wire `forbes2023`
  into a vignette" (lineage: M54 Out). Plan-gate decisions: new standalone vignette; full applied
  example; ACs pinned to fixture values; data.R typo fixed in-milestone.
- 2026-07-16: /milestone-implement start — Forbes (2023) title verified via Crossref API
  ("bass-ackward method", singular; confirms the T2 fix). Branch m62-forbes2023-vignette.

## Decisions
<!-- owner: implement / review · append-only -->

## Review
<!-- owner: review · exclusive -->
