# M097: Cross-check coverage and the matrix behind each edge

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, IP2
- **Resolves:** —
- **Surface tier:** user-facing — it changes the shipped `ackwards()` help and the test suite
- **Branch/PR:** m097-cross-check-coverage

## Goal

The record of which correlation matrix each edge is exact for, and of which fits the algebra-vs-scores cross-check certifies, matches the code.

## Scope

**In:**
- New algebra-vs-scores tests on the PCA, EFA, and ESEM engines. They cover `cor = "spearman"` with scores from column ranks (plain and tied data), and `missing = "listwise"` with scores from the complete rows.
- A rewrite of DESIGN's first Known limitations entry, the §5.4 closing sentence, and the §5.3 Spearman sentence.
- The `ackwards()` and `compute_edges()` help state that an edge is exact for the stored `x$r`. They name the ESEM settings in which lavaan fits another matrix.
- A DESIGN sweep that names the matrix at each edge "exact" claim, IP1 included, with a D-entry for the IP1 wording.
- A NEWS entry.
- This absorbs two candidate rows. One is the cross-check missing paths (M93 review F3, M95 review R1). The other is DESIGN's unnamed "exact" claims (M95 review R9).

**Out:**
- README, DESCRIPTION, the vignettes, the manuscript, and the `prune()` help keep the general wording "the correlation matrix the fit uses". The plan gate kept them, because both ESEM settings warn at fit time.
- A user-facing route for sample-realized edges stays the existing `[low]` candidate row (D-038).
- D-004's text stays as written, because DECISIONS.md is append-only.

## Acceptance criteria

- [ ] AC1: Tests in `tests/testthat/test-compute_edges.R` cover `engine = "pca"`, `"efa"`, and `"esem"` with `k_max = 4`, in three cases. Case 1 is `cor = "spearman"` on `sim16`, and case 2 is `cor = "spearman"` on `round(sim16)`. In both, the scores route gets the column ranks (`apply(<data>, 2, rank)`). Case 3 is `missing = "listwise"` on `sim16` with planted missing values, and the scores route gets its complete rows. The edges come from `compute_edges(x$levels, x$r, pairs = "all")`. In each case, the `edge_method = "algebra"` edges agree within 1e-10 with the `edge_method = "scores"` edges on every level pair. The tests pass.
- [ ] AC2: The first entry of the "Known limitations" section of `cairn/DESIGN.md` lists the settings outside the algebra-vs-scores cross-check. These are `cor = "polychoric"` on any engine, and `missing = "fiml"` on PCA or EFA (the `psych::corFiml()` matrix) and on ESEM (lavaan's saturated-model matrix). They also include `missing = "pairwise"` on data with missing values on any engine, and correlation-matrix input. The entry names the `tests/testthat/test-compute_edges.R` tests that cover two more settings. These are complete-data `cor = "spearman"` fits (by column ranks) and `missing = "listwise"` fits (by complete rows). The entry says that a setting it neither lists nor names as tested is untested. Spearman with listwise deletion on data with missing values is its example. The closing sentence of §5.4 points to that entry. §5.3 has a paragraph "Where the `scores` route runs". It says that scores from the raw items do not reproduce a Spearman `R`, but scores from the column ranks do.
- [ ] AC3: Take the opening paragraph of the `ackwards()` help (`R/ackwards.R`, and so `man/ackwards.Rd`). It says that each edge is exact for the correlation matrix stored as `x$r`. It names the ESEM settings in which the correlation matrix lavaan fits differs from `x$r`. These are `cor = "spearman"`, and `cor = "pearson"` with `missing = "pairwise"` on data with missing values. Under that second setting, ML and MLR fit the complete rows, and WLSMV and ULSMV fit pairwise covariances. The `compute_edges()` help (`R/compute_edges.R`) names the same settings.
- [ ] AC4: Take each line that `grep -n -i -E 'exact|true correlation' cairn/DESIGN.md` returns. Every sentence in it that calls the edge algebra, an edge, or a score correlation exact or true names the correlation matrix it holds for. Text quoted in a correction note as superseded wording is excluded. A D-entry in `cairn/DECISIONS.md` records the IP1 wording change.
- [ ] AC5: `Rscript tools/dod-gate.R` exits 0 on the branch head.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T5

## Tasks

- [x] T1: Add the AC1 tests to `tests/testthat/test-compute_edges.R`. Plant the missing values with a fixed seed. Suppress the once-per-session ESEM Spearman warning (`R/ackwards.R:663`) and the ordinal advisory. Do not route these condition-wrapped fits through `cached()`. Show each case red once with its planted defect (unranked data for Spearman, the unreduced rows for listwise), and log the measured gap.
- [ ] T2: Rewrite DESIGN's first Known limitations entry, the §5.4 closing sentence, and the §5.3 Spearman sentence (AC2). Name the T1 tests by title.
- [ ] T3: Edit the opening paragraph of the `ackwards()` roxygen and the `compute_edges()` roxygen (AC3). Derive each ESEM claim from `R/ackwards.R:798` and `R/engine_esem.R:434`, and measure it in R before you write it. Run `devtools::document()` and `Rscript tools/check-prose.R` on both files. Re-read every split sentence against the code (M86 lesson). Add a NEWS entry for the help change.
- [ ] T4: Run the AC4 search and record each hit with its classification as a ledger in the work log. Name the matrix at each edge claim, IP1 included. Append the D-entry for the IP1 wording, and follow D-031's procedure for an IP change.
- [ ] T5: Run `Rscript tools/dod-gate.R` and fix what it reports.

## Work log

- 2026-10-01: created by /milestone-plan. Criteria audit (full mode, fresh Opus reader) gave 9 findings, all fixed at the gate. AC1 gained tied data, all pairs, and listwise. AC2 names the tests by file and adds the §5.3 Spearman sentence. AC3 gained a third ESEM setting (WLSMV or ULSMV, pairwise). AC4 holds over its search's hits only, excludes quoted superseded text, and its ledger moved to T4.
- 2026-10-01: plan gate chose new rank-score and complete-row tests over listing Spearman and listwise as uncovered. The measured agreement was about 1e-15 on all three engines. Falsified by a Spearman or listwise fit whose algebra and rank or complete-row scores disagree beyond 1e-10.
- 2026-10-01: plan gate chose the `ackwards()` and `compute_edges()` help over rewording all 15 sites the search finds. The general wording holds except in the ESEM settings AC3 names, and those warn at fit time. Falsified by a user who reads the general wording and misreads an ESEM Spearman or pairwise fit.
- 2026-10-01: plan gate absorbed the DESIGN "exact" row (M95 review R9) over keeping it separate, because it states the same fact in the same file.
- 2026-10-01: T1 done. Three tests (Spearman, tied Spearman, listwise) on PCA, EFA, and ESEM over all 6 level pairs agree to at most 2.4e-15. Red once with planted defects: raw data or unreduced rows failed all 18 pair checks per test, gaps 0.009 to 0.016. Suite 812 tests, 0 failed.

## Decisions

## Review
