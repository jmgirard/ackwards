# M097: Cross-check coverage and the matrix behind each edge

- **Status:** review
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

- [x] AC1: Tests in `tests/testthat/test-compute_edges.R` cover `engine = "pca"`, `"efa"`, and `"esem"` with `k_max = 4`, in three cases. Case 1 is `cor = "spearman"` on `sim16`, and case 2 is `cor = "spearman"` on `round(sim16)`. In both, the scores route gets the column ranks (`apply(<data>, 2, rank)`). Case 3 is `missing = "listwise"` on `sim16` with planted missing values, and the scores route gets its complete rows. The edges come from `compute_edges(x$levels, x$r, pairs = "all")`. In each case, the `edge_method = "algebra"` edges agree within 1e-10 with the `edge_method = "scores"` edges on every level pair. The tests pass.
- [x] AC2: The first entry of the "Known limitations" section of `cairn/DESIGN.md` lists the settings outside the algebra-vs-scores cross-check. These are `cor = "polychoric"` on any engine, and `missing = "fiml"` on PCA or EFA (the `psych::corFiml()` matrix) and on ESEM (lavaan's saturated-model matrix). They also include `missing = "pairwise"` on data with missing values on any engine, and correlation-matrix input. The entry names the `tests/testthat/test-compute_edges.R` tests that cover two more settings. These are complete-data `cor = "spearman"` fits (by column ranks) and `missing = "listwise"` fits (by complete rows). The entry says that a setting it neither lists nor names as tested is untested. Spearman with listwise deletion on data with missing values is its example. The closing sentence of §5.4 points to that entry. §5.3 has a paragraph "Where the `scores` route runs". It says that scores from the raw items do not reproduce a Spearman `R`, but scores from the column ranks do.
- [x] AC3: Take the opening paragraph of the `ackwards()` help (`R/ackwards.R`, and so `man/ackwards.Rd`). It says that each edge is exact for the correlation matrix stored as `x$r`. It names every ESEM setting in which the correlation matrix lavaan fits differs from `x$r`. The code at `R/ackwards.R:798` and `R/engine_esem.R:434` determines these settings, and a sweep of `cor` × `missing` × the four accepted estimators, on complete data and on data with missing values, confirms them (lavaan 0.7.2, logged in the work log). These are `cor = "spearman"` with `missing = "pairwise"` or `"listwise"`, and `cor = "pearson"` with `missing = "pairwise"` on data with missing values. Under that second setting, ML and MLR fit the complete rows, and WLSMV fits pairwise covariances. The `compute_edges()` help (`R/compute_edges.R`) names the same settings.
- [x] AC4: Take each line that `grep -n -i -E 'exact|true correlation' cairn/DESIGN.md` returns. Every sentence in it that calls the edge algebra, an edge, or a score correlation exact or true names the correlation matrix it holds for. Text quoted in a correction note as superseded wording is excluded. A D-entry in `cairn/DECISIONS.md` records the IP1 wording change.
- [x] AC5: `Rscript tools/dod-gate.R` exits 0 on the branch head.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T5

## Tasks

- [x] T1: Add the AC1 tests to `tests/testthat/test-compute_edges.R`. Plant the missing values with a fixed seed. Suppress the once-per-session ESEM Spearman warning (`R/ackwards.R:663`) and the ordinal advisory. Do not route these condition-wrapped fits through `cached()`. Show each case red once with its planted defect (unranked data for Spearman, the unreduced rows for listwise), and log the measured gap.
- [x] T2: Rewrite DESIGN's first Known limitations entry, the §5.4 closing sentence, and the §5.3 Spearman sentence (AC2). Name the T1 tests by title.
- [x] T3: Edit the opening paragraph of the `ackwards()` roxygen and the `compute_edges()` roxygen (AC3). Derive each ESEM claim from `R/ackwards.R:798` and `R/engine_esem.R:434`, and measure it in R before you write it. Run `devtools::document()` and `Rscript tools/check-prose.R` on both files. Re-read every split sentence against the code (M86 lesson). Add a NEWS entry for the help change.
- [x] T4: Run the AC4 search and record each hit with its classification as a ledger in the work log. Name the matrix at each edge claim, IP1 included. Append the D-entry for the IP1 wording, and follow D-031's procedure for an IP change.
- [x] T5: Run `Rscript tools/dod-gate.R` and fix what it reports.

## Work log

- 2026-10-01: created by /milestone-plan. Criteria audit (full mode, fresh Opus reader) gave 9 findings, all fixed at the gate. AC1 gained tied data, all pairs, and listwise. AC2 names the tests by file and adds the §5.3 Spearman sentence. AC3 gained a third ESEM setting (WLSMV or ULSMV, pairwise). AC4 holds over its search's hits only, excludes quoted superseded text, and its ledger moved to T4.
- 2026-10-01: plan gate chose new rank-score and complete-row tests over listing Spearman and listwise as uncovered. The measured agreement was about 1e-15 on all three engines. Falsified by a Spearman or listwise fit whose algebra and rank or complete-row scores disagree beyond 1e-10.
- 2026-10-01: plan gate chose the `ackwards()` and `compute_edges()` help over rewording all 15 sites the search finds. The general wording holds except in the ESEM settings AC3 names, and those warn at fit time. Falsified by a user who reads the general wording and misreads an ESEM Spearman or pairwise fit.
- 2026-10-01: plan gate absorbed the DESIGN "exact" row (M95 review R9) over keeping it separate, because it states the same fact in the same file.
- 2026-10-01: T1 done. Three tests (Spearman, tied Spearman, listwise) on PCA, EFA, and ESEM over all 6 level pairs agree to at most 2.4e-15. Red once with planted defects: raw data or unreduced rows failed all 18 pair checks per test, gaps 0.009 to 0.016. Suite 812 tests, 0 failed.
- 2026-10-01: T2 done. The first Known limitations entry lists the four uncovered settings and names the three T1 tests by title. §5.4 closes with a pointer to it. §5.3 says rank scores reproduce a Spearman `R`, and narrows "under missing data" to pairwise or FIML.
- 2026-10-01: T3 finding. ESEM ULSMV on continuous items errors at level 1 under `missing = "pairwise"` (lavaan 0.7.2, `available.cases`), so AC3's ULSMV clause cannot be derived. Logged as a `[high]` hotfix candidate row. AC3 amendment pending at a mini gate.
- 2026-10-01: T4 ledger, AC4 search on DESIGN.md (16 lines). Fixed to name `x$r`: 119 (IP1), 186 (§4), 193 (pca row), 221 (§5 intro, was "the matrix the fit uses", marked corrected), 458 (§9 engine), 459 (§9 rotation, "W'RW identity is exact"), 466 (§9 redundancy_phi, "true correlation" and "algebra-exactness"). Already named: 233 (§5.1, "implied by `R`"). Quoted superseded wording, excluded: 150, 223. Not an edge claim: 65, 146, 529, 547, 580, 664. D-039 appended for the IP1 wording.
- re-audit: AC3 (full) — 4 findings on the ULSMV-dropped wording. A: under `cor = "spearman"` with `missing = "fiml"`, `x$r` is lavaan's Pearson FIML matrix, so "Spearman" overstates the differing set. B: under FIML lavaan fits the likelihood, not a matrix. C: the help must not imply ULSMV is covered. D: the list reads as complete but names no enumerating procedure. A and D taken into the gate wording, B and C not needed.
- 2026-10-01: ESEM grid sweep (`cor` 3 × `missing` 3 × `estimator` 4 × complete or missing data, sim16, k_max 2). lavaan's matrix differs from `x$r` only under Spearman with pairwise or listwise (0.03 to 0.046), and Pearson pairwise on missing data (ML and MLR fit the complete rows, 0.029, and WLSMV fits pairwise covariances, 0.0051). Spearman with FIML stores lavaan's Pearson FIML matrix as `x$r`. Logged as a second `[high]` hotfix candidate row.
- re-audit: AC3 (full) — 3 findings on the revised wording (Spearman pairwise or listwise, Pearson pairwise on missing data, WLSMV only). F1: "as a sweep finds" rests completeness on an uncommitted session sweep, so ground it in the code paths with the sweep as confirmation. F2: D-039's Context names Spearman without the `missing` qualifier. F3 (adjacent): the fit-time Spearman warning says "Pearson-ML" though WLSMV is not ML. This is the second AC3 re-audit, so further AC3 wording goes to the user.
- 2026-10-01: amendment, user-selected at the mini gate. AC3 drops the ULSMV clause, limits Spearman to `missing = "pairwise"` or `"listwise"`, and grounds the list in the two code sites with the sweep as confirmation (narrowing, F1 applied). D-039's Context and the DESIGN §5 correction note gained the same qualifier, D-039 edited in place on the unmerged branch at the user's choice.
- 2026-10-01: T3 done. The `ackwards()` opening paragraph and the `compute_edges()` help name `x$r` and the two ESEM settings, each clause measured in the grid sweep. NEWS entry added. `document()` rewrote both Rd files. `check-prose.R` clean on both files and NEWS.md, and `--code-unchanged master` reports no code line changed.
- 2026-10-01: T5 done. `Rscript tools/dod-gate.R` exit 0 on 34eae8d: prose clean, check 0/0/0, coverage 100%, style and lint clean, pkgdown index complete.
- claim audit: 20 claims read, 2 corrected — NEWS.md, R/ackwards.R, R/compute_edges.R, man/ackwards.Rd, man/compute_edges.Rd, tests/testthat/test-compute_edges.R
- 2026-10-01: the two claim-audit fixes (NEWS credits the ranking to the tests, and the help says why WLSMV's pairwise covariances differ from `x$r`) re-read once by the same reader, both hold. `document()` rerun, `check-prose.R` and `--code-unchanged master` clean. Status set to review.
- 2026-10-01: review in progress. AC1 to AC4 verified with fresh evidence and ticked. The DoD gate (AC5) and three reviewers are running.

## Decisions

## Review

Evidence gathered 2026-10-01 on the branch head 6be9482. The branch already contained `origin/master`, so no merge was needed.

- AC1: `devtools::test(filter = "compute_edges")` ran 23 tests with 0 failed and 0 errors. The three new tests passed: Spearman (24 expectations), tied Spearman (25), and listwise (25). Each test fits `pca`, `efa`, and `esem` at `k_max = 4` and checks all 6 level pairs from `compute_edges(x$levels, x$r, pairs = "all")` at the 1e-10 bound. The scores route gets `apply(<data>, 2, rank)` for both Spearman cases and the complete rows for listwise. The tied case asserts duplicate values in column 1, and the listwise case asserts that rows were dropped.
- AC2: I read `cairn/DESIGN.md` lines 620 to 636. The first Known limitations entry lists four settings. These are polychoric `R` on any engine, FIML on PCA or EFA and on ESEM, pairwise deletion on data with missing values, and matrix input. It names the file and the three test titles, which match the `test_that()` titles exactly. It says that a setting it neither lists nor names is untested. Its example is Spearman plus listwise on missing data. §5.4 (line 335) ends with a pointer to that entry. The §5.3 paragraph "Where the `scores` route runs" is at lines 317 to 327. It says that raw-item scores do not reproduce a Spearman `R`, but column-rank scores do.
- AC3: The opening paragraph of the `ackwards()` roxygen says that each edge is exact for the matrix stored as `x$r`. It names Spearman with pairwise or listwise, and Pearson pairwise on data with missing values. In the second setting it says that ML and MLR fit the complete rows and WLSMV fits pairwise covariances. The `compute_edges()` roxygen names the same settings. I read the two code sites. `R/ackwards.R:806` builds `x$r` with `stats::cor(method = cor, use = "pairwise.complete.obs")`, and `R/engine_esem.R:434` maps pairwise to `available.cases` for WLSMV and ULSMV and to listwise otherwise. A fresh sweep ran on lavaan 0.7.2 over `cor` 2 × `missing` 3 × `estimator` 4 × complete or missing `sim16`, at `k_max = 2`. Polychoric was left out, because lavaan computes that `R` itself. The sweep compared lavaan's sample (or FIML h1) correlations with `x$r`. Every setting agreed to 2.2e-16 or better except the named ones. Spearman pairwise or listwise differed by 0.03 to 0.046. Pearson pairwise on missing data differed by 0.029 under ML and MLR and by 0.0051 under WLSMV. ULSMV with pairwise failed at level 1, as the hotfix candidate row records. Spearman with FIML gave 0, because `x$r` there is lavaan's Pearson FIML matrix, which is the second hotfix candidate row.
- AC4: The search returned 16 lines, the same set as the T4 ledger. I read each line in context. Lines 119 (IP1), 186, 193, 221, 458, 459, and 466 name `x$r`. They do so in every sentence that calls the algebra, an edge, or a score correlation exact or true. Line 233 says the correlations "implied by `R`" are exact, so it names the matrix. Lines 150 and 223 hold quoted superseded wording. Lines 65, 146, 529, 547, 580, and 664 use "exact" about other things, such as Forbes reproduction, `augment()`, layout, and seeds. D-039 in `cairn/DECISIONS.md` records the IP1 wording change and cites D-031's procedure.
- AC5: `Rscript tools/dod-gate.R` exited 0 on the code of 6be9482. The only later commit, 1abc3d8, changes this tracking file and nothing else. The gate reported vignette freshness, ledger anchors, CI path filters, and prose clean. It reported check 0 errors, 0 warnings, 0 notes, coverage 100%, style and lint clean, and a complete pkgdown index.
- Consistency gate: `cairn_validate.py` exited 0, with 16 work-log format warnings that are all in the M84 file. `cairn_impact.py --changed` listed the IP1, IP2, and IP6 citing lines. Only IP1's text changed, and its rule is unchanged (D-039), so no citing line needs a change. `devtools::document()` left no diff. README is unchanged on the branch, NEWS.md has an entry with no milestone number, and the branch adds no top-level file. The pkgdown index and `check()` results are under AC5.

Independent review: three fresh reviewers (Opus diff, Sonnet blame history, Sonnet prior reviews). The GitHub thread probe returned no comments. No finding shows a criterion failing. The findings are merged across lenses below, with the proposed disposition of each.

- R1 (Opus 1). The `@param missing` listwise bullet in `ackwards()` says the fit and the edges "are all consistent". The new opening paragraph on the same page says ESEM Spearman listwise fits another matrix. Proposed: fix now.
- R2 (prior reviews 1). The `rotation = "varimax"` bullet in `ackwards()` says the algebra is exact for "the correlation matrix the fit uses". Proposed: fix now, name `x$r`.
- R3 (blame 1 and 2). DESIGN still says "the fit uses" or "the fit's `R`" at the §5.1 legend (line 237), §5.3 (322), §11 (484), and Known limitations (639). Line 233 calls the edges exact for that `R`. Proposed: fix now, name `x$r`.
- R4 (blame 3, Opus 8). The sweep edited text inside the M76 and M43 correction notes (lines 459 and 466) with no marker. Proposed: fix now, add a "corrected M097" marker to each.
- R5 (Opus 2). DESIGN line 466 says `|r|` on `x$r` is "the true correlation" of the components. That overclaims for a polychoric, pairwise, or FIML `x$r`. Proposed: fix now, say the correlation that `x$r` implies.
- R6 (Opus 4). On ESEM the new tests check edges against `x$r`, not against the matrix lavaan fits. The Known limitations entry does not say so. Proposed: fix now, one sentence.
- R7 (Opus 3). The Known limitations ESEM pairwise entry covers ML and MLR only. WLSMV on continuous items fits pairwise covariances. Proposed: fix now, one sentence.
- R8 (Opus 10, prior reviews 3, blame 6). The `compute_edges()` help ends with "passes a pooled `R` instead", after a stranded line break. Proposed: fix now, say "in place of `x$r`".
- R9 (Opus 5, blame 5). The `@param missing` pairwise bullet calls WLSMV and ULSMV "(ordinal)" and says `$meta` documents the inconsistency. This text predates the branch. Proposed: follow-up, absorbed into the ULSMV hotfix candidate row, which edits this bullet.
- R10 (prior reviews 2, blame 2). README, the vignettes, the manuscript, `?prune`, and `?tidy` keep "the matrix the fit uses". Proposed: reject, because the plan's Scope Out kept them on purpose, and both ESEM settings warn at fit time.
- R11 (blame 1). §5.3 says an edge is "never a sample-realized score correlation", but rank scores reproduce a Spearman edge. Proposed: reject, because rank scores are not the factor scores a user computes. R3 renames the matrix in that sentence.
- R12 (blame 4). The rows at DESIGN lines 186, 193, 458, and 466 gained `x$r` with no correction marker. Proposed: reject, because naming the matrix narrows a true claim, and D-039 records the sweep.
- R13 (Opus 6). "Pearson covariances of the raw data" does not say which rows lavaan uses. Proposed: reject, because the sentence contrasts ranked and unranked data, and the next sentence covers the rows.
- R14 (Opus 9). §5.3 does not mention matrix input. Proposed: reject, because matrix input has no scores, and Known limitations lists it.
- R15 (Opus 10). AC5 unticked and the milestone file uncommitted. Reject: stale, both done by 003e432.
- R16 (Opus 7). ULSMV pairwise fails, and ESEM Spearman with FIML stores a Pearson matrix. Noted: both are `[high]` hotfix candidate rows already.
