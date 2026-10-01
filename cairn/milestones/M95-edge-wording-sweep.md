# M95: Edge wording sweep across shipped docs and DESIGN §5

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, IP2, GP3
- **Resolves:** —
- **Surface tier:** user-facing — DESCRIPTION, help pages, vignettes, README, and NEWS ship in the package
- **Branch/PR:** m095-edge-wording-sweep

## Goal

Each edge description that the search below finds names the matrix the edge is exact for and credits Waller (2007) only with his own results.

## Scope

**In:** Wording fixes at the sites the M94 review listed (R5, R6, R7, P6, P7, R9) and at every other site the search below finds. DESIGN §5 details that predate M93 (the §5.3 pseudocode, the §5.1 scoring and Waller text, the §5.2 comments). The `compute_edges()` roxygen header. A NEWS entry. Re-precomputed vignettes and a rebuilt README where their sources change.

The swept files D are these:

- `DESCRIPTION` and `README.Rmd`
- the roxygen text (`#'` lines) of `R/*.R`
- `vignettes/*.Rmd.orig` and `vignettes/ackwards-interpret.Rmd`
- `manuscript/manuscript.qmd`
- the development section of `NEWS.md` (above `# ackwards 0.2.0`)

A shipped edge is an edge that an exported function returns. The search Q works paragraph by paragraph, so a phrase split by a line break still matches:

```
perl -00 -ne 'print "$ARGV\n$_\n" if /score[\s#'"'"']+correlation|correlations?[\s#'"'"']+(between|among)[\s#'"'"']+(the[\s#'"'"']+)?(factor[\s#'"'"']+|component[\s#'"'"']+)?scores|factor-score[\s#'"'"']+correlation|algebraic[\s#'"'"']+equivalent|between-level[\s#'"'"']+correlation|closed[\s#'"'"'-]+form|exact/i' <files>
```

**Out:** The console note in `R/print.R` ("Cross-level edges are descriptive score correlations") stays as is. GP3 asks `print()` to keep that wording, and the plan gate kept it. Released NEWS sections are a record of past releases and stay unedited. The cross-check-paths Known limitations row stays a `[low]` candidate. No R code behavior changes.

## Acceptance criteria

- [ ] AC1: Take the lines that `grep -n -i -E 'materiali[sz]'` returns over D at the review commit. No such line describes materialized scores as a route by which a shipped edge is computed. (DESCRIPTION:15-17 now says edges come "with exact linear algebra (Waller, 2007) or from materialized scores". D-038 makes the scores route internal only.)
- [ ] AC2: Take each line of D and of `cairn/DESIGN.md` that `grep -n -i waller` returns at the review commit, less reference-list entries. Read it with its sentence and paragraph. Each one credits Waller (2007) only with results that `cairn/references/waller2007.md` records he wrote. Two examples are the principal-components result in transformation-matrix form (waller2007, Eq. 14, p. 749) and its oblique form (§3, p. 749). None credits him with the `W′RW` form for general linear weights, or with the general identity.
- [ ] AC3: Take the paragraphs that search Q returns over D and over `cairn/DESIGN.md` §5 at the review commit. For `R/*.R`, only roxygen text counts. Some of these paragraphs describe a between-level edge. They say what an edge is, call an edge exact, or equate it with a named quantity. Each such description says that the edge is exact for (or implied by) the correlation matrix the fit uses. As an alternative, it points to a passage in the same document that says so. A paragraph that uses a matched phrase for another purpose needs no qualifier. An example is the redundancy prose that compares a pair's `|r|` to a threshold.
- [ ] AC4: In the engines vignette, the PCA paragraph credits Waller (2007) with the components result in his transformation-matrix terms. Every prose mention of the algebra in that vignette writes `W′RW` with the prime that the intro vignette uses. The PCA paragraph uses "exact" only in the sense "exact for the correlation matrix the fit uses".
- [ ] AC5: The `compute_edges()` help page is `man/compute_edges.Rd`, built from its roxygen. It says the algebra is exact for the fit's R under the linear scoring of PCA, EFA, and ESEM. Its account of the conditions that send a pair to the scores branch matches the function's branch conditions at the review commit. It adds that no shipped caller reaches that branch.
- [ ] AC6: `cairn/DESIGN.md` §5 matches `R/compute_edges.R` at the review commit in five respects. The §5.3 pseudocode signature lists exactly the arguments and defaults of `formals(compute_edges)`. Its body shows the abort that `edge_method = "algebra"` raises on a pair that fails the algebra conditions. Its body has no sign-alignment step. §5.1 names ten Berge as the EFA and ESEM scoring, with regression as the fallback. The §5.2 `NULL if !linear` comments are removed or marked as reached by no engine.
- [ ] AC7: The shipped generated files match their sources. Each regenerated vignette `.Rmd` carries the md5 stamp of its edited `.Rmd.orig`. It differs from master only in the edited prose and that stamp. `man/` and `README.md` match what `devtools::document()` and `devtools::build_readme()` produce from the edited sources. `manuscript/manuscript.qmd` renders without error.

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2, T3, T4, T5
- AC3 → T1, T2, T3, T4, T5
- AC4 → T3
- AC5 → T2
- AC6 → T5
- AC7 → T2, T3, T4, T7

## Tasks

- [x] T1: Run the AC1 and AC2 greps and search Q over D and DESIGN §5. Record each hit as an edge description or another use, with a count per file, in one work-log line.
- [x] T2: Edit `DESCRIPTION` (Description field, lines 14-17), `README.Rmd` (near lines 108-112), and the roxygen hits: `R/ackwards.R` lines 6-8, 13, and 21, the `R/compute_edges.R` header (lines 1-16), and any other edge description T1 found. Run `devtools::document()` and `devtools::build_readme()`.
- [x] T3: Edit the vignette sources. Known sites are intro 55 and 140-142, engines 83-87 and 460 and 540, ordinal 384, and girard 273, plus T1's other hits. Run `Rscript vignettes/precompute.R`, revert run noise line by line (M75, M87), and diff each `.Rmd` against master with the stamp line removed (M94).
- [ ] T4: Edit the manuscript sites (146-148, 304, 407, and T1's other hits). Render the manuscript.
- [ ] T5: Edit DESIGN §5: the §5.1 Waller and scoring text (line 237-240), the §5.2 comments (253-254), and the §5.3 pseudocode signature and body (260-283).
- [ ] T6: Re-read each rewritten claim against `cairn/references/waller2007.md` and the R source (M64, M67, M86). Run `Rscript tools/check-prose.R` on every edited doc file (M85). Add a NEWS entry in the development section.
- [ ] T7: Run `Rscript tools/dod-gate.R` (it includes the vignette-freshness check).

## Work log

- 2026-10-01: created by /milestone-plan. It absorbs two candidate rows. One is the M94-review edge-wording sweep (R5, R6, R7, P6, P7). The other is the stale DESIGN §5 details (M93 review F5, F8, F9, and M94 review R9).
- 2026-10-01: criteria audit (full mode, fresh Opus reader) returned 12 findings and no IP or D-entry conflict. All were fixed before the gate. The fixes are a paragraph-mode search with a narrowed Goal, DESIGN §5 in AC3, and a decidable "exact" test. Also sentence context for AC2, a whitelist from the Waller note, and the NEWS dev section in D. Also a defined "shipped edge", a named sense for AC4, AC6's count and fallback, and deliverable-only wording for AC7.
- 2026-10-01: plan gate chose one milestone for shipped docs plus DESIGN §5 over a shipped-docs-only milestone. The §5 row asks to ride the next §5 edit, and it shares the Waller slip. Falsified by §5 edits that push the plan-owned body past its cap or that need a design decision.
- 2026-10-01: plan gate chose to leave the print note unchanged over adding a qualifier. GP3 asks print() to keep "score correlations", and the change adds code and snapshot churn to a docs-only milestone. Falsified by a user who reads the print note as a sample-realized correlation.
- 2026-10-01: plan gate kept the cross-check-paths row separate over folding in its docs half. That row poses a document-or-test choice of its own. Falsified by a T5 read of §5.4 that shows the omission makes a §5 sentence false.
- 2026-10-01: plan chose a paragraph-mode `perl -00` search over the per-line grep. The audit found five edge descriptions split across line breaks that the grep missed. Falsified by a missed site whose phrase is in the pattern but spans a paragraph break.
- 2026-10-01: T1 done. Search Q found 90 paragraphs over D and DESIGN §5, 29 of them edge descriptions (E). Per file, hits/E: DESCRIPTION 1/1, README.Rmd 1/1, manuscript 14/7, NEWS dev 4/1, DESIGN §5 3/3, vignettes engines 6/4, forbes 10/1, forbes2023 4/0, girard 3/1, interpret 1/0, intro 7/3, ordinal 1/1, suggest-k 2/0, visualization 4/0. Roxygen: ackwards 6/2, compute_edges 2/2, prune 3/1, tidy 1/1, and 0 E among augment 3, autoplot 2, boot_edges 2, comparability 1, data 4, label_template 1, layout 2, predict 1, suggest_k 1. AC1 found 16 lines in D, one a shipped-edge route (DESCRIPTION). AC2 found 14 Waller lines, 4 of them reference entries. The ledger is under Decisions.
- 2026-10-01: T2 done. DESCRIPTION drops the materialized-scores route and credits Waller with the components result only. README.Rmd, `ackwards()`, `prune()`, and `tidy()` roxygen name the fit's correlation matrix. The `compute_edges()` header states the scores-branch conditions from the code and says no shipped caller reaches it. `document()` and `build_readme()` changed only those paragraphs. check-prose and its code-unchanged guard pass.
- 2026-10-01: T3 done. Edited engines, intro, ordinal, girard, and forbes `.Rmd.orig`. The engines PCA paragraph credits Waller with the transformation-matrix result, and its prose writes `W′RW`. Re-ran precompute.R. Reverted suggest-k, visualization, three PNGs, 8 timing lines, and two gt table ids (run noise, M75, M87). Each edited `.Rmd` differs from master, stamp line removed, only in the edited prose. Freshness check and code-unchanged guard pass.

## Decisions

- T1 ledger for search Q, by master line number. A hit not listed here is another use. The E paragraphs to fix are DESCRIPTION 16-19 and README.Rmd 107. In the vignettes they are engines 82, 277, 460, and 539 (prime only), girard 267, intro 53 and 137, ordinal 372, and forbes 531 (it equates the PCA `|r|` with the components' correlation). In the manuscript they are 135, 308, and 404. In DESIGN §5 they are the lead paragraph and §5.1. In roxygen they are ackwards.R 3 and 10, compute_edges.R 1-16, the prune.R `redundancy_phi` entry, and the tidy.R `what = "edges"` entry. The E paragraphs already compliant are intro 276, the manuscript abstract, manuscript 151, 304, and 382, NEWS dev 3, and the DESIGN §5.3 note. The other uses are the redundancy and artifact prose, captions that only encode `|r|`, fidelity claims, the forbes 20 and forbes2023 20 accounts of Goldberg's and Forbes's methods, and each "exactly" that means "precisely".

## Review
