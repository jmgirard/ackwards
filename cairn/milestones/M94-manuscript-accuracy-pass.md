# M94: Manuscript accuracy pass for four older claims

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M93
- **Driving RR:** —
- **Principles touched:** IP1, IP2, GP3
- **Resolves:** —
- **Surface tier:** user-facing — the manuscript is written for journal readers, and the two vignettes ship with the package
- **Branch/PR:** m094-manuscript-accuracy-pass

## Goal

The manuscript and two vignettes describe the edge algebra, Waller's (2007) contribution, and the AI-use disclosure accurately.

## Scope

**In:**

- Qualify each claim that the edge algebra is exact, in the manuscript abstract and methods passage and in the intro and engines vignettes.
- Replace the manuscript's fallback-to-scores sentence: every edge comes from the algebra, and the scores cross-check covers only complete-data linear engines.
- Limit every Waller (2007) credit to principal components, and present the weight-matrix form as the covariance of linear composites, with no citation and no novelty claim.
- Make both disclosure passages give AI-assisted development from June 2026 to submission, and name Opus, Sonnet, and Fable.

**Out:**

- `R/ackwards.R:13` and the engines vignette's "holds exactly for components" credit, which are correct.
- NEWS.md entries, which are history.
- DESIGN.md's scores-route text → M93.
- Waller's §4 point that estimated factor scores only approximate the factor correlations. This milestone corrects attributions, not the package's choice of score-composite edges (D-007).

## Acceptance criteria

- [ ] AC1: Four paragraphs call the edge algebra exact: the abstract (`manuscript/manuscript.qmd` ~37), the methods passage (~162), the intro vignette's "Between-level edges" paragraph (`vignettes/ackwards-intro.Rmd.orig` ~277), and the engines vignette's PCA paragraph (`vignettes/ackwards-engines.Rmd.orig` ~84-86). Each, read in full, either drops the exactness claim or says two things. The edges are exact for the correlation matrix supplied, and scores computed from the observed items do not in general reproduce them when that matrix is polychoric, Spearman, pairwise, or FIML. No other paragraph of the three files that contains a match of `grep -n -i 'exact'` says that the edges equal correlations of scores computed from the observed items.
- [ ] AC2: No paragraph of the three files that contains a match of `grep -n -i -E 'materiali|falls? back|nonlinear'` says that the package computes edges from materialized scores. The manuscript's methods passage (~158-162) and conclusion passage (~382-384) each say that every edge the package reports comes from the closed-form algebra. Where either passage mentions the cross-check against materialized scores, it limits the check to complete-data linear engines and excludes the polychoric and FIML paths, as the first `cairn/DESIGN.md` "Known limitations" entry does.
- [ ] AC3: In the three files, every paragraph that contains a match of `grep -n -i 'waller'`, read in full with reference-list entries aside, credits Waller (2007) with no more than the closed form for principal components and its oblique form for rotated components (waller2007 (pp. 748-749), Eqs. 9-14 and §3). Where a paragraph states the weight-matrix form `W'RW` for any linear scoring weights, it presents that form as following from the covariance of linear composites, attributes it to no source, and claims no novelty for the package.
- [ ] AC4: The front-matter disclosure and the "Use of generative AI" section both describe the period of AI-assisted development as running from June 2026 to submission. Both name the Claude model families Opus, Sonnet, and Fable.
- [ ] AC5: `quarto render manuscript.qmd`, run in `manuscript/`, exits 0 and writes `manuscript.pdf` and `manuscript.docx`. The count of em dashes in `manuscript/manuscript.qmd` (`grep -o '—' manuscript/manuscript.qmd | wc -l`) is no higher than on `master`.
- [ ] AC6: `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` exits 0 on the branch head.

## Coverage

- AC1 → T1, T2, T4, T5
- AC2 → T1, T2, T5
- AC3 → T1, T2, T4, T5
- AC4 → T3
- AC5 → T6
- AC6 → T4, T6

## Tasks

- [x] T1: Rewrite the manuscript methods passage (`manuscript/manuscript.qmd` ~150-170). Limit the Waller credit (~154-156) to components and present `W'RW` as the covariance of linear composites. Replace the fallback sentence (~160-161) and qualify "exact" (~162).
- [x] T2: Fix the Waller credit, exactness, and cross-check scope in the abstract (~37), the introduction (~103), the results text (~295), and the conclusion (~382-384).
- [ ] T3: Rewrite the front-matter disclosure (~19) and the "Use of generative AI" section (~465) to give AI-assisted development from June 2026 to submission and to name Opus, Sonnet, and Fable.
- [ ] T4: Fix the intro vignette (`vignettes/ackwards-intro.Rmd.orig` ~276-278) and the engines vignette (`vignettes/ackwards-engines.Rmd.orig` ~84-86). Re-run `Rscript vignettes/precompute.R`, then revert noise in untouched vignettes and noise lines in the two edited ones (LESSONS M61, M75, M87).
- [ ] T5: Read each rewritten paragraph in full in its final position (LESSONS M92). Check every new claim against the waller2007 page images and the R source (LESSONS M67, M86). Run the AC1-AC3 searches on the branch head.
- [ ] T6: Run `quarto render manuscript.qmd` in `manuscript/`, count em dashes against `master`, and run `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R`.

## Work log

- 2026-10-01: created by /milestone-plan, promoted from the candidate row "Manuscript accuracy pass, four older claims from the M92 review" (added 2026-09-30).
- 2026-10-01: criteria audit (full mode, user-facing tier) returned 11 findings. Fixed before the gate: the exactness caveat was too narrow and unprovable, AC1 passed if the word was deleted, AC2 missed the conclusion and engines-vignette lines, AC2 overstated where scores are used, AC3 checked lines instead of paragraphs, and the gate script's code guard was run apart. Re-reading rewritten text in place moved to T5. The gate settled the other four.
- 2026-10-01: the owner asked at the gate whether the Waller attribution is wrong. The waller2007 page images (pp. 745-752) confirm a principal-components-only result in transformation-matrix form (Eq. 14, §3), and §4 sets the factor model aside.
- 2026-10-01: plan gate chose plain wording with no citation for the general `W'RW` form over adding a textbook citation, because the identity is elementary and a citation needs a verified reference note; falsified by an editor or reviewer asking for a source.
- 2026-10-01: plan gate chose "June 2026 to submission" over "June to October 2026" so that the period never goes stale; falsified by a journal that requires a closed date range.
- 2026-10-01: plan gate chose to include the two vignette passages over a candidate row, so the same fault does not keep shipping in package docs; falsified by a regeneration that cannot stay scoped to the two passages.
- 2026-10-01: plan chose to depend on M93 over planning independently, because the new methods wording contradicts IP1 and IP2 until M93 lands; falsified by M93 being dropped.
- 2026-10-01: implement started on branch m094-manuscript-accuracy-pass. The question gate was skipped because the plan gate settled every open choice.
- 2026-10-01: T1 done. The methods passage credits Waller with the components result only. It derives `W'RW` from the covariance of linear composites and says every edge comes from the algebra. It limits the scores check to complete Pearson data on the linear engines and qualifies "exact". The Waller wording was read against the PDF text (abstract, pp. 748-749).
- 2026-10-01: T2 done. The abstract and results text no longer cite Waller or call the algebra exact. The introduction limits his credit to principal components. The conclusion says every edge comes from the closed form and limits the scores check as T1 does.

## Decisions

## Review
