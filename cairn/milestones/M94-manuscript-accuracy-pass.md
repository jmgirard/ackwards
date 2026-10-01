# M94: Manuscript accuracy pass for four older claims

- **Status:** review
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

- [x] AC1: Four paragraphs call the edge algebra exact: the abstract (`manuscript/manuscript.qmd` ~37), the methods passage (~162), the intro vignette's "Between-level edges" paragraph (`vignettes/ackwards-intro.Rmd.orig` ~277), and the engines vignette's PCA paragraph (`vignettes/ackwards-engines.Rmd.orig` ~84-86). Each, read in full, either drops the exactness claim or says two things. The edges are exact for the correlation matrix supplied, and scores computed from the observed items do not in general reproduce them when that matrix is polychoric, Spearman, pairwise, or FIML. No other paragraph of the three files that contains a match of `grep -n -i 'exact'` says that the edges equal correlations of scores computed from the observed items.
- [x] AC2: No paragraph of the three files that contains a match of `grep -n -i -E 'materiali|falls? back|nonlinear'` says that the package computes edges from materialized scores. The manuscript's methods passage (~158-162) and conclusion passage (~382-384) each say that every edge the package reports comes from the closed-form algebra. Where either passage mentions the cross-check against materialized scores, it limits the check to complete-data linear engines and excludes the polychoric and FIML paths, as the first `cairn/DESIGN.md` "Known limitations" entry does.
- [x] AC3: In the three files, every paragraph that contains a match of `grep -n -i 'waller'`, read in full with reference-list entries aside, credits Waller (2007) with no more than the closed form for principal components and its oblique form for rotated components (waller2007 (pp. 748-749), Eqs. 9-14 and §3). Where a paragraph states the weight-matrix form `W'RW` for any linear scoring weights, it presents that form as following from the covariance of linear composites, attributes it to no source, and claims no novelty for the package.
- [x] AC4: The front-matter disclosure and the "Use of generative AI" section both describe the period of AI-assisted development as running from June 2026 to submission. Both name the Claude model families Opus, Sonnet, and Fable.
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
- [x] T3: Rewrite the front-matter disclosure (~19) and the "Use of generative AI" section (~465) to give AI-assisted development from June 2026 to submission and to name Opus, Sonnet, and Fable.
- [x] T4: Fix the intro vignette (`vignettes/ackwards-intro.Rmd.orig` ~276-278) and the engines vignette (`vignettes/ackwards-engines.Rmd.orig` ~84-86). Re-run `Rscript vignettes/precompute.R`, then revert noise in untouched vignettes and noise lines in the two edited ones (LESSONS M61, M75, M87).
- [x] T5: Read each rewritten paragraph in full in its final position (LESSONS M92). Check every new claim against the waller2007 page images and the R source (LESSONS M67, M86). Run the AC1-AC3 searches on the branch head.
- [x] T6: Run `quarto render manuscript.qmd` in `manuscript/`, count em dashes against `master`, and run `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R`.

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
- 2026-10-01: T3 done. Both disclosure passages give June 2026 to submission and name Opus, Sonnet, and Fable.
- 2026-10-01: T4 done. Both vignette passages qualify "exact". The intro passage limits Waller to components and derives `W'RW` from linear composites. precompute.R ran clean. Untouched vignettes and assets were reverted, and the two edited `.Rmd` files were rebuilt from master plus the prose and stamp lines. Noise counts are 0, and the freshness check passes.
- 2026-10-01: T5 done. Every rewritten paragraph was re-read in place. The waller2007 page images (pp. 748-749) confirm "component transformation matrix" and the oblique form for rotated components. The R source confirms that all three engines set `linear = TRUE` and that `cor` accepts spearman. Every scores-agreement test (PCA, EFA, ESEM, oblique) uses complete Pearson data. The AC1-AC3 searches leave no failing paragraph. Intro vignette line 142 ("identity is exact for any fixed linear scoring") stays, because it does not say edges equal observed-score correlations.
- 2026-10-01: T6 done. `quarto render manuscript.qmd` exited 0 and rewrote both outputs. The em-dash count is 5 on the branch and on master. The first gate run failed on its prose check, because a 39-word intro vignette sentence was over the 30-word limit. The sentence was split and the stamp was updated, and the second run of `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` passed (check 0/0/0, coverage 100%).
- claim audit: 17 claims read, 3 corrected — manuscript/manuscript.qmd, vignettes/ackwards-intro.Rmd(.orig), vignettes/ackwards-engines.Rmd(.orig)
- 2026-10-01: claim-audit fixes. "The linear engines" became "all three engines" in the methods and conclusion. The intro covariance became `W_a' R W_b`. The engines vignette dropped "therefore" before its exactness claim. The same reader re-read all three, and they hold. The kept engines sentence that credits Waller with `W'RW` is out of scope by plan, and the reader marks it unclear (Waller uses transformation matrices). It is raised for review.
- 2026-10-01: after the fixes, `quarto render` exited 0 and rewrote both outputs. The em-dash count is 5 on the branch and on master. `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` passed (check 0/0/0, coverage 100%). Status set to review.

## Decisions

## Review

Evidence gathered 2026-10-01 at head 7c9d79a. Master has not moved since the branch was cut.

- AC1 evidence: the full paragraphs that match `exact` were printed from all three files. The abstract no longer contains "exact". The methods paragraph (ms:152) says the algebra is exact for the matrix given, and that observed-item scores do not in general reproduce the edges under polychoric, Spearman, pairwise, or FIML matrices. The intro paragraph (intro.orig:276) and the engines PCA paragraph (engines.orig:82) say the same. The other matches (ms:404, intro.orig:91, 133, 137, 333, 360, engines.orig:22, 125, 538) do not say that edges equal correlations of observed-item scores. ms:404 calls edges "correlations among factor scores, not parameters of a fitted model", which contrasts them with model parameters and makes no claim about observed-item scores. Pass.
- AC2 evidence: the paragraphs that match `materiali|falls? back|nonlinear` are ms:152, ms:384, intro.orig:276, and engines.orig:82. None says the package computes edges from materialized scores. ms:152's "direct route materializes" describes the general route, cited to Grice. The methods passage (ms:152) and the conclusion (ms:384) each say every edge the package reports comes from the closed-form algebra. Both limit the scores check to all three engines (all linear) with complete data and Pearson correlations, and both exclude the polychoric and FIML paths. Pass.
- AC3 evidence: the paragraphs that match `waller`, reference entries aside (intro.orig:444, engines.orig:658), are ms:102, ms:152, intro.orig:276, and engines.orig:82. ms:102 and intro.orig:276 credit a closed form for principal components. ms:152 credits the components result via transformation matrices and the oblique form for rotated components, which the page images confirm (pp. 748-749). ms:152 and intro.orig:276 state `W'RW` for any linear weights as following from the covariance of linear composites, with no source and no novelty claim. engines.orig:82 keeps "Waller (2007) showed that the between-level algebra (`W'RW`) holds exactly for components". That credit is limited to components and states no general form, and Scope Out lists it as correct. Pass.
- AC4 evidence: the front matter (ms:19-21) reads "From June 2026 to submission, the author used large language models (Anthropic's Claude Opus, Sonnet, and Fable models ...)". The "Use of generative AI" section (ms:479-482) reads "AI-assisted development ... ran from June 2026 to submission ... with models from its Opus, Sonnet, and Fable families." Pass.
