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

- [x] AC1: Four paragraphs call the edge algebra exact: the abstract (`manuscript/manuscript.qmd` ~37), the methods passage (~162), the intro vignette's "Between-level edges" paragraph (`vignettes/ackwards-intro.Rmd.orig` ~277), and the engines vignette's PCA paragraph (`vignettes/ackwards-engines.Rmd.orig` ~84-86). Each, read in full, either drops the exactness claim or says two things. The edges are exact for the correlation matrix supplied, and scores computed from the observed items do not in general reproduce them when that matrix is polychoric, Spearman, pairwise, or FIML. No other paragraph of the three files that contains a match of `grep -n -i 'exact'` says that the edges equal correlations of scores computed from the observed items.
- [x] AC2: No paragraph of the three files that contains a match of `grep -n -i -E 'materiali|falls? back|nonlinear'` says that the package computes edges from materialized scores. The manuscript's methods passage (~158-162) and conclusion passage (~382-384) each say that every edge the package reports comes from the closed-form algebra. Where either passage mentions the cross-check against materialized scores, it limits the check to complete-data linear engines and excludes the polychoric and FIML paths, as the first `cairn/DESIGN.md` "Known limitations" entry does.
- [x] AC3: In the three files, every paragraph that contains a match of `grep -n -i 'waller'`, read in full with reference-list entries aside, credits Waller (2007) with no more than the closed form for principal components and its oblique form for rotated components (waller2007 (pp. 748-749), Eqs. 9-14 and §3). Where a paragraph states the weight-matrix form `W'RW` for any linear scoring weights, it presents that form as following from the covariance of linear composites, attributes it to no source, and claims no novelty for the package.
- [x] AC4: The front-matter disclosure and the "Use of generative AI" section both describe the period of AI-assisted development as running from June 2026 to submission. Both name the Claude model families Opus, Sonnet, and Fable.
- [x] AC5: `quarto render manuscript.qmd`, run in `manuscript/`, exits 0 and writes `manuscript.pdf` and `manuscript.docx`. The count of em dashes in `manuscript/manuscript.qmd` (`grep -o '—' manuscript/manuscript.qmd | wc -l`) is no higher than on `master`.
- [x] AC6: `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` exits 0 on the branch head.

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
- [ ] T7: Add a NEWS.md entry for the intro and engines vignette corrections (review gate return 1).
- [x] T8: Fix review findings R1-R4, R8, R10, R11 in the manuscript and intro vignette, then re-run the T6 render, count, and gate.

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
- 2026-10-01: review consistency gate failed. NEWS.md has no entry for the user-visible vignette corrections (profile consistency-gate, NEWS line). Defect return 1. Status set back to in-progress. AC1-AC6 evidence and 16 ranked reviewer findings (R1-R16) are recorded in the Review section for triage at the next gate.
- 2026-10-01: implement resumed on the return. At the question gate the user chose to add the NEWS entry and fix R1-R4, R8, R10, and R11 now. R5 (the engines Waller line) stays out of scope, with a candidate row to be added at review. Minor amendment: added T7 (NEWS) and T8 (fixes).
- 2026-10-01: T8 done. R1: the Discussion now says "the PCA and EFA engines" (boot_edges and comparability both exclude ESEM). R2: the matrix caveat moved before the test-suite and rotation sentences, so "What the rotation changes" follows the orthogonality sentence again. R3: "standardized items" in the manuscript and intro vignette. R4: the conclusion names what was checked. R8: both passages list polychoric, Spearman, pairwise, and FIML as uncovered. R10: the front matter restores "During the preparation of this work", drops the repeated "models", and reads "No AI system is an author". R11: `W_a′RW_b` with the prime character, reflowed, and "full-information" hyphenated. The intro stamp was updated, and the freshness and prose checks are clean. The em-dash count is 5.

## Decisions

## Review

Evidence gathered 2026-10-01 at head 7c9d79a. Master has not moved since the branch was cut.

- AC1 evidence: the full paragraphs that match `exact` were printed from all three files. The abstract no longer contains "exact". The methods paragraph (ms:152) says the algebra is exact for the matrix given, and that observed-item scores do not in general reproduce the edges under polychoric, Spearman, pairwise, or FIML matrices. The intro paragraph (intro.orig:276) and the engines PCA paragraph (engines.orig:82) say the same. The other matches (ms:404, intro.orig:91, 133, 137, 333, 360, engines.orig:22, 125, 538) do not say that edges equal correlations of observed-item scores. ms:404 calls edges "correlations among factor scores, not parameters of a fitted model", which contrasts them with model parameters and makes no claim about observed-item scores. Pass.
- AC2 evidence: the paragraphs that match `materiali|falls? back|nonlinear` are ms:152, ms:384, intro.orig:276, and engines.orig:82. None says the package computes edges from materialized scores. ms:152's "direct route materializes" describes the general route, cited to Grice. The methods passage (ms:152) and the conclusion (ms:384) each say every edge the package reports comes from the closed-form algebra. Both limit the scores check to all three engines (all linear) with complete data and Pearson correlations, and both exclude the polychoric and FIML paths. Pass.
- AC3 evidence: the paragraphs that match `waller`, reference entries aside (intro.orig:444, engines.orig:658), are ms:102, ms:152, intro.orig:276, and engines.orig:82. ms:102 and intro.orig:276 credit a closed form for principal components. ms:152 credits the components result via transformation matrices and the oblique form for rotated components, which the page images confirm (pp. 748-749). ms:152 and intro.orig:276 state `W'RW` for any linear weights as following from the covariance of linear composites, with no source and no novelty claim. engines.orig:82 keeps "Waller (2007) showed that the between-level algebra (`W'RW`) holds exactly for components". That credit is limited to components and states no general form, and Scope Out lists it as correct. Pass.
- AC4 evidence: the front matter (ms:19-21) reads "From June 2026 to submission, the author used large language models (Anthropic's Claude Opus, Sonnet, and Fable models ...)". The "Use of generative AI" section (ms:479-482) reads "AI-assisted development ... ran from June 2026 to submission ... with models from its Opus, Sonnet, and Fable families." Pass.
- AC5 evidence: `quarto render manuscript.qmd`, run in `manuscript/` at 09:24:16, exited 0. It wrote `manuscript.docx` (09:24:39) and `manuscript.pdf` (09:24:43). The em-dash count is 5 on the branch and 5 on `master`. Pass.
- AC6 evidence: `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` exited 0. The run started at 7c9d79a, and the only commit made during it (f0968dd) touched this file alone. Vignette freshness, prose, and code-unchanged were clean. `check()` gave 0 errors, 0 warnings, 0 notes, and coverage was 100%. Style, lint, and the pkgdown index were clean. Pass.

**Consistency gate (2026-10-01).** `cairn_validate.py` passed (16 advisory work-log WARNs, all in M84). `devtools::document()` gave no diff. README was untouched. pkgdown and `check()` were clean through the gate above. **Fail:** the profile requires a NEWS.md entry for this milestone's user-visible changes. The intro and engines vignette corrections ship with the package, and the development section of NEWS.md has no entry for them. Earlier vignette rewrites got entries (NEWS.md:69-85). Defect return 1.

**Independent review (2026-10-01).** Three lenses ran: Opus diff-bug, Sonnet blame-history, and Sonnet prior-review (PR-comment probe empty, archive and LESSONS evidence used). Findings are merged across lenses and ranked most severe first. The dispositions are proposals, and the user triages them at the next approval gate.

- R1 (diff-bug 1, prior-review 1): ms:459 says comparability and bootstrap CIs exist "for the linear (PCA and EFA) engines". That now contradicts "All three engines score linearly" (ms:164, 393). Proposed: fix now ("for the PCA and EFA engines").
- R2 (diff-bug 2): ms:170-175. The new matrix caveat sits between "Nothing in it uses the orthogonality of a rotation." and "What the rotation changes ...", so the contrast breaks (LESSONS M92). Proposed: fix now by moving the caveat before the rotation sentences.
- R3 (diff-bug 5): ms:158-163 and intro.orig:278-280 call scores "linear composites of the items". With R a correlation matrix, `W_a' R W_b` is the covariance of composites of standardized items. Proposed: fix now ("standardized items").
- R4 (diff-bug 6): ms:399-401. "inherit a computation that has already been checked" follows the sentence saying the check misses the polychoric path, which the paper's own Big Five example uses. Proposed: fix now by naming what was checked.
- R5 (diff-bug 3, prior-review 4, blame 1, claim audit): engines.orig:83, "Waller (2007) showed that the between-level algebra (`W'RW`) holds exactly for components". Waller writes transformation matrices, not `W'RW`. Scope Out keeps this sentence. Proposed: follow-up, or a gated scope amendment if the user wants it fixed here.
- R6 (blame 2): DESCRIPTION:15-17 says edges are "computed with exact linear algebra (Waller, 2007) or from materialized scores". That offers the internal scores route as shipped (D-038) and credits Waller with the general form. Not in scope. Proposed: follow-up candidate row.
- R7 (diff-bug 4, prior-review 2 and 3): unedited text calls edges score correlations without qualification: ms:147-148, ms:407, intro.orig:55-56, intro.orig:140 ("identity is exact for any fixed linear scoring"), R/ackwards.R:21, and ordinal.orig:384 ("come from the exact algebra"). Proposed: follow-up candidate row (repo-wide sweep, LESSONS M76).
- R8 (diff-bug 7): ms:167-169 lists two exclusions (polychoric, FIML), while the next sentences name four bases. Proposed: fix now (drop the sentence, because "only with complete data and Pearson correlations" already covers it).
- R9 (diff-bug 8, blame 2 DESIGN part): R/compute_edges.R:11-16 roxygen says "exact for PCA and EFA" and describes a nonlinear fallback. DESIGN.md:237 says Waller "derived this for orthogonal components". Proposed: follow-up, absorbed into the existing `[low]` DESIGN §5 row.
- R10 (diff-bug 12, blame 5): ms:19-24 repeats "models" and keeps the singular "The AI system is not an author". It also drops the policy-style opening "During the preparation of this work". Proposed: fix now (wording only).
- R11 (diff-bug 10, 11): notation and style. intro.orig:279 uses `W_a' R W_b` (ASCII) where line 141 uses `W′RW`, and the line runs to 88 characters. ms:169 has "full information" unhyphenated. Proposed: fix now.
- R12 (diff-bug 9): intro.orig:336-338, "For pure PCA on Pearson correlations ... agree exactly", also holds for any linear weights. Pre-existing and not wrong. Proposed: reject.
- R13 (blame 3): a hard line break inside a ms:170s sentence. Markdown treats it as a soft break, so nothing changes. Proposed: reject.
- R14 (blame 4): the scores-check scope is narrower than the M93 F3 candidate row. The manuscript's positive wording stays accurate. Proposed: reject (noted).
- R15 (diff-bug 13): the `claim audit:` work-log line has no date and contains an em dash. That is the fixed shape the implement skill mandates. Proposed: reject.
- R16 (diff-bug 14): the model list is incomplete if any subagent ran on Haiku. cairn never uses Haiku. Proposed: reject.
