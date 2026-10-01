# M92: Manuscript wording for the oblique rotation option

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP3
- **Resolves:** —
- **Surface tier:** user-facing, because journal reviewers and readers outside the repo read the manuscript
- **Branch/PR:** m092-manuscript-oblique-wording

## Goal

The manuscript describes the shipped `rotation` option: varimax by default, oblique rotations as an option, and edge algebra that holds without orthogonality.

## Scope

**In:** the two passages of `manuscript/manuscript.qmd` that state a rotation restriction. One is the methods passage on computing between-level correlations (lines 162 to 166 at plan time). The other is the scope paragraph of the Discussion (lines 412 to 416). Wording follows the algebra verdict of RR01, D-034, D-036, and the `rotation` entry of `?ackwards`. D-034 holds that oblique edges stay interpretable total correlations, and it withdrew the confound rationale. A source comment marks the Discussion as author-owned. The owner chose this milestone at the 2026-09-30 plan gate.

**Out:** package code, roxygen, and vignettes, which M90 already updated. New manuscript content on oblique results, such as an oblique worked example, is not queued. The numbers that the manuscript computes at render time stay as they are, because the varimax default did not change.

## Acceptance criteria

- [x] AC1: The methods passage states that the `W'RW` algebra is exact for any fixed linear scoring weights, oblique rotations included. It no longer says that an orthogonal rotation is what makes the algebra exact. It states that varimax [@kaiser1958] is the default and that oblique rotations are an option.
- [x] AC2: The scope paragraph no longer says that rotation is varimax only or that correlated factors confound the cross-level signal. It no longer says that the closed-form algebra depends on orthonormality. It states that varimax is the default. It also states three facts about an oblique rotation. Each edge is a total correlation, `tidy()` reports the partialled coefficient `beta` beside it, and primary parents still follow `r`.
- [x] AC3: Search `manuscript/manuscript.qmd` case-insensitively for `oblique|orthogon|orthonormal|varimax|rotat|correlated|confound`. No sentence that contains a match says that the package offers only orthogonal rotation, or that the edge algebra needs an orthogonal rotation. No such sentence says that correlated factors confound the cross-level signal or make the edges uninterpretable. A sentence can say that an oblique edge includes overlap through correlated factors at the same level.
- [x] AC4: The rewritten passages contain no em dash, typed as `—` or as `---`. Every citation key that they use resolves in `manuscript/references.bib`.
- [x] AC5: With the package and a TeX toolchain installed, `quarto render manuscript/manuscript.qmd` exits 0. It writes both `manuscript.pdf` and `manuscript.docx`.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T1, T2, T3
- AC5 → T3

## Tasks

- [x] T1: Read the RR01 verdict (`cairn/reviews/archive/RR01-oblique-algebra-claim.md`), D-034, D-036, and the `rotation` entry of `?ackwards` (R/ackwards.R:169). Then rewrite the methods passage (manuscript.qmd:162). Keep the `@waller2007` and `@kaiser1958` citations.
- [x] T2: Rewrite the scope paragraph (manuscript.qmd:412). Keep its sentences on the sequential hierarchy and on Schmid-Leiman unchanged.
- [x] T3: Run the AC3 search and record each matching sentence with its disposition in the work log. Look for em dashes and unresolved citation keys in the rewritten passages. Render to PDF and docx.

## Work log

- 2026-09-30: created by /milestone-plan, from the candidate row added at the M90 review (finding F2).
- 2026-09-30: criteria audit ran in full mode (fresh-context Opus reader) and returned 7 findings, all fixed at the plan gate. The scope now cites D-034 for the withdrawn confound rationale, because RR01 had kept that sentence. AC3 separates the forbidden claim from the allowed overlap statement, and its search gained `correlated|confound`. AC1 says "fixed linear". AC2 adds that primary parents still follow `r`. AC4 counts `---`, and AC5 states its toolchain.
- 2026-09-30: plan gate chose a planned milestone over a later direct docs commit, at the owner's selection. The rewrite then gets a criteria audit and a review. Falsified by a review that finds no defect at all, which shows that a direct commit was enough.
- 2026-09-30: implement started on branch m092-manuscript-oblique-wording. No question gate, because the plan left no choice open.
- 2026-09-30: T1 done. The methods passage now says the closed form is exact for any fixed linear scoring weights and cites Waller's oblique form (section 3, per `references/waller2007.md`). It names varimax as the default and oblique rotations as an option.
- 2026-09-30: T2 done. The scope paragraph now names varimax as the default, calls an oblique edge a total correlation, and points to `beta` in `tidy()` and to primary parents following `r` (read against `R/tidy.R` and the fit-time advisory in `R/ackwards.R`). The sequential-hierarchy and Schmid-Leiman sentences are unchanged.
- 2026-09-30: T3 AC3 search hit 9 lines in 5 sentences, all in the rewritten passages. Methods: "exact for any fixed linear scoring weights, oblique rotations included" (allowed, states the algebra holds), "Nothing in it uses the orthogonality" (allowed), "choice of rotation therefore changes" (allowed), "varimax rotation at every level" (allowed, default), "Oblique rotations ... available as an option" (allowed). Scope: "Varimax is the default ... oblique rotations are an option" (allowed), "each edge is a total correlation, which includes overlap through correlated factors" (allowed overlap statement). No match elsewhere in the file.
- 2026-09-30: T3 checks. No `—` or `---` in either rewritten passage. Keys `waller2007`, `kaiser1958`, `grice2001` each occur once in `references.bib`. `quarto render manuscript/manuscript.qmd` exited 0 after `devtools::install()` of the branch (the package was not installed), and wrote `manuscript.pdf` and `manuscript.docx`. One render warning: no `rsvg-convert` for the docx ORCID icon, unrelated to this change.
- 2026-09-30: correction to the T3 search line. The hits fall in 7 sentences, not 5, as the line's own list of 7 shows.
- 2026-09-30: claim audit: 11 claims read, 2 corrected — manuscript/manuscript.qmd. Varimax applies "at every level with two or more factors" (level 1 is not rotated), and Waller gives "the oblique form for rotated components". The same reader re-read both, and both hold.
- 2026-09-30: after the audit fixes, the AC3 search still hits only the two rewritten passages, the methods passage has no `—` or `---`, and `quarto render` exited 0 again with both outputs written. No R code or roxygen changed, so the profile's `devtools::test()` step does not apply. Status set to review.

## Decisions

## Review

Fresh evidence, 2026-09-30, on branch head caca670. The branch contains `origin/master` 0324d51, and no PR exists yet.

- AC1: read manuscript.qmd:162-169 from `git diff master...HEAD`. It says "The closed form is exact for any fixed linear scoring weights, oblique rotations included" and "Nothing in it uses the orthogonality of a rotation". The old "Orthogonal (varimax) rotation ... is what makes the algebra exact" and "$T' = T^{-1}$" sentences are gone. It names varimax [@kaiser1958] as the default and oblique rotations as an option. Pass.
- AC2: read manuscript.qmd:415-420. The removed lines held "orthogonal (varimax) only", "would confound the very cross-level signal", and "depends on the orthonormality". The new text says varimax is the default and an oblique edge is a total correlation. It says `tidy()` reports `beta` beside each edge and primary parents follow `r`. Pass.
- AC3: `grep -n -i -E 'oblique|orthogon|orthonormal|varimax|rotat|correlated|confound'` hits lines 162-168 and 415-418 only, all inside the rewritten passages. No hit says the package offers only orthogonal rotation, that the algebra needs it, or that correlated factors confound or make edges uninterpretable. Line 417-418 is the allowed overlap statement. Pass.
- AC4: grep for `—|---` over lines 155-175 and 410-425 finds none. Keys used in the rewritten passages (`waller2007`, `kaiser1958`) each match once in `manuscript/references.bib`, as do the neighboring keys `forbes2023`, `goldberg2006`, `schmid1957`, `williams2025`, `yung1999`. Pass.
- AC5: after deleting old outputs, `quarto render manuscript.qmd` exited 0 and wrote `manuscript.pdf` (96,384 bytes) and `manuscript.docx` (308,836 bytes). The PDF text contains "oblique" 3 times. One warning: no `rsvg-convert` for the docx ORCID icon, unrelated. The installed package is 0.2.0.9000, and the branch changes no package code. Pass.

Consistency gate, 2026-09-30:

- `cairn_validate.py` exited 0. All checks passed, with 16 work-log format warnings, all in M84.
- No principle text changed, so `cairn_impact.py` does not apply.
- `devtools::document()` left no diff. `devtools::check()` gave 0 errors, 0 warnings, 0 notes. `pkgdown::check_pkgdown()` found no problems.
- No NEWS entry is needed, because `manuscript/` is in `.Rbuildignore` and is not part of the package. No README or new top-level file changed.

Independent review, 2026-09-30, full three-lens fan-out (user-facing tier). The prior-review lens found one applicable lesson (no em dashes, M74), not regressed, and zero findings. No finding shows a criterion failing, so none triggers a return. The diff lens (D) and the blame lens (B) ranked these findings. Dispositions wait for the gate.

- D1: the scope paragraph omits the cost that D-036 accepted. Under oblique rotation a primary parent can be a factor that only correlates with the real parent. The `prune()` thresholds were calibrated under varimax. The fit-time message says both (R/ackwards.R:1088-1092).
- D2: the manuscript never says why varimax is the default or that `r` equals `beta` under varimax.
- D3: "partialled coefficient `beta`" is not defined.
- D4 (same as B2): the scope paragraph's flow broke. "Scope is bounded" now leads into a widening, and "The resulting hierarchy" follows the `beta` sentence.
- D5: "therefore changes what an edge means" does not follow from the sentences before it.
- D6: "Under an oblique rotation, each edge is a total correlation" implies varimax edges are something else.
- D7: "correlated factors at the same level" does not say it is the parent level.
- D8: "exact" needs a qualifier for a polychoric R, where no observed scores reproduce the edges.
- D9 (pre-existing): manuscript.qmd:160-162 says the package falls back to scores for nonlinear scoring and cross-checks the routes. No shipped engine is nonlinear, and the cross-check runs only in tests.
- D10 (same as B4, pre-existing): manuscript.qmd:154 credits Waller with the general weight-matrix form, which `references/waller2007.md` calls ours.
- D11: "available as an option" omits that some oblique rotations need the suggested package GPArotation.
- D12 (same as B5, pre-existing): the AI-use disclosure (manuscript.qmd:20 and 456) says June to July 2026.
- B1: the Discussion source comment still says "author-owned stub", which is stale since M74.
- B3: `beta` can be `NA` for a level whose score correlation cannot be inverted, so "beside each edge" is broad.
- B6: "keeps the factors within a level uncorrelated" holds for factors, not always for regression scores. The "X, not Y" contrast in line 165 is a form an earlier style pass removed.
