# M92: Manuscript wording for the oblique rotation option

- **Status:** in-progress
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

- [ ] AC1: The methods passage states that the `W'RW` algebra is exact for any fixed linear scoring weights, oblique rotations included. It no longer says that an orthogonal rotation is what makes the algebra exact. It states that varimax [@kaiser1958] is the default and that oblique rotations are an option.
- [ ] AC2: The scope paragraph no longer says that rotation is varimax only or that correlated factors confound the cross-level signal. It no longer says that the closed-form algebra depends on orthonormality. It states that varimax is the default. It also states three facts about an oblique rotation. Each edge is a total correlation, `tidy()` reports the partialled coefficient `beta` beside it, and primary parents still follow `r`.
- [ ] AC3: Search `manuscript/manuscript.qmd` case-insensitively for `oblique|orthogon|orthonormal|varimax|rotat|correlated|confound`. No sentence that contains a match says that the package offers only orthogonal rotation, or that the edge algebra needs an orthogonal rotation. No such sentence says that correlated factors confound the cross-level signal or make the edges uninterpretable. A sentence can say that an oblique edge includes overlap through correlated factors at the same level.
- [ ] AC4: The rewritten passages contain no em dash, typed as `—` or as `---`. Every citation key that they use resolves in `manuscript/references.bib`.
- [ ] AC5: With the package and a TeX toolchain installed, `quarto render manuscript/manuscript.qmd` exits 0. It writes both `manuscript.pdf` and `manuscript.docx`.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T1, T2, T3
- AC5 → T3

## Tasks

- [x] T1: Read the RR01 verdict (`cairn/reviews/archive/RR01-oblique-algebra-claim.md`), D-034, D-036, and the `rotation` entry of `?ackwards` (R/ackwards.R:169). Then rewrite the methods passage (manuscript.qmd:162). Keep the `@waller2007` and `@kaiser1958` citations.
- [ ] T2: Rewrite the scope paragraph (manuscript.qmd:412). Keep its sentences on the sequential hierarchy and on Schmid-Leiman unchanged.
- [ ] T3: Run the AC3 search and record each matching sentence with its disposition in the work log. Look for em dashes and unresolved citation keys in the rewritten passages. Render to PDF and docx.

## Work log

- 2026-09-30: created by /milestone-plan, from the candidate row added at the M90 review (finding F2).
- 2026-09-30: criteria audit ran in full mode (fresh-context Opus reader) and returned 7 findings, all fixed at the plan gate. The scope now cites D-034 for the withdrawn confound rationale, because RR01 had kept that sentence. AC3 separates the forbidden claim from the allowed overlap statement, and its search gained `correlated|confound`. AC1 says "fixed linear". AC2 adds that primary parents still follow `r`. AC4 counts `---`, and AC5 states its toolchain.
- 2026-09-30: plan gate chose a planned milestone over a later direct docs commit, at the owner's selection. The rewrite then gets a criteria audit and a review. Falsified by a review that finds no defect at all, which shows that a direct commit was enough.
- 2026-09-30: implement started on branch m092-manuscript-oblique-wording. No question gate, because the plan left no choice open.
- 2026-09-30: T1 done. The methods passage now says the closed form is exact for any fixed linear scoring weights and cites Waller's oblique form (section 3, per `references/waller2007.md`). It names varimax as the default and oblique rotations as an option.

## Decisions

## Review
