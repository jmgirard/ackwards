# M67: Re-verify the 9 single-source reference pages against their shelf PDFs

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** IP1, IP2, IP8
- **Branch/PR:** `m67-reference-verification` · https://github.com/jmgirard/ackwards/pull/71

## Goal

Read every standing fact on the nine single-source `cairn/references/` pages back
against its shelf PDF, so each page's extraction status states a dated verification
instead of an unverified first pass.

## Scope

**In:** the nine pages with a 1:1 shelf source — `lorenzoseva2006`, `waller2007`,
`tenberge1999`, `goldberg2006`, `asparouhov2009`, `everett1983`, `goldberg1990`,
`saucier1997`, `saucier2005`. Verification target is the **standing facts** as
tracking-rules defines them: extracted values, printed formulas, page/table/equation
anchors, verbatim wordings. Discrepancies are corrected in place and marked.

**Out:**
- The three collapsed synthesis pages (`applications`, `background`, `rotation-and-k`)
  and their 17 member sources → M68.
- The repo's own interpretive prose on each page ("why this matters to us",
  "relation to our implementation") — it is repo voice, not a claim about the
  source, and nothing in the PDF can settle it.
- `forbes2023.md` — already verified by the M44 + M53 fidelity suite.
- Fixing any wrong fact that reached user-facing text (roxygen, `man/`, vignettes,
  NEWS) → `/hotfix`, enumerated by T8. This milestone corrects the reference page
  and hands off; it does not grow to match what verification happens to find.

## Acceptance criteria

- [ ] All nine pages carry an extraction status naming a verification verb
      (`verified` / `checked against` / `read against` / `read directly`) and its own
      `— observed YYYY-MM-DD`; no page still reads `unverified — first pass`.
- [ ] `cairn_validate`'s `references staleness` advisory drops from 12 to exactly 3,
      those 3 being the collapsed pages M68 owns.
- [ ] Every standing fact on each of the nine pages is checked against its shelf PDF;
      each page's work-log line names what was checked and what, if anything, was wrong.
- [ ] Each corrected fact is grepped across `R/`, `man/`, `vignettes/`, and `NEWS.md`;
      any propagation is enumerated with `file:line` and routed to `/hotfix`. A
      milestone that corrects a page silently leaving wrong user-facing text fails.
- [ ] Prior verification is absorbed, not re-derived: `goldberg2006`'s issue number
      (M63, Crossref) and `saucier1997`'s DOI + fn-14 evidence (M61) are cited in
      their status lines rather than re-checked from scratch.
- [ ] Docs-only — no files outside `cairn/` modified; `cairn_validate` all checks PASS.

## Coverage

- AC1 → T1, T2, T3, T4, T5, T6, T7
- AC2 → T9
- AC3 → T1, T2, T3, T4, T5, T6, T7
- AC4 → T8
- AC5 → T4, T7
- AC6 → T9

## Tasks

- [x] T1: `lorenzoseva2006.md` — highest stakes; the numbers backing the package's
      `.95` `redundancy_phi` default (56 practitioners, 448 loading-column pairs,
      φ ∈ {.62….97}, r = .974, the .85–.94 "fair" band, φ > .95 "equal", their Eq. 1).
- [x] T2: `waller2007.md` — Eq. 14, Eqs. 9–10, Eqs. 6–7, the §3 oblique extension, the
      §4 estimated-scores caveat, and Appendix A's `BASS(R, maxP, Print)` signature and
      column-sum-positive sign convention (backs IP1/IP2).
- [x] T3: `tenberge1999.md` — the LCP definition and closed-form solution behind the
      `tenBerge` scoring default; three page/table anchors.
- [x] T4: `goldberg2006.md` — six anchors incl. the stopping criterion and the ≥ .30
      display threshold; absorb M63's Crossref-verified issue number rather than re-derive.
- [x] T5: `asparouhov2009.md` — the ESEM definition and the WLSMV ordinal basis claim.
- [x] T6: `everett1983.md` + `goldberg1990.md` — the `comparability()` precedent pair;
      check the ⚠ mis-attribution block on `goldberg1990` still states the source correctly.
- [x] T7: `saucier1997.md` + `saucier2005.md` — the Goldberg-lab split-half pair (fn 14's
      PA-overextraction note; the explicit .90 replication gate); absorb M61's DOI evidence.
- [x] T8: grep every corrected fact across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      enumerate propagation with `file:line`; open `/hotfix` if any is user-facing.
- [x] T9: re-run `cairn_validate`; confirm staleness 12 → 3 and all checks PASS; commit.

## Work log

- 2026-07-19: created by /milestone-plan.
- 2026-07-19: implement started on `m67-reference-verification`; no question gate (plan settled depth/split/defect handling; nothing open).
- 2026-07-19: review opened draft PR #71; master had not moved since the branch was cut (0 commits behind), no merge needed. AC evidence gathering in progress; dod-gate and three review lenses running.
- 2026-07-19: T8 propagation sweep — all nine corrected facts grepped across `R/`, `man/`, `vignettes/`, `NEWS.md`, `README`: **zero hits**, every correction confined to the reference pages, so no `/hotfix` is owed. Spot-checked how the package actually cites Saucier et al. (2005): the .90 benchmark, Everett's 81%-shared-variance rationale, the single-split precedent, and `n_splits = 10` owned as our own choice are all accurate against the verified sources.
- 2026-07-19: T9 — `cairn_validate` all checks PASS, `references staleness` 12 → 3 (exactly the three collapsed pages M68 owns); branch diff is `cairn/`-only, so the r-package `verify` slot (devtools::test/document) is not triggered — no R code, roxygen, or generated file changed. Status → review.
- 2026-07-19: T7 `saucier1997` + `saucier2005` verified. saucier1997 (pp. 1296–1305): the 500 descriptors, N = 700/201, PCA+varimax 2–10, the "single estimate of factor reliability" averaging, maximize-magnitude matching, the variables-split design, the p. 1304 Everett quotation (ellipsis correct), the "well below .70" drop-off after 7 (5 for dispositions), and fn 14's PA quotation all exact; one correction ("nested" selections — paper never says nested, and containment fails between the two disposition selections). saucier2005: all three quoted passages verbatim, 3,302 adjectives and the emic 6-factor claim exact — but **two errors in the summary line**: "N ≈ 991 + 201" (no "201" anywhere in the paper; Samples are 991/429/538 — the 201 is saucier1997's acquaintance N, cross-contaminated via the shared ingest commit `254e023`) and "PCA + varimax" ("varimax"/"oblique"/"promax" occur zero times; the paper says only "1 unrotated and 2 to 10 rotated factors"). Also retired that page's undated second "verified against the PDF" note in favour of the provenance block.
- 2026-07-19: T6a `everett1983` verified against pp. 197–204 — the F₁ₜ = S₁·Vₜ procedure, comparability-vs-congruence, the ≥ .90 / 26° / 81% passage, sub-population splits, and all six Table 2 values (.98/.97/.96 at k=3; .99/.94/.82/.04 at k=4) exact; no corrections. Added the p. 204 greedy-with-removal matching procedure — source-backed, and the bijection D-022 already uses.
- 2026-07-19: T6b `goldberg1990` verified against pp. 1216–1224 — Study 1/2/3 counts, the 10 method combinations, .950–.996, 30-of-3,750, the verbatim 13-factor invariance quote, and Study 2's .86–.94 (mean .91) congruence all exact. ⚠ block's core negative claim re-confirmed (no split-halves anywhere; Studies 2–3 use independent samples). Two corrections: its "docs currently over-attribute" was stale — all four flagged sites were repointed by M61/M63 and none survives (grep-verified); and "established the Big Five label-set" contradicted p. 1217, which calls I–IV traditional.
- 2026-07-19: T5 `asparouhov2009` verified against pp. 397–401 — m² restrictions, Jennrich-extending rotated SEs + fit tests, and Eq. 3's probit/underlying-normal formulation with limited-information WLS (Muthén 1984) all exact. One correction: the Browne (2001, p. 113) quotation had dropped "the examination of model" with the elision unmarked; restored verbatim with the p. 398 block-quote location. No propagation into user-facing text. Issue number not printed in source.
- 2026-07-19: T4 `goldberg2006` verified against pp. 347–358 — stopping criterion, one-large-correlational-analysis step, adjoining-levels path coefficients, FUPC, `4/1` size labels, "virtually identical", the ≥ .30 threshold, all six figure descriptions, the §5 sequential/part-whole caveats, and the Waller footnote all exact; no corrections. M63's Crossref issue number absorbed by citation (not printed in the source). Also cross-confirmed Ashton/Lee/Goldberg 2004 and Yung/Thissen/McLeod 1999 in the reference list (the latter feeds M68's `background.md`).
- 2026-07-19: T3 `tenberge1999` verified against pp. 311–318 — Eq. 3, Eqs. 5/7/9, all three loss functions, Thm 1/2/3, Cor. 1, and the orthogonality-coincidence conclusion all exact; no corrections. Checked the one suspicious anchor: Eq. 9 first appears in Lemma 1 (singular Ψ) but Thm 1 is precisely about it, so the page's "(Eq. 9, Thm 1)" is correct. Issue number not printed in source.
- 2026-07-19: T2 `waller2007` verified against pp. 745–752 — Eq. 14, S = [I|0], Eqs. 9–10, §3 oblique form, §4 caveat (Guttman 1955; McDonald & Mulaik 1979), Appendix A signature + column-sum-positive sign convention all exact. Two corrections: the "Eqs. 6–7" anchor over-attributed the rotated `W_i = QΛ^{-1/2}T_i` (Eq. 7 is unrotated `W`, Eq. 6 is `P`); dropped the unsourced "first public R implementation" priority claim (M64 lesson). Noted the issue number is not printed in the source.
- 2026-07-19: T1 `lorenzoseva2006` verified against pp. 57–61 — all 14 standing facts confirmed (Eq. 1, four invariance properties, prior thresholds .80/.90/MacCallum bands, 56 judges, φ set, r = .974, .85–.94 fair, >.95 equal, judge-disparity and label/variable-count findings). No errors. Clarified 48 distinct pairs vs 448 evaluations — both printed in the source.

## Decisions

## Review
