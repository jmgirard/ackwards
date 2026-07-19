# M67: Re-verify the 9 single-source reference pages against their shelf PDFs

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** IP1, IP2, IP8
- **Branch/PR:** `m67-reference-verification`

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
- [ ] T3: `tenberge1999.md` — the LCP definition and closed-form solution behind the
      `tenBerge` scoring default; three page/table anchors.
- [ ] T4: `goldberg2006.md` — six anchors incl. the stopping criterion and the ≥ .30
      display threshold; absorb M63's Crossref-verified issue number rather than re-derive.
- [ ] T5: `asparouhov2009.md` — the ESEM definition and the WLSMV ordinal basis claim.
- [ ] T6: `everett1983.md` + `goldberg1990.md` — the `comparability()` precedent pair;
      check the ⚠ mis-attribution block on `goldberg1990` still states the source correctly.
- [ ] T7: `saucier1997.md` + `saucier2005.md` — the Goldberg-lab split-half pair (fn 14's
      PA-overextraction note; the explicit .90 replication gate); absorb M61's DOI evidence.
- [ ] T8: grep every corrected fact across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      enumerate propagation with `file:line`; open `/hotfix` if any is user-facing.
- [ ] T9: re-run `cairn_validate`; confirm staleness 12 → 3 and all checks PASS; commit.

## Work log

- 2026-07-19: created by /milestone-plan.
- 2026-07-19: implement started on `m67-reference-verification`; no question gate (plan settled depth/split/defect handling; nothing open).
- 2026-07-19: T2 `waller2007` verified against pp. 745–752 — Eq. 14, S = [I|0], Eqs. 9–10, §3 oblique form, §4 caveat (Guttman 1955; McDonald & Mulaik 1979), Appendix A signature + column-sum-positive sign convention all exact. Two corrections: the "Eqs. 6–7" anchor over-attributed the rotated `W_i = QΛ^{-1/2}T_i` (Eq. 7 is unrotated `W`, Eq. 6 is `P`); dropped the unsourced "first public R implementation" priority claim (M64 lesson). Noted the issue number is not printed in the source.
- 2026-07-19: T1 `lorenzoseva2006` verified against pp. 57–61 — all 14 standing facts confirmed (Eq. 1, four invariance properties, prior thresholds .80/.90/MacCallum bands, 56 judges, φ set, r = .974, .85–.94 fair, >.95 equal, judge-disparity and label/variable-count findings). No errors. Clarified 48 distinct pairs vs 448 evaluations — both printed in the source.

## Decisions

## Review
