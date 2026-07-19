# M68: Re-verify the 3 collapsed synthesis pages against their member sources

- **Status:** review
- **Priority:** normal
- **Depends on:** M67
- **Branch:** `m68-reference-verification-collapsed-pages`
- **Principles touched:** IP8

## Goal

Read every standing fact on the three collapsed `cairn/references/` synthesis pages
back against the shelf PDFs of their 17 member sources, clearing the last of the
`references staleness` advisory.

## Scope

**In:** `rotation-and-k.md` (6 members: crawford1970, browne2001a, horn1965,
velicer1976, lim2019, achim2021), `applications.md` (6 members: kim2015, markon2005,
wright2014a, forbush2018, forbush2024, cowan2024), `background.md` (5 members:
kotov2017, schmid1957, yung1999, saucier1996, widiger2019). All 17 have shelf PDFs.
Same bar as M67 — standing facts only, corrections marked in place.

**Out:**
- The nine single-source pages → M67 (this milestone depends on it).
- Interpretive prose framing each collapsed set (the italic scope note, the
  "precedent for" clauses) — repo voice, not a source claim.
- User-facing propagation of any corrected fact → `/hotfix`, enumerated by T4.
- Splitting any collapsed page back into per-source notes. These pages were
  deliberately collapsed; if verification argues one member deserves its own page,
  that is a candidate row, not a scope grab.

## Acceptance criteria

- [x] All three pages carry an extraction status naming a verification verb and its
      own `— observed YYYY-MM-DD`.
- [x] `cairn_validate`'s `references staleness` advisory reaches 0.
- [x] Every member entry's standing facts are checked against that member's shelf PDF;
      each page's work-log line names the members checked and any correction made.
- [x] M61's already-verified `lim2019` / `achim2021` entries are absorbed by citation,
      not re-derived — including its finding that a vignette clause overclaimed
      `lim2019` (best among 13 PA *variants*, not all criteria).
- [x] Each corrected fact is grepped across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      any propagation is enumerated with `file:line` and routed to `/hotfix`.
- [x] Docs-only — no files outside `cairn/` modified; `cairn_validate` all checks PASS.

## Coverage

- AC1 → T1, T2, T3
- AC2 → T5
- AC3 → T1, T2, T3
- AC4 → T1
- AC5 → T4
- AC6 → T5

## Tasks

- [x] T1: `rotation-and-k.md` — the varimax-default and k-advice supports; check the
      CF-varimax κ = 1/p claim, Horn's PA and Velicer's MAP descriptions, and absorb
      M61's verified DOIs for `lim2019` / `achim2021`.
- [x] T2: `applications.md` — six published-application entries; verify each citation's
      standing facts and the "precedent for" claims that name a specific finding.
- [x] T3: `background.md` — five framing entries; check the Schmid–Leiman and
      Yung et al. descriptions especially, since DESIGN cites them for what is
      deliberately out of scope, and the `saucier1996` "ruled out as comparability
      source" note against what that paper actually says.
- [x] T4: grep every corrected fact across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      enumerate propagation with `file:line`; open `/hotfix` if any is user-facing.
- [x] T5: re-run `cairn_validate`; confirm staleness 0 and all checks PASS; commit.

## Work log

- 2026-07-19: created by /milestone-plan.
- 2026-07-19: status → in-progress; branch `m68-reference-verification-collapsed-pages` cut from master.
- 2026-07-19: review pass — PR #72 opened; all 6 ACs verified with fresh evidence; `cairn_validate` exit 0; dod-gate PASSED; full CI matrix green. Independent review: blame-history and prior-PR lenses clean; diff-bug lens returned 9 findings, **6 of them introduced by this milestone's own correction pass** (scores 96/93/92/90/88/88) — all fixed on-branch: the wright2014a ⚠ block was inverted (the hierarchy used *orthogonal* Geomin, p. 49, and is a support for DESIGN §9 rather than a departure), markon2005's "those terms do not appear" was false (the paper has a section headed "The Big Two, Three, and Four."), and four page anchors were wrong. **Supersedes the T1 line below:** browne2001a's Table 1 is p. **118**, not p. 119. Three sub-threshold findings logged, not actioned (F7/68 DESIGN.md:418, F8/75, F9/45).
- 2026-07-19: T5 — `cairn_validate` all checks PASS, `references staleness` **0** (12 → 3 at M67 → 0 now); 80 dangling-id-token advisories are pre-existing and untouched by this diff. First run reported staleness 2: `Extraction:` must begin its own line for the parser to see it, and two pages had it mid-line after `per entry.` — fixed. Diff is `cairn/`-only (AC6). `Rscript tools/dod-gate.R` PASSED (check 0/0/0, coverage 100%, style/lint clean, pkgdown index complete); status → review.
- 2026-07-19: T4 propagation sweep — every corrected fact grepped across `R/`, `man/`, `vignettes/`, `NEWS.md`. Zero user-facing propagation, so no `/hotfix` is owed. One internal code comment is wrong (`R/engine_esem.R:5-6`, WLSMV attributed to forbush2024) — enumerated and dispositioned in Decisions. Confirmed clean: the intro vignette's `CF(κ = 1/p)` claim for kim2015 (exact, p. 1067) and `comparability()`'s citations (Everett 1983 / saucier1997 / saucier2005 — correctly not saucier1996).
- 2026-07-19: T3 `background.md` verified — **no corrections**; all five entries confirmed, including the two claims about goldberg2006 (p. 348's "received wisdom" quote and p. 356–357's "sequential" concession) and saucier1996's "ruled out as comparability source" negative (re-confirmed: "split"/"comparab" occur zero times in it). kotov2017's 40-author count verified exactly and its level definitions quoted from the paper. Sharpened saucier1996's `(.95–.84)` shorthand into the actual non-monotone series, and flagged its congruence correlations as not comparability coefficients. Noted kotov2017's shelf PDF is online-first (no journal pagination).
- 2026-07-19: T2 `applications.md` verified — kim2015 and forbush2018 accurate to the digit (kim2015 p. 1067 recorded as the verbatim source for DESIGN §9's orthogonal-rotation rationale and κ = 1/p); **three corrections** — markon2005's "goldberg2006 §5 singles it out" (Goldberg p. 356 names Saucier 1997; Markon is one of four "recent reports"), forbush2024's conflation of two results (92.4%/58.7% own variance vs the 0.88–334× / 1.95–80.8× comparisons, the low end below parity), cowan2024's "pure Goldberg recipe" (minres EFA + regression scores, not PCA + component scores). Added wright2014a's oblique geomin rotation and forbush2024's ULSMV estimator. One propagation found (`R/engine_esem.R:5-6`) — see Decisions.
- 2026-07-19: T1 `rotation-and-k.md` verified — crawford1970/browne2001a/horn1965/velicer1976 accurate (anchors added: CF Eq. 7 p. 323, orthomax equivalence pp. 324–326, Browne Table 1 p. 119 for κ = 1/p, MAP Eq. 9 p. 323); lim2019 abstract-verbatim, M61's variant-scoping absorbed as a ⚠ note; **achim2021 corrected** — the ~17% is the share of noise-factor simulations (16.8%, p. 70), not of "correct" retentions. Added crawford1970 p. 331's explicit varimax-when-k-unknown recommendation. Noted lim2019's shelf PDF is online-first (no journal pagination). No propagation into `R/`, `man/`, `vignettes/`, `NEWS.md`.

## Decisions

- **2026-07-19 (T2/T4) — one propagation found; not a hotfix.**
  `R/engine_esem.R:5-6` says ordinal data "uses WLSMV estimation (Kim & Eaton
  2015; **Forbush et al. 2024**)". Forbush et al. (2024, p. 633) used a
  mean-and-variance-adjusted **unweighted** least squares estimator (ULSMV),
  not WLSMV; kim2015 is the correct and only WLSMV precedent. **Disposition:**
  the line is an internal code comment — not rendered in any help page,
  vignette, or NEWS entry — so it is not user-visible and the hotfix bar
  (tracking-rules: "user-visible bug") is not met. It is also barred from this
  milestone by AC6 (docs-only, nothing outside `cairn/`). Routed as a
  **trivial** comment fix for a direct commit to master after this merge.
  The neighbouring `R/ackwards.R:14-15` varimax claim was checked and stands:
  forbush2024 used an orthogonal Crawford-Ferguson rotation, so "matches" is
  fair, though the paper does not state κ.

## Review

Reviewed 2026-07-19. PR: https://github.com/jmgirard/ackwards/pull/72

### Acceptance criteria — fresh evidence

- **AC1 — extraction status with a verification verb + own `— observed` stamp.**
  All three provenance blocks read `Extraction: verified 2026-07-19 (M68) …
  — observed 2026-07-19`, confirmed by extracting each block and matching both
  fields. Two mechanical fixes were needed to get here: `Extraction:` must
  begin its own line for the parser to see it (two pages had it mid-line), and
  rotation-and-k's stamp was hard-wrapped between `observed` and its date —
  both corrected on-branch.
- **AC2 — `references staleness` reaches 0.** `cairn_validate` reports
  `OK    references staleness` (was WARN (3) at M67's close, WARN (12) before
  it). Exit 0.
- **AC3 — every member's standing facts checked against its shelf PDF.** All
  17 members are named in their page and were read against the PDF in
  `cairn/references/sources/`: 6 in rotation-and-k, 6 in applications, 5 in
  background. Three work-log lines (T1/T2/T3) name the members checked and
  every correction made.
- **AC4 — M61 absorbed by citation, not re-derived.** M61's verified DOIs are
  carried forward unchanged (not re-derived), and its vignette-overclaim
  finding is cited at `rotation-and-k.md:106-108` as a standing ⚠ note
  scoping "performs best" to the 13 PA variants. Note the boundary: the
  achim2021 **17% figure was outside M61's scope** (M61 verified DOIs and the
  ±1-range claim), so AC3 required checking it — which produced this
  milestone's achim2021 correction.
- **AC5 — corrected facts grepped across `R/`, `man/`, `vignettes/`,
  `NEWS.md`.** Sweep run; **zero user-facing propagation**, so no `/hotfix` is
  owed. One non-user-facing hit enumerated with `file:line` and dispositioned
  in Decisions (`R/engine_esem.R:5-6`, WLSMV misattributed to forbush2024 —
  an internal code comment, routed to a trivial post-merge commit).
  Confirmed clean: the intro vignette's `CF(κ = 1/p)` claim for kim2015
  (exact, p. 1067) and `comparability()`'s citations (Everett 1983 /
  saucier1997 / saucier2005 — correctly not saucier1996).
- **AC6 — docs-only; `cairn_validate` all checks PASS.**
  `git diff --name-only master...HEAD | grep -v '^cairn/'` returns nothing.
  `cairn_validate` exit 0, all checks PASS (80 `dangling id tokens`
  advisories are pre-existing and untouched by this diff).

### Consistency gate

Universal: `cairn_validate` exit 0 — `scaffold present` PASS, `coverage
complete` PASS, `references staleness` OK. No IP/GP principle changed, so
`cairn_impact` was skipped.

### Independent review — 3 lenses + scorer

- **[S] blame-history:** no findings. Confirmed the diff does not undo M61's
  vignette-overclaim fix (`man/suggest_k.Rd` still correctly scoped), M63's
  Goldberg issue-number fix, or M67's asparouhov2009 elision / goldberg1990
  corrections; every deleted paragraph is replaced by an informational
  superset.
- **[S] prior-PR-comments:** no findings. PRs #71 and #62 carry no substantive
  GitHub review comments (Codecov bot only) — this repo reviews locally, so
  the lens no-ops. It cross-checked the diff against M67's two defect *classes*
  (unmarked elision, unsourced priority claim) and found neither reintroduced.
- **[O] diff-bug:** 9 findings. **The correction pass introduced 6 of them** —
  the same failure mode M67 hit, and the reason this lens exists.

**Actioned (score ≥ 80) — all 6 fixed on-branch:**

- **F1 (96) — the wright2014a ⚠ block was inverted.** It asserted the paper
  used *oblique* geomin for the hierarchy and told readers "do not cite it as
  support for" DESIGN §9's orthogonal-only stance. The oblique Geomin is from
  the *conjoint five-factor* EFA (p. 47); the **hierarchical** analysis used
  **orthogonal** Geomin (p. 49) for DESIGN §9's own reason, reached
  independently. Rewritten into a two-analyses/two-rotations block that cites
  it as support. Verified by page render.
- **F2 (93) — markon2005 "those terms do not appear in it" was false.** The
  paper has a section headed "The Big Two, Three, and Four." and explicitly
  relates Big Two → Big Three → Big Four (15/12/7 occurrences). My original
  grep was defective: the PDF text extraction strips spaces, so `Big Two`
  matched nothing. Rewritten to credit the labels to the paper; the true
  zero-occurrence claim ("bass-ackwards") retained.
- **F3 (92) — browne2001a Table 1 is p. 118, not p. 119**, wrong in all three
  places, and it anchors DESIGN §9's κ = 1/p indexing. Verified by page render.
- **F4 (90) — forbush2024 "Method detail (p. 633)" is p. 632**; only the
  ULSMV/"crawfer" sentence is on p. 633.
- **F5 (88) — cowan2024's orthogonal-rotation rationale is p. 7, not p. 8.**
  Verified by page render.
- **F6 (88) — widiger2019's bass-ackwards mention is p. 348, not p. 347.**
  Verified by page render.

**Below threshold — logged, not actioned (IP3):**

- **F7 (68) — `cairn/DESIGN.md:418` repeats the forbush2024/WLSMV error.**
  Scored down as arguably outside M68's stated scope (AC5's grep list named
  `R/`, `man/`, `vignettes/`, `NEWS.md`, not `cairn/`). The defect itself is
  real and independently confirmed. Folded into the same post-merge trivial
  follow-up as `R/engine_esem.R:5-6` (see Decisions) rather than fixed here,
  and surfaced to the user at the approval gate.
- **F8 (75) — kotov2017 "top-down definitions" framing.** The paper presents
  its levels "from narrowest to broadest"; the page lists them broadest-first.
  The four quoted definitions are verbatim-exact, so this is presentation
  framing, not a wrong fact.
- **F9 (45) — forbush2018's bare "Verified exact."** is an undated
  repo-state claim, but the provenance block already carries the dated stamp
  covering all six entries, and sibling pages use the same unstamped idiom.

Toolchain (`r-package` profile `consistency-gate` slot): `devtools::document()`
produced no diff; README.md in sync with README.Rmd; `pkgdown::check_pkgdown()`
reference index complete; no NEWS entry owed (zero non-`cairn/` files changed,
so no user-visible change); `Rscript tools/dod-gate.R` **PASSED** — check
0 errors / 0 warnings / 0 notes, coverage 100%, style and lint clean. Re-run
after the final on-branch edit so the evidence matches the merged tree.
