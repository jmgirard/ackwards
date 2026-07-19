# M69: Author + verify the 8 secondary-methods source notes against their shelf PDFs

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —
- **Branch/PR:** `m69-reference-verification-secondary-sources` · https://github.com/jmgirard/ackwards/pull/73

## Goal

Author a verified `cairn/references/` source note for each of the eight
newly-added secondary-methods PDFs, pinning the value or attribution each backs
in the code against its source and recording any drift for downstream routing.

## Scope

**In:** Eight new source notes — `hu1999`, `kaiser1974`, `ruscio2012a`,
`maccallum1999`, `revelle1979`, `xia2019`, `zhang2020`, `forbes2021` — each
authored from `templates/source-note.md`, read against its shelf PDF, carrying a
full citation, a `**Provenance.**` block with a dated M69 verification, the
specific backing value/attribution quoted with a page/table anchor, and a
`Traces to` list naming the exact code sites. Eight new `INDEX.md` lines.
Discrepancies (bad attribution, wrong page anchor, wrong user-visible number)
recorded **in the page** and enumerated in the work log with a routing
disposition.

**Out:**
- User-visible propagation of any drift found (a wrong printed number, roxygen,
  or vignette text) → `/hotfix`, with the regression test there. No user-facing
  text is corrected inside M69 (gate decision, M67/M68 precedent).
- Internal code-comment corrections → a trivial post-merge follow-up commit,
  not this milestone.
- Executable guard tests asserting the constants match the source → declined
  this milestone (documentary-only, gate decision); a `candidate` if wanted later.
- Re-verifying the existing 12 pages — M67/M68 already did that.

## Acceptance criteria

- [x] All eight source notes exist under `cairn/references/`, each authored from
      `source-note.md` with the full citation as printed and a `**Provenance.**`
      block whose `Extraction:` line records a dated M69 verification against the
      shelf PDF; the specific value/attribution each backs is quoted verbatim
      with its page/table anchor.
- [x] Each page's `Traces to` names the exact code sites it backs (the file:line
      map in the Tasks below), and each backing claim is confirmed against the
      source or its discrepancy recorded in-page.
- [x] The two user-visible constant sets are each checked against their source
      with a page anchor and confirmed-or-recorded: Hu & Bentler CFI/TLI `.95`,
      RMSEA `.06`, SRMR `.08` (`R/tidy.R:283-286`) against `hu1999`; the Kaiser
      KMO bands `.90/.80/.70/.60/.50` (`R/factorability.R:185`) against `kaiser1974`.
- [x] `INDEX.md` carries one filename-first line per new page; `cairn_validate`
      exits 0 with its references check clean and the `references staleness`
      advisory showing the eight new pages as verified (not first-pass), with no
      regression on the existing 12.
- [x] Every discrepancy found is enumerated in the work log with its routing
      disposition (hotfix / trivial / none); no user-facing text is corrected
      inside M69.
- [x] Diff is `cairn/`-only (no `R/`, `tests/`, `man/`, `NAMESPACE`, or
      `NEWS.md` change); the DoD gate stays green (no code touched).

## Coverage

- AC1 → T1, T2, T3, T4, T5
- AC2 → T1, T2, T3, T4, T5
- AC3 → T1
- AC4 → T6, T8
- AC5 → T7
- AC6 → T8

## Tasks

- [x] **T1** — Author + verify `hu1999` and `kaiser1974`, the two user-visible
      constant backers. Read the fit-index cutoffs against Hu & Bentler (1999)
      and quote each verbatim with a page anchor (CFI/TLI `.95`, RMSEA `.06`,
      SRMR `.08`); read the KMO band labels against Kaiser (1974) (`.90/.80/.70/
      .60/.50`). `Traces to`: `R/tidy.R:283-286`, `R/tidy.R:275`,
      `R/autoplot.R:720`, `R/autoplot.R:792`; `R/factorability.R:178-193`,
      `R/factorability.R:349`.
- [x] **T2** — Author + verify `ruscio2012a` and `revelle1979`, the two
      `suggest_k` method backers. Confirm the Comparison Data algorithm as we
      describe it and the VSS-1/VSS-2 criterion. `Traces to`: `R/suggest_k.R:400`,
      `R/suggest_k.R:156`; `R/suggest_k.R:20`, `R/suggest_k.R:366`,
      `cairn/DESIGN.md:378`.
- [x] **T3** — Author + verify `maccallum1999` — the sample-size-depends-on-
      communalities-not-a-fixed-ratio claim printed by `factorability()`. Confirm
      the claim wording against the source. `Traces to`: `R/factorability.R:26`,
      `R/factorability.R:349`.
- [x] **T4** — Author + verify `xia2019` and `zhang2020`, the two fit-index
      caveat backers. Confirm the "CFI/TLI/RMSEA badly optimistic under WLSMV"
      claim and the "approximate under the two-step FIML→EFA route regardless of
      N" claim. `Traces to`: `R/engine_esem.R:269`, `R/tidy.R:44`;
      `R/ackwards.R:58`, `R/ackwards.R:161`, `cairn/DECISIONS.md:151`.
- [x] **T5** — Author + verify `forbes2021` — dataset provenance for the
      N = 3,175 Australian general-population sample. Confirm the N and sample
      description. `Traces to`: `R/data.R:139`, `R/data.R:157`,
      `vignettes/ackwards-forbes2023.Rmd.orig:24`.
- [x] **T6** — Add eight filename-first lines to `INDEX.md` (a new
      "Secondary methods & diagnostics backers" section is fine), and cross-link
      related pages with `[[…]]` where natural (e.g. `hu1999` ↔ the ESEM notes).
- [x] **T7** — Drift ledger: enumerate in the work log every discrepancy found
      across T1–T5 and its routing disposition (hotfix / trivial / none). Make no
      user-facing correction here.
- [x] **T8** — Gate: `cairn_validate` exit 0 with references check clean and the
      staleness advisory showing the eight new pages verified; confirm the diff is
      `cairn/`-only; run the DoD gate as a code-untouched safety check.

## Work log

- 2026-07-19: created by /milestone-plan. Shape set at the question gate — pure M67/M68 docs-only mold (drift routes to /hotfix), documentary-only (no guard tests).
- 2026-07-19: T1 — hu1999 + kaiser1974 authored, verified against rendered source pages. hu1999 p.27: CFI/TLI .95, RMSEA .06, SRMR .08 all match `.fit_cutoffs()` (`R/tidy.R:283-286`); noted the paper's two-index rule pairs .95 with SRMR .09. kaiser1974 p.35: six KMO bands + cutoffs match `.kmo_band()`; DRIFT — code prints "marvellous" (Commonwealth) vs Kaiser's "marvelous" (user-visible label → drift ledger, /hotfix owner-call).
- 2026-07-19: T2 — ruscio2012a + revelle1979 authored, verified. ruscio p.282 abstract + pp.285-286: CD generate-and-increase logic faithful; two user-visible wording notes — source term is RMSR (7×, never RMSE) vs roxygen "RMSE" (EFAtools' term); roxygen title "of a known" vs paper "of Known" (→ drift ledger, low-severity). revelle p.403/405-406: VSS complexity-1/2 matches roxygen + DESIGN:378, no drift. DOIs not printed on hu1999/kaiser1974/revelle scans — marked registered-not-source-verified (ruscio DOI was printed).
- 2026-07-19: T3 — maccallum1999 authored, verified against p.84 abstract. The "required N depends on communalities and overdetermination, not a fixed ratio" claim (roxygen + printout footnote) is faithful; citation matches `R/factorability.R:52-53`. No drift.
- 2026-07-19: T4 — xia2019 + zhang2020 authored, verified against abstracts (both source-printed DOIs). xia p.409: DWLS/ULS optimistic vs ML "regardless of scaling" — supports the direction of the code caveat, but the code's naive-vs-scaled framing (engine_esem:269 comment, tidy.R:44 roxygen) is a nuance note, not clearly false → ledger low-severity. zhang p.219: FIML AFIs distorted at population level → "approximate regardless of N" faithful; scope-note added (they study FIML-in-SEM, code applies to its two-step route). Fixed the D-entry ref (D-020, not D-016).
- 2026-07-19: T5 — forbes2021 authored, verified against p.139 header+abstract (source-printed DOI). N=3,175 Australian general-population + 18 DSM disorders match `R/data.R:139-140`; noted the US clinical n=1,775 sample is correctly not the one shipped. No drift.
- 2026-07-19: T6 — new INDEX.md section "Secondary methods & diagnostics backers" (8 filename-first lines); cross-links [[hu1999]] added to xia2019/zhang2020, [[forbes2023]] in forbes2021. cairn_validate exit 0: references index<->disk PASS, staleness OK (all 8 verified). Fixed the M60/M68 trap — `Extraction:` must begin its own physical line (7 of 8 initially reported "no extraction status").
- 2026-07-19: T7 — drift ledger consolidated in Decisions below (4 discrepancies, all wording/attribution; no wrong user-visible numbers). Per the M69 gate none corrected here — routed to /hotfix (owner's call) or trivial follow-up.
- 2026-07-19: T8 — diff confirmed cairn/-only (R package byte-unchanged → last green DoD carries); cairn_validate exit 0. Status → review.
- 2026-07-19: review — PR #73; DoD gate re-run PASSED (check 0/0/0, coverage 100%, style/lint/pkgdown clean). 3 independent lenses: blame-history + diff-bug clean; prior-PR-comments caught 2 quote-accuracy defects in revelle1979.md (fabricated "[as the estimate]" ending p.403; dropped leading "The optimal" p.406) — both confirmed by rendering the source pages and FIXED. Diff-bug lens false-negatived these (trusted flattened text) — recorded. All 6 ACs verified.

## Decisions

**M69 drift ledger** (2026-07-19). Four discrepancies found; all wording/
attribution, none a wrong user-visible number. Per the M69 gate, none corrected
in M69 — recorded here and in each source page, routed out:

- **kaiser1974** — code prints `"marvellous"` where Kaiser (p. 35) prints
  `"marvelous"`. User-visible band label → `/hotfix` (owner's call; possibly
  intentional Commonwealth spelling).
- **ruscio2012a** — roxygen says `RMSE` where the source says **RMSR** (7×,
  never RMSE); the delegated `EFAtools::CD()` wrapper uses "RMSE". User-visible
  roxygen → `/hotfix`, low-severity.
- **ruscio2012a** — `@references` title reads "…comparison data **of a** known
  factorial structure"; the paper's title is "…of Known Factorial Structure".
  User-visible → `/hotfix`, low-severity.
- **xia2019** — the naive-vs-scaled framing (`R/engine_esem.R:269` internal
  comment; `R/tidy.R:44` roxygen) is a nuance, not false: Xia & Yang's optimism
  is DWLS/ULS-vs-ML "regardless of scaling", so "only the scaled are defensible"
  is a lavaan-reporting choice they do not endorse as sufficient. Internal
  comment → trivial follow-up if reworded; roxygen owner's call.

No drift: **hu1999, revelle1979, maccallum1999, zhang2020, forbes2021** — the
code values/claims match their sources.

## Review

**Reviewed 2026-07-19. PR #73.**

### Acceptance-criterion evidence

- **AC1** (8 notes, citation + dated-M69 Provenance + quoted value w/ anchor):
  all 8 files present under `cairn/references/`; `cairn_validate` `references
  staleness` = OK (all 8 read as `verified 2026-07-19 (M69)`, none first-pass).
  Each carries a full citation and a page/table-anchored value.
- **AC2** (Traces to + confirmed-or-recorded): each page has a `Traces to` list
  of concrete file:lines; drift recorded in-page for the 4 discrepant sources,
  "no drift" stated for the other 4.
- **AC3** (user-visible constants vs source): Hu & Bentler p. 27 rendered —
  CFI/TLI .95, SRMR .08, RMSEA .06 match `.fit_cutoffs()` (`R/tidy.R:283-286`);
  Kaiser p. 35 rendered — six KMO bands + cutoffs match `.kmo_band()`
  (`R/factorability.R:185-197`), one spelling variant recorded.
- **AC4** (INDEX + validate): INDEX.md gained 8 filename-first lines in a new
  "Secondary methods & diagnostics backers" section; `cairn_validate` exit 0
  (`references index<->disk` PASS, `references staleness` OK, `weight caps` PASS).
- **AC5** (drift ledger + no in-milestone fix): 4 discrepancies enumerated in
  the work log and the Decisions drift ledger with routing dispositions; diff
  touches no user-facing text (proven cairn/-only, below).
- **AC6** (cairn/-only + DoD green): `git diff --name-only master...HEAD` = all
  under `cairn/` (R package byte-unchanged). DoD-gate PASSED — check 0 err/0
  warn/0 note, coverage 100%, style/lint clean, pkgdown index complete.

### Consistency gate

- `cairn_validate` exit 0; all structural checks PASS. Advisories only: 81
  pre-existing dangling-id tokens (entombed milestone refs, not from this diff).
- No principle change (Principles touched: —) → `cairn_impact` skipped.

### Independent review

Three fresh-context lenses (ref-based git; shared tree). Two of the three cleared
cleanly; the prior-PR lens caught two real quote-accuracy defects, both fixed.

- **[O] diff-bug (Opus):** no findings — independently re-verified all 8 page
  anchors, DOIs, code traces, drift ledger, and INDEX against the shelf PDFs.
  (Note: it cleared revelle's quotes as "verbatim" — a false negative; it trusted
  flattened text and missed the two defects the prior-PR lens caught. Recorded
  honestly, not silently.)
- **[S] blame-history (Sonnet):** no findings — no contradiction of
  DECISIONS/DESIGN, no M67/M68 lesson repeated, no prior work undone, INDEX
  conventions correct.
- **[S] prior-PR-comments (Sonnet):** 2 findings on `revelle1979.md`, both
  regressing the M67 "unmarked-elision / a quotation is a new claim" lesson.

**Findings triage** (both adjudicated by rendering the source page directly —
primary-source ground truth, stronger than a confidence score; both ≫80, fixed):

1. `revelle1979.md` abstract quote (was) "…is taken [as the estimate]." — p. 403
   (rendered) actually reads "…is taken **as being the optimal number of factors
   to extract**." Fabricated bracketed ending under a "verbatim" label. **Fixed.**
2. `revelle1979.md` complexity quote dropped the leading "The optimal" and began
   mid-sentence with lowercase "the number", unmarked. p. 406 (rendered) reads
   "**The optimal** number of interpretable factors (of complexity v) is the
   number of factors, k, which maximizes VSSvk." **Fixed** (re-anchored p. 406).

No sub-80 findings to log. Both fixes re-verified against the rendered PDF pages.
