# M69: Author + verify the 8 secondary-methods source notes against their shelf PDFs

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —
- **Branch/PR:** `m69-reference-verification-secondary-sources`

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

- [ ] All eight source notes exist under `cairn/references/`, each authored from
      `source-note.md` with the full citation as printed and a `**Provenance.**`
      block whose `Extraction:` line records a dated M69 verification against the
      shelf PDF; the specific value/attribution each backs is quoted verbatim
      with its page/table anchor.
- [ ] Each page's `Traces to` names the exact code sites it backs (the file:line
      map in the Tasks below), and each backing claim is confirmed against the
      source or its discrepancy recorded in-page.
- [ ] The two user-visible constant sets are each checked against their source
      with a page anchor and confirmed-or-recorded: Hu & Bentler CFI/TLI `.95`,
      RMSEA `.06`, SRMR `.08` (`R/tidy.R:283-286`) against `hu1999`; the Kaiser
      KMO bands `.90/.80/.70/.60/.50` (`R/factorability.R:185`) against `kaiser1974`.
- [ ] `INDEX.md` carries one filename-first line per new page; `cairn_validate`
      exits 0 with its references check clean and the `references staleness`
      advisory showing the eight new pages as verified (not first-pass), with no
      regression on the existing 12.
- [ ] Every discrepancy found is enumerated in the work log with its routing
      disposition (hotfix / trivial / none); no user-facing text is corrected
      inside M69.
- [ ] Diff is `cairn/`-only (no `R/`, `tests/`, `man/`, `NAMESPACE`, or
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
- [ ] **T4** — Author + verify `xia2019` and `zhang2020`, the two fit-index
      caveat backers. Confirm the "CFI/TLI/RMSEA badly optimistic under WLSMV"
      claim and the "approximate under the two-step FIML→EFA route regardless of
      N" claim. `Traces to`: `R/engine_esem.R:269`, `R/tidy.R:44`;
      `R/ackwards.R:58`, `R/ackwards.R:161`, `cairn/DECISIONS.md:151`.
- [ ] **T5** — Author + verify `forbes2021` — dataset provenance for the
      N = 3,175 Australian general-population sample. Confirm the N and sample
      description. `Traces to`: `R/data.R:139`, `R/data.R:157`,
      `vignettes/ackwards-forbes2023.Rmd.orig:24`.
- [ ] **T6** — Add eight filename-first lines to `INDEX.md` (a new
      "Secondary methods & diagnostics backers" section is fine), and cross-link
      related pages with `[[…]]` where natural (e.g. `hu1999` ↔ the ESEM notes).
- [ ] **T7** — Drift ledger: enumerate in the work log every discrepancy found
      across T1–T5 and its routing disposition (hotfix / trivial / none). Make no
      user-facing correction here.
- [ ] **T8** — Gate: `cairn_validate` exit 0 with references check clean and the
      staleness advisory showing the eight new pages verified; confirm the diff is
      `cairn/`-only; run the DoD gate as a code-untouched safety check.

## Work log

- 2026-07-19: created by /milestone-plan. Shape set at the question gate — pure
  M67/M68 docs-only mold (drift routes to /hotfix), documentary-only (no guard tests).
- 2026-07-19: T1 — hu1999 + kaiser1974 authored, verified against rendered source pages. hu1999 p.27: CFI/TLI .95, RMSEA .06, SRMR .08 all match `.fit_cutoffs()` (`R/tidy.R:283-286`); noted the paper's two-index rule pairs .95 with SRMR .09. kaiser1974 p.35: six KMO bands + cutoffs match `.kmo_band()`; DRIFT — code prints "marvellous" (Commonwealth) vs Kaiser's "marvelous" (user-visible label → drift ledger, /hotfix owner-call).
- 2026-07-19: T2 — ruscio2012a + revelle1979 authored, verified. ruscio p.282 abstract + pp.285-286: CD generate-and-increase logic faithful; two user-visible wording notes — source term is RMSR (7×, never RMSE) vs roxygen "RMSE" (EFAtools' term); roxygen title "of a known" vs paper "of Known" (→ drift ledger, low-severity). revelle p.403/405-406: VSS complexity-1/2 matches roxygen + DESIGN:378, no drift. DOIs not printed on hu1999/kaiser1974/revelle scans — marked registered-not-source-verified (ruscio DOI was printed).
- 2026-07-19: T3 — maccallum1999 authored, verified against p.84 abstract. The "required N depends on communalities and overdetermination, not a fixed ratio" claim (roxygen + printout footnote) is faithful; citation matches `R/factorability.R:52-53`. No drift.

## Decisions

## Review
