# M61: Enrich suggest_k() docs with verified k-selection citations

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —
- **Branch/PR:** m61-suggest-k-citations

## Goal

Source `suggest_k()`'s two currently-uncited claims — the consensus-range stance and
"PA-PC overextracts" — with the three verified-at-ingest citations (Lim & Jahng 2019,
Achim 2021, Saucier 1997), in both the roxygen help page and the suggest-k vignette.

## Scope

**In:** roxygen edits to `R/suggest_k.R` (two prose sentences + three `@references`
entries); the same citation pass on `vignettes/ackwards-suggest-k.Rmd.orig` (PA-PC
"typical behavior" passage, consensus-range statement, vignette References) with
regenerated precomputed `.Rmd`; a one-line NEWS entry in the 0.1.1 section.

**Out:** any behavior or argument change to `suggest_k()` (none wanted); citation
passes on other help pages/vignettes (add a candidate row if a gap is found);
`comparability()` citations (shipped in hotfix #61).

## Acceptance criteria

- [ ] AC1 — `?suggest_k` "Interpreting the output" states the ±1-advisory-range stance
      citing Lim & Jahng (2019) paired with Achim (2021) (the pair is mandatory per
      `rotation-and-k.md`), and "A note on overextraction" cites Saucier (1997, fn 14)
      alongside Forbes (2023); all three appear in `@references` with titles
      and DOIs matching the ingested notes (`rotation-and-k.md`, `saucier1997.md`).
- [ ] AC2 — the suggest-k vignette's PA-PC passage and consensus-range discussion carry
      the same citations; its References section gains the three entries; the
      precomputed `.Rmd` is regenerated from the `.orig` via `vignettes/precompute.R`
      (no stale-output ship).
- [ ] AC3 — NEWS.md 0.1.1 section records the enrichment in one line (no milestone
      numbers in user-facing text).
- [ ] AC4 — profile verify clean: `devtools::document()` run with `man/` committed and
      no residual diff; `TESTTHAT_CPUS=8 devtools::test()` all green.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T1, T2, T4

## Tasks

- [x] T1 — `R/suggest_k.R`: add the ±1-range sentence + Lim & Jahng/Achim pair to
      `@section Interpreting the output` (~line 35); add the Saucier fn-14 in-the-wild
      case to `@section A note on overextraction` (~line 118); add the three
      `@references` entries (alphabetical, `\doi{}` from the reference notes). Run
      `devtools::document()`; commit `man/` with it.
- [x] T2 — `vignettes/ackwards-suggest-k.Rmd.orig`: cite Saucier (1997) in the PA-PC
      "Typical behavior" passage (~line 76) and Lim & Jahng/Achim at the consensus-range
      statement (~line 61 or ~line 220); add the three entries to `## References`
      (~line 489). Re-run `Rscript vignettes/precompute.R`; commit the regenerated
      `.Rmd` (+ assets if changed).
- [x] T3 — NEWS.md: one-line bullet in the 0.1.1 section (extends the existing
      citation-lineage bullet's neighborhood).
- [x] T4 — `TESTTHAT_CPUS=8` suite run clean; spot-render `?suggest_k` to confirm the
      references block formats.

## Work log

- 2026-07-16: created by /milestone-plan (promotes the 2026-07-16 candidate row from
  the comparability-citations hotfix #61 follow-up; vignette + NEWS scope added at the
  plan gate).
- 2026-07-16: T1 done — ±1-range sentence (Lim & Jahng + Achim) in Interpreting-the-output,
  Saucier fn-14 case in overextraction note, 3 @references entries; document() clean,
  suite 2300 pass / 0 fail.
- 2026-07-16: T2 done — same citations in the vignette (consensus passage, PA-PC
  typical-behavior, References); precompute.R re-run; restored 4 unrelated vignettes +
  3 PNGs whose diffs were pure run noise (cli timings, gt element IDs).
- 2026-07-16: T3 done — NEWS 0.1.1 bullet after the comparability citation-lineage entry.
- 2026-07-16: T4 done — suite green at T1 (2300 pass, 0 fail; no R code changed since);
  Rd2txt spot-render confirms all three citations + DOI macros format. Status → review.

## Decisions

## Review
