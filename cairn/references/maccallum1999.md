# maccallum1999 — the "required N depends on communalities, not a fixed ratio" claim in `factorability()`

**Full citation.** MacCallum, R. C., Widaman, K. F., Zhang, S., & Hong, S.
(1999). Sample size in factor analysis. *Psychological Methods, 4*(1), 84–99.
Registered DOI 10.1037/1082-989X.4.1.84 *(not printed on the scan — vol/issue/
pages are from the source p. 84 header; the DOI is the publisher-registered
identifier, not source-verified)*.

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/maccallum1999.pdf` (gitignored). Pagination: journal
pages (84–99). Extraction: verified 2026-07-19 (M69) — citation and the
sample-size claim read directly against the source (p. 84 abstract).

## Why this is a primary source

The source `factorability()` cites when it warns the reader **not** to treat the
N:p rules of thumb as settled — required sample size depends on communalities
and factor overdetermination, not on any fixed ratio.

## Extracted values

Abstract (p. 84), verbatim: "A fundamental misconception about this issue is
that the minimum sample size, or the minimum ratio of sample size to the number
of variables, is invariant across studies. In fact, necessary sample size is
dependent on several aspects of any given study, including the level of
communality of the variables and the level of overdetermination of the factors."
The abstract closes: "Results demonstrate the lack of validity of common rules
of thumb…" (p. 84).

The package's two paraphrases both hold:

- roxygen (`R/factorability.R:26-28`): "the required `N` in particular depends
  on communalities and factor overdetermination far more than on any fixed
  ratio" — faithful.
- printout footnote (`R/factorability.R:349-352`): "required N depends on
  communalities and factor overdetermination (MacCallum et al. 1999)" —
  faithful.

Citation (authors, year, venue, 4(1), 84–99) matches the roxygen `@references`
entry at `R/factorability.R:52-53` exactly. No drift.

## Traces to

- `R/factorability.R:26-28` — the roxygen "required N depends on communalities…"
  caveat.
- `R/factorability.R:52-53` — the `@references` citation entry.
- `R/factorability.R:349-352` — the printout footnote citing MacCallum et al.
  1999.

## Open questions

- No code drift found. Nothing routed to `/hotfix` from this source. — observed
  2026-07-19.
