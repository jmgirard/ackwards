# kaiser1974 — the KMO sampling-adequacy band labels printed by `factorability()`

**Full citation.** Kaiser, H. F. (1974). An index of factorial simplicity.
*Psychometrika, 39*(1), 31–36. Registered DOI 10.1007/BF02291575 *(not printed
on the scan — vol/issue/pages are from the source header; the DOI is the
publisher-registered identifier, not source-verified)*.

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/kaiser1974.pdf` (gitignored). Pagination: journal
pages (31–36). Extraction: verified 2026-07-19 (M69) — the six-band evaluation
and its cutoffs were read directly against the source (p. 35, rendered) and the
citation header against p. 31.

## Why this is a primary source

The source of the **KMO band labels** that `factorability()` prints for the
Kaiser-Meyer-Olkin overall and per-item measures of sampling adequacy. Kaiser
introduces the index as the *Index of Factorial Simplicity* (IFS); the same
statistic is now universally called the KMO measure.

## Extracted values

The verbatim six-band evaluation of the index (p. 35), which Kaiser prefaces
"Subjective reflection, based upon Table 1 and primarily observing IFSs for a
substantial number of factor analyses from the real world, suggests the
following evaluation of levels of our index of factorial simplicity:"

- "in the .90s, marvelous" — p. 35.
- "in the .80s, meritorious" — p. 35.
- "in the .70s, middling" — p. 35.
- "in the .60s, mediocre" — p. 35.
- "in the .50s, miserable" — p. 35.
- "below .50, unacceptable" — p. 35.

The cutoffs (`.90 / .80 / .70 / .60 / .50`) match the code's `.kmo_band()`
thresholds exactly. The labels match the printed labels with **one spelling
variant**: the code prints "marvellous" (Commonwealth double-l) where Kaiser
prints "marvelous" (single-l). Kaiser frames the bands as "subjective
reflection", i.e. a rule of thumb — exactly how `factorability()` presents them.

## Traces to

- `R/factorability.R:178` — the `.kmo_band()` comment attributing the bands to
  Kaiser 1974.
- `R/factorability.R:185-197` — the band cutoffs and labels themselves.
- `R/factorability.R:349` — the printout footnote ("KMO bands (Kaiser 1974) …").
- `R/factorability.R:9`, `R/factorability.R:21` — roxygen naming the KMO measure
  and framing every cutoff as a convention, not a verdict.

## Open questions

- The "marvellous" vs. "marvelous" spelling deviation is a user-visible printed
  label; recorded in the M69 drift ledger. Per the M69 gate, user-visible
  propagation routes to `/hotfix` (owner's call whether the Commonwealth
  spelling is intentional house style) — not corrected in M69. — observed
  2026-07-19.
