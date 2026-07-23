# grice2001 — factor-score indeterminacy; why the scoring method is a deliberate choice (backs the tenBerge scoring default)

**Provenance.** Ingested 2026-07-23 by M70 from
`cairn/references/sources/grice2001.pdf` (gitignored).
Pagination: journal pages (430–450).
Extraction: verified 2026-07-23 against the rendered source — the title, abstract, and the three-researchers illustration (p. 430) were read directly — observed 2026-07-23.

**Citation.** Grice, J. W. (2001). Computing and evaluating factor scores.
*Psychological Methods, 6*(4), 430–450. DOI `10.1037/1082-989X.6.4.430`
(printed on source).

**Role.** The reference treatment of **factor-score indeterminacy** — a
common-factor model does not determine a unique set of factor scores, so
different scoring methods yield different (sometimes contradictory) scores. This
is the backdrop for `ackwards`' `tenBerge` scoring default (DESIGN §9 `scores`
row; D-007) and for the `redundancy_phi` auto-resolve, whose rationale turns on
EFA/ESEM factor-score *indeterminacy* vs. PCA component *determinacy*
(DESIGN §9 `redundancy_phi` row).

## Extracted values

- Scope — "This article reviews the history and nature of factor score indeterminacy", abstract, p. 430.
- Recommendation — "factor score indeterminacy should be routinely assessed and reported as part of any exploratory factor analysis and that factor scores should be thoroughly evaluated before they are reported or used in subsequent statistical analyses", abstract, p. 430.
- The consequence the method choice carries — in the three-researchers illustration, "All three researchers will likely arrive at a different set of factor scores, possibly yielding widely discrepant rankings of the individuals along the extracted factors", p. 430.

## Traces to

- `cairn/DESIGN.md` §9 `scores (method)` and `redundancy_phi` rows — the tenBerge default and the PCA-determinate vs EFA-indeterminate φ-filter rationale (M70 adds this citation).
- `R/ackwards.R` `@references` — the scoring citation block (M70).

## Open questions

- The note relies on the abstract + p. 430 framing; Grice's specific indeterminacy indices and program outputs (interior pages) were not transcribed — not needed for the default's rationale — observed 2026-07-23.
