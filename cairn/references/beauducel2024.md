# beauducel2024 — the determinacy vs inter-factor-correlation trade-off (backs the tenBerge / correlation-preserving default)

**Provenance.** Ingested 2026-07-23 by M70 from
`cairn/references/sources/beauducel2024.pdf` (gitignored).
Pagination: journal pages (289–313).
Extraction: verified 2026-07-23 against the rendered source — the title and abstract (p. 289) were read directly; author name and DOI cross-checked against Crossref — observed 2026-07-23.

**Citation.** Beauducel, A., Hilger, N., & Kuhl, T. (2024). The trade-off between
factor score determinacy and the preservation of inter-factor correlations.
*Educational and Psychological Measurement, 84*(2), 289–313.
DOI `10.1177/00131644231171137` (printed on source). The copyright line reads
"© The Author(s) 2023" (online-ahead 2023; issue dated 2024). Third author is
recorded as "Kuhl" (not "Kühl") on both the source byline and in Crossref.

**Role.** The sharpest current formalization of the trade-off `ackwards`' default
scoring sits inside. `ackwards` uses **tenBerge correlation-preserving** factor
scores (DESIGN §9 `scores` row; D-007) precisely because the between-level
correlations are the method's whole output — but correlation-preserving scores
buy that unbiased inter-correlation at the cost of *determinacy* relative to
regression scores. This paper states that trade-off exactly, and so also backs
the `redundancy_phi` rationale (factor-score indeterminacy is why EFA/ESEM adds a
congruence guard that determinate PCA components do not need).

## Extracted values

- Regression scores maximize determinacy but bias correlations — "Regression factor score predictors have the maximum factor score determinacy, that is, the maximum correlation with the corresponding factor, but they do not have the same inter-correlations as the factors", abstract, p. 289.
- The cost of correlation-preservation — "correlation-preserving factor score predictors have smaller correlations with the corresponding factors (factor score determinacy) than regression factor score predictors", abstract, p. 289.
- The trade-off stated — "higher factor score determinacy goes along with bias of the inter-correlations and unbiased inter-correlations go along with lower factor score determinacy", abstract, p. 289.

## Traces to

- `cairn/DESIGN.md` §9 `scores (method)` and `redundancy_phi` rows — the tenBerge default and the indeterminacy-driven φ-filter rationale (M70 adds this citation beside Lorenzo-Seva & ten Berge 2006).
- `R/ackwards.R` `@references` — the scoring citation block (M70).

## Open questions

- The note relies on the abstract; the Appendix syntax (for computing correlation-preserving predictors from regression predictors) and the simulation conditions were not transcribed — not needed for the default's rationale — observed 2026-07-23.
