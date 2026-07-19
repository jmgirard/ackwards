# zhang2020 — why FIML-path EFA fit indices are flagged as approximate

**Full citation.** Zhang, X., & Savalei, V. (2020). Examining the effect of
missing data on RMSEA and CFI under normal theory full-information maximum
likelihood. *Structural Equation Modeling: A Multidisciplinary Journal, 27*(2),
219–239. https://doi.org/10.1080/10705511.2019.1642111 *(DOI printed on source
p. 1; published online 5 Sep 2019)*.

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/zhang2020.pdf` (gitignored). Pagination: journal pages
(219–239).
Extraction: verified 2026-07-19 (M69) — citation and the core finding
read directly against the source (p. 219 abstract).

## Why this is a primary source

The source cited for the caveat that EFA fit indices ([[hu1999]] RMSEA/CFI)
computed on the FIML missing-data path (`missing = "fiml"`) are *approximate* — the value the package
emits via cli and the D-020 rationale in DECISIONS.

## Extracted values

Abstract finding, verbatim (p. 219): "it is not commonly known that approximate
fit indices (AFIs) can be distorted, relative to their complete data
counterparts, when FIML is used to handle missing data. In this article, we show
that [the] two most popular AFIs, the root-mean-square error of approximation
(RMSEA) and the comparative fit index (CFI), often approach different population
values under FIML estimation when missing data are present."

- The distortion is at the **population** level (the AFIs "approach different
  population values"), i.e. a bias in the estimand — so it does **not** shrink
  with sample size. This is exactly what the package's "approximate … regardless
  of N (Zhang & Savalei, 2020)" claim asserts. Faithful.

**Scope note (standing fact).** Zhang & Savalei study FIML fit indices *within
SEM directly*. The package applies the citation to its **two-step**
`corFiml()` → normal-theory EFA route (a FIML correlation matrix fed into a
normal-theory EFA), which carries its own additional layer of approximation. The
citation correctly backs the general point — FIML-based RMSEA/CFI are
distorted/approximate under missing data, independent of N — which is what
motivates flagging the two-step EFA indices. No overreach beyond that.

## Traces to

- `R/ackwards.R:54-59` — the roxygen "those [EFA fit indices] are *approximate*
  … regardless of N (Zhang & Savalei, 2020)".
- `R/ackwards.R:161` — the `@references` citation entry.
- `cairn/DECISIONS.md:151` (D-020 consequences) — "EFA fit indices are
  approximate under the two-step route regardless of N (Zhang & Savalei 2020)".

## Open questions

- No code drift found — the "approximate regardless of N" attribution is faithful
  to the population-level distortion the paper derives. Nothing routed to
  `/hotfix`. — observed 2026-07-19.
