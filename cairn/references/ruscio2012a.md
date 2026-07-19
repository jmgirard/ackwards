# ruscio2012a — the Comparison Data (CD) method offered by `suggest_k()`

**Full citation.** Ruscio, J., & Roche, B. (2012). Determining the number of
factors to retain in an exploratory factor analysis using comparison data of
known factorial structure. *Psychological Assessment, 24*(2), 282–292.
https://doi.org/10.1037/a0025697

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/ruscio2012a.pdf` (gitignored). Pagination: journal
pages (282–292).
Extraction: verified 2026-07-19 (M69) — the citation and the
method's generate-and-increase logic read directly against the source (p. 282
abstract; pp. 285–286 method), and the RMSR terminology count taken from the
full-text extraction.

## Why this is a primary source

The source of the **Comparison Data** number-of-factors method `suggest_k()`
offers as `"cd"` (delegated to `EFAtools::CD()`). The roxygen and the vignette
describe the algorithm and cite Ruscio & Roche (2012).

## Extracted values

Method summary, verbatim from the abstract (p. 282): "Comparison data (CD) with
known factorial structure are first generated using 1 factor, and then the
number of factors is increased until the reproduction of the observed
eigenvalues fails to improve significantly." The decision statistic is the
**root-mean-square residual (RMSR)** between observed and reproduced eigenvalue
profiles, with a nonparametric significance test at each step (pp. 285–286).

The package's roxygen description — "Resamples from the observed item
distributions to generate comparison eigenvalue profiles; retains factors until
adding one no longer improves RMSE beyond chance" — is faithful in substance.
Two wording notes:

- **RMSR vs. RMSE.** The source's term is **RMSR** (root-mean-square residual;
  7 occurrences, never "RMSE"). The roxygen says "RMSE", following the term
  `EFAtools::CD()` uses in the wrapper the package actually calls. Same
  quantity, different label.
- **Title wording.** The roxygen `@references` entry reads "…comparison data
  **of a** known factorial structure"; the paper's printed title is "…Comparison
  Data **of Known** Factorial Structure" (no "a").

## Traces to

- `R/suggest_k.R:24-29` — the `"cd"` method roxygen description.
- `R/suggest_k.R:156-158` — the `@references` citation entry.
- `R/suggest_k.R:400` — the CD implementation branch (`# Comparison Data
  (Ruscio & Roche 2012)`).
- `vignettes/ackwards-suggest-k.Rmd.orig:171`, `:518` — the vignette prose and
  reference list.

## Open questions

- Both wording notes above (RMSR→RMSE, "of a known"→"of Known") are user-visible
  (rendered `@references`/roxygen); recorded in the M69 drift ledger as
  low-severity. Per the M69 gate, any correction routes to `/hotfix` (owner's
  call — the RMSE term matches the delegated `EFAtools` wrapper), not corrected
  in M69. — observed 2026-07-19.
