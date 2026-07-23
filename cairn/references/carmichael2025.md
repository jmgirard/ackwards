# carmichael2025 — HiTOP-TBI: Forbes' extended bass-ackwards applied to traumatic brain injury (citation precedent)

**Provenance.** Ingested 2026-07-23 by M71 from
`cairn/references/sources/carmichael2025.pdf` (gitignored). Pagination: journal
pages (714–730); the masthead prints the locus as "42:714–730 (April 2025)".
Extraction: verified 2026-07-23 against the rendered source — the title/abstract (p. 714) and the "Extended bass-ackwards modeling" methods section (p. 721) were read directly; DOI printed on the source — observed 2026-07-23.

**Citation.** Carmichael, J., Ponsford, J., Gould, K. R., Tiego, J., Forbes,
M. K., Kotov, R., Fornito, A., & Spitz, G. (2025). A transdiagnostic,
hierarchical taxonomy of psychopathology following traumatic brain injury
(HiTOP-TBI). *Journal of Neurotrauma, 42*, 714–730.
https://doi.org/10.1089/neu.2024.0006

**Role.** Manuscript/vignette **citation precedent**: an applied use of the
bass-ackwards method in a new clinical domain (psychopathology after
moderate-severe traumatic brain injury), nothing in `R/` traces here. Notable
because it uses **Forbes' extended bass-ackwards** ([[forbes2023]]) — with
Miriam K. Forbes herself as a co-author — so it doubles as an in-the-wild
witness that the forbes2023 redundancy contract is applied as `ackwards`
implements it: redundant factors removed at `r ≥ |0.90|` conjoined with Tucker
congruence `r_c > 0.95` (the `redundancy_r` = .9 + φ > .95 conjunction of
D-017). Method is **bass-ackwards** (top-down PCA-style hierarchy on the HiTOP
scale set), not an adjacent method.

## Extracted values

- The approach — "Using a top-down, exploratory latent variable approach (bass-ackwards modeling), we subsequently constructed a hierarchical model of psychopathological dimensions tailored to TBI", abstract, p. 714.
- Construction input — "Bass-ackwards analysis was applied to construct a hierarchical model from the 57 original and customized psychopathology scales", p. 721.
- Stopping choice — "Both parallel and Velicer's minimum average partial analyses suggested extracting a maximum of 10 factors. However, the 8-, 9-, and 10-factor models each included a factor that was not well-defined, possessing fewer than three primary indicators. Thus, we present the hierarchical model up to seven factors", p. 721.
- Redundancy pruning (the forbes2023 contract in use) — "Following Forbes' extended bass-ackwards approach, we removed strings of statistically redundant factors (i.e., r ≥ |0.90| and r_c > 0.95)", p. 721.
- Sample — 410 individuals with moderate-severe TBI (abstract, p. 714).
- Result — the HiTOP-TBI model "encompassed broad internalizing and externalizing spectra, splitting into seven narrower dimensions" (abstract, p. 714).

## Traces to

Nothing in the repo reads this page yet. It is a **citation precedent** for the
BRM manuscript / vignettes (applied bass-ackwards in a TBI sample) and an
independent corroboration of the D-017 redundancy criteria; it backs no code
value or oracle.

## Open questions

- The full 57-scale loading tables and the level-by-level edge correlations were not transcribed — the headline method, stopping choice, and redundancy criteria (above) are recorded; the per-level detail lives in the paper's Table/Figure set (Fig. 3, pp. 722–723) and Supplementary materials — observed 2026-07-23.
