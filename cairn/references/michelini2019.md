# michelini2019 — bass-ackwards on child + adult psychopathology in the ABCD study (citation precedent)

**Provenance.** Ingested 2026-07-23 by M71 from
`cairn/references/sources/michelini2019.pdf` (gitignored). *Translational
Psychiatry* is an article-number journal — the article carries no sequential
journal page numbers (article 9:261), so anchors below are the PDF/manuscript
pages the journal itself prints ("Page N of 15" in the running head).
Extraction: verified 2026-07-23 against the rendered source — the title/abstract (PDF p. 1) and the "Statistical analysis" methods column (PDF p. 3) were read directly; DOI printed on the source — observed 2026-07-23.

**Citation.** Michelini, G., Barch, D. M., Tian, Y., Watson, D., Klein, D. N.,
& Kotov, R. (2019). Delineating and validating higher-order dimensions of
psychopathology in the Adolescent Brain Cognitive Development (ABCD) study.
*Translational Psychiatry, 9*, 261. https://doi.org/10.1038/s41398-019-0593-4
Open Access, CC-BY 4.0.

**Role.** Manuscript/vignette **citation precedent**: bass-ackwards applied to
child/adolescent psychopathology (N = 9987 nine- and ten-year-olds plus their
parents, ABCD study), yielding a general 'p' factor over five specific
dimensions. Method is **Goldberg's bass-ackwards** run on **EFA** factor
scores (PCA extraction, geomin rotation) — an EFA-engine precedent (like
[[cowan2024]]), *not* Goldberg's original PCA-with-component-scores recipe, and
its rotation is **oblique geomin**, differing from `ackwards`' orthogonal
default (D-002). Its adjacent-level edge threshold (≥ 0.65) is a concrete
published cut worth citing beside `prune()`'s `redundancy_r`. States the
strongest available "why bass-ackwards" rationale for the EFA case.

## Extracted values

- Extraction + rotation — "we used exploratory factor analysis (EFA) to empirically extract (with principal component analysis) and rotate (with geomin) factor solutions with an increasing number of factors", PDF p. 3.
- The bass-ackwards step — "we correlated factor scores on adjacent levels of the hierarchy to describe transitions between levels using Goldberg's bass-ackwards hierarchical method", PDF p. 3.
- Edge threshold — "The paths between levels in the hierarchical model reflect correlations ≥0.65 between the factor scores", PDF p. 3.
- Why bass-ackwards — "The bass-ackwards approach was chosen to be consistent with previous studies … and because, to our knowledge, it is the only method that allows for the delineation of multiple hierarchical levels from factors derived through EFA", PDF p. 3.
- Stopping — "The maximum number of factors to extract was determined with parallel analyses (extraction was stopped when eigenvalues fell within the 95% confidence interval of eigenvalues from simulated data)", PDF p. 3.
- Result — "A hierarchical structure with a general psychopathology ('p') factor at the apex and five specific factors (internalizing, somatoform, detachment, neurodevelopmental, and externalizing) emerged in children" (abstract, PDF p. 1).

## Traces to

Nothing in the repo reads this page yet. It is a **citation precedent** for the
BRM manuscript / vignettes (bass-ackwards in a large child+adult sample; the
EFA-engine and oblique-rotation variant of the method); it backs no code value
or oracle.

## Open questions

- The validator-correlation and familial-aggregation tables (the paper's incremental-validity results) were not transcribed — the method and edge-threshold facts above are what the repo relies on — observed 2026-07-23.
