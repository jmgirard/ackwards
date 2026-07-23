# forbes2025 — Forbes' DSM-5 reconstruction: extended bass-ackward continuation + fidelity drift-watch

**Provenance.** Ingested 2026-07-23 by M71 from
`cairn/references/sources/forbes2025.pdf` (gitignored). Pagination: journal
pages (462–488); PDF page N → journal p. 461+N (PDF p. 1 = p. 462).
Extraction: verified 2026-07-23 against the rendered source — the title/abstract (p. 462), the "Data analysis" methods (p. 467) with the analytic-pipeline figure (Fig. 1, p. 468), and the higher-order Results stopping detail (p. 468) were read directly; every drift-watch quote below cross-checked against the extracted text; DOI printed on the source — observed 2026-07-23.

**Citation.** Forbes, M. K., Baillie, A., Batterham, P. J., Calear, A., Kotov,
R., Krueger, R. F., Markon, K. E., Mewton, L., Pellicano, E., Roberts, M.,
Rodriguez-Seijas, C., Sunderland, M., Watson, D., Watts, A. L., Wright,
A. G. C., & Clark, L. A. (2025). Reconstructing psychopathology: A data-driven
reorganization of the symptoms in the *Diagnostic and Statistical Manual of
Mental Disorders*. *Clinical Psychological Science, 13*(3), 462–488.
https://doi.org/10.1177/21677026241268345 (© The Author(s) 2024.)

**Role.** Manuscript/vignette **citation precedent** and a **Forbes-authored
extended-bass-ackward continuation** of [[forbes2023]] — Miriam K. Forbes is
first author. It applies "hierarchical principal components analysis (hPCA;
i.e., an extended bass-ackward approach; Forbes, 2023b)" to a data-driven
reconstruction of the DSM-5 symptom set (N = 14,762; split primary n = 11,762 /
hold-out n = 3,000), cross-validated against Ward's hierarchical clustering.
Method is **extended bass-ackward** (the forbes2023 method), so it doubles as a
**drift-watch** against the forbes2023 fidelity contract (see below). Nothing in
`R/` traces here.

## Extracted values

- The method — "using both hierarchical principal components analysis (hPCA; i.e., an extended bass-ackward approach; Forbes, 2023b) and hierarchical clustering (i.e., Ward's hierarchical agglomerative clustering; J. H. Ward, 1963)", p. 467.
- Rotation — "We used oblique (oblimin rotation) principal components analysis (PCA) rather than factor analysis to maximize computational economy given the complex hierarchy with many levels and components expected", p. 467.
- Extraction stopping — "The number of components extracted was guided by parallel analysis, but we required at least three variables to have a unique primary loading ≥ |.4| (i.e., with all cross-loadings < |.4|) on each dimension so future work can operationalize the constructs as latent variables", p. 467.
- Redundancy pruning — "After redundant components were removed from each solution (correlations > .9 and congruence coefficients > .95 for all variables in a chain), the hierarchical clustering solution in each sample was examined for convergence/divergence with the hPCA structure to determine which principal components in the latter were possible statistical artifacts (see also Forbes, 2023b)", p. 467.
- Stopping in practice — "Parallel analysis indicated a maximum of 30 components to be extracted, but only the first six components had at least three variables with a unique primary loading ≥ |.4|, so we extracted one to six components", p. 468.

## Forbes fidelity-contract drift-watch

Whether this paper's extended-BA/hPCA procedure diverges from the [[forbes2023]]
fidelity contract (`ackwards`' O1/O2 reproduction target; redundancy chased at
`r ≥ .9` conjoined with Tucker φ > .95 per D-017):

- **Redundancy pruning — no divergence.** forbes2025 uses "correlations > .9 and congruence coefficients > .95 for all variables in a chain" (p. 467) — the same conjunctive criterion as forbes2023, i.e. `prune("redundant")` with `redundancy_phi = 0.95` and the direct/chain criterion (D-017) — observed 2026-07-23.
- **Stopping/extraction — a study-design heuristic `ackwards` does not adopt as a default.** It caps extraction by a component well-definedness rule (parallel analysis, then "at least three variables … unique primary loading ≥ |.4|"; p. 467–468). `ackwards` requires an explicit `k_max` and retains every level 1..k_max for inspection (D-013; [[tong2025]]); it neither auto-selects nor imposes this well-definedness cap. Not a fidelity-contract divergence (the contract is the redundancy chase, which matches); it is a downstream analytic choice — observed 2026-07-23.
- **Rotation — oblique oblimin, versus `ackwards`' orthogonal-only default (D-002).** A departure at the rotation layer, shared with [[partsch2022]]'s oblique Promax and wright2014a's conjoint-step oblique Geomin (`applications.md`); already catalogued as a Goldberg/Forbes-vs-ackwards departure in [[source-departures]] — observed 2026-07-23.

## Traces to

Nothing in the repo reads this page yet. It is a **citation precedent** (a
Forbes-authored applied continuation of the method) and a drift-watch witness;
it backs no code value or oracle.

## Open questions

- The full component-loading tables and level-by-level edges live in the OSF supplement (https://osf.io/urpqt/ per p. 467); only the method, rotation, stopping, and redundancy facts the drift-watch needs were transcribed — observed 2026-07-23.
