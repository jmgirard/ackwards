# forbes2025a — Forbes' youth hierarchy: extended bass-ackward continuation + a deliberate fidelity-contract divergence

**Provenance.** Ingested 2026-07-23 by M71 from
`cairn/references/sources/forbes2025a.pdf` (gitignored). Pagination: journal
pages (278–300); PDF page N → journal p. 277+N (PDF p. 1 = p. 278).
Extraction: verified 2026-07-23 against the rendered source — the title/abstract (p. 278) and the "Bass-ackward models" methods + Results (pp. 287–288) were read directly; the divergence and stopping quotes below cross-checked against the extracted text; DOI printed on the source — observed 2026-07-23.

**Citation.** Forbes, M. K., Watts, A. L., Twose, M., Barrett, A., Hudson,
J. L., Lyneham, H. J., McLellan, L., Newton, N. C., Sicouri, G., Chapman, C.,
McKinnon, A., Rapee, R. M., Slade, T., Teesson, M., Markon, K., & Sunderland,
M. (2025). A hierarchical model of the symptom-level structure of
psychopathology in youth. *Clinical Psychological Science, 13*(2), 278–300.
https://doi.org/10.1177/21677026241257852 (© The Author(s) 2024.)

**Role.** Manuscript/vignette **citation precedent** and a **Forbes-authored
extended-bass-ackward continuation** of [[forbes2023]] — Miriam K. Forbes is
first author. It builds a symptom-level hierarchy of youth psychopathology
(children/adolescents, most aged 11–17; N = 18,290) using "extended
bass-ackward" cross-validated against hierarchical clustering, yielding 15
narrow dimensions under four broad domains. Method is **extended bass-ackward**
(the forbes2023 method), so it doubles as a **drift-watch** — and it is the one
M71 source that records a *deliberate* departure from the forbes2023 redundancy
criterion (below). Nothing in `R/` traces here.

## Extracted values

- The method — "We used extended bass-ackward … The traditional bass-ackward method (Goldberg, 2006) … The extended bass-ackward approach modifies this traditional method by removing redundant components and probable artifacts from the bass-ackward solution and looking at the component correlations among all remaining components to elucidate a simpler yet more complete hierarchy", pp. 286–287.
- Extraction stopping — "Parallel analysis suggested extracting up to 44 components, but only the first 19 components had at least three unique primary indicators ≥ .4 for each component. We therefore carried forward the first 19 levels of the hierarchy", p. 287.
- Top-of-hierarchy choice — low convergence in the upper levels "led us to take Level 7 as the top of the hierarchy in the extended bass-ackward approach", p. 288.
- Artifact removal — "We removed all redundant and artifactual components as per the extended bass-ackward approach (Forbes, 2023) … (n = 14 [8.2%] of the 171 total components)", p. 288.
- Result — "included 15 narrow dimensions nested under four broad dimensions of (a) internalizing, (b) externalizing, (c) eating pathology, and (d) uncontrollable worry, obsessions, and compulsions" (abstract, p. 278).

## Forbes fidelity-contract drift-watch

Whether this paper's extended-BA procedure diverges from the [[forbes2023]]
fidelity contract (`ackwards`' O1/O2 reproduction target; redundancy chased at
`r ≥ .9` conjoined with Tucker φ > .95 per D-017):

- **Redundancy pruning — a deliberate, stated divergence.** The authors dropped the Tucker-φ conjunction: "Forbes (2023) suggested using both r ≥ .9 and a congruence coefficient > .95 to establish redundancy between components, but we used a less stringent cutoff here—requiring only r ≥ .9—because the traditional bass-ackward solution included an intractable number of components (n = 171). This less stringent cutoff is consistent with Goldberg's (2006) proposed approach for identifying components that perpetuate between levels of the hierarchy" (p. 287) — observed 2026-07-23.
  - **Disposition — no capability gap; already parameterized in `ackwards`.** This exact looser cutoff is `prune("redundant")` with the φ filter off (`redundancy_phi = NA`, or the PCA auto-default of no-φ; DESIGN §9), whereas forbes2023's stricter conjunction is `redundancy_phi = 0.95`. `ackwards` supports both, so this is a *witness that the `redundancy_phi` choice matters in practice*, not a missing feature — **no new ROADMAP candidate is warranted** (creating one would be a no-op against existing capability) — observed 2026-07-23.
- **Stopping/extraction — same well-definedness heuristic as [[forbes2025]], which `ackwards` does not adopt as a default** (parallel analysis capped by "at least three unique primary indicators ≥ .4"; p. 287). `ackwards` requires explicit `k_max` and retains all levels (D-013; [[tong2025]]). Not a fidelity-contract divergence — observed 2026-07-23.

## Traces to

Nothing in the repo reads this page yet. It is a **citation precedent** (a
Forbes-authored applied continuation) and the drift-watch witness for the
`redundancy_phi`-off cutoff; it backs no code value or oracle.

## Open questions

- The full component set and level structure are in the OSF supplement (https://osf.io/5ptsn/ per p. 287) and Table S3; only the method, stopping, and the divergence facts the drift-watch needs were transcribed — observed 2026-07-23.
