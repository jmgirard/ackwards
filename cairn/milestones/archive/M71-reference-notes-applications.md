# M71: Author + verify the 5 application source notes as citation precedents + Forbes drift-watch

**Status:** done (2026-07-23, PR #76 https://github.com/jmgirard/ackwards/pull/76)

**Goal:** Author a verified `cairn/references/` source note for each of the 5 application PDFs (carmichael2025, michelini2019, partsch2022, forbes2025, forbes2025a) as manuscript/vignette citation precedents, the two Forbes notes doubling as extended-bass-ackward drift-watch against the forbes2023 fidelity contract.

**Outcome:** Five standalone `references/` notes shipped (docs-only, `cairn/`-only), each verified against its rendered shelf PDF with page-anchored facts and a dated `Extraction:` line. Domains + methods: carmichael2025 (TBI; Forbes' extended BA, redundancy r≥.9 & φ>.95, p.721), michelini2019 (ABCD child/adolescent; Goldberg's BA on EFA/geomin scores, edge ≥.65), partsch2022 (VIA character strengths; Goldberg's BA via Waller-2007 code, oblique Promax — bass-ackwards confirmed p.832), forbes2025 (DSM-5; hPCA/extended-BA, redundancy matches the contract), forbes2025a (youth; extended-BA). Added a new `## Applications (one note per source)` INDEX subsection (5 filename-first lines). The two Forbes notes carry a fidelity drift-watch: forbes2025 matches the contract; forbes2025a records a deliberate divergence (dropped the Tucker-φ conjunction, r≥.9 only, p.287).

**Decisions:** MD-1 — the forbes2025a redundancy divergence gets no ROADMAP candidate row (the looser r≥.9-only cutoff is already expressible via `prune()`'s `redundancy_phi`; no capability gap). Owner accepted at the merge gate.

**Review:** 3-lens fan-out. 4 findings, all ground-truth-confirmed and fixed on-branch: F1 (100) forbes2025 DOI wrong (`…268545`→`…268345`); F4 (95) michelini "oblique geomin" unsourced → reframed as rotation-family departure, orientation unconfirmed (regressed the M68 orthogonal-geomin lesson); F2/F3 (85) two dangling wikilinks to collapse-only citekeys (cowan2024, wright2014a) → plain text. Blame-history + prior-PR lenses otherwise clean.
