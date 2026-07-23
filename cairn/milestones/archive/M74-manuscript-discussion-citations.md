# M74: Draft the manuscript Discussion + enrich the seeded sections with method-backer citations

**Status:** done (2026-07-23, PR #78 https://github.com/jmgirard/ackwards/pull/78)

**Goal:** Replace the manuscript's `[AUTHOR TO DRAFT]` Discussion stub with a complete citation-backed first draft, and backfill the seeded Method/Package sections with the method-backer citations their uncited claims are owed.

**Outcome:** `manuscript/manuscript.qmd` `# Discussion` is now a five-paragraph draft (four beats: what the package enables over bespoke scripts; interpretive cautions; scope/limitations with the descriptive-not-confirmatory + sequential-not-hierarchical contrast vs. Schmid–Leiman/higher-order; concrete future directions + availability). Seeded sections gained five citations-only backfills (`kaiser1958` varimax, `asparouhov2009` ESEM, `grice2001` factor scores, `revelle1979`+`ruscio2012a` suggest_k criteria) with no technical wording change. `references.bib` gained 7 keys, each DOI verified live against Crossref. Renders clean to PDF + docx. Diff is `manuscript/` + `cairn/` only (both `.Rbuildignore`d) — no package code touched.

**Decisions:** none (the citekey→note provenance map recorded in the live file was verification evidence, not a choice).

**Review:** AC1–AC6 all PASS on fresh evidence; consistency gate clean (cairn_validate exit 0; R-toolchain vacuous — build package untouched; pkgdown OK). Three-lens fan-out + scorer: Finding B (em dashes reverting the `1629d26` style pass, score 93) actioned/fixed; A (opening over-claim vs. Statement of Need, 75) + D (Goldberg verbatim-attribution, 42) fixed anyway; C (CRAN-availability while 0.1.1 pending, 62) rejected for M74 (mirrors unchanged Package-section claim; owner release-tail task) and surfaced to author. Merged local-green per repo non-release CI override.
