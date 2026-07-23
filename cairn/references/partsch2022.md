# partsch2022 — bass-ackwards on the 24 VIA character strengths via Waller's code (citation precedent)

**Provenance.** Ingested 2026-07-23 by M71 from
`cairn/references/sources/partsch2022.pdf` (gitignored). Pagination: journal
pages (825–845). The PDF/journal map is PDF page N → journal p. 824+N (PDF p. 1
= p. 825).
Extraction: verified 2026-07-23 against the rendered source — the title/abstract (p. 825) and the "Step 2: Bass-ackwards analysis" methods section (p. 832) were read directly; DOI printed on the source — observed 2026-07-23.

**Citation.** Partsch, M. V., Bluemke, M., & Lechner, C. M. (2022). Revisiting
the hierarchical structure of the 24 VIA character strengths: Three global
dimensions may suffice to capture their essence. *European Journal of
Personality, 36*(5), 825–845. https://doi.org/10.1177/08902070211017760
(© The Author(s) 2021; received 1 July 2020, accepted 21 April 2021.)

**Role.** Manuscript/vignette **citation precedent**: bass-ackwards applied in
a **personality / character-strengths** domain (the 24 VIA strengths measured
with the IPIP-VIA-R in two samples, total N ≈ 2,000 from Germany and the UK),
concluding three global dimensions (*positivity, dependability, mastery*).
Method is confirmed **Goldberg's (2006) bass-ackwards** (AC3): a top-down PCA
hierarchy built with **Waller's (2007) R code** ([[waller2007]], the same
algebraic basis `ackwards` wraps) — but modified to apply **oblique** Promax
rotation, differing from `ackwards`' orthogonal-only default (D-002). Edges are
adjacent-level component correlations. A clean, code-transparent bass-ackwards
precedent outside psychopathology.

## Extracted values

- The bass-ackwards step (AC3 anchor) — "Step 2: Bass-ackwards analysis. To unfold the solutions-hierarchy of the VIA strengths, we conducted Bass-ackwards analyses (Goldberg, 2006) both in Germany and in the UK", p. 832.
- Procedure — "Using the unit-weighted scale scores for each of the 24 character strengths as input, we conducted PCAs extracting an increasing number of components—first one component, then two components, and so forth. We then computed correlations between components of adjacent levels", p. 832.
- Code + rotation — "For our Bass-ackwards analyses, we used the R code provided by Waller (2007). We modified Waller's function to apply to obliquely rotated principal components (Promax rotation, m = 4)", p. 832.
- Tucker's Phi bands used for cross-country replicability — "we interpreted values of Tucker's Phi of Φ ≥ .95 as essentially equivalent, Φ ≥ .90 as highly similar, and Φ ≥ .85 as fairly similar", p. 832 (cf. [[lorenzoseva2006]]).
- Spelling — the paper writes "Bass-ackwards" (medial hyphen, no medial *s*) throughout.
- Result — "three global dimensions suffice to capture the essence of character strengths: Level III recovered more than 50% of the total variation of the 24 character strengths" (abstract, p. 825).

## Traces to

Nothing in the repo reads this page yet. It is a **citation precedent** for the
BRM manuscript / vignettes (bass-ackwards outside psychopathology; a
Waller-2007-code implementation, and an oblique-rotation variant); it backs no
code value or oracle.

## Open questions

- The loading matrices and the level-by-level Tucker's Phi congruency tables (the cross-country replicability evidence) were not transcribed — the method confirmation and rotation/code facts above are what the repo relies on; the note earlier rested on secondary sources because the PDF was access-blocked, now resolved against the shelf PDF — observed 2026-07-23.
