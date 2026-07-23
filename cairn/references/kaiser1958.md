# kaiser1958 — origin of the varimax rotation criterion (backs the varimax default)

**Provenance.** Ingested 2026-07-23 by M70 from
`cairn/references/sources/kaiser1958.pdf` (gitignored).
Pagination: journal pages (article begins p. 187; range 187–200).
Extraction: verified 2026-07-23 against the rendered source — the title page and abstract (p. 187) were read directly. (The 1958 source prints no DOI, predating DOI registration; the registered `10.1007/BF02289233` is carried from the DOI registry, as the Citation block records.) — observed 2026-07-23.

**Citation.** Kaiser, H. F. (1958). The varimax criterion for analytic rotation
in factor analysis. *Psychometrika, 23*(3), 187–200. The source prints no DOI;
commonly cited as `10.1007/BF02289233`.

**Role.** The founding paper of the varimax rotation criterion — the rotation
`ackwards` hardcodes at every level (DESIGN §9 `rotation` row; D-002). It is
here to anchor the *varimax* name and its original justification: an
**analytic** (objective, reproducible) alternative to subjective graphical
rotation. DESIGN §9 already derives varimax as the CF(κ = 1/p) special case
(Crawford & Ferguson 1970; Browne 2001); this note supplies the criterion's own
origin.

## Extracted values

- Title as printed — "THE VARIMAX CRITERION FOR ANALYTIC ROTATION IN FACTOR ANALYSIS", p. 187.
- The case for analytic rotation — "An analytic criterion for rotation is defined. The scientific advantage of analytic criteria over subjective (graphical) rotational procedures is discussed", abstract, p. 187.
- varimax vs simple structure — "the normal varimax solution probably coincides closely to the application of the principle of simple structure", abstract, p. 187.
- The ultimate criterion — "it is proposed that the ultimate criterion of a rotational procedure is factorial invariance, not simple structure", abstract, p. 187.
- Name provenance — "Dr. John Caffrey suggested the name *varimax*", acknowledgment footnote, p. 187.

## Traces to

- `cairn/DESIGN.md` §9 `rotation` row — the varimax default's rationale (M70 adds this citation beside Crawford & Ferguson 1970 / Browne 2001).
- `R/ackwards.R` `@references` — the rotation citation block (M70).

## Open questions

- Only the title page/abstract (p. 187) was read directly; the raw-vs-normal varimax derivation on the interior pages was not transcribed — not needed, since the note backs the default's origin, not its formula — observed 2026-07-23.
