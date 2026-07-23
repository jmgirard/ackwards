# williams2025 — SEM vs factor-score hierarchies for criterion validity (backs the factor-score approach + its categorical-indicator caveat)

**Provenance.** Ingested 2026-07-23 by M70 from
`cairn/references/sources/williams2025.pdf` (gitignored; the shelf file was
renamed from `williams2025a.pdf` at M70).
Pagination: journal pages (128–145).
Extraction: verified 2026-07-23 against the rendered source — the title and abstract (p. 128) were read directly; DOI confirmed via Crossref — observed 2026-07-23.

**Citation.** Williams, A. L., Conway, C. C., Olino, T. M., Revelle, W.,
Zinbarg, R. E., & HiTOP Utility Workgroup. (2025). Testing criterion validity in
hierarchical models of psychopathology: Comparison of latent-variable and
factor-score approaches. *Clinical Psychological Science, 13*(1), 128–145.
DOI `10.1177/21677026231225414` (printed on source; Crossref-confirmed).

**Role.** A simulation study comparing **latent-variable (SEM)** and
**factor-score** operationalizations of hierarchical psychopathology dimensions
for criterion-validity estimation — the factor-score side being the family
`ackwards` computes (tenBerge scores → between-level edges). It backs (a) the
legitimacy of the factor-score approach for continuous indicators and (b) the
caveat that factor scores are *less* accurate than SEM for dichotomous /
categorical indicators — echoing DESIGN §9's ordinal-ESEM handling (polychoric +
WLSMV) and GP3 descriptive honesty about what factor-score edges do and don't
support.

## Extracted values

- Continuous indicators — "In models based on continuously distributed psychopathology indicators (e.g., symptom composites), SEM and factor-score methods both tended to yield unbiased estimates of criterion validity coefficients", abstract, p. 128.
- Dichotomous indicators — "for models based on dichotomous indicators (e.g., categorical diagnoses), SEM led to more accurate estimates than factor scores in most cases", abstract, p. 128.
- Reusable tooling — the authors "provide an R function (https://osf.io/u3j5d/) that investigators can use to apply the approaches studied here in real-world data sets", abstract, p. 128.

## Traces to

- `cairn/DESIGN.md` §9 `scores (method)` row and the ordinal/ESEM handling — the factor-score approach and its categorical-indicator caveat (M70 adds this citation).
- `R/ackwards.R` `@references` — the scoring citation block (M70).

## Open questions

- The note relies on the abstract's headline findings; the specific simulation conditions and the magnitude of the factor-score disadvantage under dichotomous indicators were not transcribed — observed 2026-07-23.
