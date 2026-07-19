# revelle1979 — the Very Simple Structure (VSS) criterion offered by `suggest_k()`

**Full citation.** Revelle, W., & Rocklin, T. (1979). Very simple structure: An
alternative procedure for estimating the optimal number of interpretable
factors. *Multivariate Behavioral Research, 14*(4), 403–414.
*(The 1979 scan prints no DOI; the roxygen `@references` entry likewise carries
none.)*

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/revelle1979.pdf` (gitignored). Pagination: journal
pages (403–414).
Extraction: verified 2026-07-19 (M69) — the title, citation,
and the complexity-v VSS criterion read directly against the source (p. 403
title/abstract; pp. 405–406 criterion definition).

## Why this is a primary source

The source of the **VSS** number-of-factors method `suggest_k()` offers as
`"vss"` — reported at complexities 1 and 2 (VSS-1, VSS-2), delegated to
`psych::vss()`. The roxygen and the DESIGN k-advice table cite Revelle &
Rocklin (1979).

## Extracted values

- Method (abstract, p. 403), verbatim: "The new method evaluates the magnitude
  of the Very Simple Structure index of goodness of fit for factor solutions of
  increasing rank. The number of factors which maximizes the VSS criterion is
  taken [as the estimate]."
- Complexity parameter (p. 405–406): the criterion is defined "For a Very Simple
  Structure solution of factor complexity v" as `VSS_vk`, and "the number of
  interpretable factors (of complexity v) is the number of factors, k, which
  maximizes VSS_vk." The two commonly-used cases are complexity **v = 1**
  (VSS-1) and **v = 2** (VSS-2) — exactly what `psych::vss()` returns and what
  `suggest_k()` reports.

The package's roxygen — "Very Simple Structure fit at complexities 1 and 2
(VSS-1 and VSS-2). Finds the k maximising the fit of a very simple loading
structure" — matches the source. No drift found.

## Traces to

- `R/suggest_k.R:20-23` — the `"vss"` method roxygen description.
- `R/suggest_k.R:152-154` — the `@references` citation entry.
- `R/suggest_k.R:366`, `R/suggest_k.R:395-396` — the VSS-1/VSS-2 computation
  (shared `psych::vss()` call; `which.max` over `vss1_vals`/`vss2_vals`).
- `cairn/DESIGN.md:378` — the VSS-1/VSS-2 row of the k-advice table.

## Open questions

- No code drift found. Nothing routed to `/hotfix` from this source. — observed
  2026-07-19.
