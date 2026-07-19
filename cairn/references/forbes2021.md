# forbes2021 — provenance for the N = 3,175 AMH sample behind the `forbes2023` dataset

**Full citation.** Forbes, M. K., Sunderland, M., Rapee, R. M., … Krueger, R. F.
(2021). A detailed hierarchical model of psychopathology: From individual
symptoms up to the general factor. *Clinical Psychological Science, 9*(2),
139–168. https://doi.org/10.1177/2167702620954799 *(DOI printed
on source p. 139; many-author paper — ellipsis elides the middle authors between
Rapee and Krueger)*.

**Provenance.** Ingested 2026-07-19 by M69 from
`cairn/references/sources/forbes2021.pdf` (gitignored). Pagination: journal
pages (139–168). Extraction: verified 2026-07-19 (M69) — citation, sample size,
and disorder count read directly against the source (p. 139 header + abstract).

## Why this is a primary source

The origin study for the **Assessing Mental Health (AMH)** sample whose Spearman
correlation matrix ships as the `forbes2023` dataset. `?forbes2023` cites Forbes
et al. (2021) as the source of the N = 3,175 Australian general-population data;
the matrix itself was published in Forbes (2023, the fidelity oracle — see
[[forbes2023]]).

## Extracted values

Abstract (p. 139), verbatim: the study explored the hierarchical structure of
"symptoms spanning 18 DSM disorders in two large samples—one from the general
population in Australia (n = 3,175) and the other a treatment-seeking clinical
sample from the United States (n = 1,775)."

- **N = 3,175**, Australian general population — matches `R/data.R:139` exactly.
- **18 DSM disorders** — matches "spanning symptoms of 18 DSM disorders"
  (`R/data.R:139-140`).
- The second sample (US clinical, n = 1,775) is **not** the one the package
  uses; `forbes2023` is the Australian AMH matrix only. The docs correctly name
  only the Australian sample.

No drift found.

## Traces to

- `R/data.R:139-140` — "the Australian general-population sample (N = 3,175) of
  Forbes et al. (2021) -- spanning symptoms of 18 DSM disorders".
- `R/data.R:157` — the `@source` note "The underlying data are from the Assessing
  Mental Health study (Forbes et al., 2021)".
- `vignettes/ackwards-forbes2023.Rmd.orig:24` — the vignette's "N = 3,175;
  Forbes et al., 2021" attribution.

## Open questions

- No code drift found. The N, the country, and the 18-disorder scope all match
  the source. Nothing routed to `/hotfix`. — observed 2026-07-19.
