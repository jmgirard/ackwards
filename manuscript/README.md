# Manuscript: *ackwards* (Behavior Research Methods)

A reproducible APA-7 manuscript describing the **ackwards** R package, targeted at
*Behavior Research Methods*. Every number, table, and figure is computed at render
time from the installed package — nothing is hardcoded.

This directory is excluded from the R package build (`^manuscript$` in
`.Rbuildignore`), so it never ships to CRAN.

## Contents

| File | Purpose |
|------|---------|
| `manuscript.qmd` | The manuscript source (Quarto, apaquarto format). |
| `references.bib` | BibTeX bibliography (primary sources + software). |
| `_extensions/` | Vendored [apaquarto](https://github.com/wjschne/apaquarto) extension (committed, so no `quarto add` / network is needed to render). |

Built outputs (`manuscript.pdf`, `manuscript.docx`, `manuscript_files/`,
`.quarto/`) are git-ignored; regenerate them with the render command below.

## Requirements

- **Quarto** ≥ 1.4 (developed with 1.9.38)
- **R** ≥ 4.1 (developed with R 4.6.1)
- A **LaTeX** toolchain for PDF output — [TinyTeX](https://yihui.org/tinytex/)
  is sufficient (`tinytex::install_tinytex()`); apaquarto uses XeLaTeX.
- The **ackwards** package installed and loadable (`library(ackwards)`),
  along with its `ggplot2` suggest (used by `autoplot()`):

  ```r
  # from the repository root, install the current source:
  install.packages("remotes")
  remotes::install_local(".", dependencies = TRUE)   # or install from CRAN
  install.packages("ggplot2")
  ```

## Render

From this directory:

```sh
quarto render manuscript.qmd            # both PDF and Word
quarto render manuscript.qmd --to apaquarto-pdf
quarto render manuscript.qmd --to apaquarto-docx
```

## Reproducibility record

Last verified render:

- R version 4.6.1 (2026-06-24)
- ackwards 0.1.1
- Quarto 1.9.38
- apaquarto 5.0.18

The applied example (`ackwards(forbes2023, k_max = 10, pairs = "all")`) reproduces
Forbes's (2023) reference implementation; this is guarded by the package test
`tests/testthat/test-forbes-fidelity.R` (agreement to ~1e-14 across all 45
between-level pairs; OSF `pcwm8`).

## Status

The **Introduction** is a complete first draft awaiting the author's revision;
the **Discussion** remains a structured stub (marked `[AUTHOR TO DRAFT]`)
awaiting scholarly prose. The Abstract, Method, package description, worked
examples, tables, and figures are complete.
