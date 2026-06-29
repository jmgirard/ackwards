#' Big Five Inventory -- 25-item IPIP example dataset
#'
#' A 1 000-row subset of the 25 IPIP Big Five personality-marker items from
#' Revelle's `psych` package (`psych::bfi`). Included so that package examples
#' and vignettes run without reaching into `psych`'s namespace directly.
#'
#' @format A data frame with 1 000 rows and 25 integer columns (Likert responses
#'   scored 1--6, with some `NA`s reflecting genuine missing values in the
#'   original survey). The 25 items span the Big Five personality domains:
#'
#'   | Columns | Domain |
#'   |---------|--------|
#'   | `A1`--`A5` | Agreeableness |
#'   | `C1`--`C5` | Conscientiousness |
#'   | `E1`--`E5` | Extraversion |
#'   | `N1`--`N5` | Neuroticism |
#'   | `O1`--`O5` | Openness |
#'
#' @details
#' Derived from `psych::bfi[, 1:25]` by sampling 1 000 rows with
#' `set.seed(42)`. Missing values are preserved as in the original data.
#' The full dataset (`psych::bfi`) contains responses from 2 800 participants
#' collected via the SAPA project.
#'
#' To regenerate this dataset, run `source("data-raw/bfi25.R")` from the
#' package root.
#'
#' @source
#' Revelle, W. (2026). *psych: Procedures for Psychological, Psychometric, and
#'   Personality Research*. Northwestern University, Evanston, Illinois.
#'   R package. <https://CRAN.R-project.org/package=psych>
#'
#' Data collected via the SAPA (Synthetic Aperture Personality Assessment)
#' project (<https://www.sapa-project.org/>). Items are public-domain IPIP
#' (International Personality Item Pool) markers developed by Goldberg (1999):
#'
#' Goldberg, L. R. (1999). A broad-bandwidth, public domain, personality
#'   inventory measuring the lower-level facets of several five-factor models.
#'   In I. Mervielde, I. Deary, F. De Fruyt, & F. Ostendorf (Eds.),
#'   *Personality Psychology in Europe*, Vol. 7, pp. 7--28. Tilburg University Press.
#'
#' @examples
#' dim(bfi25)
#' head(bfi25)
"bfi25"
