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

#' Simulated continuous bass-ackwards teaching example (16 items, known hierarchy)
#'
#' A 1 000-row, fully continuous dataset simulated from a population model
#' with a known 1 -> 2 -> 4 bass-ackwards hierarchy, for showcasing the
#' default `cor = "pearson"` extraction path without the ordinal-detection
#' warning that `bfi25`'s Likert items trigger.
#'
#' @format A data frame with 1 000 rows and 16 numeric columns (`i1`--`i16`,
#'   continuous, no missing values).
#'
#' @details
#' **Population model.** `Sigma = Lambda %*% Phi %*% t(Lambda) + Psi`, an
#' oblique common-factor model with 4 true group factors, sampled via a
#' Cholesky factorization (base R; no `MASS` dependency):
#'
#' | Factor | Items | Metatrait |
#' |---|---|---|
#' | `f1` | `i1`--`i4`   | 1 (with `f2`) |
#' | `f2` | `i5`--`i8`   | 1 (with `f1`) |
#' | `f3` | `i9`--`i12`  | 2 (with `f4`) |
#' | `f4` | `i13`--`i16` | 2 (with `f3`) |
#'
#' All items load `0.75` on their true factor (no cross-loadings). Factor
#' correlations: `0.45` within a metatrait (`f1`-`f2`, `f3`-`f4`), `0.15`
#' between metatraits. Uniquenesses are `1 - communality` (uniform `0.4375`
#' by the symmetry of the design above).
#'
#' **Ground-truth hierarchy** (verified against `ackwards(engine = "efa")`):
#' `k=1` recovers a single general factor across all 16 items; `k=2` splits
#' along the metatrait line (`i1`-`i8` vs. `i9`-`i16`); `k=4` recovers the 4
#' true group factors exactly; `suggest_k()` reaches a 6-criteria consensus
#' of `k = 4`.
#'
#' **Idealized by design.** The planted signal is strong and clean, so all six
#' `suggest_k()` criteria converge on `k = 4` -- deliberately the *easy* case,
#' for building intuition about what recovering a known hierarchy looks like.
#' Real data rarely agree this cleanly: on `bfi25` the same criteria span
#' `k = 4`--`6`. The two datasets are complementary teaching foils -- `sim16`
#' for "watch the method recover a structure we planted," `bfi25` for
#' "reason about a hierarchy when the criteria disagree." Present `sim16`'s
#' consensus as the ideal, not the norm.
#'
#' **Deliberate overextraction artifact at k=5.** The population has exactly
#' 4 factors, so requesting a 5th finds no real dimension: EFA produces an
#' orphan factor with zero primary-loading items. With
#' `prune(x, "artifact")` (default `min_items = 3`, `orphan_r = 0.5`), that
#' factor is flagged both `few_items` and `orphan`. Because the true (non-
#' splitting) factors persist essentially unchanged from `k=3` onward, their
#' parent-child score correlations approach 1 and are flagged by
#' `prune(x, "redundant")` (`|r| >= .9` and, by default, Tucker's phi `>=
#' .95`). This is a textbook overextraction artifact, included so the
#' Forbes/redundancy examples have a guaranteed finding to teach against
#' (unlike `bfi25`, which does not reliably trigger one).
#'
#' To regenerate this dataset, run `source("data-raw/sim16.R")` from the
#' package root (`set.seed(42)`).
#'
#' @source Simulated; see `data-raw/sim16.R` for the full generative model.
#'
#' @examples
#' dim(sim16)
#' head(sim16)
"sim16"
