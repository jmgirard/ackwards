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
#' Each item column carries its public-domain IPIP stem (Goldberg, 1999) as a
#' `label` attribute, so `ackwards()` captures it at fit time and `top_items()`
#' prints the wording as `code: label` (e.g. `E4: Make friends easily`) with
#' no setup. These are plain attributes: base row-subsetting (e.g.
#' `na.omit(bfi25)`, `bfi25[rows, ]`) drops them, as base R does for any
#' non-`labelled`-class vector, so fit on `bfi25` **directly** -- its `NA`s are
#' handled by the `missing` argument of `ackwards()` -- to keep the labels.
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
#' attr(bfi25$E4, "label") # public-domain IPIP item stem
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
#' true group factors exactly; all six `suggest_k()` recommendations (five
#' criteria -- VSS reports at complexities 1 and 2) reach a consensus of
#' `k = 4`.
#'
#' **Idealized by design.** The planted signal is strong and clean, so all six
#' `suggest_k()` recommendations converge on `k = 4` -- deliberately the *easy* case,
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
#' `prune(x, "redundant")` (`|r| >= .9` and, under the EFA auto-default,
#' Tucker's phi `> .95`). This is a textbook overextraction artifact, included so the
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

#' Assessing Mental Health symptom correlation matrix (Forbes 2023 applied example)
#'
#' The 155 x 155 Spearman correlation matrix among 155 mental-health symptom
#' variables that forms the applied example in Forbes (2023). It is a real,
#' deep hierarchy: `ackwards(forbes2023, k_max = 10)` unfolds a general factor of
#' psychopathology at the top down to 10 fine-grained components, the worked
#' example that motivates the Forbes extension (`pairs = "all"`, `prune()`).
#'
#' Where `sim16` and `bfi25` are teaching foils, `forbes2023` is a
#' fidelity/reproduction dataset: a large, messy, published case bundled so the
#' package can reproduce the exact applied example analyzed in Forbes's paper.
#'
#' @format A 155 x 155 numeric matrix of Spearman correlations: symmetric, unit
#'   diagonal, correlations in roughly `0.01`--`0.94`. Row and column names are
#'   the 155 symptom-variable labels (e.g. `Impulsivity`, `Blurting`).
#'
#' @details
#' The correlations come from the Assessing Mental Health (AMH) study -- the
#' Australian general-population sample (N = 3,175) of Forbes et al. (2021) --
#' spanning symptoms of 18 DSM disorders. Being a correlation matrix, it carries
#' no per-variable sample size; supply `n_obs = 3175` to `ackwards()` if you want
#' EFA/ESEM fit statistics scaled to the original sample.
#'
#' This matrix reproduces Forbes's published results exactly: the package
#' regression test `test-forbes-fidelity.R` runs `ackwards()` on this exported
#' `forbes2023` and matches her reference implementation's between-level
#' correlations to `1.3e-14` across all 45 level-pairs at `k_max = 10`.
#'
#' To regenerate this dataset, run `source("data-raw/forbes2023.R")` from the
#' package root (downloads the source CSV from OSF).
#'
#' @source
#' Forbes's OSF project for the 2023 paper (file `corSpearman_AMH.csv`),
#' <https://osf.io/pcwm8/>, redistributed here under its Creative Commons
#' Attribution 4.0 International (CC-BY 4.0) license (see the package's
#' `LICENSE.note`). The underlying data are from the Assessing Mental Health
#' study (Forbes et al., 2021).
#'
#' @references
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' Forbes, M. K., Sunderland, M., Rapee, R. M., Batterham, P. J., Calear, A. L.,
#'   Carragher, N., Ruggero, C., Zimmerman, M., Baillie, A. J., Lynch, S. J.,
#'   Mewton, L., Slade, T., & Krueger, R. F. (2021). A detailed hierarchical
#'   model of psychopathology: From individual symptoms up to the general factor
#'   of psychopathology. *Clinical Psychological Science*, 9(2), 139--168.
#'   \doi{10.1177/2167702620954799}
#'
#' @seealso `vignette("ackwards-forbes2023")` for a full reproduction of the
#'   Forbes (2023) applied example on this dataset, and
#'   `vignette("ackwards-forbes")` for the method extension itself.
#'
#' @examples
#' dim(forbes2023)
#' forbes2023[1:3, 1:3]
#' \donttest{
#' # The Forbes (2023) applied example: a 10-level hierarchy from 155 symptoms.
#' x <- ackwards(forbes2023, k_max = 10, pairs = "all")
#' x
#' }
"forbes2023"
