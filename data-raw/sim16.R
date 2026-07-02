## data-raw/sim16.R
##
## Creates the sim16 example dataset: a simulated, fully continuous 16-item
## bass-ackwards teaching example with a known ground-truth hierarchy.
##
## Population model: Sigma = Lambda %*% Phi %*% t(Lambda) + Psi (an oblique
## common-factor model), sampled via a Cholesky factorization -- no MASS
## dependency, base R only.
##
##   4 true group factors (f1-f4), 4 items each, loading = 0.75:
##     f1 = i1-i4     f2 = i5-i8     f3 = i9-i12    f4 = i13-i16
##   Factor correlations (Phi):
##     within metatrait (f1-f2, f3-f4)  = 0.45  -- {f1,f2} = metatrait 1 (i1-i8)
##     between metatrait (f1/f2-f3/f4)  = 0.15  -- {f3,f4} = metatrait 2 (i9-i16)
##   Psi = 1 - communality (uniform 0.4375 by symmetry of the design above).
##
## This yields a clean 1 -> 2 -> 4 bass-ackwards hierarchy (verified against
## ackwards(engine = "efa")): k=1 general factor, k=2 splits along the
## metatrait line, k=4 recovers the 4 true group factors exactly. Because
## the population has exactly 4 factors, k=5 has no real 5th dimension to
## find: EFA produces an orphan factor with zero primary-loading items,
## signalled by prune(x, "artifact") (few_items + orphan), while the
## persisting (non-splitting) factors from k=3 onward form redundant chains
## (|r| >= .9, and Tucker's phi > .95 under the EFA auto-default) under
## prune(x, "redundant") -- a textbook overextraction artefact, deliberately
## exercised so the Forbes vignette has a guaranteed finding to teach against
## (unlike bfi25, where it does not reliably trigger). See R/data.R for the
## full ground-truth documentation.
##
## To regenerate: source this script from the package root via
##   source("data-raw/sim16.R")

set.seed(42)

loading <- 0.75
n_items <- 16L
group_factors <- c("f1", "f2", "f3", "f4")

Lambda <- matrix(0, nrow = n_items, ncol = 4L)
dimnames(Lambda) <- list(paste0("i", seq_len(n_items)), group_factors)
Lambda[1:4, "f1"] <- loading
Lambda[5:8, "f2"] <- loading
Lambda[9:12, "f3"] <- loading
Lambda[13:16, "f4"] <- loading

within_metatrait <- 0.45
between_metatrait <- 0.15

Phi <- diag(4L)
dimnames(Phi) <- list(group_factors, group_factors)
Phi["f1", "f2"] <- Phi["f2", "f1"] <- within_metatrait
Phi["f3", "f4"] <- Phi["f4", "f3"] <- within_metatrait
for (a in c("f1", "f2")) {
  for (b in c("f3", "f4")) {
    Phi[a, b] <- Phi[b, a] <- between_metatrait
  }
}

Sigma_common <- Lambda %*% Phi %*% t(Lambda)
uniquenesses <- 1 - diag(Sigma_common)
stopifnot(all(uniquenesses > 0)) # positive uniquenesses => valid population model
R_pop <- Sigma_common + diag(uniquenesses)

n_obs <- 1000L
Z <- matrix(rnorm(n_obs * n_items), nrow = n_obs, ncol = n_items)
X <- Z %*% chol(R_pop)
colnames(X) <- colnames(R_pop)
sim16 <- as.data.frame(X)
rownames(sim16) <- NULL

usethis::use_data(sim16, overwrite = TRUE, compress = "xz")
