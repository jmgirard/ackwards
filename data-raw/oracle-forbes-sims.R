# Generate the simulation fidelity fixture
#   tests/testthat/fixtures/forbes2023_sims.rds                            (M57)
#
# Forbes's (2023) three simulation studies, from her OSF project
# https://osf.io/pcwm8/ (CC-BY 4.0). Each study defines a factor structure
# (fx loadings + optional Phi factor correlations), simulates n = 5000 cases
# under set.seed(123) via psych::sim.structure(), and takes the Spearman
# correlation matrix. The expected between-level values (comp_corr, cong,
# corr_chase) are then computed with HER reference implementation
# ("ExtendedBassAckwards functions with annotation.R"), so no Forbes code or
# network runs at test time -- only ackwards() runs (see test-forbes-fidelity.R).
#
# WHY THE MATRICES ARE REGENERATED HERE (M57)
# -------------------------------------------------------------------------
# The pre-M57 fixture carried the three matrices but no committed generator, and
# its exact random realizations proved unreproducible (no fx/Phi/seed/n/method
# combination reproduces them -- even n = 1e5 Spearman stays ~0.036 off). This
# script replaces those murky artifacts with the DETERMINISTIC realizations of
# Forbes's exact recipe: given a fixed R + psych, set.seed(123) makes each matrix
# bit-reproducible (two runs agree to 0). The fidelity contract is unaffected --
# ackwards reproduces her ExtendedBassAckwards on the regenerated matrices to
# ~4e-15, and the redundancy topology is preserved (the prune expectations in
# test-forbes-fidelity.R still hold). See M57 decision log.
#
# The fx / Phi specifications below are transcribed verbatim from her OSF script
# "R script for simulations and applied example_R2.R" (guid ztngp); the md5 pin
# guards that that script has not changed under us. num.comp = 4 for all three
# (her Simulation-3 DGP has 5 factors but she analyses it at 4 levels, k = 4).
#
# Re-run after any change to the DGP or expected-value definitions:
#   Rscript data-raw/oracle-forbes-sims.R
# Requires network access to osf.io. `psych` must be installed.

# Forbes's functions call fa.sort()/pca() unqualified, so psych must be attached.
suppressPackageStartupMessages(library(psych))

## --- OSF sources (project pcwm8, CC-BY 4.0; guids resolved via the OSF API) ---
osf <- list(
  functions = list(
    guid = "7jfkw", name = "ExtendedBassAckwards functions with annotation.R",
    url = "https://osf.io/download/7jfkw/", md5 = "a3e85df897d2a4a4310b9c45dcb068d9"
  ),
  sim_script = list(
    guid = "ztngp", name = "R script for simulations and applied example_R2.R",
    url = "https://osf.io/download/ztngp/", md5 = "8de5ec351992e8fb11c0b8d7015f3329"
  )
)

tmp <- tempfile(fileext = "_sims")
dir.create(tmp)
fun_path <- file.path(tmp, "ExtendedBassAckwards.R")
scr_path <- file.path(tmp, "sim_script.R")
utils::download.file(osf$functions$url, fun_path, mode = "wb", quiet = TRUE)
utils::download.file(osf$sim_script$url, scr_path, mode = "wb", quiet = TRUE)

## Integrity guards: pin the exact published files. The functions file is sourced
## (it defines the reference implementation); the sim script is pinned only to
## guarantee the transcribed fx / Phi below still match her published DGP.
stopifnot(
  unname(tools::md5sum(fun_path)) == osf$functions$md5,
  unname(tools::md5sum(scr_path)) == osf$sim_script$md5
)
source(fun_path, local = TRUE) # defines ExtendedBassAckwards(), FindRedundantComp()

## --- The three DGPs (transcribed from the md5-pinned sim script) --------------
## Sim 1: a hierarchy that unfolds over multiple levels (Phi links factors 2-3).
fx1 <- matrix(c(
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0,   0,
  0.4, 0,   0.6, 0,
  0.4, 0,   0.6, 0,
  0.4, 0,   0.6, 0,
  0.4, 0,   0.6, 0,
  0.3, 0,   0,   0.8,
  0.3, 0,   0,   0.8,
  0.8, 0,   0,   0.3,
  0.8, 0,   0,   0.3
), 12, 4, byrow = TRUE)
Phi1 <- matrix(c(
  1, 0,   0,   0,
  0, 1,   0.6, 0,
  0, 0.6, 1,   0,
  0, 0,   0,   1
), 4, 4, byrow = TRUE)

## Sim 2: item-level cross-loadings a cluster analysis cannot capture (Phi 0.3).
fx2 <- matrix(c(
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0,   0,
  0.4, 0.6, 0.3, 0,
  0.4, 0.3, 0.6, 0,
  0.4, 0,   0.6, 0,
  0.4, 0,   0.6, 0,
  0.4, 0,   0.6, 0,
  0.3, 0,   0,   0.8,
  0.3, 0,   0,   0.8,
  0.8, 0,   0,   0.3,
  0.8, 0,   0,   0.3
), 12, 4, byrow = TRUE)
Phi2 <- matrix(c(
  1, 0,   0,   0,
  0, 1,   0.3, 0,
  0, 0.3, 1,   0,
  0, 0,   0,   1
), 4, 4, byrow = TRUE)

## Sim 3: a two-tier hierarchy; five orthogonal group factors (no Phi).
fx3 <- matrix(c(
  0.6, 0.4, 0,   0,   0,
  0.6, 0.4, 0,   0,   0,
  0.6, 0.4, 0,   0,   0,
  0.6, 0,   0.4, 0,   0,
  0.6, 0,   0.4, 0,   0,
  0.6, 0,   0.4, 0,   0,
  0.6, 0,   0,   0.4, 0,
  0.6, 0,   0,   0.4, 0,
  0.6, 0,   0,   0.4, 0,
  0.6, 0,   0,   0,   0.4,
  0.6, 0,   0,   0,   0.4,
  0.6, 0,   0,   0,   0.4
), 12, 5, byrow = TRUE)

## One simulation -> Spearman matrix + expected values from her reference impl.
## set.seed(123) immediately before each draw (as her script does) makes the
## matrix deterministic; num.comp = 4 and last component "d4" match her k = 4.
K <- 4L
one_sim <- function(fx, Phi = NULL) {
  set.seed(123)
  obs <- if (is.null(Phi)) {
    sim.structure(fx = fx, n = 5000)$observed
  } else {
    sim.structure(fx = fx, Phi = Phi, n = 5000)$observed
  }
  R <- cor(obs, use = "pairwise.complete.obs", method = "spearman")
  fb <- ExtendedBassAckwards(R, num.comp = K, fm = "pca")
  chase <- unlist(FindRedundantComp(fb$comp.corr, fb$cong, "d4")$corr.chase)
  stopifnot(
    dim(R) == c(12L, 12L),
    isSymmetric(unname(R), tol = 1e-8),
    all(abs(diag(R) - 1) < 1e-8),
    length(fb$comp.corr) == choose(K, 2), # 6 level-pairs
    length(fb$cong) == choose(K, 2),
    length(chase) == sum(2:K) # b1..d4 = 9 components (a1 excluded)
  )
  list(R = R, comp_corr = fb$comp.corr, cong = fb$cong, corr_chase = chase)
}

sims <- list(
  sim1 = one_sim(fx1, Phi1),
  sim2 = one_sim(fx2, Phi2),
  sim3 = one_sim(fx3) # orthogonal: no Phi
)

## Structured, top-level provenance (uniform target for the M57 guard test).
attr(sims, "provenance") <- list(
  source = paste0(
    "Forbes, M. K. (2023). Improving hierarchical models of individual ",
    "differences: An extension of Goldberg's bass-ackward method. ",
    "Psychological Methods. doi:10.1037/met0000546"
  ),
  osf = "https://osf.io/pcwm8/",
  license = "CC-BY 4.0 International (https://creativecommons.org/licenses/by/4.0/)",
  files = vapply(osf, function(f) sprintf("%s (guid %s, md5 %s) <%s>", f$name, f$guid, f$md5, f$url), character(1)),
  recipe = paste0(
    "Per study: set.seed(123); psych::sim.structure(fx[, Phi], n = 5000)$observed; ",
    "Spearman correlation matrix. Expected comp_corr/cong/corr_chase from Forbes's ",
    "ExtendedBassAckwards(R, num.comp = 4, fm = 'pca'). fx/Phi transcribed from the ",
    "md5-pinned OSF sim script (guid ztngp)."
  ),
  note = paste0(
    "Matrices are the deterministic realizations of Forbes's recipe (M57): the ",
    "pre-M57 fixture's draws were unreproducible and were replaced. Only ackwards() ",
    "runs at test time."
  ),
  generator = "data-raw/oracle-forbes-sims.R",
  generated = as.character(Sys.Date()),
  R_version = R.version.string,
  psych = as.character(utils::packageVersion("psych"))
)

out <- file.path("tests", "testthat", "fixtures", "forbes2023_sims.rds")
saveRDS(sims, out, compress = "xz")
message(sprintf("Wrote %s (%.0f KB)", out, file.size(out) / 1024))
