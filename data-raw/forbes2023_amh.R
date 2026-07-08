# Generate tests/testthat/fixtures/forbes2023_amh.rds  (M53)
#
# Forbes's (2023) 155-variable "Assessing Mental Health" (AMH) applied example,
# from her OSF project https://osf.io/pcwm8/ (CC-BY 4.0 International). This
# generator downloads the published Spearman correlation matrix and her own
# reference implementation, computes the between-level correlations and loading
# congruences with HER code, and stores them alongside the matrix so the test
# suite runs only ackwards() at test time (no vendored Forbes code, no network).
#
# The AMH matrix is redistributed here under CC-BY 4.0 (see LICENSE.note);
# the simulation fixture (forbes2023_sims.rds) needs no such note because those
# matrices are seed-regenerated, not Forbes data.
#
# Re-run after any change to the expected-value definitions:
#   Rscript data-raw/forbes2023_amh.R
# Requires network access to osf.io. `psych` must be installed.

# Forbes's functions call fa.sort() unqualified, so psych must be attached.
suppressPackageStartupMessages(library(psych))

## --- OSF sources (guids resolved via the OSF API; project pcwm8, CC-BY 4.0) ---
osf <- list(
  matrix = list(
    guid = "s9bjz", name = "corSpearman_AMH.csv",
    url = "https://osf.io/download/s9bjz/"
  ),
  functions = list(
    guid = "7jfkw", name = "ExtendedBassAckwards functions with annotation.R",
    url = "https://osf.io/download/7jfkw/"
  )
)

tmp <- tempfile(fileext = "_amh")
dir.create(tmp)
csv_path <- file.path(tmp, "corSpearman_AMH.csv")
fun_path <- file.path(tmp, "ExtendedBassAckwards.R")
utils::download.file(osf$matrix$url, csv_path, mode = "wb", quiet = TRUE)
utils::download.file(osf$functions$url, fun_path, mode = "wb", quiet = TRUE)

## --- The 155x155 Spearman matrix (row/col names carried through) ---
R <- as.matrix(utils::read.csv(csv_path, row.names = 1, check.names = FALSE))
colnames(R) <- rownames(R)
stopifnot(dim(R) == c(155L, 155L), isSymmetric(unname(R)), all(diag(R) == 1))

## --- Forbes's own reference implementation -> expected values ---
source(fun_path, local = TRUE) # defines ExtendedBassAckwards()
K <- 10L # her applied example uses k = 10
fb <- ExtendedBassAckwards(R, num.comp = K, fm = "pca")
stopifnot(length(fb$comp.corr) == 45L, length(fb$cong) == 45L) # choose(10, 2)

## Her redundancy chase (ChaseCorrPaths) for every component b1..j10: the
## direct/skip-level rule the paper uses. Stored as her raw "X--Y" strings
## ("X--null" = no chase), mirroring the simulation fixture's corr_chase.
corr_chase <- unlist(FindRedundantComp(fb$comp.corr, fb$cong, last_component = "j10")$corr.chase)
stopifnot(length(corr_chase) == sum(2:K)) # b1..j10 = 54 components (a1 excluded)

## comp.corr / cong enumerate pairs as: for c in 2..K, for i in 1..(c-1).
## comp.corr[[idx]] = t(W_i) R W_c  (unstandardized, unaligned) -> |.| == our edges.
## cong[[idx]]      = psych::factor.congruence (rounded to 2 dp) for the same pair.
amh <- list(
  R          = R,
  comp_corr  = fb$comp.corr,
  cong       = fb$cong,
  corr_chase = corr_chase,
  k_max      = K
)

attr(amh, "provenance") <- list(
  source = paste0(
    "Forbes, M. K. (2023). Improving hierarchical models of individual ",
    "differences: An extension of Goldberg's bass-ackward method. ",
    "Psychological Methods. doi:10.1037/met0000546"
  ),
  osf = "https://osf.io/pcwm8/",
  license = "CC-BY 4.0 International (https://creativecommons.org/licenses/by/4.0/)",
  files = vapply(osf, function(f) sprintf("%s (guid %s) <%s>", f$name, f$guid, f$url), character(1)),
  note = paste0(
    "R = published 155-variable AMH Spearman matrix; comp_corr/cong computed ",
    "with Forbes's ExtendedBassAckwards reference implementation. ",
    "Only ackwards() runs at test time."
  ),
  generated = as.character(Sys.Date()),
  R_version = R.version.string,
  psych = as.character(utils::packageVersion("psych"))
)

## Wrap in a named list mirroring forbes2023_sims.rds (single applied example).
fixture <- list(amh = amh)

out <- file.path("tests", "testthat", "fixtures", "forbes2023_amh.rds")
saveRDS(fixture, out, compress = "xz")
message(sprintf("Wrote %s (%.0f KB)", out, file.size(out) / 1024))
