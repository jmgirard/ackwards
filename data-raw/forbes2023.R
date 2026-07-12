# Generate BOTH the exported dataset data/forbes2023.rda AND the fidelity fixture
# tests/testthat/fixtures/forbes2023_amh.rds  (M54; supersedes data-raw/forbes2023_amh.R)
#
# Forbes's (2023) 155-variable "Assessing Mental Health" (AMH) applied example,
# from her OSF project https://osf.io/pcwm8/ (CC-BY 4.0 International). One
# md5-pinned download feeds both artifacts, so the shipped matrix (`forbes2023`)
# and the fixture's expected values can never drift apart:
#
#   * data/forbes2023.rda                   -- the exported user dataset (matrix only)
#   * tests/testthat/fixtures/forbes2023_amh.rds -- expected comp_corr/cong/corr_chase
#       computed with HER reference implementation (no matrix; the test reads the
#       matrix from the exported forbes2023). No vendored Forbes code, no network at
#       test time -- only ackwards() runs.
#
# The AMH matrix is redistributed under CC-BY 4.0 (see LICENSE.note; Forbes is
# listed as data copyright holder `cph` in DESCRIPTION). The project was briefly
# CC-BY-NC in early July 2026; Forbes switched it to CC-BY on learning of the
# NonCommercial implications, so the matrix can be bundled (see legacy MILESTONES M53).
#
# Re-run after any change to the expected-value definitions:
#   Rscript data-raw/forbes2023.R
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

## Integrity guard: pin the exact published file that both artifacts are built
## from. If OSF ever re-publishes the matrix this stops silently, so the shipped
## forbes2023 and the fixture's expected values are regenerated together or not at all.
stopifnot(
  unname(tools::md5sum(csv_path)) == "c1dd9eca009c2738c268487179d43e87"
)

## --- The 155x155 Spearman matrix (row/col names carried through) ---
forbes2023 <- as.matrix(utils::read.csv(csv_path, row.names = 1, check.names = FALSE))
colnames(forbes2023) <- rownames(forbes2023)
stopifnot(
  dim(forbes2023) == c(155L, 155L),
  isSymmetric(unname(forbes2023), tol = 1e-8),
  all(abs(diag(forbes2023) - 1) < 1e-8)
)

## (1) The exported dataset -----------------------------------------------------
usethis::use_data(forbes2023, overwrite = TRUE, compress = "xz")

## (2) The fidelity fixture (expected values only; matrix comes from forbes2023) -
## Forbes's own reference implementation -> expected between-level values.
source(fun_path, local = TRUE) # defines ExtendedBassAckwards()
K <- 10L # her applied example uses k = 10
fb <- ExtendedBassAckwards(forbes2023, num.comp = K, fm = "pca")
stopifnot(length(fb$comp.corr) == 45L, length(fb$cong) == 45L) # choose(10, 2)

## Her redundancy chase (ChaseCorrPaths) for every component b1..j10: the
## direct/skip-level rule the paper uses. Stored as her raw "X--Y" strings
## ("X--null" = no chase), mirroring the simulation fixture's corr_chase.
corr_chase <- unlist(FindRedundantComp(fb$comp.corr, fb$cong, last_component = "j10")$corr.chase)
stopifnot(length(corr_chase) == sum(2:K)) # b1..j10 = 54 components (a1 excluded)

## comp.corr / cong enumerate pairs as: for c in 2..K, for i in 1..(c-1).
## comp.corr[[idx]] = t(W_i) R W_c  (unstandardized, unaligned) -> |.| == our edges.
## cong[[idx]]      = psych::factor.congruence (rounded to 2 dp) for the same pair.
## No `R` field: the test reads the matrix from the exported `forbes2023` (M54).
amh <- list(
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
  matrix_md5 = "c1dd9eca009c2738c268487179d43e87",
  files = vapply(osf, function(f) sprintf("%s (guid %s) <%s>", f$name, f$guid, f$url), character(1)),
  note = paste0(
    "Expected comp_corr/cong/corr_chase computed with Forbes's ",
    "ExtendedBassAckwards reference implementation from the same md5-pinned ",
    "matrix now exported as data/forbes2023.rda. Only ackwards() runs at test time."
  ),
  generator = "data-raw/forbes2023.R",
  generated = as.character(Sys.Date()),
  R_version = R.version.string,
  psych = as.character(utils::packageVersion("psych"))
)

## Wrap in a named list mirroring forbes2023_sims.rds (single applied example).
## The structured provenance is also attached at the TOP level so the M57 guard
## test (test-oracle-provenance.R) finds a uniform `generator` + `source` contract
## on every fixture.
fixture <- list(amh = amh)
attr(fixture, "provenance") <- attr(amh, "provenance")

out <- file.path("tests", "testthat", "fixtures", "forbes2023_amh.rds")
saveRDS(fixture, out, compress = "xz")
message(sprintf(
  "Wrote data/forbes2023.rda and %s (%.0f KB)", out, file.size(out) / 1024
))
