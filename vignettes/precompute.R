# Precompute expensive vignettes for CRAN.
#
# The heaviest vignettes refit `ackwards()` on the full bfi25 polychoric basis
# (and run parallel analysis, split-half comparability, and bootstrap edges) on
# every build. CRAN flagged the resulting re-build time (~317s) as too costly to
# check regularly. To fix this without dumbing down the examples, those vignettes
# are authored as `*.Rmd.orig` and knitted *here*, ahead of time, into plain
# `*.Rmd` with all results and figures baked in. R then builds those `*.Rmd`
# with pandoc only (near-zero cost) -- the R code never re-runs on CRAN.
#
# This is the rOpenSci precompute pattern:
#   https://ropensci.org/blog/2019/12/08/precompute-vignettes/
#
# Vignettes that stay live (built normally on CRAN) have no `.Rmd.orig`:
#   - ackwards-interpret.Rmd (bfi25 subset; needs IPIP labels; already ~2s)
#
# sim16 was considered for the live set but degrades the two candidates: the
# girard vignette's comparability lesson needs realistically messy data, and the
# visualization reference needs above-threshold negative edges to demo the
# sign-encoding aesthetics -- sim16's clean planted structure has neither.
#
# USAGE: run from the package root after editing any `*.Rmd.orig` AND before
# every release (the generated `*.Rmd` + `figure/` must be committed):
#
#   Rscript vignettes/precompute.R
#
# Requires the Suggests stack (lavaan, ggplot2, gt, EFAtools) so the ESEM,
# plotting, and CD chunks bake in real output.

stopifnot(file.exists("DESCRIPTION")) # must run from package root

devtools::load_all(".", quiet = TRUE)

old_wd <- setwd("vignettes")
on.exit(setwd(old_wd), add = TRUE)

orig_files <- list.files(pattern = "\\.Rmd\\.orig$")
if (length(orig_files) == 0L) {
  stop("No *.Rmd.orig files found in vignettes/.")
}

for (f in orig_files) {
  base <- sub("\\.Rmd\\.orig$", "", f)
  message("Precomputing ", f, " -> ", base, ".Rmd")
  # Unique figure prefix per vignette so the shared assets/ dir never collides
  # on common chunk labels (e.g. "setup"). The directory is deliberately NOT
  # named "figure"/"figures": R CMD check flags those as leftover knitr output.
  knitr::opts_chunk$set(fig.path = file.path("assets", paste0(base, "-")))
  knitr::knit(f, output = paste0(base, ".Rmd"), quiet = TRUE)
}

message("Done. Commit the generated *.Rmd and vignettes/assets/.")
