## data-raw/bfi25.R
##
## Creates the bfi25 example dataset: the first 25 item columns of psych::bfi,
## subsampled to 1000 rows with a fixed seed.
##
## Source: W. Revelle's {psych} package (>= 2.0), itself drawn from the SAPA
## (Synthetic Aperture Personality Assessment) project. The 25 items measure the
## Big Five personality domains via public-domain IPIP (International Personality
## Item Pool) markers developed by Goldberg (1999).
##
## To regenerate: source this script from the package root via
##   source("data-raw/bfi25.R")

set.seed(42)

bfi_full <- psych::bfi[, 1:25] # first 25 IPIP items
row_idx <- sample(nrow(bfi_full), size = 1000L) # random 1000-row subset
bfi25 <- bfi_full[row_idx, , drop = FALSE]
rownames(bfi25) <- NULL # tidy up rownames

usethis::use_data(bfi25, overwrite = TRUE, compress = "xz")
