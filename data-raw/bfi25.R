## data-raw/bfi25.R
##
## Creates the bfi25 example dataset: the first 25 item columns of psych::bfi,
## subsampled to 1000 rows with a fixed seed, with each item column carrying a
## human-readable "label" attribute (its IPIP item stem).
##
## Source: W. Revelle's {psych} package (>= 2.0), itself drawn from the SAPA
## (Synthetic Aperture Personality Assessment) project. The 25 items measure the
## Big Five personality domains via public-domain IPIP (International Personality
## Item Pool) markers developed by Goldberg (1999).
##
## Item labels: the 25 IPIP marker stems below are in the public domain
## (ipip.ori.org; Goldberg, 1999). They are sourced here directly as public-
## domain IPIP items -- the text is necessarily identical to psych::bfi.dictionary
## because both draw the same items from the same public pool, but this is not a
## copy of that object. Trailing periods are stripped so the stems read cleanly
## as variable labels (e.g. top_items() prints "Make friends easily (E4)").
##
## To regenerate: source this script from the package root via
##   source("data-raw/bfi25.R")

set.seed(42)

bfi_full <- psych::bfi[, 1:25] # first 25 IPIP items
row_idx <- sample(nrow(bfi_full), size = 1000L) # random 1000-row subset
bfi25 <- bfi_full[row_idx, , drop = FALSE]
rownames(bfi25) <- NULL # tidy up rownames

# Public-domain IPIP marker stems (Goldberg, 1999; ipip.ori.org), keyed by the
# item code. Order matches the 25 columns of bfi25.
ipip_labels <- c(
  A1 = "Am indifferent to the feelings of others",
  A2 = "Inquire about others' well-being",
  A3 = "Know how to comfort others",
  A4 = "Love children",
  A5 = "Make people feel at ease",
  C1 = "Am exacting in my work",
  C2 = "Continue until everything is perfect",
  C3 = "Do things according to a plan",
  C4 = "Do things in a half-way manner",
  C5 = "Waste my time",
  E1 = "Don't talk a lot",
  E2 = "Find it difficult to approach others",
  E3 = "Know how to captivate people",
  E4 = "Make friends easily",
  E5 = "Take charge",
  N1 = "Get angry easily",
  N2 = "Get irritated easily",
  N3 = "Have frequent mood swings",
  N4 = "Often feel blue",
  N5 = "Panic easily",
  O1 = "Am full of ideas",
  O2 = "Avoid difficult reading material",
  O3 = "Carry the conversation to a higher level",
  O4 = "Spend time reflecting on things",
  O5 = "Will not probe deeply into a subject"
)

stopifnot(identical(names(ipip_labels), names(bfi25)))
for (col in names(bfi25)) {
  attr(bfi25[[col]], "label") <- unname(ipip_labels[[col]])
}

usethis::use_data(bfi25, overwrite = TRUE, compress = "xz")
