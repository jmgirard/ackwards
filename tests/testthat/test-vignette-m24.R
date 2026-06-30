# Guard tests for vignette idioms introduced in M24:
# wide + delta comparison tables in engines/ordinal/forbes vignettes.

# Helper: primary loading per item at a given level
.primary_loading <- function(x, level_k, items) {
  d <- tidy(x, what = "loadings")
  d <- d[d$level == level_k & d$item %in% items, ]
  do.call(rbind, lapply(split(d, d$item), function(r) r[which.max(abs(r$loading)), ]))
}

test_that("engines vignette: loadings wide-pivot produces expected columns and shape", {
  # Guards loadings-compare chunk in ackwards-engines.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_pca <- ackwards(bfi, k_max = 3, cor = "polychoric")
  x_efa <- ackwards(bfi, k_max = 3, engine = "efa", cor = "polychoric")

  anchors <- c("N1", "N2", "E1", "E2", "C1", "C2")
  pca_p <- .primary_loading(x_pca, 3, anchors)
  efa_p <- .primary_loading(x_efa, 3, anchors)

  # Alignment assertion must pass
  pca_s <- pca_p[order(pca_p$item), ]
  efa_s <- efa_p[order(efa_p$item), ]
  expect_equal(pca_s$factor, efa_s$factor)
  expect_equal(sign(pca_s$loading), sign(efa_s$loading))

  wide <- merge(
    pca_p[, c("factor", "item", "loading")],
    efa_p[, c("item", "loading")],
    by = "item", suffixes = c("_pca", "_efa")
  )
  wide$delta <- abs(wide$loading_efa) - abs(wide$loading_pca)

  expect_named(wide, c("item", "factor", "loading_pca", "loading_efa", "delta"),
    ignore.order = TRUE
  )
  expect_equal(nrow(wide), length(anchors))
  expect_true(all(wide$delta < 0),
    info = "EFA |loading| should be uniformly smaller than PCA |loading| for these anchor items"
  )
})

test_that("engines vignette: edges wide-pivot produces expected columns and no NA mismatch", {
  # Guards edges-compare chunk in ackwards-engines.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_pca <- ackwards(bfi, k_max = 3, cor = "polychoric")
  x_efa <- ackwards(bfi, k_max = 3, engine = "efa", cor = "polychoric")

  pca_e <- tidy(x_pca, what = "edges", primary_only = TRUE)
  efa_e <- tidy(x_efa, what = "edges", primary_only = TRUE)

  wide_e <- merge(
    pca_e[, c("from", "to", "r")],
    efa_e[, c("from", "to", "r")],
    by = c("from", "to"), all = TRUE, suffixes = c("_pca", "_efa")
  )
  wide_e$delta <- wide_e$r_efa - wide_e$r_pca

  expect_named(wide_e, c("from", "to", "r_pca", "r_efa", "delta"), ignore.order = TRUE)
  expect_true(nrow(wide_e) > 0L)
  # On well-structured data, primary parents should agree between engines
  expect_false(anyNA(wide_e),
    info = "Expect no engine disagreement on primary parents for BFI at k=3"
  )
})

test_that("ordinal vignette: loadings wide-pivot produces expected columns and positive deltas", {
  # Guards load-compare chunk in ackwards-ordinal.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_pearson <- suppressWarnings(ackwards(bfi, k_max = 5))
  x_poly <- ackwards(bfi, k_max = 5, cor = "polychoric")

  n_items <- c("N1", "N2", "N3", "N4", "N5")
  pear_p <- .primary_loading(x_pearson, 5, n_items)
  poly_p <- .primary_loading(x_poly, 5, n_items)

  # Alignment assertion must pass
  pear_s <- pear_p[order(pear_p$item), ]
  poly_s <- poly_p[order(poly_p$item), ]
  expect_equal(pear_s$factor, poly_s$factor)
  expect_equal(sign(pear_s$loading), sign(poly_s$loading))

  wide <- merge(
    pear_p[, c("factor", "item", "loading")],
    poly_p[, c("item", "loading")],
    by = "item", suffixes = c("_pearson", "_poly")
  )
  wide$delta <- wide$loading_poly - wide$loading_pearson

  expect_named(wide, c("item", "factor", "loading_pearson", "loading_poly", "delta"),
    ignore.order = TRUE
  )
  expect_equal(nrow(wide), length(n_items))
  expect_true(all(wide$delta > 0),
    info = "Polychoric loadings should be uniformly larger for Neuroticism items"
  )
})

test_that("ordinal vignette: edges wide-pivot produces expected columns and no NA mismatch", {
  # Guards edge-compare chunk in ackwards-ordinal.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_pearson <- suppressWarnings(ackwards(bfi, k_max = 5))
  x_poly <- ackwards(bfi, k_max = 5, cor = "polychoric")

  pear_e <- tidy(x_pearson, what = "edges", primary_only = TRUE)
  poly_e <- tidy(x_poly, what = "edges", primary_only = TRUE)

  wide_e <- merge(
    pear_e[, c("from", "to", "r")],
    poly_e[, c("from", "to", "r")],
    by = c("from", "to"), all = TRUE, suffixes = c("_pearson", "_poly")
  )
  wide_e$delta <- wide_e$r_poly - wide_e$r_pearson

  expect_named(wide_e, c("from", "to", "r_pearson", "r_poly", "delta"), ignore.order = TRUE)
  expect_true(nrow(wide_e) > 0L)
  expect_false(anyNA(wide_e),
    info = "Expect no basis disagreement on primary parents for BFI at k=5"
  )
})

test_that("forbes vignette: prune-nodes table idiom returns expected columns and n_redundant", {
  # Guards prune-nodes chunk in ackwards-forbes.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_prune <- ackwards(bfi,
    k_max = 5, cor = "polychoric",
    pairs = "all", prune = "redundant"
  )

  nodes <- tidy(x_prune, what = "nodes")
  n_redundant <- sum(nodes$pruned)

  expect_named(nodes, c("id", "level", "pruned", "prune_reason"), ignore.order = TRUE)
  expect_true(n_redundant > 0L, info = "Expect some redundant nodes for BFI at k=5")
  expect_true(all(nodes$prune_reason[nodes$pruned] == "redundant"),
    info = "Flagged nodes should carry prune_reason = 'redundant'"
  )
  expect_true(all(is.na(nodes$prune_reason[!nodes$pruned])),
    info = "Unflagged nodes should have NA prune_reason"
  )
})

test_that("forbes vignette: skip-edge inline-R values are computed correctly", {
  # Guards the edge-counts and top-skip-edge inline R in ackwards-forbes.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_adj <- ackwards(bfi, k_max = 5, cor = "polychoric")
  x_all <- ackwards(bfi, k_max = 5, cor = "polychoric", pairs = "all")

  n_adj <- nrow(tidy(x_adj, what = "edges"))
  n_all <- nrow(tidy(x_all, what = "edges"))
  expect_true(n_adj > 0L)
  expect_true(n_all > n_adj)

  edges <- tidy(x_all, what = "edges", sort = "strength")
  skip <- edges[
    abs(edges$level_to - edges$level_from) > 1 & abs(edges$r) >= 0.5,
    c("from", "to", "level_from", "level_to", "r")
  ]
  expect_true(nrow(skip) > 0L)

  top_r <- round(skip$r[1], 2)
  expect_true(abs(top_r) >= 0.5)
  expect_true(abs(skip$level_to[1] - skip$level_from[1]) > 1)
})
