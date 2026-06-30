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
  wide_e$delta <- abs(wide_e$r_efa) - abs(wide_e$r_pca)

  expect_named(wide_e, c("from", "to", "r_pca", "r_efa", "delta"), ignore.order = TRUE)
  expect_true(nrow(wide_e) > 0L)
  expect_true(all(is.finite(wide_e$delta)))
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
  wide$delta <- abs(wide$loading_poly) - abs(wide$loading_pearson)

  expect_named(wide, c("item", "factor", "loading_pearson", "loading_poly", "delta"),
    ignore.order = TRUE
  )
  expect_equal(nrow(wide), length(n_items))
  expect_true(all(wide$delta > 0),
    info = "Polychoric |loading| should be uniformly larger for Neuroticism items"
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
  wide_e$delta <- abs(wide_e$r_poly) - abs(wide_e$r_pearson)

  expect_named(wide_e, c("from", "to", "r_pearson", "r_poly", "delta"), ignore.order = TRUE)
  expect_true(nrow(wide_e) > 0L)
  expect_true(all(is.finite(wide_e$delta)))
  expect_false(anyNA(wide_e),
    info = "Expect no basis disagreement on primary parents for BFI at k=5"
  )

  # Magnitude delta keeps the "positive = stronger connection" reading correct
  # even for the negatively-signed primary edge (m2f2 -> m3f2), which strengthens
  # under the polychoric basis. A signed delta would render that row negative.
  neg_edge <- wide_e[wide_e$r_pearson < 0 & wide_e$r_poly < 0, ]
  if (nrow(neg_edge) > 0L) {
    stronger <- abs(neg_edge$r_poly) > abs(neg_edge$r_pearson)
    expect_true(all(neg_edge$delta[stronger] > 0),
      info = "A strengthening negative edge must show a positive magnitude delta"
    )
  }
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

  # Inline-derivation helpers used in the prose (flagged_ids + lvl_summary)
  flagged_ids <- nodes$id[nodes$pruned]
  expect_equal(length(flagged_ids), n_redundant)

  lvl_summary <- paste(
    vapply(
      sort(unique(nodes$level[nodes$pruned])),
      function(L) {
        ids_L <- nodes$id[nodes$pruned & nodes$level == L]
        if (length(ids_L) == sum(nodes$level == L)) {
          paste0("the entire k = ", L, " level")
        } else {
          paste0(length(ids_L), " factor at k = ", L)
        }
      },
      character(1)
    ),
    collapse = ", "
  )
  expect_type(lvl_summary, "character")
  expect_true(nzchar(lvl_summary))
  # On BFI k=5 the whole k=4 level is flagged, so the phrase must appear
  expect_match(lvl_summary, "the entire k = 4 level", fixed = TRUE)
})

test_that("forbes vignette: prune-artefact table idiom returns expected columns", {
  # Guards the prune-artefact-nodes chunk in ackwards-forbes.Rmd
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_art <- ackwards(bfi,
    k_max = 5, cor = "polychoric",
    pairs = "all", prune = "artefact"
  )

  art_nodes <- tidy(x_art, what = "nodes")
  n_artefact <- sum(art_nodes$pruned)

  expect_named(art_nodes, c("id", "level", "pruned", "prune_reason"), ignore.order = TRUE)
  expect_type(n_artefact, "integer")
  # Artefact is never auto-flagged (DESIGN §14.21); BFI flags nothing
  expect_equal(n_artefact, 0L)
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

test_that("vignette table idiom: knitr::kable fallback renders without gt", {
  # Guards the `else { knitr::kable(...) }` branch of every comparison chunk.
  # gt is installed in dev/CI so that branch never runs during vignette builds;
  # this exercises the fallback idiom directly so it cannot rot unnoticed.
  skip_if_not_installed("psych")
  bfi <- na.omit(bfi25)
  x_pca <- ackwards(bfi, k_max = 3, cor = "polychoric")
  x_efa <- ackwards(bfi, k_max = 3, engine = "efa", cor = "polychoric")

  pca_e <- tidy(x_pca, what = "edges", primary_only = TRUE)
  efa_e <- tidy(x_efa, what = "edges", primary_only = TRUE)
  wide_e <- merge(pca_e[, c("from", "to", "r")], efa_e[, c("from", "to", "r")],
    by = c("from", "to"), all = TRUE, suffixes = c("_pca", "_efa")
  )
  wide_e$delta <- abs(wide_e$r_efa) - abs(wide_e$r_pca)

  out <- knitr::kable(
    wide_e,
    digits = 2, row.names = FALSE,
    col.names = c("From", "To", "PCA", "EFA", "Δ (|EFA| - |PCA|)")
  )
  expect_s3_class(out, "knitr_kable")
  expect_true(length(out) > 0L)
})

test_that("vignette edges idiom: disagreeing primary parents surface as NA in the merge", {
  # Guards the documented behaviour behind the edge-table footnotes: when two
  # engines/bases pick a different primary parent for the same child, the
  # all = TRUE merge yields NA rather than silently dropping the row. BFI never
  # triggers this, so demonstrate it on synthetic edge frames.
  a <- data.frame(from = c("m1f1", "m1f1"), to = c("m2f1", "m2f2"), r = c(0.9, 0.4))
  b <- data.frame(from = c("m1f1", "m1f2"), to = c("m2f1", "m2f2"), r = c(0.8, 0.5))

  wide <- merge(a, b,
    by = c("from", "to"), all = TRUE, suffixes = c("_a", "_b")
  )
  wide$delta <- abs(wide$r_b) - abs(wide$r_a)

  # m2f2 has from = m1f1 in `a` but from = m1f2 in `b`: two rows, each half-NA
  disagreement <- wide[wide$to == "m2f2", ]
  expect_equal(nrow(disagreement), 2L)
  expect_true(anyNA(disagreement$r_a))
  expect_true(anyNA(disagreement$r_b))
})
