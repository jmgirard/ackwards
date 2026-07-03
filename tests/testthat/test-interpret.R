test_that("top_items returns correct S3 class and structure", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  out <- top_items(x)
  expect_s3_class(out, "top_items")
  expect_named(out, c(
    "data", "levels_shown", "cut", "n", "sort", "by",
    "item_labels", "factor_labels", "engine", "k_max"
  ))
  expect_true(is.data.frame(out$data))
  # No labels on this fixture -> no label column, default grouping is by factor
  expect_named(out$data, c("level", "factor", "item", "loading"))
  expect_identical(out$by, "factor")
  expect_null(out$item_labels)
  expect_null(out$factor_labels)
})

test_that("top_items cut filters correctly", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  out <- top_items(x, cut = 0.3)
  expect_true(all(abs(out$data$loading) >= 0.3))

  out_high <- top_items(x, cut = 0.99)
  expect_equal(nrow(out_high$data), 0L)

  out_zero <- top_items(x, cut = 0.0)
  # With cut = 0 we should get all items across all factors (may duplicate items
  # that load on multiple factors)
  expect_true(nrow(out_zero$data) > 0L)
})

test_that("top_items level argument subsets levels", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  out1 <- top_items(x, level = 1)
  expect_equal(unique(out1$data$level), 1L)
  expect_equal(out1$levels_shown, 1L)

  out23 <- top_items(x, level = c(2, 3))
  expect_setequal(unique(out23$data$level), c(2L, 3L))
})

test_that("top_items errors on invalid level", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))
  expect_error(top_items(x, level = 99), class = "rlang_error")
})

test_that("top_items n caps items per factor", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  out <- top_items(x, cut = 0, n = 2)
  counts <- tapply(out$data$item, out$data$factor, length)
  expect_true(all(counts <= 2L))
})

test_that("top_items sort = FALSE preserves item order", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 2, engine = "pca"))

  out_unsorted <- top_items(x, cut = 0, sort = FALSE)
  full_loadings <- tidy(x, what = "loadings")

  for (lev in unique(out_unsorted$data$level)) {
    for (fac in unique(out_unsorted$data$factor[out_unsorted$data$level == lev])) {
      items_unsorted <- out_unsorted$data$item[
        out_unsorted$data$level == lev & out_unsorted$data$factor == fac
      ]
      # tidy() preserves loading-matrix row order; unsorted must match that order
      items_tidy_order <- full_loadings$item[
        full_loadings$level == lev & full_loadings$factor == fac
      ]
      expected_order <- items_tidy_order[items_tidy_order %in% items_unsorted]
      expect_identical(items_unsorted, expected_order)
    }
  }
})

test_that("top_items works with engine = 'efa'", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "efa"))

  out <- top_items(x)
  expect_s3_class(out, "top_items")
  expect_identical(out$engine, "efa")
  # Filtering and structure work identically to PCA
  expect_true(all(abs(out$data$loading) >= 0.3))
  expect_named(out$data, c("level", "factor", "item", "loading"))
})

test_that("top_items data matches tidy(what = 'loadings')", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  full_loadings <- tidy(x, what = "loadings")
  out <- top_items(x, cut = 0, sort = FALSE)

  # Every row in top_items$data should match a row in tidy loadings
  for (i in seq_len(nrow(out$data))) {
    row <- out$data[i, ]
    match_row <- full_loadings[
      full_loadings$level == row$level &
        full_loadings$factor == row$factor &
        full_loadings$item == row$item, ,
      drop = FALSE
    ]
    expect_equal(nrow(match_row), 1L)
    expect_equal(row$loading, match_row$loading, tolerance = 1e-10)
  }
})

test_that("top_items errors on non-ackwards input", {
  expect_error(top_items(list()), class = "rlang_error")
  expect_error(top_items("string"), class = "rlang_error")
})

test_that("top_items errors on bad cut/n/sort", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 2, engine = "pca"))

  expect_error(top_items(x, cut = -0.1), class = "rlang_error")
  expect_error(top_items(x, cut = 1.5), class = "rlang_error")
  expect_error(top_items(x, n = 0), class = "rlang_error")
  expect_error(top_items(x, sort = "yes"), class = "rlang_error")
})

# ── by = "item" grouping + variable labels (M36) ──────────────────────────────

test_that("top_items(by = 'item') keeps the same rows but inverts grouping", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  by_fac <- top_items(x, level = 3, cut = 0.25)
  by_item <- top_items(x, level = 3, cut = 0.25, by = "item")
  expect_identical(by_item$by, "item")
  # Same underlying (level, factor, item, loading) set, just reordered
  key <- function(d) {
    paste(d$level, d$factor, d$item, round(d$loading, 8))
  }
  expect_setequal(key(by_fac$data), key(by_item$data))
})

test_that("top_items(by = 'item', n = ) caps factors per item", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  out <- top_items(x, cut = 0, n = 2, by = "item")
  counts <- tapply(out$data$factor, list(out$data$level, out$data$item), length)
  expect_true(all(counts[!is.na(counts)] <= 2L))
})

test_that("top_items exposes captured variable labels", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- .make_esem_data()
  attr(d$x1, "label") <- "First indicator"
  attr(d$x4, "label") <- "Fourth indicator"
  x <- cached(ackwards(d, k_max = 3, engine = "pca"))

  out <- top_items(x, cut = 0)
  expect_true("label" %in% names(out$data))
  expect_identical(out$item_labels, c(x1 = "First indicator", x4 = "Fourth indicator"))
  # Labelled items carry their text; unlabelled ones are NA in the column
  expect_identical(out$data$label[out$data$item == "x1"][1L], "First indicator")
  expect_true(all(is.na(out$data$label[out$data$item == "x2"])))
})

test_that("top_items(show_labels = FALSE) drops the label column and text", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- .make_esem_data()
  attr(d$x1, "label") <- "First indicator"
  x <- cached(ackwards(d, k_max = 3, engine = "pca"))

  out <- top_items(x, cut = 0, show_labels = FALSE)
  expect_false("label" %in% names(out$data))
  expect_null(out$item_labels)
})

test_that(".format_item_label formats 'id: label' with bare-id fallback", {
  fmt <- ackwards:::.format_item_label
  labs <- c(x1 = "First indicator")
  expect_identical(fmt("x1", labs), "x1: First indicator")
  expect_identical(fmt("x2", labs), "x2") # no label -> bare id
  expect_identical(fmt("x1", NULL), "x1") # no labels at all -> bare id
})

test_that("top_items errors on bad by/show_labels", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 2, engine = "pca"))

  expect_error(top_items(x, by = "column"), class = "rlang_error")
  expect_error(top_items(x, show_labels = "yes"), class = "rlang_error")
})

test_that("print.top_items(by = 'item') and label display run without error", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- .make_esem_data()
  attr(d$x1, "label") <- "First indicator"
  x <- cached(ackwards(d, k_max = 3, engine = "pca"))

  expect_no_error(print(top_items(x, level = 3, cut = 0.25, by = "item")))
  expect_no_error(print(top_items(x, level = 3, cut = 0.25))) # id: label path
})

test_that("print.top_items renders 'code: label' in header (by=item) and body (by=factor)", {
  skip_if_not_installed("psych")
  set.seed(1)
  d <- .make_esem_data()
  attr(d$x1, "label") <- "First indicator"
  x <- cached(ackwards(d, k_max = 3, engine = "pca"))

  # by = "item": the item is the group HEADER, so the label appears there.
  by_item <- cli::cli_fmt(print(top_items(x, level = 2, cut = 0, by = "item")))
  expect_true(any(grepl("x1: First indicator", by_item, fixed = TRUE)))

  # by = "factor": the item is a BODY entry under each factor.
  by_factor <- cli::cli_fmt(print(top_items(x, level = 2, cut = 0)))
  expect_true(any(grepl("x1: First indicator", by_factor, fixed = TRUE)))

  # show_labels = FALSE prints the bare code, never the label text.
  bare <- cli::cli_fmt(print(top_items(x, level = 2, cut = 0, show_labels = FALSE)))
  expect_false(any(grepl("First indicator", bare, fixed = TRUE)))
  expect_true(any(grepl("x1", bare, fixed = TRUE)))
})

test_that("print.top_items runs without error", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))
  out <- top_items(x)
  expect_no_error(print(out))
  # returns invisibly
  expect_identical(print(out), out)
})

test_that("print.top_items handles zero-row result gracefully", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 2, engine = "pca"))
  out <- top_items(x, cut = 0.99)
  expect_no_error(print(out))
})

test_that("interpret vignette edges idiom: adjacent primary filter and column select", {
  # Guards the pattern used in ackwards-interpret.Rmd (edges-weak chunk):
  #   edges[edges$is_primary & edges$level_to == edges$level_from + 1, c("from","to","r")]
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  edges <- tidy(x, what = "edges")
  primary <- edges[edges$is_primary & edges$level_to == edges$level_from + 1, ]
  expect_true(nrow(primary) > 0L)
  expect_true(all(primary$level_to == primary$level_from + 1L))
  expect_true(all(primary$is_primary))

  # Column-select exactly as the vignette does
  sub <- primary[order(abs(primary$r)), c("from", "to", "r")]
  expect_named(sub, c("from", "to", "r"))
  expect_true(nrow(sub) > 0L)
})

test_that("interpret vignette top_items idioms run", {
  # Guards the top_items() signatures shown in ackwards-interpret.Rmd.
  skip_if_not_installed("psych")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  expect_no_error(top_items(x, level = 3))
  expect_no_error(top_items(x, level = 3, cut = 0.45))
  expect_no_error(top_items(x, level = 3, cut = 0.3, n = 4))
  expect_no_error(top_items(x, level = 2, cut = 0.25)) # cross-loadings idiom

  ti <- top_items(x, level = 3, cut = 0.3)
  expect_true(is.data.frame(ti$data)) # $data access shown in vignette
})

test_that("print.top_items skips levels with no items (nrow(df_k) == 0)", {
  # Covers the `if (nrow(df_k) == 0L) next` branch. Rather than rely on a cut
  # that may or may not empty a level on a given dataset (the old version
  # skipped on bfi25), build a real top_items object and deterministically
  # drop one shown level's rows: df stays non-empty (level 3 remains), so the
  # early return doesn't fire, but k = 1 hits `next`.
  skip_if_not_installed("psych")
  x <- cached(ackwards(bfi25[, 1:6], k_max = 3))
  ti <- top_items(x, level = c(1L, 3L), cut = 0.3)
  ti$data <- ti$data[ti$data$level != 1L, , drop = FALSE]
  # Preconditions for the target branch:
  expect_true(1L %in% ti$levels_shown) # level 1 still requested...
  expect_false(any(ti$data$level == 1L)) # ...but has no rows -> `next`
  expect_gt(nrow(ti$data), 0L) # df non-empty -> early return skipped
  expect_no_error(print(ti))
  expect_invisible(print(ti))
})

test_that("interpret/visualization vignette idiom: label_template feeds autoplot cleanly", {
  # Guards the round-trip autoplot(x, node_labels = label_template(x, style)).
  # Every ID label_template emits must resolve to a real factor, so no
  # "match no factor ID" warning should fire for any style.
  skip_if_not_installed("psych")
  skip_if_not_installed("ggplot2")
  set.seed(1)
  x <- cached(ackwards(.make_esem_data(), k_max = 3, engine = "pca"))

  for (style in c("id", "forbes", "blank")) {
    labs <- suppressMessages(label_template(x, style = style))
    p <- expect_no_warning(ggplot2::autoplot(x, node_labels = labs))
    expect_s3_class(p, "ggplot")
  }

  # Blank-slate idiom: fill specific IDs, pass back (vignette label-blank chunk)
  labs <- suppressMessages(label_template(x, style = "blank"))
  labs["m3f1"] <- "Alpha"
  p <- expect_no_warning(ggplot2::autoplot(x, node_labels = labs))
  expect_s3_class(p, "ggplot")
})
