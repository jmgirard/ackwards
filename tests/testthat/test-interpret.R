test_that("top_items returns correct S3 class and structure", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- top_items(x)
  expect_s3_class(out, "top_items")
  expect_named(out, c("data", "levels_shown", "cut", "n", "sort", "engine", "k_max"))
  expect_true(is.data.frame(out$data))
  expect_named(out$data, c("level", "factor", "item", "loading"))
})

test_that("top_items cut filters correctly", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

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
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out1 <- top_items(x, level = 1)
  expect_equal(unique(out1$data$level), 1L)
  expect_equal(out1$levels_shown, 1L)

  out23 <- top_items(x, level = c(2, 3))
  expect_setequal(unique(out23$data$level), c(2L, 3L))
})

test_that("top_items errors on invalid level", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")
  expect_error(top_items(x, level = 99), class = "rlang_error")
})

test_that("top_items n caps items per factor", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- top_items(x, cut = 0, n = 2)
  counts <- tapply(out$data$item, out$data$factor, length)
  expect_true(all(counts <= 2L))
})

test_that("top_items sort = FALSE preserves item order", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")

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
  skip_if_not_installed("GPArotation")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "efa")

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
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

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
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")

  expect_error(top_items(x, cut = -0.1), class = "rlang_error")
  expect_error(top_items(x, cut = 1.5), class = "rlang_error")
  expect_error(top_items(x, n = 0), class = "rlang_error")
  expect_error(top_items(x, sort = "yes"), class = "rlang_error")
})

test_that("print.top_items runs without error", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")
  out <- top_items(x)
  expect_no_error(print(out))
  # returns invisibly
  expect_identical(print(out), out)
})

test_that("print.top_items handles zero-row result gracefully", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")
  out <- top_items(x, cut = 0.99)
  expect_no_error(print(out))
})

test_that("interpret vignette edges idiom: adjacent primary filter and column select", {
  # Guards the pattern used in ackwards-interpret.Rmd (edges-weak chunk):
  #   edges[edges$is_primary & edges$level_to == edges$level_from + 1, c("from","to","r")]
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

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
