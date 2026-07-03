test_that("set_factor_labels / factor_labels round-trip and validate", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))

  # No labels by default
  expect_null(factor_labels(x))

  # Set some labels; getter round-trips
  x2 <- set_factor_labels(x, c(m4f1 = "Alpha", m4f2 = "Beta"))
  expect_identical(factor_labels(x2), c(m4f1 = "Alpha", m4f2 = "Beta"))
  expect_s3_class(x2, "ackwards")

  # Original object is untouched (functional API, copy-on-modify)
  expect_null(factor_labels(x))

  # Repeated calls merge/update
  x3 <- set_factor_labels(x2, c(m4f1 = "Alpha2", m1f1 = "General"))
  expect_identical(
    factor_labels(x3),
    c(m4f1 = "Alpha2", m4f2 = "Beta", m1f1 = "General")
  )

  # NA / "" removes just that one
  x4 <- set_factor_labels(x3, c(m4f2 = NA))
  expect_false("m4f2" %in% names(factor_labels(x4)))
  expect_true("m4f1" %in% names(factor_labels(x4)))
  x5 <- set_factor_labels(x4, c(m4f1 = ""))
  expect_false("m4f1" %in% names(factor_labels(x5)))

  # NULL clears all
  x6 <- set_factor_labels(x3, NULL)
  expect_null(factor_labels(x6))
})

test_that("set_factor_labels errors on unknown IDs and bad input", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))

  expect_error(set_factor_labels(x, c(m9f9 = "Nope")), "no factor ID")
  expect_error(set_factor_labels(x, c("Unnamed")), "named")
  expect_error(set_factor_labels(x, 1:3), "character")
  expect_error(set_factor_labels(42, c(m1f1 = "x")), "ackwards")
  expect_error(factor_labels(42), "ackwards")
})

test_that("labels survive prune() and ride in meta", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))
  x <- set_factor_labels(x, c(m4f1 = "Alpha"))
  xp <- prune(x, rules = "redundant")
  expect_identical(factor_labels(xp), c(m4f1 = "Alpha"))
})

test_that("tidy() adds label columns only when labels are set", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))

  # Unlabeled: no label columns at all
  expect_false("factor_label" %in% names(tidy(x, what = "loadings")))
  expect_false("factor_label" %in% names(tidy(x, what = "variance")))
  expect_false("from_label" %in% names(tidy(x, what = "edges")))

  x <- set_factor_labels(x, c(m4f1 = "Alpha", m2f1 = "Broad"))

  ld <- tidy(x, what = "loadings")
  expect_true("factor_label" %in% names(ld))
  expect_equal(unique(ld$factor_label[ld$factor == "m4f1"]), "Alpha")
  expect_true(all(is.na(ld$factor_label[ld$factor == "m4f2"])))

  vr <- tidy(x, what = "variance")
  expect_true("factor_label" %in% names(vr))
  expect_equal(vr$factor_label[vr$factor == "m2f1"], "Broad")

  ed <- tidy(x, what = "edges")
  expect_true(all(c("from_label", "to_label") %in% names(ed)))
  # m2f1 labeled -> appears as a label on whichever endpoint it is
  expect_true(any(ed$from_label == "Broad", na.rm = TRUE) ||
    any(ed$to_label == "Broad", na.rm = TRUE))
})

test_that("tidy(scores) gets a factor_label column when set", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca", keep_scores = TRUE))
  x <- set_factor_labels(x, c(m4f1 = "Alpha"))
  sc <- tidy(x, what = "scores")
  expect_true("factor_label" %in% names(sc))
  expect_equal(unique(sc$factor_label[sc$factor == "m4f1"]), "Alpha")
})

test_that("print / summary show label (id) when set", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))
  x <- set_factor_labels(x, c(m4f1 = "Alpha"))

  print_out <- cli::cli_fmt(print(x))
  expect_true(any(grepl("Factor labels", print_out)))

  sum_out <- cli::cli_fmt(print(summary(x)))
  expect_true(any(grepl("Alpha \\(m4f1\\)", sum_out)))
})

test_that("top_items(by = 'factor') headers show label (id)", {
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))
  x <- set_factor_labels(x, c(m4f1 = "Alpha"))
  out <- cli::cli_fmt(print(top_items(x, level = 4)))
  expect_true(any(grepl("Alpha \\(m4f1\\)", out)))
})

test_that("autoplot uses stored labels as baseline, call-time overrides", {
  skip_if_not_installed("ggplot2")
  x <- cached(ackwards(sim16, k_max = 4, engine = "pca"))
  x <- set_factor_labels(x, c(m4f1 = "Alpha", m4f2 = "Beta"))

  p <- ggplot2::autoplot(x)
  lab_layer <- p$layers[[which(vapply(
    p$layers, function(l) inherits(l$geom, "GeomText"), logical(1)
  ))[[1]]]]
  # Build the node label data the plot draws from
  lay <- ba_layout(x)$nodes
  # Stored labels appear in the layout-derived text via autoplot's baseline
  built <- ggplot2::ggplot_build(p)
  txt <- unlist(lapply(built$data, function(d) if ("label" %in% names(d)) d$label))
  expect_true("Alpha" %in% txt)
  expect_true("Beta" %in% txt)

  # Call-time node_labels override wins for that node
  p2 <- ggplot2::autoplot(x, node_labels = c(m4f1 = "Override"))
  built2 <- ggplot2::ggplot_build(p2)
  txt2 <- unlist(lapply(built2$data, function(d) if ("label" %in% names(d)) d$label))
  expect_true("Override" %in% txt2)
  expect_false("Alpha" %in% txt2)
  expect_true("Beta" %in% txt2) # untouched node keeps stored label
})
