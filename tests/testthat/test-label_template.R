test_that("label_template returns named character vector of correct length", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- label_template(x)
  expect_type(out, "character")
  # 1 + 2 + 3 = 6 factors total
  expect_length(out, 6L)
  expect_named(out)
})

test_that("label_template names match ba_layout node ids", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  layout_ids <- ba_layout(x)$nodes$id
  out <- label_template(x)
  expect_identical(names(out), layout_ids)
})

test_that("label_template style='id' is a round-trip no-op for node_labels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  labs <- label_template(x, style = "id")
  # values equal names
  expect_identical(unname(labs), names(labs))

  # Passing to autoplot should produce same node labels as default (no-op)
  lay_default <- ba_layout(x)$nodes$label
  lay_labeled <- ba_layout(x)$nodes
  lay_labeled$label[match(names(labs), lay_labeled$id)] <- labs[names(labs)]
  expect_identical(lay_default, lay_labeled$label)
})

test_that("label_template style='forbes' assigns letter+index values", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- label_template(x, style = "forbes")
  # Level 1: one factor -> "A1"
  expect_identical(out[["m1f1"]], "A1")
  # Level 2: two factors -> "B1", "B2"
  expect_true("B1" %in% out)
  expect_true("B2" %in% out)
  # Level 3: three factors -> "C1", "C2", "C3"
  expect_true("C1" %in% out)
  expect_true("C2" %in% out)
  expect_true("C3" %in% out)
  # All values are non-empty
  expect_true(all(nchar(out) > 0L))
})

test_that("label_template style='blank' returns all empty strings", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- label_template(x, style = "blank")
  expect_true(all(out == ""))
  expect_length(out, 6L)
})

test_that("label_template errors on non-ackwards input", {
  expect_error(label_template(list()), class = "rlang_error")
})

test_that("label_template errors on invalid style", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")
  expect_error(label_template(x, style = "bad"), class = "error")
})

test_that("label_template returns invisibly", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")
  # Should print but return the vector invisibly
  out <- withVisible(label_template(x))
  expect_false(out$visible)
  expect_type(out$value, "character")
})

test_that("label_template includes level-1 factor in all styles", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 2, engine = "pca")
  # Level 1 always has exactly one factor
  out_id <- label_template(x, style = "id")
  expect_true("m1f1" %in% names(out_id))
  expect_identical(out_id[["m1f1"]], "m1f1")

  out_forbes <- label_template(x, style = "forbes")
  expect_identical(out_forbes[["m1f1"]], "A1")

  out_blank <- label_template(x, style = "blank")
  expect_identical(out_blank[["m1f1"]], "")
})

test_that("label_template forbes uses LETTERS correctly for deeper hierarchies", {
  skip_if_not_installed("psych")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "pca")

  out <- label_template(x, style = "forbes")
  # Check level prefixes: names starting with m1 -> "A", m2 -> "B", m3 -> "C"
  m1_vals <- out[grep("^m1", names(out))]
  m2_vals <- out[grep("^m2", names(out))]
  m3_vals <- out[grep("^m3", names(out))]
  expect_true(all(grepl("^A", m1_vals)))
  expect_true(all(grepl("^B", m2_vals)))
  expect_true(all(grepl("^C", m3_vals)))
})

test_that("label_template works with engine = 'efa'", {
  skip_if_not_installed("psych")
  skip_if_not_installed("GPArotation")
  set.seed(1)
  x <- ackwards(.make_esem_data(), k_max = 3, engine = "efa")

  out <- label_template(x, style = "id")
  expect_type(out, "character")
  expect_length(out, 6L) # 1 + 2 + 3 factors
  expect_identical(names(out), ba_layout(x)$nodes$id)

  out_forbes <- label_template(x, style = "forbes")
  expect_identical(out_forbes[["m1f1"]], "A1")
})
