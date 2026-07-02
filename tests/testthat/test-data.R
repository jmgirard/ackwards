test_that("bfi25 has the expected shape and structure", {
  expect_equal(dim(bfi25), c(1000L, 25L))
  expect_equal(
    colnames(bfi25),
    c(
      paste0("A", 1:5), paste0("C", 1:5), paste0("E", 1:5),
      paste0("N", 1:5), paste0("O", 1:5)
    )
  )
  expect_true(all(vapply(bfi25, is.integer, logical(1))))
  expect_true(anyNA(bfi25))
  expect_true(all(bfi25 >= 1L | is.na(bfi25)))
  expect_true(all(bfi25 <= 6L | is.na(bfi25)))
})

test_that("bfi25 ships public-domain IPIP labels that flow into top_items()", {
  skip_if_not_installed("psych")

  # Every column carries a non-empty scalar "label" attribute.
  labs <- vapply(
    bfi25,
    function(col) attr(col, "label", exact = TRUE) %||% NA_character_,
    character(1)
  )
  expect_equal(sum(!is.na(labs)), 25L)
  expect_true(all(nzchar(labs)))
  expect_identical(unname(labs[["E4"]]), "Make friends easily")

  # Fit on the raw dataset (labels survive; NAs handled by `missing`) -> the
  # labels are captured and surfaced by top_items() with zero user setup.
  x <- suppressWarnings(ackwards(bfi25, k_max = 3))
  expect_length(x$meta$item_labels, 25L)
  ti <- top_items(x, level = 3, cut = 0.4)
  expect_true("label" %in% names(ti$data))
  expect_true(any(ti$data$label == "Make friends easily", na.rm = TRUE))
})

test_that("sim16 has the expected shape and structure", {
  expect_equal(dim(sim16), c(1000L, 16L))
  expect_equal(colnames(sim16), paste0("i", 1:16))
  expect_true(all(vapply(sim16, is.double, logical(1))))
  expect_false(anyNA(sim16))
  # Continuous, not ordinal-looking: far more than 7 distinct values per column.
  expect_true(all(vapply(sim16, function(v) length(unique(v)) > 7L, logical(1))))
})

test_that("sim16 triggers no ordinal-detection warning under pca or efa", {
  skip_if_not_installed("psych")
  x_pca <- suppressMessages(ackwards(sim16, k_max = 4, engine = "pca"))
  x_efa <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa"))
  expect_false(x_pca$meta$ordinal_warned)
  expect_false(x_efa$meta$ordinal_warned)
})

test_that("sim16 recovers its known 1 -> 2 -> 4 hierarchy under efa", {
  skip_if_not_installed("psych")
  x <- suppressMessages(ackwards(sim16, k_max = 4, engine = "efa"))

  .primary_groups <- function(level) {
    load_mat <- x$levels[[as.character(level)]]$loadings
    unname(apply(abs(load_mat), 1, which.max))
  }

  # k = 1: a single general factor spans all 16 items.
  expect_equal(length(unique(.primary_groups(1))), 1L)

  # k = 2: splits exactly along the metatrait line (i1-i8 vs. i9-i16).
  g2 <- .primary_groups(2)
  expect_equal(length(unique(g2[1:8])), 1L)
  expect_equal(length(unique(g2[9:16])), 1L)
  expect_false(g2[1] == g2[9])

  # k = 4: recovers the 4 true group factors (i1-4, i5-8, i9-12, i13-16) exactly.
  g4 <- .primary_groups(4)
  true_groups <- rep(1:4, each = 4)
  expect_equal(length(unique(g4)), 4L)
  for (grp in 1:4) {
    idx <- which(true_groups == grp)
    expect_equal(length(unique(g4[idx])), 1L)
  }
})

test_that("suggest_k(sim16) converges on the true k = 4 (idealized case)", {
  skip_if_not_installed("psych")
  skip_if_not_installed("EFAtools")
  sk <- suggest_k(sim16, n_iter = 5L, seed = 1L)

  # sim16 is deliberately the *idealized* teaching case: a strong, clean planted
  # signal on which the criteria converge (contrast bfi25, where the same
  # criteria span k = 4-6). We assert that convergence robustly rather than
  # pinning all six criteria to 4 exactly -- PA-PC/PA-FA/CD resample, so their
  # exact value can wobble +/-1 across platforms/RNG even under a fixed seed.
  expect_equal(sk$k_map, 4L) # MAP is deterministic -> a stable anchor on k = 4

  ks <- c(
    sk$k_parallel_pc, sk$k_parallel_fa, sk$k_map,
    sk$k_vss1, sk$k_vss2, sk$k_cd
  )
  # The majority of criteria (the "consensus") lands on the true k = 4.
  expect_equal(as.integer(names(which.max(table(ks)))), 4L)
})

test_that("sim16 at k_max = 5 guarantees a redundant chain and an artifact signal", {
  skip_if_not_installed("psych")
  x <- suppressMessages(
    ackwards(sim16, k_max = 5, engine = "efa") |>
      prune(c("redundant", "artifact"))
  )

  # Redundancy: the true (non-splitting) factors persist across levels, so at
  # least one adjacent-level chain is flagged (|r| >= .9 and phi >= .95).
  expect_gte(sum(x$prune$nodes$pruned), 1L)
  expect_gte(nrow(subset(x$prune$phi, abs_phi >= 0.95)), 1L)

  # Artifact: k = 5 has no real 5th dimension, so the extra factor is an
  # orphan with too few primary-loading items.
  structural <- x$prune$structural
  expect_gte(sum(structural$few_items | structural$orphan), 1L)

  # Stronger: the artefact factor is a primary loader for *zero* items (the
  # doc's "orphan with zero primary-loading items"), and that same factor is
  # exactly the one flagged both few_items and orphan.
  load5 <- x$levels[["5"]]$loadings
  primary5 <- colnames(load5)[apply(abs(load5), 1, which.max)]
  counts5 <- table(factor(primary5, levels = colnames(load5)))
  expect_true(any(counts5 == 0L))
  orphan_ids <- names(counts5)[counts5 == 0L]
  flagged_ids <- subset(structural, few_items & orphan)$id
  expect_true(all(orphan_ids %in% flagged_ids))
})

test_that("data-raw/sim16.R regenerates sim16 identically (fixed-seed reproducibility)", {
  script_path <- testthat::test_path("..", "..", "data-raw", "sim16.R")
  skip_if(!file.exists(script_path), "data-raw/sim16.R not reachable from this test context")

  # Drop the trailing usethis::use_data() call (writes to disk; not what this
  # test is checking) without referencing usethis:: literally -- R CMD check's
  # unstated-dependency scan flags any pkg:: token in tests/, even in a string.
  exprs <- parse(script_path)
  is_use_data_call <- function(e) {
    is.call(e) && grepl("use_data", deparse(e[[1]])[1], fixed = TRUE)
  }
  exprs <- Filter(Negate(is_use_data_call), exprs)

  env <- new.env(parent = globalenv())
  for (e in exprs) eval(e, envir = env)

  expect_equal(env$sim16, sim16, ignore_attr = TRUE)
})
