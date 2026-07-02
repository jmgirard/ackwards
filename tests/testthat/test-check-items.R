# check_items(): pre-analysis item screen + the shared internal screen that
# ackwards() runs (error on constant, warn on near-constant).

.mk_items <- function(n = 300, seed = 1) {
  set.seed(seed)
  mk <- function(p) sample(seq_along(p), n, replace = TRUE, prob = p)
  d <- as.data.frame(lapply(1:5, function(i) mk(rep(0.2, 5))))
  names(d) <- paste0("i", 1:5)
  d
}

test_that("check_items() returns a per-item report with the expected shape", {
  d <- .mk_items()
  out <- check_items(d, cor = "polychoric")
  expect_s3_class(out, "check_items")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(d))
  expect_named(
    out,
    c("item", "n_valid", "pct_missing", "n_distinct", "min_count", "top_prop", "flag")
  )
  expect_identical(attr(out, "basis"), "polychoric")
  # clean five-category items -> nothing flagged
  expect_true(all(out$flag == "ok"))
})

test_that("check_items() flags constant, near-constant, and sparse items", {
  n <- 300
  d <- .mk_items(n = n)
  d$const <- rep(3L, n) # no variance
  d$nearconst <- c(rep(2L, n - 2L), 3L, 4L) # 298,1,1 -> dominant category
  d$sparse <- c(rep(1L, 80), rep(2L, 80), rep(3L, 80), rep(4L, 57), rep(5L, 3))

  out <- check_items(d, cor = "polychoric")
  fl <- setNames(out$flag, out$item)
  expect_identical(unname(fl["const"]), "constant")
  expect_identical(unname(fl["nearconst"]), "near-constant")
  expect_identical(unname(fl["sparse"]), "sparse category")
})

test_that("check_items() sparse flag is polychoric-specific", {
  n <- 300
  d <- .mk_items(n = n)
  d$sparse <- c(rep(1L, 80), rep(2L, 80), rep(3L, 80), rep(4L, 57), rep(5L, 3))
  # Under pearson, a sparse category is not a problem -> not flagged.
  out_p <- check_items(d, cor = "pearson")
  expect_identical(out_p$flag[out_p$item == "sparse"], "ok")
  # A constant item is still flagged on any basis.
  d$const <- rep(2L, n)
  out_p2 <- check_items(d, cor = "pearson")
  expect_identical(out_p2$flag[out_p2$item == "const"], "constant")
})

test_that("check_items() flags high missingness", {
  n <- 300
  d <- .mk_items(n = n)
  d$holey <- d$i1
  d$holey[1:100] <- NA # 33% missing
  out <- check_items(d, cor = "pearson")
  expect_identical(out$flag[out$item == "holey"], "high missing")
  expect_gt(out$pct_missing[out$item == "holey"], 0.3)
})

test_that("check_items() print runs for clean and flagged inputs", {
  d <- .mk_items()
  expect_no_error(print(check_items(d, cor = "pearson")))
  # Flag both a constant (red) and a near-constant (yellow) item so both print
  # branches are exercised.
  d$const <- rep(1L, nrow(d))
  d$nc <- c(rep(2L, nrow(d) - 2L), 3L, 4L)
  expect_no_error(print(check_items(d)))
})

test_that("check_items() validates its input", {
  expect_error(check_items(1:10), "data frame or numeric matrix")
  expect_error(check_items(data.frame()), "no columns")
})

test_that("ackwards() errors on a constant item, naming it", {
  skip_if_not_installed("psych")
  n <- 300
  d <- .mk_items(n = n)
  d$dead <- rep(3L, n)
  expect_error(
    suppressWarnings(ackwards(d, k_max = 3, engine = "efa", cor = "polychoric", correct = 0)),
    "dead"
  )
  expect_error(
    suppressWarnings(ackwards(d, k_max = 3, engine = "efa", cor = "polychoric", correct = 0)),
    "no variance"
  )
})

test_that("ackwards() warns on a near-constant item but still builds", {
  skip_if_not_installed("psych")
  n <- 300
  d <- .mk_items(n = n)
  d$nc <- c(rep(2L, n - 2L), 3L, 4L)
  rlang::reset_warning_verbosity("ackwards_degenerate_items")
  # A near-constant item also drives the matrix non-positive-definite (smoothing)
  # and psych emits its own notes; capture all warnings so none leak, and assert
  # the near-constant advisory is among them.
  w <- testthat::capture_warnings(
    x <- suppressMessages(
      ackwards(d, k_max = 3, engine = "efa", cor = "polychoric", correct = 0)
    )
  )
  expect_true(any(grepl("near-constant", w)))
  expect_s3_class(x, "ackwards")
})

test_that("ackwards() raises no item warning on clean data", {
  skip_if_not_installed("psych")
  expect_no_warning(
    suppressMessages(ackwards(sim16, k_max = 3, engine = "efa")),
    message = "near-constant|no variance"
  )
})
