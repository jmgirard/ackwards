# tests/testthat/test-predict.R
# M45: predict.ackwards() -- out-of-sample scoring front door

test_that("predict() is identical to augment(append = FALSE)", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  train <- d[1:500, ]
  test <- d[501:nrow(d), ]
  suppressWarnings(x <- ackwards(train, k_max = 3))

  expect_identical(
    predict(x, test),
    augment(x, data = test, append = FALSE)
  )
  expect_identical(
    predict(x, test, scaling = "sample"),
    augment(x, data = test, append = FALSE, scaling = "sample")
  )
})

test_that("predict() returns one row per newdata row, score columns only", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d[1:400, ], k_max = 3))
  sc <- predict(x, d[401:450, ])
  expect_s3_class(sc, "data.frame")
  expect_equal(nrow(sc), 50L)
  expect_named(sc, paste0(".", unlist(lapply(x$levels, `[[`, "labels"))))
})

test_that("predict() requires newdata and inherits augment's validation", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 2))
  expect_error(predict(x), "newdata")
  expect_error(predict(x, NULL), "newdata")
  expect_error(predict(x, d[, 1:3]), "missing")
  # cor-matrix objects carry no fit-time moments: default scaling errors,
  # scaling = "sample" works
  suppressWarnings(suppressMessages(xr <- ackwards(cor(d), k_max = 2)))
  expect_error(predict(xr, d), "fit-time item means")
  expect_no_error(predict(xr, d, scaling = "sample"))
})

test_that("predict(scaling='fit') is engine-agnostic: EFA value-level check", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  train <- d[1:500, ]
  test <- d[501:550, ]
  suppressWarnings(x <- ackwards(train, k_max = 3, engine = "efa"))
  sc <- predict(x, test)

  mu <- colMeans(as.matrix(train))
  sg <- apply(as.matrix(train), 2, stats::sd)
  Z <- sweep(sweep(as.matrix(test), 2, mu, "-"), 2, sg, "/")
  W <- x$levels[["3"]]$scoring$weights
  S <- sweep(Z %*% W, 2, sqrt(x$levels[["3"]]$scoring$score_var), "/")
  for (j in 1:3) {
    expect_equal(sc[[paste0(".", colnames(W)[j])]], unname(S[, j]),
      tolerance = 1e-12
    )
  }
})

test_that("predict(scaling='fit') on a polychoric-basis object: metric-consistent + warns", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)[1:400, 1:10]
  suppressWarnings(x <- ackwards(d, k_max = 2, cor = "polychoric"))

  # The model-implied-SD caveat still fires for materialized scores.
  rlang::reset_warning_verbosity("ackwards_nonpearson_scores")
  expect_warning(full <- predict(x, d), "model-implied SDs")

  # Fit-moment scaling keeps subsets on the same metric as the full set.
  sub <- suppressWarnings(predict(x, d[1:25, ]))
  expect_equal(sub, full[1:25, ], ignore_attr = TRUE)
})

test_that("predict() works on a truncated hierarchy (converged levels only)", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # 6 variables: lavaan::efa() fits k <= 3 only; k_max = 5 triggers truncation
  suppressWarnings(x <- ackwards(d, k_max = 5, engine = "esem"))
  expect_equal(x$k_max, 3L)
  sc <- predict(x, d[1:10, ])
  expect_equal(ncol(sc), 1L + 2L + 3L)
  expect_equal(nrow(sc), 10L)
  expect_false(anyNA(sc))
})

test_that("predict() propagates NA rows in place", {
  skip_if_not_installed("psych")
  d <- na.omit(ackwards::bfi25)
  suppressWarnings(x <- ackwards(d, k_max = 2))
  test <- d[1:10, ]
  test[3, 1] <- NA
  sc <- suppressWarnings(predict(x, test))
  expect_true(all(is.na(sc[3, ])))
  expect_false(anyNA(sc[-3, ]))
})
