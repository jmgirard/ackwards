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
