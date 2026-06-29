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
