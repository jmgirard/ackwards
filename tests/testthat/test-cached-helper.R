# Unit tests for the suite-wide cached() fit memo in helper-data.R (M48).
# The memo's key semantics are load-bearing for the whole suite's validity:
# a false cache hit would silently make unrelated tests assert against the
# wrong object. These tests use cheap local functions, never package fits.

test_that("cached() evaluates once per key and returns the memoized value", {
  cnt <- new.env()
  cnt$n <- 0L
  f <- function(x) {
    cnt$n <- cnt$n + 1L
    sum(x)
  }
  d <- c(1, 2)
  expect_identical(cached(f(d)), 3)
  expect_identical(cached(f(d)), 3)
  expect_identical(cnt$n, 1L)
})

test_that("cached() distinguishes same call text on different data", {
  cnt <- new.env()
  cnt$n <- 0L
  f <- function(x) {
    cnt$n <- cnt$n + 1L
    sum(x)
  }
  d <- c(1, 2)
  expect_identical(cached(f(d)), 3)
  d <- c(5, 5) # same deparse text `f(d)`, different value
  expect_identical(cached(f(d)), 10)
  expect_identical(cnt$n, 2L)
})

test_that("cached() distinguishes same-named functions with different bodies", {
  f <- function(x) sum(x)
  d <- c(1, 2)
  expect_identical(cached(f(d)), 3)
  f <- function(x) prod(x) # same call text and same `d`, different function
  expect_identical(cached(f(d)), 2)
})

test_that("mutating a cached() result does not poison the memo", {
  f <- function(x) list(v = x)
  d <- 1:3
  a <- cached(f(d))
  a$v <- 99L
  expect_identical(cached(f(d))$v, 1:3)
})

test_that("cached() muffles fit-time conditions identically on miss and hit", {
  f <- function(x) {
    warning("fit-time warning")
    message("fit-time message")
    x
  }
  d <- 7L
  expect_no_condition(first <- cached(f(d))) # miss: evaluated, muffled
  expect_no_condition(second <- cached(f(d))) # hit: nothing evaluated
  expect_identical(first, 7L)
  expect_identical(second, 7L)
})
