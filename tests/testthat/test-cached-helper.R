# Unit tests for the suite-wide cached() fit memo in helper-data.R (M48).
# The memo's key semantics are load-bearing for the whole suite's validity:
# a false cache hit would silently make unrelated tests assert against the
# wrong object. These tests use cheap local functions, never package fits.

# Each test probes cache-key semantics in isolation, so it starts from an empty
# .fit_cache. Two sibling tests that define a textually identical local `f(d)`
# on the same data hash to the *same* key once source references are dropped (as
# under R CMD check, and on the R-devel CI runner), so without this reset the
# earlier test's entry is a false hit in the later one (observed as `cnt$n == 1`
# where 2 was expected). Clearing only drops memoized fits that other files
# would recompute (a perf, not correctness, concern) and is a no-op under the
# parallel runner, where every file has its own process and cache.
clean_fit_cache <- function() {
  rm(list = ls(.fit_cache, all.names = TRUE), envir = .fit_cache)
}

test_that("cached() evaluates once per key and returns the memoized value", {
  clean_fit_cache()
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
  clean_fit_cache()
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
  clean_fit_cache()
  f <- function(x) sum(x)
  d <- c(1, 2)
  expect_identical(cached(f(d)), 3)
  f <- function(x) prod(x) # same call text and same `d`, different function
  expect_identical(cached(f(d)), 2)
})

test_that("mutating a cached() result does not poison the memo", {
  clean_fit_cache()
  f <- function(x) list(v = x)
  d <- 1:3
  a <- cached(f(d))
  a$v <- 99L
  expect_identical(cached(f(d))$v, 1:3)
})

test_that("cached() muffles fit-time conditions identically on miss and hit", {
  clean_fit_cache()
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
