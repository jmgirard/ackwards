test_that(".format_r() strips leading zero and pads trailing zeros", {
  fmt <- ackwards:::.format_r

  expect_equal(fmt(0.23, 2L), ".23")
  expect_equal(fmt(-0.23, 2L), "-.23")
  expect_equal(fmt(0.3, 2L), ".30")
  expect_equal(fmt(-0.3, 2L), "-.30")
  expect_equal(fmt(0, 2L), ".00")
  expect_equal(fmt(1, 2L), "1.00")
  expect_equal(fmt(-1, 2L), "-1.00")
})

test_that(".format_r() does not produce '-.00' when magnitude rounds to zero", {
  fmt <- ackwards:::.format_r

  # -0.003 rounds to .00 at 2 digits; sign should be suppressed
  expect_equal(fmt(-0.003, 2L), ".00")
  # Exact zero is also clean
  expect_equal(fmt(-0, 2L), ".00")
})

test_that(".format_r() respects digits argument", {
  fmt <- ackwards:::.format_r

  expect_equal(fmt(0.234, 3L), ".234")
  expect_equal(fmt(0.2, 1L), ".2")
})

test_that(".format_r() is vectorised", {
  fmt <- ackwards:::.format_r

  result <- fmt(c(0.23, -0.45, 1, 0), 2L)
  expect_equal(result, c(".23", "-.45", "1.00", ".00"))
})

test_that("inst/CITATION has two entries with correct years and no unknown year", {
  skip_if_not_installed("ackwards")
  cites <- citation("ackwards")
  expect_length(cites, 2)
  # Goldberg (2006) method paper
  expect_equal(cites[[1]]$year, "2006")
  # Package entry: year is set, version note is present, no "????" fallback
  pkg_year <- cites[[2]]$year
  expect_false(is.null(pkg_year))
  expect_false(is.na(pkg_year))
  expect_match(pkg_year, "^[0-9]{4}$")
  expect_match(cites[[2]]$note, "R package version")
})
