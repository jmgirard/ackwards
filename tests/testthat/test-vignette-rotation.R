# Guard tests for the "Orthogonal or oblique rotation" section of
# ackwards-engines.Rmd: EFA with varimax and with oblimin on the polychoric
# basis of na.omit(bfi25), k_max = 5. Each test pins findings that the section's
# prose states. The helpers below repeat the section's chunk code.

.bfi_rot <- na.omit(bfi25)

.rot_fits <- function() {
  skip_if_not_installed("GPArotation")
  list(
    var = cached(ackwards(.bfi_rot, k_max = 5, engine = "efa", cor = "polychoric")),
    obl = cached(ackwards(.bfi_rot,
      k_max = 5, engine = "efa", cor = "polychoric",
      rotation = "oblimin", seed = 1
    ))
  )
}

# Item-by-factor loading matrix at level k (the vignette's loading_matrix()).
.rot_loading_matrix <- function(x, k) {
  d <- tidy(x, what = "loadings")
  unclass(stats::xtabs(loading ~ item + factor, data = d[d$level == k, ]))
}

# Match each oblimin factor to the varimax factor with the largest absolute
# Tucker congruence of loadings, level by level (the vignette's matches table).
.rot_matches <- function(x_var, x_obl) {
  do.call(rbind, lapply(1:5, function(k) {
    cong <- psych::factor.congruence(
      .rot_loading_matrix(x_obl, k), .rot_loading_matrix(x_var, k),
      digits = 3
    )
    best <- apply(abs(cong), 1, which.max)
    data.frame(
      level = k,
      oblimin = rownames(cong),
      varimax = colnames(cong)[best],
      congruence = abs(cong[cbind(seq_along(best), best)])
    )
  }))
}

# Secondary edges flagged above_cut (the vignette's above_cut_secondary()).
.rot_secondary <- function(x) {
  e <- tidy(x)
  e[e$above_cut & !e$is_primary, c("from", "to", "r", "beta")]
}

test_that("rotation section (1): the largest oblimin factor correlation is m5f3-m5f4 at .34", {
  fits <- .rot_fits()
  fc <- tidy(fits$obl, what = "factor_cor")
  fc <- fc[order(-abs(fc$cor)), ]
  # The three largest, in order: the prose names the first and calls the
  # other two "the next two", citing them for the secondary edges.
  expect_identical(
    paste(fc$factor_a, fc$factor_b, sep = "-")[1:3],
    c("m5f3-m5f4", "m3f1-m3f3", "m4f1-m4f3")
  )
  expect_identical(fc$level[1], 5L)
  expect_identical(round(fc$cor[1:3], 2), c(0.34, 0.33, 0.30))
})

test_that("rotation section (2): one-to-one match, smallest congruence .94, only level 5 swaps IDs", {
  fits <- .rot_fits()
  m <- .rot_matches(fits$var, fits$obl)
  # One-to-one within every level.
  expect_false(anyDuplicated(m$varimax) > 0L)
  expect_identical(nrow(m), 15L)
  # Same IDs everywhere except m5f2 and m5f3, which trade places.
  swapped <- m[m$oblimin != m$varimax, c("oblimin", "varimax")]
  expect_identical(swapped$oblimin, c("m5f2", "m5f3"))
  expect_identical(swapped$varimax, c("m5f3", "m5f2"))
  # Smallest congruence per level (levels 2 to 5). The prose states the
  # level-5 value, .94.
  min_by_level <- as.vector(tapply(m$congruence, m$level, min))
  expect_identical(round(min_by_level[2:5], 2), c(0.99, 0.98, 0.98, 0.94))
  expect_true(m$oblimin[which.min(m$congruence)] %in% c("m5f2", "m5f3"))
})

test_that("rotation section (3): after the match, the primary-parent trees are equal", {
  fits <- .rot_fits()
  m <- .rot_matches(fits$var, fits$obl)
  to_var <- stats::setNames(m$varimax, m$oblimin)
  p_var <- tidy(fits$var, primary_only = TRUE)[, c("from", "to")]
  p_obl <- tidy(fits$obl, primary_only = TRUE)[, c("from", "to")]
  p_obl_as_var <- data.frame(
    from = unname(to_var[p_obl$from]),
    to = unname(to_var[p_obl$to])
  )
  key <- function(d) sort(paste(d$from, d$to, sep = "->"))
  expect_identical(key(p_obl_as_var), key(p_var))
  # The raw oblimin edges differ only where level 5 swaps IDs.
  raw_only <- setdiff(key(p_obl), key(p_var))
  expect_identical(raw_only, c("m4f1->m5f3", "m4f3->m5f2"))
})

test_that("rotation section (4): above-cut secondary edges, none for varimax, three for oblimin", {
  fits <- .rot_fits()
  expect_identical(nrow(.rot_secondary(fits$var)), 0L)

  s <- .rot_secondary(fits$obl)
  expect_identical(
    paste(s$from, s$to, sep = "->"),
    c("m3f1->m4f3", "m4f1->m5f2", "m4f3->m5f4")
  )
  expect_true(all(s$r >= 0.3 & s$r < 0.35))
  # Every beta is positive and smaller than its r.
  expect_true(all(s$beta > 0 & s$beta < s$r))
  beta <- stats::setNames(s$beta, paste(s$from, s$to, sep = "->"))
  # Pinned at the precision the prose states them.
  expect_identical(round(beta[["m3f1->m4f3"]], 3), 0.016)
  expect_identical(signif(beta[["m4f1->m5f2"]], 2), 0.0015)
  expect_identical(round(beta[["m4f3->m5f4"]], 2), 0.16)
})
