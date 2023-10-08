test_that("wba results match Waller (2006)", {
  res <- wba(cor(skiers), nfactors = 3)

  # Correlations level 1 to 2
  r12 <- round(res$correlations[[2]], 3)
  expect_equal(r12[[1]],  0.904)
  expect_equal(r12[[2]], -0.427)

  # Correlations level 2 to 3
  r23 <- round(res$correlations[[3]], 3)
  expect_equal(r23[1, 1],  1.000)
  expect_equal(r23[2, 1],  0.000)
  expect_equal(r23[3, 1], -0.021)
  expect_equal(r23[1, 2],  0.000)
  expect_equal(r23[2, 2],  1.000)
  expect_equal(r23[3, 2],  0.003)

  # Loadings level 1 (not in paper)
  l1 <- round(res$loadings[[1]], 3)
  expect_equal(l1[1, 1], -0.500)
  expect_equal(l1[2, 1],  0.357)
  expect_equal(l1[3, 1],  0.891)
  expect_equal(l1[4, 1],  0.919)

  # Loadings level 2
  l2 <- round(res$loadings[[2]], 3)
  expect_equal(l2[1, 1], -0.087)
  expect_equal(l2[2, 1], -0.072)
  expect_equal(l2[3, 1],  0.997)
  expect_equal(l2[4, 1],  0.998)
  expect_equal(l2[1, 2],  0.988)
  expect_equal(l2[2, 2], -0.989)
  expect_equal(l2[3, 2],  0.026)
  expect_equal(l2[4, 2], -0.040)

  # Loadings level 3
  l3 <- round(res$loadings[[3]], 3)
  expect_equal(l3[1, 1], -0.084)
  expect_equal(l3[2, 1], -0.070)
  expect_equal(l3[3, 1],  0.998)
  expect_equal(l3[4, 1],  0.997)
  expect_equal(l3[1, 2],  0.987)
  expect_equal(l3[2, 2], -0.989)
  expect_equal(l3[3, 2],  0.026)
  expect_equal(l3[4, 2], -0.040)
  expect_equal(l3[1, 3],  0.134)
  expect_equal(l3[2, 3],  0.129)
  expect_equal(l3[3, 3],  0.033)
  expect_equal(l3[4, 3], -0.054)
})
