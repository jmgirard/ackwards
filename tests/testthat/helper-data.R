# Shared simulated datasets loaded by testthat before every test file.

.make_esem_data <- function(seed = 42, n = 200) {
  set.seed(seed)
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  data.frame(
    x1 = f1 + 0.4 * rnorm(n), x2 = f1 + 0.4 * rnorm(n),
    x3 = f1 + 0.4 * rnorm(n), x4 = f2 + 0.4 * rnorm(n),
    x5 = f2 + 0.4 * rnorm(n), x6 = f2 + 0.4 * rnorm(n)
  )
}

.make_ordinal_data <- function(seed = 42, n = 300) {
  set.seed(seed)
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  d <- data.frame(
    x1 = f1 + 0.4 * rnorm(n), x2 = f1 + 0.4 * rnorm(n),
    x3 = f1 + 0.4 * rnorm(n), x4 = f2 + 0.4 * rnorm(n),
    x5 = f2 + 0.4 * rnorm(n), x6 = f2 + 0.4 * rnorm(n)
  )
  breaks <- c(-Inf, -1, -0.5, 0.5, 1, Inf)
  as.data.frame(lapply(d, function(x) as.integer(cut(x, breaks))))
}
