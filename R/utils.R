# Generic utilities -------------------------------------------------------

# Power operator for matrices
"%^%" <- function(S, power) with(eigen(S), vectors %*% (values ^ power * t(vectors)))

# Calculate ten berge factor scores
#   From seminr package
calc_ten_berge <- function(X, Lambda, Phi, i.means, i.sds) {
  if (any(is.na(X))) {
    # if any missing, impute using person average
    p.means <- rowMeans(X, na.rm = TRUE)
    missings <- which(is.na(X), arr.ind = TRUE)
    X[is.na(X)] <- p.means[missings[, 1]]
    X <- scale(X)
  } else {
    X <- t((t(X) - i.means) / i.sds)
  }
  R <- stats::cor(X, use = "pairwise")
  R.sqrt.i <- R %^% -0.5
  Phi.sqrt <- Phi %^% 0.5
  L <- Lambda %*% Phi.sqrt
  C <- R.sqrt.i %*% L %*% ((t(L) %*% chol2inv(chol(R)) %*% L) %^% -0.5)
  W <- R.sqrt.i %*% C %*% Phi.sqrt
  colnames(W) <- colnames(Lambda)
  rownames(W) <- rownames(Lambda)
  scores <- X %*% W
  colnames(scores) <- colnames(Lambda)
  scores
}

make_seq_names <- function(level) {
  stopifnot(level <= 26)
  paste0(LETTERS[[level]], 1:level)
}

validate_nfactors <- function(nfactors) {
  stopifnot(length(nfactors) == 1)
  stopifnot(is.integerish(nfactors))
  stopifnot(nfactors > 0)
}

is.integerish <- function(x) {
  is.numeric(x) & (x %% 1 == 0)
}

# lavaan utilities --------------------------------------------------------

# Adapted from seminr package
ten_berge_lavaan <- function(fit) {
  x <- lavaan::lavInspect(fit, what = "data")
  m <- colMeans(x)
  s <- sqrt(diag(lavaan::lavInspect(fit, what = "sampstat")$cov))
  lambda <- lavaan::lavInspect(fit, what = "std.lv")$lambda
  phi <- matrix(lavaan::lavInspect(fit, what = "cor.lv"), ncol(lambda))
  calc_ten_berge(x, lambda, phi, m, s)
}

make_lavaan_efa_syntax <- function(nfactors, varnames) {
  paste0(
    paste('efa("efa")*f', 1:nfactors, sep = "", collapse = " + "),
    " =~ ",
    paste(varnames, collapse = " + ")
  )
}
