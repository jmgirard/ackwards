#' Waller's Bass-Ackwards Analysis using PCA
#'
#' @param r A square correlation matrix.
#' @param nfactors The maximum number of factors (i.e., PCs) to estimate).
#' @param eps The tolerance for stopping varimax; the relative change in the sum
#'   of singular values (default = `1e-15`).
#'
#' @return A list object containing basic results
#' @export
#'
#' @examples
#' wba(cor(skiers), nfactors = 3)
wba <- function(r, nfactors, eps = 1e-15) {

  # Get or create variable names
  varnames <- rownames(r, do.NULL = FALSE, prefix = "x")

  # Count variables
  nvar <- ncol(r)

  # Spectral decomposition of correlation matrix
  sdr <- eigen(r)

  # Extract raw eigenvalues
  eval <- sdr$values

  # Extract raw eigenvectors
  evec <- sdr$vectors

  # Flip the signs of eigenvectors that sum negative
  evec <- evec %*% diag(sign(colSums(evec)))

  # TODO: Give this line a description (???)
  P <- evec %*% diag(sqrt(eval))

  # Preallocate list objects
  correlations <- vector(mode = "list", length = nfactors)
  rotations <- vector(mode = "list", length = nfactors)
  loadings <- vector(mode = "list", length = nfactors)

  # Loop through number of factors to get loadings and rotations
  for (i in 1:nfactors) {
    if (i == 1) {
      # If one factor, just take unrotated solution from P[, i]
      loadings[[i]] <- P[, i, drop = FALSE]
      rotations[[i]] <- diag(1)
    } else {
      # If more than one factor, get varimax solution on P[, 1:i]
      vout <- stats::varimax(P[, 1:i], normalize = TRUE, eps = eps)
      loadings[[i]] <- vout$loadings[1:nvar, ]
      rotations[[i]] <- vout$rotmat
    }
    rownames(loadings[[i]]) <- varnames
    colnames(loadings[[i]]) <- make_seq_names(i)
  }

  # Loop through number of factors to get cross-level correlations
  for (i in 2:nfactors) {
    # TODO: Give this line a description (???)
    S <- cbind(diag(i - 1), matrix(0, i - 1, 1))

    # Base correlation on rotations and S
    correlations[[i]] <- t(rotations[[i - 1]]) %*% S %*% rotations[[i]]

    # Add sequential names
    rownames(correlations[[i]]) <- colnames(loadings[[i - 1]])
    colnames(correlations[[i]]) <- colnames(loadings[[i]])

    # Transpose for pyramid output
    correlations[[i]] <- t(correlations[[i]])
  }

  # Create and return list object
  list(
    correlations = correlations,
    loadings = loadings,
    rotations = rotations,
    details = list(eigenvalues = eval, eigenvectors = evec)
  )
}
