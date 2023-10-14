#' Waller's Bass-Ackwards Analysis using PCA
#'
#' @param r A square correlation matrix.
#' @param nfactors The maximum number of factors (i.e., PCs) to estimate).
#' @param rotate A string containing the method for PCA rotation: either varimax
#'   (orthogonal) or promax (oblique). (default = `"varimax"`).
#'
#' @return A list object containing basic results
#' @export
#'
#' @examples
#' wba(cor(skiers), nfactors = 3)
wba <- function(r, nfactors, rotate = "varimax") {

  validate_nfactors(nfactors)
  rotate <- match.arg(rotate, choices = c("varimax", "promax"))

  # Get or create variable names
  varnames <- rownames(r, do.NULL = FALSE, prefix = "x")

  # Count variables
  nvar <- ncol(r)

  # Spectral decomposition of correlation matrix
  sdr <- eigen(r)

  # Extract raw eigenvalues
  eval <- sdr$values
  eval[eval < 0] <- 0 # TODO: Check this is legit

  # Extract raw eigenvectors (W)
  evec <- sdr$vectors

  # Flip the signs of eigenvectors that sum negative
  evec <- evec %*% diag(sign(colSums(evec)))

  # Calculate the unrotated PCA loadings
  pcl <- evec %*% diag(sqrt(eval))

  # Preallocate list objects
  correlations <- rep(list(vector("list", length = nfactors)), times = nfactors)
  rotations <- vector(mode = "list", length = nfactors)
  loadings <- vector(mode = "list", length = nfactors)

  # Loop through number of factors to get loadings and rotations
  for (i in 1:nfactors) {
    if (i == 1) {
      # If one factor, just take unrotated solution from pcl[, i]
      loadings[[i]] <- pcl[, i, drop = FALSE]
      correlations[[i]][[i]] <- matrix(1, 1, 1)
    } else {
      # If more than one factor, get rotated solution on pcl[, 1:i]
      if (rotate == "varimax") {
        vout <- stats::varimax(pcl[, 1:i], normalize = TRUE, eps = 1e-15)
        correlations[[i]][[i]] <- diag(i)
      } else if (rotate == "promax") {
        vout <- psych::Promax(pcl[, 1:i])
        correlations[[i]][[i]] <- vout$Phi
      }
      loadings[[i]] <- vout$loadings[1:nvar, ]
      rotations[[i]] <- vout$rotmat
    }
    rownames(loadings[[i]]) <- varnames
    colnames(loadings[[i]]) <- make_seq_names(i)
    rownames(correlations[[i]][[i]]) <- make_seq_names(i)
    colnames(correlations[[i]][[i]]) <- make_seq_names(i)
  }

  # Loop through number of factors to get cross-level correlations
  for (i in 1:nfactors) {
    for (j in 1:nfactors) {
      if (i != j) {
      correlations[[i]][[j]] <-
        get_clr_waller(i, j, rotations, loadings, rotate)
      }
    }
  }

  # Create and return list object
  new_bar(
    correlations = correlations,
    loadings = loadings,
    details = list(
      method = "Waller",
      engine = "eigen",
      rotation = rotate,
      eigenvalues = eval,
      eigenvectors = evec,
      rotations = rotations
    )
  )
}

# Calculate cross-level correlations using Wallers' methods
get_clr_waller <- function(i, j, rotations, loadings, rotate) {

  lo <- min(i, j)
  hi <- max(i, j)

  if (i != j) {
    if (lo == 1) {
      # Pull correlations from rotation matrix
      out <- rotations[[hi]][1, , drop = FALSE]
    } else {
      S <- cbind(diag(lo), matrix(0, lo, 1))
      if (rotate == "varimax") {
        # Calculate correlations for orthogonal factors
        out <- t(rotations[[lo]]) %*% S %*% rotations[[hi]]
      } else if (rotate == "promax") {
        # Calculate correlations for oblique factors
        out <- solve(rotations[[lo]]) %*% S %*% solve(t(rotations[[hi]]))
      }
    }
  }

  rownames(out) <- colnames(loadings[[lo]])
  colnames(out) <- colnames(loadings[[hi]])

  if (j > i) {
    out <- t(out)
  }

  out
}
