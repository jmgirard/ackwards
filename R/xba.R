#' Extended Bass-Ackwards Analysis using psych
#'
#' @param r A correlation matrix (either as a matrix or data frame).
#' @param nfactors The maximum number of factors to estimate.
#' @param fm The factoring method from psych (default = `"minres"`).
#' @param rotate The factor rotation method from psych (default = `"varimax"`).
#' @param scores The factor scoring method from psych (default = `"tenBerge"`).
#' @param ... Additional arguments to pass on to `psych::fa()`.
#'
#' @return A list object with the basic results
#' @export
#'
#' @examples
#' xba_psych(forbes, 10)
xba_psych <- function(r, nfactors, fm = "minres", rotate = "varimax",
                      scores = "tenBerge", ...) {
  # Assertions
  stopifnot(is.matrix(r) || is.data.frame(r))
  stopifnot(ncol(r) == nrow(r))
  stopifnot(is.numeric(nfactors))
  stopifnot(length(nfactors) == 1)
  stopifnot(is.finite(nfactors))
  stopifnot(nfactors >= 1)
  stopifnot(floor(nfactors) == ceiling(nfactors))
  fm <- match.arg(
    fm,
    choices = c(
      "minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank",
      "old.min", "alpha", "pca"
    )
  )
  rotate <- match.arg(
    rotate,
    choices = c(
      "none", "varimax", "quartimax", "bentlerT", "equamax", "varimin",
      "geominT", "bifactor", "Promax", "promax", "oblimin", "simplimax",
      "bentlerQ", "geominQ", "biquartimin", "cluster"
    )
  )

  # Coerce r to a matrix if needed
  if (is.matrix(r) == FALSE) {
    r <- as.matrix(r)
  }

  # Preallocate vectors
  fits <- vector(mode = "list", length = nfactors)
  corrs <- vector(mode = "list", length = nfactors)
  congs <- vector(mode = "list", length = nfactors)

  # Loop through number of factors from 1 to nfactors
  for (i in seq_along(fits)) {
    # Fit a factor analytic model to r with i factors and save results in fits
    if (fm == "pca") {
      fits[[i]] <- psych::pca(
        r = r,
        nfactors = i,
        rotate = rotate,
        ...
      )
    } else {
      fits[[i]] <- psych::fa(
        r = r,
        nfactors = i,
        fm = fm,
        rotate = rotate,
        scores = scores,
        ...
      )
    }

    # Rename the factors in loadings and weights as A1, B1:B2, C1:C3, etc.
    fnames <- paste0(LETTERS[[i]], 1:i)
    colnames(fits[[i]]$loadings) <- fnames
    colnames(fits[[i]]$weights) <- fnames

    if (i > 1) {
      # Preallocate inner lists for cross-level correlations and congruences
      corrs[[i]] <- vector(mode = "list", length = i - 1)
      congs[[i]] <- vector(mode = "list", length = i - 1)

      # Loop through all levels lowers than i
      for (j in 1:(i - 1)) {
        # Calculate factor correlations between levels i and j
        corrs[[i]][[j]] <- t(fits[[j]]$weights) %*% r %*% fits[[i]]$weights
        # Calculate factor congruences between levels i and j
        congs[[i]][[j]] <- psych::factor.congruence(
          fits[[j]]$loadings,
          fits[[i]]$loadings
        )
      }
    }
  }

  # Save results to output list
  out <- list(
    models = fits,
    correlations = corrs,
    congruences = congs
  )

  # Return output list
  out
}
