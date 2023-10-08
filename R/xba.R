#' Extended Bass-Ackwards Analysis using psych
#'
#' @param data A data frame or matrix containing the raw data.
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
#' xba_psych(r = forbes, 10)
xba_psych <- function(data = NULL, r = NULL, nfactors, fm = "minres",
                      rotate = "varimax", scores = "tenBerge", ...) {
  # Assertions
  stopifnot(xor(is.null(data), is.null(r)))
  stopifnot(is.null(data) || is.matrix(data) || is.data.frame(data))
  stopifnot(is.null(r) || is.matrix(r) || is.data.frame(r))
  stopifnot(is.null(r) || ncol(r) == nrow(r))
  stopifnot(is.numeric(nfactors))
  stopifnot(length(nfactors) == 1)
  stopifnot(is.finite(nfactors))
  stopifnot(nfactors >= 1)
  stopifnot(floor(nfactors) == ceiling(nfactors))
  stopifnot(is.character(rotate))
  stopifnot(length(rotate) == 1)
  fm <- match.arg(
    fm,
    choices = c(
      "minres", "uls", "ols", "wls", "gls", "pa", "ml", "minchi", "minrank",
      "old.min", "alpha", "pca"
    )
  )

  # Coerce r to a matrix if needed
  if (!is.null(r) && is.matrix(r) == FALSE) {
    r <- as.matrix(r)
  }

  # Preallocate vectors
  fits <- vector(mode = "list", length = nfactors)
  loadings <- vector(mode = "list", length = nfactors)
  fscores <- vector(mode = "list", length = nfactors)
  weights <- vector(mode = "list", length = nfactors)
  corrs <- vector(mode = "list", length = nfactors)
  congs <- vector(mode = "list", length = nfactors)

  # Loop through number of factors from 1 to nfactors
  for (i in seq_along(fits)) {

    # Fit a factor analytic model to r with i factors and save results in fits
    if (fm == "pca") {
      if (!is.null(r)) {
        fits[[i]] <- psych::pca(
          r = r,
          nfactors = i,
          rotate = rotate,
          ...
        )
      } else {
        fits[[i]] <- psych::pca(
          r = data,
          nfactors = i,
          rotate = rotate,
          ...
        )
      }
    } else {
      if (!is.null(r)) {
        fits[[i]] <- psych::fa(
          r = r,
          nfactors = i,
          fm = fm,
          rotate = rotate,
          scores = scores,
          ...
        )
      } else {
        fits[[i]] <- psych::fa(
          r = data,
          nfactors = i,
          fm = fm,
          rotate = rotate,
          scores = scores,
          ...
        )
      }
    }

    # Extract the loadings and factor scores (or beta weights)
    loadings[[i]] <- fits[[i]]$loadings
    weights[[i]] <- fits[[i]]$weights
    if (!is.null(data)) {
      fscores[[i]] <- fits[[i]]$scores
    }

    # Rename the factors in loadings and weights as A1, B1:B2, C1:C3, etc.
    fnames <- make_seq_names(i)
    colnames(loadings[[i]]) <- fnames
    colnames(weights[[i]]) <- fnames
    if (!is.null(data)) {
      colnames(fscores[[i]]) <- fnames
    }

    if (i > 1) {
      # Preallocate inner lists for cross-level correlations and congruences
      corrs[[i]] <- vector(mode = "list", length = i - 1)
      congs[[i]] <- vector(mode = "list", length = i - 1)

      # Loop through all levels lowers than i
      for (j in 1:(i - 1)) {
        # Calculate factor correlations between levels i and j
        if (!is.null(r)) {
          corrs[[i]][[j]] <- t(weights[[j]]) %*% r %*% weights[[i]]
        } else {
          corrs[[i]][[j]] <- stats::cor(fscores[[j]], fscores[[i]])
        }

        # Calculate factor congruences between levels i and j
        congs[[i]][[j]] <- psych::fa.congruence(
          loadings[[j]],
          loadings[[i]]
        )
      }
    }
  }

  # Save results to output list
  out <- list(
    models = fits,
    loadings = loadings,
    weights = weights,
    fscores = fscores,
    correlations = corrs,
    congruences = congs
  )

  # Return output list
  out
}

#' Extended Bass-Ackwards Analysis using lavaan
#'
#' @param data A data frame or matrix containing the raw data.
#' @param nfactors The maximum number of factors to estimate.
#' @param rotate The factor rotation method from lavaan (default = `"varimax"`).
#' @param ... Additional arguments to pass on to `lavaan::cfa()`.
#'
#' @return A list object with the basic results
#' @export
#'
#' @examples
#' data("PoliticalDemocracy", package = "lavaan")
#' xba_lavaan(PoliticalDemocracy, nfactors = 3, rotate = "varimax")
xba_lavaan <- function(data, nfactors, rotate = "varimax", ...) {

  # Assertions
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.numeric(nfactors))
  stopifnot(length(nfactors) == 1)
  stopifnot(is.finite(nfactors))
  stopifnot(nfactors >= 1)
  stopifnot(floor(nfactors) == ceiling(nfactors))
  stopifnot(is.character(rotate))
  stopifnot(length(rotate) == 1)

  # Preallocate vectors
  fits <- vector(mode = "list", length = nfactors)
  loadings <- vector(mode = "list", length = nfactors)
  fscores <- vector(mode = "list", length = nfactors)
  corrs <- vector(mode = "list", length = nfactors)
  congs <- vector(mode = "list", length = nfactors)

  # Loop through number of factors from 1 to nfactors
  for (i in seq_along(fits)) {
    # Fit a factor analytic model to r with i factors and save results in fits
    fits[[i]] <- lavaan::cfa(
      model = make_lavaan_efa_syntax(i, colnames(data)),
      data = data,
      rotation = rotate,
      ...
    )

    # Extract loadings and factor scores
    loadings[[i]] <- lavaan::lavInspect(fits[[i]], what = "std.lv")$lambda
    fscores[[i]] <- ten_berge_lavaan(fits[[i]])

    # Rename the factors in loadings and weights as A1, B1:B2, C1:C3, etc.
    fnames <- make_seq_names(i)
    colnames(loadings[[i]]) <- fnames
    colnames(fscores[[i]]) <- fnames

    if (i > 1) {
      # Preallocate inner lists for cross-level correlations and congruences
      corrs[[i]] <- vector(mode = "list", length = i - 1)
      congs[[i]] <- vector(mode = "list", length = i - 1)

      # Loop through all levels lowers than i
      for (j in 1:(i - 1)) {
        # Calculate factor correlations between levels i and j
        corrs[[i]][[j]] <- stats::cor(fscores[[i]], fscores[[j]])
        # Calculate factor congruences between levels i and j
        congs[[i]][[j]] <- psych::fa.congruence(
          loadings[[j]],
          loadings[[i]]
        )
      }
    }
  }

  # Save results to output list
  out <- list(
    models = fits,
    loadings = loadings,
    weights = weights,
    fscores = fscores,
    correlations = corrs,
    congruences = congs
  )

  # Return output list
  out

}

