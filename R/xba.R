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

  new_bar(
    correlations = corrs,
    loadings = loadings,
    details = list(
      method = "Forbes",
      engine = "psych",
      models = fits,
      weights = weights,
      fscores = fscores,
      congruences = congs
    )
  )
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

  new_bar(
    correlations = corrs,
    loadings = loadings,
    details = list(
      method = "Forbes",
      engine = "lavaan",
      models = fits,
      fscores = fscores,
      congruences = congs
    )
  )

}




ChaseCorrPaths <- function (comp.corr, component = "levelnum") #This function is called on below, not used alone. It calculates the paths of correlations >.9
{
  chased_levels <- vector() #initializing empty objects that are populated below
  chased_to_level <- vector()
  chased_to <- list()
  sub_revcomp.corr <- list()
  for (i in (length(comp.corr):1))
  {
    if (component %in% colnames(comp.corr[[i]]))
      #if the component we're interested is in the matrix
      chased_levels[[i]] <-
        (max(comp.corr[[i]][, component]) >= .9) #tell me if the maximum component correlation for the relevant column is >=.9
  }
  revcomp.corr <-
    rev(comp.corr) #reverse order of comp.corr to work from the bottom up
  component_level <-
    (length(chased_levels[!is.na(chased_levels)]) + 1) #level of current component (calculated as number of upward comparison matrices +1)
  chased_levels <-
    rev(chased_levels) #reverses order, so looking at lowest levels of hierarchy first

  if (any(chased_levels, na.rm = TRUE))
  {
    chased_levels <-
      (which.min(chased_levels)) - 1 #counts number of consecutive true values before first false
  }
  else {
    chased_levels <- 0
  }

  chased_to_level <- (component_level - chased_levels)

  if (chased_levels == 0)
  {
    chased_to <- "null"
  } #if no trues
  else {
    #isolate block of matrices in revcomp.corr relevant to component
    #end range
    end_comp.corr <- ((component_level * (component_level - 1) / 2))
    #start range
    start_comp.corr <- (end_comp.corr - (component_level - 2))

    sub_comp.corr <-
      comp.corr[start_comp.corr:end_comp.corr] #subset of matrices

    chased_to <-
      rownames(as.data.frame(which.max(sub_comp.corr[[chased_to_level]][, component]))) #name of component chased to

  }
  if ((component == "b1") &
      (max(comp.corr[[1]][, "b1"]) >= .9))
    #need to calculate the b->a level separately
  {
    chased_to <- "a1"
  }
  if ((component == "b2") & (max(comp.corr[[1]][, "b2"]) >= .9))
  {
    chased_to <- "a1"
  }
  result <- list(component, chased_to)
  result <- paste(result, sep = " ", collapse = "--")
  return(result)
}
