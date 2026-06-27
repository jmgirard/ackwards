# PCA engine -- internal, not exported
#' @importFrom stats setNames
#
# Uses psych::pca() (varimax rotation), which is the same function
# psych::bassAckward() uses internally — ensures our PCA path matches
# psych's reference implementation within floating-point tolerance.
#
# Returns list(levels = <named list per §4 contract>, fits = <named list | NULL>)

pca_levels <- function(R, k_max, rotation, cor_type = "pearson", keep_fits = FALSE) {
  rlang::check_installed("psych", reason = "for the PCA engine")

  p <- nrow(R)
  result <- vector("list", k_max)
  names(result) <- as.character(seq_len(k_max))
  fits_list <- if (keep_fits) vector("list", k_max) else NULL
  if (keep_fits) names(fits_list) <- as.character(seq_len(k_max))

  # cfT (orthogonal CF ≈ varimax) is the only supported rotation.
  # cfQ is caught before reaching here by ackwards(); this switch is a safety net.
  psych_rotate <- switch(rotation,
    cfT = "varimax",
    cli::cli_abort(
      "rotation = {.val {rotation}} is not supported. \\
       Only {.val cfT} (orthogonal CF) is available."
    )
  )

  for (k in seq_len(k_max)) {
    if (k == 1L) {
      # k=1: no rotation (single component is already determined)
      fit <- psych::pca(R, nfactors = 1L, rotate = "none")
      L_rot <- unclass(fit$loadings)
      # Ensure positive manifold
      if (sum(L_rot) < 0) {
        L_rot <- -L_rot
        fit$weights <- -fit$weights
      }
    } else {
      fit <- psych::pca(R, nfactors = k, rotate = psych_rotate)
      L_rot <- unclass(fit$loadings)
    }

    W <- unclass(fit$weights)

    # psych::pca() already sorts by descending variance explained.
    # Label columns with our stable m{k}f{j} scheme.
    colnames(L_rot) <- make_labels(k)
    rownames(L_rot) <- rownames(R)
    colnames(W) <- make_labels(k)
    rownames(W) <- rownames(R)

    # Score variances: diag(W' R W) — NOT assumed to be 1
    score_var <- diag(crossprod(W, R %*% W))

    # Variance explained per factor and cumulative
    var_per_factor <- unname(colSums(L_rot^2) / p)
    variance <- c(
      setNames(var_per_factor, make_labels(k)),
      cumulative = sum(var_per_factor)
    )

    # Eigenvalues as the "fit" summary for PCA levels
    eig <- fit$values[seq_len(k)]
    fit_info <- setNames(eig, paste0("eigenvalue.", make_labels(k)))

    result[[as.character(k)]] <- list(
      k = k,
      loadings = L_rot,
      loadings_se = NULL, # PCA does not produce rotation-aware SEs
      variance = variance,
      fit = fit_info,
      converged = TRUE,
      factor_cor = diag(k), # orthogonal: identity
      labels = make_labels(k),
      scoring = list(
        linear    = TRUE,
        method    = "components",
        basis     = cor_type,
        weights   = W,
        score_var = score_var
      )
    )
    if (keep_fits) fits_list[[as.character(k)]] <- fit
  }

  list(levels = result, fits = fits_list)
}
