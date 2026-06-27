# PCA engine -- internal, not exported
#' @importFrom stats setNames
#
# Uses psych::pca() (varimax rotation by default), which is the same function
# psych::bassAckward() uses internally. This ensures our PCA path matches
# psych's reference implementation within floating-point tolerance.
#
# For later milestones, cfQ (oblique) will use GPArotation::cfQ directly.
#
# Returns a list indexed by k (as character), each element matching the
# DESIGN.md §4 level contract.

pca_levels <- function(R, k_max, rotation, cor_type = "pearson") {
  rlang::check_installed("psych", reason = "for the PCA engine")

  p <- nrow(R)
  result <- vector("list", k_max)
  names(result) <- as.character(seq_len(k_max))

  # Map our rotation label to psych's rotate argument.
  # "cfT" (orthogonal CF ≈ varimax) maps to "varimax" for M1 via psych::pca().
  psych_rotate <- switch(rotation,
    cfT = "varimax",
    cfQ = cli::cli_abort(
      "rotation = {.val cfQ} is not yet implemented for the PCA engine. \\
       Oblique CF rotation is planned for a future milestone."
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
  }

  result
}
