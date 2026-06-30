# EFA engine -- internal, not exported
#' @importFrom stats setNames
#
# Wraps psych::fa() for each level 1..k_max, returning the standard s.4 level
# contract. Scoring uses tenBerge weights (linear S = ZW), which keeps
# compute_edges() on the algebra path. Convergence failures truncate the
# hierarchy at the last successful level; Heywood cases warn but continue.

efa_levels <- function(R, k_max, fm, n_obs, cor = "pearson",
                       keep_fits = FALSE) {
  p <- nrow(R)
  result <- list()
  fits_list <- if (keep_fits) list() else NULL

  for (k in seq_len(k_max)) {
    rotate_k <- if (k == 1L) "none" else "varimax"

    # Run psych::fa(), intercepting convergence warnings so we can act on them
    warn_msgs <- character(0L)
    fit <- tryCatch(
      withCallingHandlers(
        psych::fa(R, nfactors = k, rotate = rotate_k, fm = fm, n.obs = n_obs),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        cli::cli_warn(
          c(
            "!" = "EFA failed at k = {k}: {conditionMessage(e)}",
            "i" = "Truncating hierarchy at level {k - 1L}."
          )
        )
        NULL
      }
    )

    if (is.null(fit)) break

    # Convergence failure detected via psych warning messages
    converge_fail <- any(grepl(
      "did not converge|failed to converge|no convergence|not converge",
      warn_msgs,
      ignore.case = TRUE
    ))
    if (converge_fail) { # nocov start
      cli::cli_warn(
        c(
          "!" = "EFA did not converge at k = {k}.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        )
      )
      break
    } # nocov end

    # Heywood case: warn but do NOT truncate (convergence is data, not an error)
    heywood <- any(fit$uniquenesses < 0, na.rm = TRUE) ||
      any(fit$communalities > 1, na.rm = TRUE)
    if (heywood) {
      cli::cli_warn(
        c(
          "!" = "Heywood case at k = {k}: communality > 1 or uniqueness < 0.",
          "i" = "Results may be unreliable. Consider reducing {.arg k} or changing {.arg fm}."
        )
      )
    }

    L_rot <- unclass(fit$loadings)
    labels_k <- make_labels(k)

    # Positive manifold anchor for k = 1 (matches PCA engine behaviour).
    # nocov: only fires when a single-factor solution loads net-negative, which
    # does not occur for positive-manifold data; the PCA analogue is excluded
    # the same way (engine_pca.R).
    flip <- (k == 1L) && (sum(L_rot) < 0)
    if (flip) L_rot <- -L_rot # nocov

    colnames(L_rot) <- labels_k
    rownames(L_rot) <- rownames(R)

    # tenBerge weights: W = R^{-1} L (L' R^{-1} L)^{-1/2}
    # These make scores exactly uncorrelated with unit variance for orthogonal
    # factors, keeping the algebra path in compute_edges() valid.
    weight_method <- "tenBerge"
    W <- tryCatch(
      .tenBerge_weights(R, L_rot),
      error = function(e) { # nocov start
        weight_method <<- "regression" # honest label on fallback (Invariant 6)
        cli::cli_warn(
          c(
            "!" = "tenBerge weights failed at k = {k}: {conditionMessage(e)}",
            "i" = "Falling back to regression (Thurstone) weights."
          )
        )
        w_fall <- unclass(fit$weights)
        if (flip) w_fall <- -w_fall
        colnames(w_fall) <- labels_k
        rownames(w_fall) <- rownames(R)
        w_fall
      } # nocov end
    )

    # Score variances: diag(W' R W); exact 1 for tenBerge (orthogonal), but
    # always compute rather than assume -- Invariant 1.
    score_var <- diag(crossprod(W, R %*% W))

    # Variance explained (sum of squared loadings / p)
    var_per_factor <- unname(colSums(L_rot^2) / p)
    variance <- c(
      setNames(var_per_factor, labels_k),
      cumulative = sum(var_per_factor)
    )

    # Fit indices -- available only when n.obs was supplied; NA otherwise.
    # fit$chi/dof/PVAL/TLI/BIC are plain scalars; fit$RMSEA is a named vector.
    fit_info <- setNames(
      c(
        if (!is.null(fit$chi)) unname(fit$chi)[[1L]] else NA_real_,
        if (!is.null(fit$dof)) unname(fit$dof)[[1L]] else NA_real_,
        if (!is.null(fit$PVAL)) unname(fit$PVAL)[[1L]] else NA_real_,
        if (!is.null(fit$RMSEA)) unname(fit$RMSEA["RMSEA"])[[1L]] else NA_real_,
        if (!is.null(fit$TLI)) unname(fit$TLI)[[1L]] else NA_real_,
        if (!is.null(fit$BIC)) unname(fit$BIC)[[1L]] else NA_real_
      ),
      c("chi", "dof", "p_value", "RMSEA", "TLI", "BIC")
    )

    result[[as.character(k)]] <- list(
      k = k,
      loadings = L_rot,
      loadings_se = NULL, # EFA (psych::fa) does not produce rotation-aware SEs
      variance = variance,
      fit = fit_info,
      converged = TRUE,
      factor_cor = diag(k), # orthogonal: I_k
      labels = labels_k,
      scoring = list(
        linear    = TRUE,
        method    = weight_method, # "tenBerge" normally; "regression" on fallback
        basis     = cor, # reflects actual R basis, not assumed "pearson"
        weights   = W,
        score_var = score_var
      )
    )
    if (keep_fits) fits_list[[as.character(k)]] <- fit
  }

  list(levels = result, fits = fits_list)
}

# Compute tenBerge factor-score weights from a correlation matrix R and a
# (rotated, sign-aligned) loading matrix L.
#
# Formula (orthogonal factors):  W = R^{-1} L (L' R^{-1} L)^{-1/2}
#
# The resulting W satisfies W'RW = I, so tenBerge scores have unit variance --
# the D standardization in compute_edges() divides by 1 but is still applied
# for numerical safety and to satisfy Invariant 1.
.tenBerge_weights <- function(R, L) {
  Ri <- solve(R) # p x p
  A <- Ri %*% L # p x k: R^{-1} L
  B <- crossprod(L, A) # k x k: L' R^{-1} L  (symmetric PD for full-rank L)

  # Matrix inverse square root of B via spectral decomposition
  eig <- eigen(B, symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps) # guard against tiny negatives
  Binvsqrt <- eig$vectors %*%
    diag(1 / sqrt(vals), nrow = length(vals)) %*%
    t(eig$vectors)

  W <- A %*% Binvsqrt
  rownames(W) <- rownames(L)
  colnames(W) <- colnames(L)
  W
}
