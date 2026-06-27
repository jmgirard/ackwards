# ESEM engine -- internal, not exported
#' @importFrom stats cov2cor
#
# Wraps lavaan::efa() for each level 1..k_max, returning the standard §4 level
# contract. For ordinal/polychoric data uses WLSMV estimation (Kim & Eaton 2015;
# Forbush et al. 2024); for continuous data uses ML. Scoring uses self-computed
# tenBerge weights from lavaan's estimated loadings and the underlying correlation
# matrix, keeping compute_edges() on the algebra path (DESIGN.md §14 item 12).
#
# Key additions over the EFA engine:
#   - loadings_se: rotation-aware delta-method SEs from lavaan's standardized solution
#   - fit: CFI, TLI, RMSEA, SRMR from lavaan::fitMeasures()
#   - Polychoric basis: for cor_type = "polychoric", ordered variables trigger WLSMV
#
# Arguments:
#   data       — numeric matrix, items in columns (raw data, NOT a correlation matrix)
#   k_max      — maximum number of factors to extract
#   rotation   — "cfT" (→ lavaan "varimax") or "cfQ" (errors; not yet implemented)
#   estimator  — "ML", "WLSMV", "ULSMV", etc.; passed directly to lavaan::efa()
#   cor_type   — "pearson", "spearman", or "polychoric"; controls ordered= argument
#   n_obs      — number of observations (for record-keeping only; lavaan uses data)
#   R_external — pre-computed correlation matrix for tenBerge weights (used for
#                continuous paths); NULL for polychoric (lavaan's R is extracted)
#
# Returns list(levels = <named list>, r_lv = <p x p matrix>)

esem_levels <- function(data, k_max, rotation, estimator, cor_type, n_obs,
                        R_external = NULL) {
  rlang::check_installed("lavaan", reason = "for the ESEM engine")
  if (!exists("efa", envir = asNamespace("lavaan"), inherits = FALSE)) {
    cli::cli_abort(
      c(
        "!" = "{.fn lavaan::efa} is not available in the installed version of lavaan.",
        "i" = "Please update lavaan: {.run install.packages('lavaan')} (>= 0.6-13 required)."
      )
    )
  }

  p          <- ncol(data)
  item_names <- colnames(data)

  # cfT ≈ varimax (CF with κ = 1/p); cfQ not yet implemented for ESEM
  lav_rotation <- switch(rotation,
    cfT = "varimax",
    cfQ = cli::cli_abort(
      "rotation = {.val cfQ} is not yet implemented for the ESEM engine. \\
       Oblique CF rotation is planned for a future milestone."
    )
  )

  # For polychoric/WLSMV: treat all columns as ordered categorical
  ordered_cols <- if (cor_type == "polychoric") item_names else NULL

  result <- list()
  r_lv   <- R_external   # set externally for continuous; extracted from lavaan for polychoric

  for (k in seq_len(k_max)) {
    # No rotation needed for a single factor
    rotate_k  <- if (k == 1L) "none" else lav_rotation
    warn_msgs <- character(0L)

    # lavaan::efa() returns an efaList; extract the single lavaan fit object.
    fit_raw <- tryCatch(
      withCallingHandlers(
        lavaan::efa(
          data      = as.data.frame(data),
          nfactors  = k,
          rotation  = rotate_k,
          estimator = estimator,
          ordered   = ordered_cols
        ),
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        cli::cli_warn(
          c(
            "!" = "ESEM failed at k = {k}: {conditionMessage(e)}",
            "i" = "Truncating hierarchy at level {k - 1L}."
          )
        )
        NULL
      }
    )

    if (is.null(fit_raw)) break
    fit <- if (inherits(fit_raw, "efaList")) fit_raw[[1L]] else fit_raw
    if (is.null(fit)) break

    # Convergence check
    converged <- isTRUE(tryCatch(
      lavaan::lavInspect(fit, "converged"),
      error = function(e) FALSE
    ))
    if (!converged) {
      cli::cli_warn(
        c(
          "!" = "ESEM did not converge at k = {k}.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        )
      )
      break
    }

    # Extract correlation matrix from lavaan once (same across levels).
    # For WLSMV + ordered: sampstat$cov IS the polychoric correlation matrix.
    # For ML + continuous: sampstat$cov is the sample covariance; cov2cor() converts.
    if (is.null(r_lv)) {
      r_lv <- tryCatch({
        sstat <- lavaan::lavInspect(fit, "sampstat")
        if (cor_type == "polychoric") sstat$cov else cov2cor(sstat$cov)
      }, error = function(e) NULL)
    }
    if (is.null(r_lv)) {
      r_lv <- stats::cor(data, method = "pearson", use = "pairwise.complete.obs")
    }

    # Standardized solution (std.all = pattern matrix on correlation scale --
    # consistent with psych::fa output and with tenBerge weights using R).
    std_sol <- tryCatch(
      lavaan::standardizedSolution(fit, type = "std.all"),
      error = function(e) NULL
    )
    if (is.null(std_sol)) {
      cli::cli_warn(
        c(
          "!" = "Could not extract standardized solution at k = {k}.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        )
      )
      break
    }

    load_rows   <- std_sol[std_sol$op == "=~", , drop = FALSE]
    factors_lav <- unique(load_rows$lhs)   # lavaan's internal factor names

    # Build p × k loading and SE matrices indexed by item_names
    L    <- matrix(0,        p, k, dimnames = list(item_names, factors_lav))
    L_se <- matrix(NA_real_, p, k, dimnames = list(item_names, factors_lav))

    for (fj in seq_along(factors_lav)) {
      frows <- load_rows[load_rows$lhs == factors_lav[fj], , drop = FALSE]
      L[frows$rhs, fj]    <- frows$est.std
      L_se[frows$rhs, fj] <- if ("se" %in% names(frows)) frows$se else NA_real_
    }

    labels_k <- make_labels(k)

    # Sort factors by descending variance explained (consistent with PCA/EFA convention)
    var_unsorted <- colSums(L^2) / p
    ord  <- order(var_unsorted, decreasing = TRUE)
    L    <- L[, ord, drop = FALSE]
    L_se <- L_se[, ord, drop = FALSE]
    colnames(L)    <- labels_k
    colnames(L_se) <- labels_k

    # Positive manifold anchor for k = 1 (matches PCA/EFA engine behaviour)
    if (k == 1L && sum(L) < 0) {
      L <- -L
      # L_se contains SEs (always non-negative); leave unchanged
    }

    # tenBerge weights from lavaan's loadings + lavaan's correlation matrix.
    # Keeps compute_edges() on the algebra path (DESIGN.md §14 item 12).
    weight_method <- "tenBerge"
    W <- tryCatch(
      .tenBerge_weights(r_lv, L),
      error = function(e) {
        weight_method <<- "regression"
        cli::cli_warn(
          c(
            "!" = "tenBerge weights failed at k = {k}: {conditionMessage(e)}",
            "i" = "Falling back to regression (Thurstone) weights."
          )
        )
        tryCatch({
          Ri    <- solve(r_lv)
          W_reg <- Ri %*% L
          colnames(W_reg) <- labels_k
          rownames(W_reg) <- item_names
          W_reg
        }, error = function(e2) NULL)
      }
    )
    if (is.null(W)) {
      cli::cli_warn(
        c(
          "!" = "Could not compute scoring weights at k = {k}.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        )
      )
      break
    }
    colnames(W) <- labels_k
    rownames(W) <- item_names

    score_var <- diag(crossprod(W, r_lv %*% W))

    # Variance explained (sum of squared standardized loadings / p)
    var_per_factor <- colSums(L^2) / p
    variance <- c(
      setNames(var_per_factor, labels_k),
      cumulative = sum(var_per_factor)
    )

    # Fit indices: chi, dof, p_value, CFI, TLI, RMSEA, SRMR
    fitmeas <- tryCatch(
      lavaan::fitMeasures(
        fit,
        c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")
      ),
      error = function(e) {
        setNames(rep(NA_real_, 7L),
                 c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
      }
    )
    fit_info <- setNames(
      as.numeric(fitmeas),
      c("chi", "dof", "p_value", "CFI", "TLI", "RMSEA", "SRMR")
    )

    # Within-level factor correlations (identity for orthogonal rotation)
    factor_cor <- tryCatch({
      Phi <- lavaan::lavInspect(fit, "cor.lv")
      if (is.matrix(Phi) && nrow(Phi) == k) {
        rownames(Phi) <- labels_k
        colnames(Phi) <- labels_k
        Phi
      } else {
        diag(k)
      }
    }, error = function(e) diag(k))

    result[[as.character(k)]] <- list(
      k           = k,
      loadings    = L,
      loadings_se = L_se,
      variance    = variance,
      fit         = fit_info,
      converged   = TRUE,
      factor_cor  = factor_cor,
      labels      = labels_k,
      scoring     = list(
        linear    = TRUE,
        method    = weight_method,
        basis     = cor_type,
        weights   = W,
        score_var = score_var
      )
    )
  }

  list(levels = result, r_lv = r_lv)
}
