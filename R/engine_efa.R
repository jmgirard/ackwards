# EFA engine -- internal, not exported
#' @importFrom stats setNames
#
# Wraps psych::fa() for each level 1..k_max, returning the standard s.4 level
# contract. Scoring uses tenBerge weights (linear S = ZW), which keeps
# compute_edges() on the algebra path. Convergence failures truncate the
# hierarchy at the last successful level; Heywood cases warn but continue.

efa_levels <- function(R, k_max, fm, n_obs, cor = "pearson",
                       keep_fits = FALSE, rotation = "varimax") {
  p <- nrow(R)
  result <- list()
  fits_list <- if (keep_fits) list() else NULL

  for (k in seq_len(k_max)) {
    rotate_k <- if (k == 1L) "none" else rotation

    # Run psych::fa(), intercepting convergence warnings so we can act on them.
    # suppressMessages() muffles psych's message-stream chatter (e.g. the
    # per-level "determinant of the smoothed correlation was zero" / "smcs < 0"
    # notes it prints on an ill-conditioned matrix) -- ackwards raises a single
    # clear near-singular warning at matrix-construction time instead.
    warn_msgs <- character(0L)
    fit <- tryCatch(
      withCallingHandlers(
        suppressMessages(
          psych::fa(R, nfactors = k, rotate = rotate_k, fm = fm, n.obs = n_obs)
        ),
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

    # An oblique rotation's failure (M90 Decisions, finding F1). psych::fa()
    # rotates from 20 random starts (n.rotations) and keeps the best, so a
    # GPArotation warning can come from a discarded start: re-raise it and
    # keep the level. A failed final step leaves no rotation matrix, and that
    # level is not the requested rotation, so it truncates (Invariant 7).
    if (k > 1L && .is_oblique(rotation)) {
      rot_msgs <- .psych_rotation_warnings(warn_msgs)
      if (is.null(fit$rot.mat)) {
        detail <- if (length(rot_msgs) > 0L) {
          rot_msgs[[1L]]
        } else {
          "psych returned no rotation matrix."
        }
        cli::cli_warn(c(
          "!" = "The {rotation} rotation failed at k = {k}: {detail}",
          "i" = "Truncating hierarchy at level {k - 1L}."
        ))
        break
      }
      if (length(rot_msgs) > 0L) {
        cli::cli_warn(c(
          "!" = "psych reported a rotation problem at k = {k}: {rot_msgs[[1L]]}",
          "i" = "psych rotates from 20 random starts and keeps the best one, \\
                 so the warning can come from a discarded start. The level is kept."
        ))
      }
    }

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

    # Within-level factor correlation. psych sets $Phi only under an oblique
    # rotation and sorts and sign-flips it with the loadings; varimax leaves it
    # absent and .engine_phi() returns the identity. The weights, the variance,
    # and factor_cor all read this one matrix. A non-finite Phi truncates.
    Phi_k <- tryCatch(
      .engine_phi(fit, k),
      error = function(e) {
        cli::cli_warn(c(
          "!" = "EFA failed at k = {k}: {conditionMessage(e)}",
          "i" = "Truncating hierarchy at level {k - 1L}."
        ))
        NULL
      }
    )
    if (is.null(Phi_k)) break

    # tenBerge weights carry Phi, so the scores reproduce the factor
    # correlation (identity under varimax: uncorrelated, unit variance). The
    # algebra path in compute_edges() stays valid either way.
    weight_method <- "tenBerge"
    W <- tryCatch(
      .tenBerge_weights(R, L_rot, Phi_k),
      error = function(e) {
        weight_method <<- "regression" # honest label on fallback (Invariant 6)
        cli::cli_warn(
          c(
            "!" = "tenBerge weights failed at k = {k}: {conditionMessage(e)}",
            "i" = "Falling back to regression (Thurstone) weights."
          )
        )
        # psych's regression weights are R^{-1} L Phi, the oblique regression
        # rule (R^{-1} L under varimax).
        w_fall <- unclass(fit$weights)
        # nocov: the k = 1 positive-manifold flip (see `flip` above), never
        # an oblique branch, because k = 1 is not rotated.
        if (flip) w_fall <- -w_fall # nocov
        colnames(w_fall) <- labels_k
        rownames(w_fall) <- rownames(R)
        w_fall
      }
    )

    # Score variances: diag(W' R W); exact 1 for tenBerge, but always compute
    # rather than assume -- Invariant 1.
    score_var <- .score_var(W, R)

    # Variance explained, diag(Phi L'L) / p (colSums(L^2) / p under varimax)
    variance <- .variance_explained(L_rot, p, labels_k, Phi_k)

    # Fit indices -- available only when n.obs was supplied; NA otherwise.
    # fit$STATISTIC/dof/PVAL/TLI/BIC are plain scalars; fit$RMSEA is a named
    # vector. `chi` is psych's $STATISTIC -- the likelihood-ratio chi-square,
    # the statistic that $PVAL, $RMSEA, and $TLI are all derived from -- so the
    # whole fit row shares one statistical framing (M42/C1). psych's $chi (the
    # residual-based *empirical* chi-square) is a different statistic; pairing
    # it with $PVAL misreports, so it is deliberately not used here.
    fit_info <- setNames(
      c(
        if (!is.null(fit$STATISTIC)) unname(fit$STATISTIC)[[1L]] else NA_real_,
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
      # psych's order is already the stored order, so the carry order is the
      # identity and the signs are unit (ackwards() applies align_signs).
      factor_cor = .label_phi(
        .carry_factor_cor(Phi_k, seq_len(k), rep(1, k)),
        labels_k
      ),
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

# Compute tenBerge factor-score weights from a correlation matrix R, a
# (rotated) pattern matrix L, and the factor correlation Phi of L's columns.
#
# Formula (ten Berge, Krijnen, Wansbeek & Shapiro 1999, Eq. 3 with Eq. 9's
# C, Thm 1), with L* = L Phi^{1/2}:
#   W = R^{-1/2} C Phi^{1/2},  C = R^{-1/2} L* (L*' R^{-1} L*)^{-1/2}
#     = R^{-1} L* (L*' R^{-1} L*)^{-1/2} Phi^{1/2}
# For a full-rank L the scores reproduce Phi: W'RW = Phi, so they have unit
# variance and correlate as the factors do. At Phi = I this is the orthogonal
# formula W = R^{-1} L (L' R^{-1} L)^{-1/2} (all correlation-preserving
# methods coincide there, the paper's Thm 3). A Phi within 1e-12 of the
# identity takes that formula unchanged: varimax, including lavaan's orthogonal
# cor.lv with its rounding-level off-diagonals, keeps its pre-oblique
# floating-point path. The D standardization in compute_edges() is still
# applied for numerical safety and to satisfy Invariant 1. If L is
# rank-deficient (two factors collinear in the R^{-1} metric -- degenerate; no
# shipped engine emits this) B = L*'R^{-1}L* is singular, W'RW is no longer
# Phi, and we warn: compute_edges() still standardizes by the *actual* score
# SDs so edges stay valid, but the factors are poorly separated. A Phi that is
# not positive definite errors, and the engines fall back to regression
# weights with a warning.
.tenBerge_weights <- function(R, L, Phi) {
  k <- ncol(L)
  dn <- dimnames(L)
  Phi <- as.matrix(Phi)
  stopifnot(nrow(Phi) == k, ncol(Phi) == k)
  oblique <- !.near_identity(Phi)
  if (oblique) {
    eig_phi <- eigen(Phi, symmetric = TRUE)
    if (min(eig_phi$values) <= 0) {
      stop("the factor correlation matrix is not positive definite.")
    }
    Phi_half <- eig_phi$vectors %*%
      diag(sqrt(eig_phi$values), nrow = k) %*%
      t(eig_phi$vectors)
    L <- L %*% Phi_half # L* = L Phi^{1/2}
  }

  Ri <- solve(R) # p x p
  A <- Ri %*% L # p x k: R^{-1} L*
  B <- crossprod(L, A) # k x k: L*' R^{-1} L*  (symmetric PD for full-rank L)

  # Matrix inverse square root of B via spectral decomposition. Clamp
  # eigenvalues below a *relative* tolerance (fp noise scales with |B|, so an
  # absolute floor misses near-zeros when |B| is large under a near-singular R).
  eig <- eigen(B, symmetric = TRUE)
  tol <- .Machine$double.eps * max(eig$values)
  n_clamp <- sum(eig$values < tol)
  if (n_clamp > 0L) {
    cli::cli_warn(
      c(
        "!" = "A factor level is near rank-deficient: {n_clamp} eigenvalue{?s} \\
               of {.code L' R^-1 L} {?is/are} numerically zero.",
        "i" = "ten Berge scores for this level neither have unit variance nor \\
               reproduce the factor correlations; \\
               {.fn compute_edges} still standardizes by the actual score SDs \\
               (edges stay valid), but the factors are poorly separated."
      ),
      .frequency = "once",
      .frequency_id = "ackwards_tenberge_rank_deficient"
    )
  }
  vals <- pmax(eig$values, tol) # guard against zero / tiny-negative eigenvalues
  Binvsqrt <- eig$vectors %*%
    diag(1 / sqrt(vals), nrow = length(vals)) %*%
    t(eig$vectors)

  W <- A %*% Binvsqrt
  if (oblique) W <- W %*% Phi_half
  dimnames(W) <- dn
  W
}
