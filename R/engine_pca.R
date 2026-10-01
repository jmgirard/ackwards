# PCA engine -- internal, not exported
#' @importFrom stats setNames
#
# Uses psych::pca() (varimax by default, or the oblique `rotation`), which is
# the same function psych::bassAckward() uses internally -- ensures our PCA
# path matches psych's reference implementation within floating-point
# tolerance. A failed oblique rotation or a non-finite component correlation
# truncates the hierarchy at the level before it (Invariant 7).
#
# Returns list(levels = <named list per s.4 contract>, fits = <named list | NULL>)

pca_levels <- function(R, k_max, cor = "pearson", keep_fits = FALSE,
                       rotation = "varimax") {
  p <- nrow(R)
  result <- list()
  fits_list <- if (keep_fits) list() else NULL

  for (k in seq_len(k_max)) {
    if (k == 1L) {
      # k=1: no rotation (single component is already determined)
      fit <- psych::pca(R, nfactors = 1L, rotate = "none")
      L_rot <- unclass(fit$loadings)
      # Ensure positive manifold
      if (sum(L_rot) < 0) { # nocov start
        L_rot <- -L_rot
        fit$weights <- -fit$weights
      } # nocov end
    } else {
      # psych::pca() rotates from one start (n.rotations = 1), so a rotation
      # warning describes the stored solution: GPArotation did not converge,
      # or psych used Promax in place of the requested rotation. That level
      # truncates (M90 Decisions, finding F1). Other warnings pass through.
      # An error inside an oblique fit (a singular rotation, say) truncates
      # too (Invariant 7). Under varimax it propagates, as before.
      rot_msgs <- character(0)
      fit <- tryCatch(
        withCallingHandlers(
          psych::pca(R, nfactors = k, rotate = rotation),
          warning = function(w) {
            if (length(.psych_rotation_warnings(conditionMessage(w))) > 0L) {
              rot_msgs <<- c(rot_msgs, conditionMessage(w))
              invokeRestart("muffleWarning")
            }
          }
        ),
        error = function(e) {
          if (!.is_oblique(rotation)) stop(e)
          cli::cli_warn(c(
            "!" = "PCA failed at k = {k}: {conditionMessage(e)}",
            "i" = "Truncating hierarchy at level {k - 1L}."
          ))
          NULL
        }
      )
      if (is.null(fit)) break
      if (length(rot_msgs) > 0L) {
        cli::cli_warn(c(
          "!" = "The {rotation} rotation failed at k = {k}: {rot_msgs[[1L]]}",
          "i" = "Truncating hierarchy at level {k - 1L}."
        ))
        break
      }
      L_rot <- unclass(fit$loadings)
    }

    # Component weights. Under an oblique rotation psych returns R^{-1} L Phi
    # (fa.stats() multiplies the pattern by Phi; oblique.scores matters only
    # for raw-data input), the exact components, whose correlation is Phi;
    # under varimax that is R^{-1} L.
    W <- unclass(fit$weights)

    # psych::pca() already sorts the columns, by colSums(L^2). Under varimax
    # that is the variance order. Under an oblique rotation the stored
    # variance is diag(Phi L'L) / p, which can fall out of descending order;
    # psych's order is kept so columns match psych and Forbes's output.
    # Label columns with our stable m{k}f{j} scheme.
    colnames(L_rot) <- make_labels(k)
    rownames(L_rot) <- rownames(R)
    colnames(W) <- make_labels(k)
    rownames(W) <- rownames(R)

    # Score variances: diag(W' R W) -- NOT assumed to be 1
    score_var <- .score_var(W, R)

    # Within-level component correlation: psych's $Phi under an oblique
    # rotation, the identity under varimax (.engine_phi()). A non-finite Phi
    # truncates.
    Phi_k <- tryCatch(
      .engine_phi(fit, k),
      error = function(e) {
        cli::cli_warn(c(
          "!" = "PCA failed at k = {k}: {conditionMessage(e)}",
          "i" = "Truncating hierarchy at level {k - 1L}."
        ))
        NULL
      }
    )
    if (is.null(Phi_k)) break

    # Variance explained per component and cumulative: diag(Phi L'L) / p,
    # colSums(L^2) / p under varimax
    variance <- .variance_explained(L_rot, p, make_labels(k), Phi_k)

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
      # psych sets $Phi only under an oblique rotation; varimax leaves it
      # absent and .engine_phi() returns the identity. psych's column sort is
      # already applied to both loadings and Phi, so the carry order is the
      # identity and the signs are unit (ackwards() applies align_signs).
      factor_cor = .label_phi(
        .carry_factor_cor(Phi_k, seq_len(k), rep(1, k)),
        make_labels(k)
      ),
      labels = make_labels(k),
      scoring = list(
        linear    = TRUE,
        method    = "components",
        basis     = cor,
        weights   = W,
        score_var = score_var
      )
    )
    if (keep_fits) fits_list[[as.character(k)]] <- fit
  }

  list(levels = result, fits = fits_list)
}
