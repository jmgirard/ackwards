# ESEM engine -- internal, not exported
#' @importFrom stats cov2cor
#
# Wraps lavaan::efa() for each level 1..k_max, returning the standard s.4 level
# contract. For ordinal/polychoric data uses WLSMV estimation (Kim & Eaton 2015;
# Forbush et al. 2024 use ULSMV); for continuous data uses ML. Scoring uses self-computed
# tenBerge weights from lavaan's estimated loadings and the underlying correlation
# matrix, keeping compute_edges() on the algebra path (DESIGN.md s.14 item 12).
#
# Key additions over the EFA engine:
#   - loadings_se: rotation-aware delta-method SEs from lavaan's standardized solution
#   - fit: CFI, TLI, RMSEA, SRMR from lavaan::fitMeasures()
#   - Polychoric basis: for cor = "polychoric", ordered variables trigger WLSMV
#
# Performance (M26): the sample statistics lavaan derives from the raw data --
# thresholds, the polychoric correlation matrix, and the (expensive) asymptotic
# weight matrix NACOV/WLS.V -- are IDENTICAL across all levels (they depend only
# on the data, estimator, and missing handling, never on nfactors/rotation). They
# are therefore computed ONCE at the anchor level (k = 1) and reused for every
# deeper level via lavaan's sample-stats slot argument (slot_sample_stats= in
# lavaan >= 0.7, slotSampleStats= before; see .esem_ss_argname()), which skips the recompute
# and yields bit-identical solutions. The per-level model fits are mutually
# independent and are dispatched through .esem_lapply(), which parallelises via
# future.apply when the user has set a future::plan() (serial otherwise).

# --- Parallel dispatch -------------------------------------------------------
# future.apply is an optional (Suggests) dependency. When present, the per-level
# fits run under the user's future::plan() (sequential by default -- no behaviour
# change unless the user opts in to a parallel plan). When absent, fall back to
# serial lapply(). future.seed = TRUE sets up reproducible per-task RNG streams
# (lavaan::efa() draws random rotation starts, so each level gets its own
# stream, and this also silences future's RNG advisory).
.esem_lapply <- function(X, FUN) {
  if (rlang::is_installed("future.apply")) {
    future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

# lavaan 0.7 renamed lavaan()'s slot arguments to snake_case (slotSampleStats=
# -> slot_sample_stats=) and its option parser rejects the old name outright.
# Detect by capability (is the new name a formal of lavaan::lavaan()?), not by
# version comparison, so a future rename fails the guard test rather than
# silently breaking. DESCRIPTION supports lavaan >= 0.6-13, so both names must
# keep working.
# lav_formals is a parameter (defaulting to the real thing) so both branches
# are testable regardless of which lavaan version is installed.
.esem_ss_argname <- function(lav_formals = names(formals(lavaan::lavaan))) {
  if ("slot_sample_stats" %in% lav_formals) {
    "slot_sample_stats"
  } else {
    "slotSampleStats"
  }
}

# lavaan's rotation argument(s) for one level. The unrotated k = 1 level passes
# the method name alone. Every rotation, varimax included, also turns on
# lavaan's rotation warnings (rotation.args$warn, FALSE by default), so that
# .esem_fit_one() can count the random starts that did not converge (M098).
# lavaan >= 0.7 takes rotation options inside `rotation`, a list whose first
# element is the method (its `rotation_args` is deprecated); earlier versions
# take `rotation.args`. Detected by capability (a `rotation_args` formal of
# lavaan::efa()), and efa_formals is a parameter so both branches are testable.
.esem_rotation_args <- function(rotation, efa_formals = names(formals(lavaan::efa))) {
  if (rotation == "none") {
    list(rotation = rotation)
  } else if ("rotation_args" %in% efa_formals) {
    list(rotation = list(rotation, warn = TRUE))
  } else {
    list(rotation = rotation, rotation.args = list(warn = TRUE))
  }
}

# --- Single-level fit + extraction -------------------------------------------
# Fits one ESEM level and returns the slim level contract (NOT the heavy lavaan
# fit, which would carry a duplicate NACOV per level and be costly to serialise
# back from a parallel worker).
#
#   ss_in    -- NULL at the anchor level: fit from raw data and harvest the
#               SampleStats slot (+ r_lv) for reuse. Non-NULL for deeper levels:
#               reuse via slotSampleStats= (skip the polychoric/NACOV recompute).
#   r_lv_in  -- pre-computed correlation matrix for tenBerge weights (continuous
#               paths); NULL means extract from the fit (polychoric / FIML).
#   rotation -- lavaan rotation for k >= 2 ("varimax", or the oblique
#               "oblimin" / "geomin"); validated upstream by .check_rotation().
#
# Returns one of:
#   list(status = "ok", k, level, bad_resid, ss, r_lv, fit_raw,
#        rot_failed, rot_starts)
#   list(status = "error",        k, error_msg)
#   list(status = "nonconverged", k)                          -- the model fit
#   list(status = "nonconverged", k, rot_failed, rot_starts)  -- the rotation
.esem_fit_one <- function(k, data_df, estimator, cor, ordered_cols, lav_missing,
                          ss_in, r_lv_in, item_names, p, keep_fit,
                          rotation = "varimax") {
  # No rotation needed for a single factor.
  rotate_k <- if (k == 1L) "none" else rotation

  efa_args <- list(
    data      = data_df,
    nfactors  = k,
    rotation  = rotate_k,
    estimator = estimator,
    ordered   = ordered_cols,
    missing   = lav_missing
  )
  # Reuse cached sample statistics for deeper levels (M26). lavaan uses the slot
  # and skips recomputing thresholds / polychorics / NACOV. The argument name is
  # lavaan-version-dependent (renamed in 0.7); see .esem_ss_argname().
  if (!is.null(ss_in)) efa_args[[.esem_ss_argname()]] <- ss_in
  # lavaan's rotation warnings are off by default (rotation.args$warn =
  # FALSE). Every rotation turns them on, so that its non-convergence reaches
  # the caller (M90 Decisions, finding F1; varimax since M098).
  rot_args <- .esem_rotation_args(rotate_k)
  efa_args[names(rot_args)] <- rot_args

  # lavaan emits fit-time warnings (e.g. non-PD vcov) that are not actionable
  # here; muffle them around the fit only (matches pre-M26 muffling behaviour).
  # They are recorded so that the rotation's failed starts can be counted.
  fit_warnings <- character(0)
  fit_raw <- tryCatch(
    withCallingHandlers(
      do.call(lavaan::efa, efa_args),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "esem_fit_error")
    }
  )
  if (inherits(fit_raw, "esem_fit_error")) {
    return(list(status = "error", k = k, error_msg = fit_raw$msg))
  }
  fit <- if (inherits(fit_raw, "efaList")) fit_raw[[1L]] else fit_raw
  if (is.null(fit)) { # nocov start
    return(list(status = "error", k = k, error_msg = "lavaan returned NULL fit"))
  } # nocov end

  # Convergence check
  converged <- isTRUE(tryCatch(
    lavaan::lavInspect(fit, "converged"),
    error = function(e) FALSE
  ))
  if (!converged) { # nocov start
    return(list(status = "nonconverged", k = k))
  } # nocov end

  # lavaan rotates from rotation.args$rstarts random starts (30 by default) and
  # keeps the one with the best criterion value, warning once per start that
  # did not converge. lavaan does not record whether the kept start converged.
  # When every start failed it cannot have, so the level counts as not
  # converged (Invariant 7). With rstarts = 0 lavaan rotates once, from the
  # identity, and that one start is every start. lavInspect(fit, "converged")
  # describes the model fit and stays TRUE when the rotation fails.
  rot_failed <- if (rotate_k == "none") {
    0L
  } else {
    sum(grepl("rotation algorithm did not converge", fit_warnings, fixed = TRUE))
  }
  rot_starts <- .esem_rotation_starts(fit)
  if (rot_failed >= rot_starts) {
    return(list(
      status = "nonconverged", k = k,
      rot_failed = rot_failed, rot_starts = rot_starts
    ))
  }

  # Improper solution / Heywood check: zero or negative residual variances.
  # lavaan clamps theta to 0 when residual variance would go negative (Heywood
  # case); flag at <= 0 to catch both the clamped-to-zero and unconstrained-
  # negative cases. Flag but do NOT truncate (Invariant 7) -- reported by the
  # caller via the returned count.
  theta <- tryCatch(
    lavaan::lavInspect(fit, "theta"),
    error = function(e) NULL
  )
  bad_resid <- if (!is.null(theta)) sum(diag(theta) <= 0) else 0L

  # Harvest the SampleStats slot at the anchor level for downstream reuse.
  harvested_ss <- if (is.null(ss_in)) fit@SampleStats else NULL

  # Extract correlation matrix for tenBerge weights.
  # Sources differ by path:
  # - polychoric/WLSMV: lavaan's sampstat$cov IS the polychoric matrix
  # - continuous pairwise/listwise: r_lv_in already set (passed through)
  # - continuous FIML: extract from FIML saturated model (h1) so the edge R
  #   reflects FIML-estimated sufficient statistics, not observed pairwise
  r_lv <- r_lv_in
  if (is.null(r_lv)) {
    if (lav_missing == "fiml") {
      r_lv <- tryCatch(
        {
          h1 <- lavaan::lavInspect(fit, "h1")
          # lavInspect("h1") returns list(cov = <matrix>, mean = <vector>)
          # for single-group models (the only case ackwards uses).
          # Multi-group would return list(g1 = list(cov=...), g2 = ...) --
          # guard against that with is.list(h1[[1L]]).
          cov_h1 <- if (is.list(h1[[1L]])) h1[[1L]]$cov else h1$cov
          cov2cor(cov_h1)
        },
        error = function(e) NULL
      )
    } else {
      r_lv <- tryCatch(
        {
          sstat <- lavaan::lavInspect(fit, "sampstat")
          # cov2cor branch: unreachable via ackwards() since r_lv_in is always
          # pre-computed for cor != "polychoric" (r_lv would not be NULL).
          if (cor == "polychoric") sstat$cov else cov2cor(sstat$cov) # nocov
        },
        error = function(e) NULL
      )
    }
  }
  if (is.null(r_lv)) { # nocov start
    r_lv <- stats::cor(data_df, method = "pearson", use = "pairwise.complete.obs")
  } # nocov end

  # Standardized solution (std.all = pattern matrix on correlation scale --
  # consistent with psych::fa output and with tenBerge weights using R).
  std_sol <- tryCatch(
    lavaan::standardizedSolution(fit, type = "std.all"),
    error = function(e) NULL
  )
  if (is.null(std_sol)) { # nocov start
    return(list(
      status = "error", k = k,
      error_msg = "could not extract standardized solution"
    ))
  } # nocov end

  load_rows <- std_sol[std_sol$op == "=~", , drop = FALSE]
  factors_lav <- unique(load_rows$lhs) # lavaan's internal factor names

  # Build p x k loading and SE matrices indexed by item_names
  L <- matrix(0, p, k, dimnames = list(item_names, factors_lav))
  L_se <- matrix(NA_real_, p, k, dimnames = list(item_names, factors_lav))

  for (fj in seq_along(factors_lav)) {
    frows <- load_rows[load_rows$lhs == factors_lav[fj], , drop = FALSE]
    L[frows$rhs, fj] <- frows$est.std
    L_se[frows$rhs, fj] <- if ("se" %in% names(frows)) frows$se else NA_real_
  }

  labels_k <- make_labels(k)

  # Within-level factor correlations (.esem_read_phi()). It is the identity
  # under orthogonal rotation, up to rounding, so a failed read takes the
  # identity there. An oblique fit whose correlation cannot be read is an
  # error (the level truncates with a warning), never a silent identity.
  phi_err <- NULL
  Phi_lav <- tryCatch(.esem_read_phi(fit, factors_lav), error = function(e) {
    phi_err <<- conditionMessage(e)
    NULL
  })
  if (is.null(Phi_lav)) {
    if (rotate_k %in% c("none", "varimax")) {
      Phi_lav <- diag(k)
    } else {
      return(list(
        status = "error", k = k,
        error_msg = paste0("could not extract the factor correlations: ", phi_err)
      ))
    }
  }

  # Sort factors by descending variance explained. This is the EFA engine's
  # order (psych::fa() sorts by the same key); PCA keeps psych::pca()'s order,
  # which matches it under varimax only (engine_pca.R).
  # `ord` is applied to the loadings here and to lavaan's factor correlation
  # (.carry_factor_cor), so Phi and factor_cor always match the column order
  # of L. .variance_key() is the variance explained times the constant p, in
  # lavaan's factor order: diag(Phi L'L) under an oblique rotation,
  # colSums(L^2) under varimax. The variance vector itself is computed once,
  # post-sort, from the same key (M60).
  ord <- order(.variance_key(L, Phi_lav), decreasing = TRUE)
  L <- L[, ord, drop = FALSE]
  L_se <- L_se[, ord, drop = FALSE]
  colnames(L) <- labels_k
  colnames(L_se) <- labels_k
  # Unit signs here; ackwards() applies the align_signs flips.
  Phi <- .carry_factor_cor(Phi_lav, ord, rep(1, k))

  # Positive manifold anchor for k = 1 (matches PCA/EFA engine behaviour)
  if (k == 1L && sum(L) < 0) {
    L <- -L
    # L_se contains SEs (always non-negative); leave unchanged
  }

  # tenBerge weights from lavaan's loadings, factor correlation, and
  # correlation matrix. Keeps compute_edges() on the algebra path (DESIGN.md
  # s.14 item 12).
  weight_method <- "tenBerge"
  W <- tryCatch(
    .tenBerge_weights(r_lv, L, Phi),
    error = function(e) {
      weight_method <<- "regression"
      # Oblique regression rule R^{-1} L Phi (R^{-1} L when Phi = I).
      tryCatch(.regression_weights(r_lv, L, Phi), error = function(e2) NULL)
    }
  )
  if (is.null(W)) {
    return(list(
      status = "error", k = k,
      error_msg = "could not compute scoring weights"
    ))
  }
  colnames(W) <- labels_k
  rownames(W) <- item_names

  score_var <- .score_var(W, r_lv)

  # Variance explained, diag(Phi L'L) / p (colSums(L^2) / p under varimax)
  variance <- .variance_explained(L, p, labels_k, Phi)

  # Fit indices: chi, dof, p_value, CFI, TLI, RMSEA, SRMR, BIC.
  # lavaan silently *omits* requested names that don't apply to the fitted
  # estimator (e.g. chisq.scaled under ML) rather than erroring, so requesting
  # the scaled variants alongside the naive ones is safe for every estimator.
  fitmeas <- tryCatch(
    lavaan::fitMeasures(
      fit,
      c(
        "chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "bic",
        "chisq.scaled", "df.scaled", "pvalue.scaled",
        "cfi.scaled", "tli.scaled", "rmsea.scaled"
      )
    ),
    error = function(e) NULL # nocov
  )
  .fm <- function(nm) {
    if (is.null(fitmeas) || !nm %in% names(fitmeas)) NA_real_ else unname(fitmeas[[nm]])
  }
  # Prefer the scaled variant for every scale-sensitive index whenever the
  # estimator produces one. lavaan only emits `<index>.scaled` names under a
  # scaled test (WLSMV/ULSMV mean-and-variance-adjusted; MLR Yuan-Bentler), so
  # the presence of the scaled name is the discriminator:
  #   * WLSMV/ULSMV -- the naive chi-square has no valid reference distribution
  #     (lavaan's own summary() labels its p-value "Unknown"/NA); the naive
  #     CFI/TLI/RMSEA are known to be badly optimistic (Xia & Yang 2019). Only
  #     the scaled test/indices are defensible here.
  #   * MLR -- the whole point of robust ML is the scaled (Yuan-Bentler) test;
  #     reporting the naive ML statistics for an MLR fit would defeat it.
  #   * ML -- no scaled names emitted, so this falls through to the naive
  #     values, which are the correct ones for ML.
  # SRMR has no scaled variant (unaffected). BIC needs a proper log-likelihood,
  # which WLSMV/ULSMV lack -- lavaan returns NA there (inapplicable, not a bug)
  # and a real value for ML/MLR. Keeping one rule for the whole row guarantees
  # chi/dof/p_value and CFI/TLI/RMSEA never mix scaled and naive framings.
  .fm_scaled <- function(base) {
    scaled <- paste0(base, ".scaled")
    if (!is.null(fitmeas) && scaled %in% names(fitmeas)) .fm(scaled) else .fm(base)
  }
  fit_info <- c(
    chi = .fm_scaled("chisq"),
    dof = .fm_scaled("df"),
    p_value = .fm_scaled("pvalue"),
    CFI = .fm_scaled("cfi"),
    TLI = .fm_scaled("tli"),
    RMSEA = .fm_scaled("rmsea"),
    SRMR = .fm("srmr"),
    BIC = .fm("bic")
  )

  factor_cor <- .label_phi(Phi, labels_k)

  level <- list(
    k = k,
    loadings = L,
    loadings_se = L_se,
    variance = variance,
    fit = fit_info,
    converged = TRUE,
    factor_cor = factor_cor,
    labels = labels_k,
    scoring = list(
      linear    = TRUE,
      method    = weight_method,
      basis     = cor,
      weights   = W,
      score_var = score_var
    )
  )

  # Some but not all starts failed (rot_failed > 0, counted above): the kept
  # start may have converged, so the caller reports the count and keeps the
  # level.
  list(
    status       = "ok",
    k            = k,
    level        = level,
    bad_resid    = bad_resid,
    ss           = harvested_ss,
    r_lv         = r_lv,
    fit_raw      = if (keep_fit) fit else NULL,
    rot_failed   = rot_failed,
    rot_starts   = rot_starts
  )
}

# The number of random starts lavaan rotated one level from: rotation.args$rstarts
# from the fit's options, and 1 when it is 0 (one start, from the identity). A
# value that cannot be read also counts as one start, so that a rotation
# failure there ends the hierarchy rather than being kept unseen.
.esem_rotation_starts <- function(fit) {
  rstarts <- tryCatch(
    lavaan::lavInspect(fit, "options")$rotation.args$rstarts,
    error = function(e) NULL
  )
  if (!is.numeric(rstarts) || length(rstarts) != 1L || is.na(rstarts)) {
    return(1L) # nocov
  }
  max(as.integer(rstarts), 1L)
}

# lavaan's within-level factor correlation (cor.lv) for one fit, with rows and
# columns in lavaan's loading-column order `factors_lav`. The rows are matched
# by lavaan's factor names, never by position, so a cor.lv without matching
# names is an error, as is a non-finite one. Returned unnamed.
.esem_read_phi <- function(fit, factors_lav) {
  Phi <- lavaan::lavInspect(fit, "cor.lv")
  if (is.null(rownames(Phi)) || !all(factors_lav %in% rownames(Phi))) {
    stop("lavaan's cor.lv does not carry the loading columns' factor names.")
  }
  Phi <- as.matrix(Phi[factors_lav, factors_lav, drop = FALSE])
  if (!all(is.finite(Phi))) stop("lavaan's cor.lv is not finite.")
  unname(Phi)
}

# --- Driver: fit all levels 1..k_max -----------------------------------------
esem_levels <- function(data, k_max, estimator, cor,
                        R_external = NULL, keep_fits = FALSE,
                        missing = "pairwise", rotation = "varimax") {
  rlang::check_installed("lavaan", reason = "for the ESEM engine")
  if (!exists("efa", envir = asNamespace("lavaan"), inherits = FALSE)) { # nocov start
    cli::cli_abort(
      c(
        "!" = "{.fn lavaan::efa} is not available in the installed version of lavaan.",
        "i" = "Please update lavaan: {.run install.packages('lavaan')} (>= 0.6-13 required)."
      )
    )
  } # nocov end

  p <- ncol(data)
  item_names <- colnames(data)
  data_df <- as.data.frame(data)

  # For polychoric/WLSMV: treat all columns as ordered categorical. On the
  # continuous path pass FALSE (not NULL): lavaan 0.7 requires an explicit
  # ordered = FALSE to allow WLSMV/ULSMV with continuous data, and FALSE is
  # equivalent to NULL (no ordered variables) on every other estimator/version
  # (verified bit-identical on unrotated fits under 0.6.21 and 0.7-2).
  ordered_cols <- if (cor == "polychoric") item_names else FALSE

  # missing= mapping to lavaan's vocabulary (constant across levels):
  #   "fiml"     -> "fiml"           (ML/MLR only; validated upstream)
  #   "pairwise" -> "available.cases" for WLSMV, and for ULSMV on the
  #                 polychoric basis (uses all rows via pairwise statistics --
  #                 MCAR-valid, honest N);
  #                 "listwise" for ML/MLR and for continuous ULSMV (lavaan
  #                 default; edge R is the separately-computed pairwise
  #                 stats::cor). lavaan 0.7.2 efa() errors for ULSMV on
  #                 continuous items under "available.cases", complete data
  #                 included.
  #   "listwise" -> "listwise"       (data already reduced upstream)
  pairwise_cases <- estimator == "WLSMV" ||
    (estimator == "ULSMV" && cor == "polychoric")
  lav_missing <- if (missing == "fiml") {
    "fiml"
  } else if (missing == "pairwise" && pairwise_cases) {
    "available.cases"
  } else {
    "listwise"
  }

  # --- Phase 1: anchor level (k = 1) -- harvest sample stats + r_lv -----------
  anchor <- .esem_fit_one(
    1L, data_df, estimator, cor, ordered_cols, lav_missing,
    ss_in = NULL, r_lv_in = R_external,
    item_names = item_names, p = p, keep_fit = keep_fits,
    rotation = rotation
  )
  if (anchor$status != "ok") { # nocov start
    # Anchor failed: nothing to build on. Warn and return empty; ackwards()
    # aborts on k_eff < 2. Mirrors the pre-M26 break-at-k=1 behaviour.
    if (anchor$status == "error") {
      cli::cli_warn(c(
        "!" = "ESEM failed at k = 1: {anchor$error_msg}",
        "i" = "No levels could be built."
      ))
    } else {
      cli::cli_warn(c(
        "!" = "ESEM did not converge at k = 1.",
        "i" = "No levels could be built."
      ))
    }
    return(list(levels = list(), r_lv = R_external, fits = NULL))
  } # nocov end
  ss <- anchor$ss
  r_lv <- anchor$r_lv

  # --- Phase 2: deeper levels 2..k_max -- independent, reuse cached stats -----
  # Dispatched through .esem_lapply (parallel under a future::plan(), else serial).
  rest <- if (k_max >= 2L) {
    .esem_lapply(seq.int(2L, k_max), function(k) {
      .esem_fit_one(
        k, data_df, estimator, cor, ordered_cols, lav_missing,
        ss_in = ss, r_lv_in = r_lv,
        item_names = item_names, p = p, keep_fit = keep_fits,
        rotation = rotation
      )
    })
  } else {
    list() # nocov  (k_max >= 2 enforced upstream)
  }

  all_res <- c(list(anchor), rest) # ordered k = 1, 2, ...

  # --- Assembly: walk in order, truncate at first failure (Invariant 7) ------
  # All cli warnings are emitted here, in the main process, in deterministic
  # level order -- the parallel workers never signal conditions themselves.
  result <- list()
  fits_list <- if (keep_fits) list() else NULL
  for (res in all_res) {
    k <- res$k
    if (res$status == "error") {
      cli::cli_warn(c(
        "!" = "ESEM failed at k = {k}: {res$error_msg}",
        "i" = "Truncating hierarchy at level {k - 1L}."
      ))
      break
    }
    if (res$status == "nonconverged") {
      if (!is.null(res$rot_starts)) {
        cli::cli_warn(c(
          "!" = "lavaan's rotation did not converge at k = {k}: \\
                 {res$rot_failed} of {res$rot_starts} random starts did not \\
                 converge, so the kept rotation did not converge.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        ))
      } else {
        cli::cli_warn(c( # nocov start
          "!" = "ESEM did not converge at k = {k}.",
          "i" = "Truncating hierarchy at level {k - 1L}."
        )) # nocov end
      }
      break
    }
    # status == "ok"
    if (res$bad_resid > 0L) {
      cli::cli_warn(
        c(
          "!" = "Improper solution at k = {k}: \\
                 {res$bad_resid} residual variance{?s} \\
                 at or below zero (Heywood case).",
          "i" = "Results may be unreliable. Consider reducing \\
                 {.arg k}, changing {.arg estimator}, \\
                 or inspecting the solution."
        )
      )
    }
    if (res$rot_failed > 0L) {
      cli::cli_warn(c(
        "!" = "lavaan reported a rotation problem at k = {k}: \\
               {res$rot_failed} of {res$rot_starts} random starts did not \\
               converge.",
        "i" = "lavaan keeps the start with the best criterion value, which \\
               can be one that did not converge. The level is kept."
      ))
    }
    if (identical(res$level$scoring$method, "regression")) {
      cli::cli_warn(c(
        "!" = "tenBerge weights failed at k = {k}; \\
               used regression (Thurstone) weights instead."
      ))
    }
    result[[as.character(k)]] <- res$level
    if (keep_fits) fits_list[[as.character(k)]] <- res$fit_raw
  }

  list(levels = result, r_lv = r_lv, fits = fits_list)
}
