#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy an ackwards object into a long data frame
#'
#' Returns structured data from an `ackwards` object in tidy format. The
#' default (`what = "edges"`) returns the graph edge list that drives diagrams.
#'
#' @param x An `ackwards` object.
#' @param what What to extract:
#'   * `"edges"` *(default)* -- one row per directed between-level edge:
#'     `from`, `to`, `level_from`, `level_to`, `r`, `is_primary`, `above_cut`.
#'   * `"loadings"` -- one row per item x factor x level:
#'     `level`, `factor`, `item`, `loading`, `se`, `ci_lower`, `ci_upper`.
#'     `se`, `ci_lower`, and `ci_upper` are populated only for
#'     `engine = "esem"` (which produces rotation-aware loading SEs); they are
#'     `NA` for PCA and EFA. The confidence level is controlled by `conf_level`.
#'   * `"variance"` -- one row per factor x level:
#'     `level`, `factor`, `proportion`, `cumulative`. Both are proportions of
#'     total item variance on a 0-1 scale (multiply by 100 for a percentage).
#'   * `"fit"` -- one row per fit statistic x level: `level`, `statistic`,
#'     `value`. For PCA objects the statistics are eigenvalues; for EFA
#'     objects they are `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC` --
#'     where `chi` is the likelihood-ratio chi-square ([psych::fa()]'s
#'     `STATISTIC`), so `chi`, `dof`, `p_value`, `RMSEA`, and `TLI` all share
#'     one statistical framing (psych's residual-based *empirical* chi-square
#'     is a different statistic and is not reported); for
#'     ESEM they are `chi`, `dof`, `p_value`, `CFI`, `TLI`, `RMSEA`, `SRMR`,
#'     `BIC`. For ESEM under a scaled-test estimator (`"WLSMV"`/`"ULSMV"` for
#'     ordinal data, `"MLR"` for continuous), the whole row -- `chi`/`dof`/
#'     `p_value` **and** `CFI`/`TLI`/`RMSEA` -- reports lavaan's
#'     mean-and-variance-adjusted ("scaled") variant, so every quantity shares
#'     one scaling. This matters most for WLSMV/ULSMV: the naive chi-square
#'     has no valid reference distribution (lavaan's own `summary()` labels
#'     its p-value "Unknown"), and the naive `CFI`/`TLI` are badly optimistic
#'     for ordinal data (Xia & Yang, 2019). `"ML"` has no scaled variant, so
#'     it reports the naive values (the correct ones for ML). `SRMR` has no
#'     scaled variant and is reported as-is. `BIC` is `NA` under WLSMV/ULSMV
#'     (no proper log-likelihood for a limited-information estimator) and
#'     populated under ML/MLR. Use `format = "wide"` for one row per
#'     **non-anchor** level (k >= 2; the saturated 1-factor anchor is dropped,
#'     matching `summary()` and `autoplot(what = "fit")`), one column per
#'     statistic. Conventional fit cutoffs (Hu & Bentler 1999) are shown as
#'     reference lines in `autoplot(what = "fit")` and inline in `summary()`,
#'     but are not returned as a pass/fail column here -- they are contested
#'     thresholds, report-only, and never gate anything (see those functions'
#'     docs). `format` is oriented to the EFA/ESEM model-fit statistics; for
#'     PCA the "statistics" are per-component eigenvalues.
#'   * `"nodes"` -- Forbes-extension pruning annotations (requires `prune != "none"`
#'     when the object was created). One row per factor across all levels:
#'     `id`, `level`, `pruned`, `prune_reason`. Returns an empty data frame with
#'     the same columns when no pruning was applied.
#'   * `"scores"` -- long-format per-observation factor scores (requires
#'     `keep_scores = TRUE` at fit time or use [augment.ackwards()] for on-the-fly
#'     computation). Columns: `obs` (row index), `level`, `factor`, `score`.
#' @param primary_only For `what = "edges"` only. When `TRUE`, returns just each
#'   factor's primary-parent edge (`is_primary == TRUE`) -- the lineage tree that
#'   the diagram draws as solid arrows. Default `FALSE` (all edges). Errors for
#'   any other value of `what`.
#' @param sort For `what = "edges"` only. One of `"none"` (default, natural order)
#'   or `"strength"` (descending `|r|`). Ignored for all other values of `what`.
#' @param format For `what = "fit"` only. One of `"long"` (default, one row per
#'   statistic x level) or `"wide"` (one row per level, one column per
#'   statistic). Errors for all other values of `what`.
#' @param conf_level For `what = "loadings"` only. Confidence level for the
#'   loading intervals; default `0.95`. The intervals are computed as
#'   `loading ± qnorm((1 + conf_level) / 2) * se` and are `NA` for engines
#'   that carry no SEs (PCA, EFA). Errors for all other values of `what`.
#' @param ... Ignored.
#'
#' @return A data frame (class `data.frame`).
#'
#' @seealso [glance.ackwards()], [print.ackwards()]
#'
#' @examples
#' x <- ackwards(bfi25, k_max = 5)
#' tidy(x) # edges in natural order
#' tidy(x, sort = "strength") # strongest edges first
#' tidy(x, primary_only = TRUE) # just the primary-parent lineage
#' tidy(x, what = "loadings")
#' tidy(x, what = "variance")
#' tidy(x, what = "fit")
#' tidy(x, what = "fit", format = "wide")
#'
#' @export
tidy.ackwards <- function(
  x,
  what = c("edges", "loadings", "variance", "fit", "nodes", "scores"),
  primary_only = FALSE,
  sort = c("none", "strength"),
  format = c("long", "wide"),
  conf_level = 0.95,
  ...
) {
  what <- match.arg(what)
  sort <- match.arg(sort)
  format <- match.arg(format)
  if (sort != "none" && what != "edges") {
    cli::cli_abort(
      "{.arg sort} is only supported for {.code what = \"edges\"}, \\
       not {.code what = \"{what}\"}."
    )
  }
  if (isTRUE(primary_only) && what != "edges") {
    cli::cli_abort(
      "{.arg primary_only} is only supported for {.code what = \"edges\"}, \\
       not {.code what = \"{what}\"}."
    )
  }
  if (format != "long" && what != "fit") {
    cli::cli_abort(
      "{.arg format} is only supported for {.code what = \"fit\"}, \\
       not {.code what = \"{what}\"}."
    )
  }
  if (!identical(conf_level, 0.95) && what != "loadings") {
    cli::cli_abort(
      "{.arg conf_level} is only supported for {.code what = \"loadings\"}, \\
       not {.code what = \"{what}\"}."
    )
  }
  out <- switch(what,
    edges    = .tidy_edges(x),
    loadings = .tidy_loadings(x, conf_level = conf_level),
    variance = .tidy_variance(x),
    fit      = .tidy_fit(x),
    nodes    = .tidy_nodes(x),
    scores   = .tidy_scores(x)
  )
  if (what == "edges") {
    if (isTRUE(primary_only)) {
      out <- out[out$is_primary, , drop = FALSE]
    }
    if (sort == "strength") {
      out <- out[order(abs(out$r), decreasing = TRUE), , drop = FALSE]
    }
    rownames(out) <- NULL
  }
  if (what == "fit" && format == "wide") {
    out <- .fit_long_to_wide(out)
  }
  out
}

.tidy_edges <- function(x) {
  x$edges$tidy
}

.tidy_nodes <- function(x) {
  if (!is.null(x$prune) && !is.null(x$prune$nodes)) {
    return(x$prune$nodes)
  }
  data.frame(
    id = character(0L),
    level = integer(0L),
    pruned = logical(0L),
    prune_reason = character(0L),
    stringsAsFactors = FALSE
  )
}

.tidy_scores <- function(x) {
  if (is.null(x$scores)) {
    cli::cli_abort(c(
      "!" = "Factor scores are not stored in this {.cls ackwards} object.",
      "i" = "Refit with {.code keep_scores = TRUE}, or use {.fn augment} to \\
             compute scores on the fly: {.code augment(x, data = your_data)}."
    ))
  }
  rows <- lapply(names(x$scores), function(ki) {
    S <- x$scores[[ki]]
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(S)), function(j) {
      data.frame(
        obs = seq_len(nrow(S)),
        level = k,
        factor = colnames(S)[j],
        score = S[, j],
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.tidy_loadings <- function(x, conf_level = 0.95) {
  z <- stats::qnorm((1 + conf_level) / 2)
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    L <- lev$loadings
    SE <- lev$loadings_se # NULL for PCA/EFA
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(L)), function(j) {
      se_col <- if (!is.null(SE)) SE[, j] else rep(NA_real_, nrow(L))
      loading_col <- L[, j]
      data.frame(
        level = k,
        factor = colnames(L)[j],
        item = rownames(L),
        loading = loading_col,
        se = se_col,
        ci_lower = loading_col - z * se_col,
        ci_upper = loading_col + z * se_col,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.tidy_fit <- function(x) {
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    k <- as.integer(ki)
    fv <- lev$fit
    data.frame(
      level = k,
      statistic = names(fv),
      value = unname(fv),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Hu & Bentler (1999) conventional thresholds — direction "hi" means higher is
# better (CFI/TLI), "lo" means lower is better (RMSEA/SRMR). Statistics not in
# this list (chi, dof, p_value, BIC, eigenvalues) have no defined threshold.
# Used internally as reference lines/inline annotations (autoplot, summary);
# not exposed as a pass/fail column in tidy() -- these thresholds are
# conventional and contested, report-only, and never gate anything.
.fit_cutoffs <- function() {
  list(
    CFI   = list(threshold = 0.95, direction = "hi"),
    TLI   = list(threshold = 0.95, direction = "hi"),
    RMSEA = list(threshold = 0.06, direction = "lo"),
    SRMR  = list(threshold = 0.08, direction = "lo")
  )
}

.fit_long_to_wide <- function(df) {
  # Exclude the anchor level (k = 1). The 1-factor anchor is the saturated
  # baseline; summary() and autoplot(what = "fit") both drop it, so the wide
  # reporting table matches them (one row per non-anchor level).
  df <- df[df$level > 1L, , drop = FALSE]
  levels <- sort(unique(df$level))
  all_stat <- unique(df$statistic) # preserve original statistic order
  rows <- lapply(levels, function(lv) {
    sub <- df[df$level == lv, , drop = FALSE]
    row_vals <- list(level = lv)
    for (stat in all_stat) {
      m <- match(stat, sub$statistic)
      row_vals[[stat]] <- if (!is.na(m)) sub$value[m] else NA_real_
    }
    as.data.frame(row_vals, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.tidy_variance <- function(x) {
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    k <- as.integer(ki)
    fac_labels <- lev$labels
    var_vals <- lev$variance[fac_labels]
    data.frame(
      level = k,
      factor = fac_labels,
      proportion = var_vals,
      cumulative = cumsum(var_vals),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Glance at an ackwards object
#'
#' Returns a one-row data frame of top-level model metadata. For EFA and ESEM
#' objects, fit indices at the deepest converged level are included. The same
#' five columns (`CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`) are present across all
#' engines; columns unavailable for a given engine or estimator are `NA`
#' (e.g., `CFI` and `SRMR` are `NA` for EFA; all five are `NA` for PCA; for
#' ESEM, `BIC` is `NA` under `estimator = "WLSMV"`/`"ULSMV"` -- these
#' limited-information estimators have no proper log-likelihood -- and
#' populated under `"ML"`/`"MLR"`). Under a scaled-test estimator
#' (`"WLSMV"`/`"ULSMV"`/`"MLR"`) the `CFI`/`TLI`/`RMSEA` reported here are the
#' scaled variants (see [tidy.ackwards()] for the rationale).
#'
#' @param x An `ackwards` object.
#' @param ... Ignored.
#'
#' @return A one-row `data.frame`.
#'
#' @seealso [tidy.ackwards()], [print.ackwards()]
#'
#' @examples
#' x <- ackwards(bfi25, k_max = 5)
#' glance(x)
#'
#' @export
glance.ackwards <- function(x, ...) {
  fc <- .glance_fit(x)
  data.frame(
    engine            = x$engine,
    rotation          = x$rotation,
    cor               = x$cor,
    k_max             = x$k_max,
    n_obs             = x$n_obs,
    deepest_converged = x$meta$deepest_converged,
    n_edges           = nrow(x$edges$tidy),
    CFI               = fc$CFI,
    TLI               = fc$TLI,
    RMSEA             = fc$RMSEA,
    SRMR              = fc$SRMR,
    BIC               = fc$BIC,
    stringsAsFactors  = FALSE
  )
}

# Extract fit indices from deepest converged non-anchor level.
# Returns a one-row data.frame with columns CFI, TLI, RMSEA, SRMR, BIC.
# Missing indices for a given engine are NA; PCA returns all NA.
.glance_fit <- function(x) {
  dc <- x$meta$deepest_converged
  na_row <- data.frame(
    CFI = NA_real_, TLI = NA_real_, RMSEA = NA_real_,
    SRMR = NA_real_, BIC = NA_real_
  )
  if (is.null(dc) || dc < 2L) {
    return(na_row)
  }
  fv <- x$levels[[as.character(dc)]]$fit
  if (is.null(fv) || length(fv) == 0L) {
    return(na_row)
  }
  data.frame(
    CFI   = if ("CFI" %in% names(fv)) fv[["CFI"]] else NA_real_,
    TLI   = if ("TLI" %in% names(fv)) fv[["TLI"]] else NA_real_,
    RMSEA = if ("RMSEA" %in% names(fv)) fv[["RMSEA"]] else NA_real_,
    SRMR  = if ("SRMR" %in% names(fv)) fv[["SRMR"]] else NA_real_,
    BIC   = if ("BIC" %in% names(fv)) fv[["BIC"]] else NA_real_
  )
}
