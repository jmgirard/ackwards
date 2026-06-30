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
#'     `level`, `factor`, `variance_pct`, `cumulative_pct`.
#'   * `"fit"` -- one row per fit index x level: `level`, `index`, `value`.
#'     For PCA objects the indices are eigenvalues; for EFA objects they are
#'     `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC`; for ESEM they are
#'     `chi`, `dof`, `p_value`, `CFI`, `TLI`, `RMSEA`, `SRMR`. Use
#'     `format = "wide"` for one row per **non-anchor** level (k >= 2; the
#'     saturated 1-factor anchor is dropped, matching `summary()` and
#'     `autoplot(what = "fit")`), one column per index. Add `cutoffs = TRUE` to
#'     append a `meets` column flagging each index against conventional
#'     thresholds (Hu & Bentler 1999: CFI/TLI >= .95, RMSEA <= .06,
#'     SRMR <= .08); indices without a defined threshold (e.g. `chi`, `BIC`,
#'     eigenvalues) return `NA` for `meets`. Thresholds are conventional and
#'     contested; they are report-only and never gate anything. `format` and
#'     `cutoffs` are oriented to the EFA/ESEM model-fit indices; for PCA the
#'     "indices" are per-component eigenvalues.
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
#'   index x level) or `"wide"` (one row per level, one column per index).
#'   Errors for all other values of `what`.
#' @param cutoffs For `what = "fit"` only. When `TRUE`, appends a logical `meets`
#'   column indicating whether each index meets its conventional threshold (Hu &
#'   Bentler 1999). `NA` when no threshold is defined for that index. Default
#'   `FALSE`. Errors for all other values of `what`.
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
#' tidy(x, what = "fit", cutoffs = TRUE)
#'
#' @export
tidy.ackwards <- function(
  x,
  what = c("edges", "loadings", "variance", "fit", "nodes", "scores"),
  primary_only = FALSE,
  sort = c("none", "strength"),
  format = c("long", "wide"),
  cutoffs = FALSE,
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
  if (!isFALSE(cutoffs) && what != "fit") {
    cli::cli_abort(
      "{.arg cutoffs} is only supported for {.code what = \"fit\"}, \\
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
  if (what == "fit") {
    if (isTRUE(cutoffs)) out <- .flag_fit(out)
    if (format == "wide") out <- .fit_long_to_wide(out)
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
      index = names(fv),
      value = unname(fv),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Hu & Bentler (1999) conventional thresholds — direction "hi" means higher is
# better (CFI/TLI), "lo" means lower is better (RMSEA/SRMR). Indices not in
# this list (chi, dof, p_value, BIC, eigenvalues) receive NA for `meets`.
.fit_cutoffs <- function() {
  list(
    CFI   = list(threshold = 0.95, direction = "hi"),
    TLI   = list(threshold = 0.95, direction = "hi"),
    RMSEA = list(threshold = 0.06, direction = "lo"),
    SRMR  = list(threshold = 0.08, direction = "lo")
  )
}

.flag_fit <- function(df) {
  cuts <- .fit_cutoffs()
  meets <- vapply(seq_len(nrow(df)), function(i) {
    cut <- cuts[[df$index[i]]]
    if (is.null(cut) || is.na(df$value[i])) {
      return(NA)
    }
    if (cut$direction == "hi") df$value[i] >= cut$threshold else df$value[i] <= cut$threshold
  }, logical(1L))
  df$meets <- meets
  df
}

.fit_long_to_wide <- function(df) {
  # Exclude the anchor level (k = 1). The 1-factor anchor is the saturated
  # baseline; summary() and autoplot(what = "fit") both drop it, so the wide
  # reporting table matches them (one row per non-anchor level).
  df <- df[df$level > 1L, , drop = FALSE]
  has_meets <- "meets" %in% names(df)
  levels <- sort(unique(df$level))
  all_idx <- unique(df$index) # preserve original index order
  rows <- lapply(levels, function(lv) {
    sub <- df[df$level == lv, , drop = FALSE]
    row_vals <- list(level = lv)
    for (idx in all_idx) {
      m <- match(idx, sub$index)
      row_vals[[idx]] <- if (!is.na(m)) sub$value[m] else NA_real_
      if (has_meets) {
        row_vals[[paste0(idx, "_meets")]] <- if (!is.na(m)) sub$meets[m] else NA
      }
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
      variance_pct = round(var_vals * 100, 2),
      cumulative_pct = round(
        cumsum(var_vals) * 100, 2
      ),
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
#' engines; columns unavailable for a given engine are `NA` (e.g., `CFI` and
#' `SRMR` are `NA` for EFA; all five are `NA` for PCA).
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
