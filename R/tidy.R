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
#' A factor is a summary variable standing in for a group of items that move
#' together, and a loading is the correlation between an item and a factor.
#'
#' @param x An `ackwards` object.
#' @param what What to extract:
#'   * `"edges"` *(default)*: one row per directed between-level edge, with
#'     columns
#'     `from`, `to`, `level_from`, `level_to`, `r`, `is_primary`, `above_cut`.
#'     If [boot_edges()] has been run on the object, four bootstrap columns are
#'     appended: `se`, `lo`, `hi` (bootstrap standard error and percentile
#'     confidence-interval endpoints), and `n_boot_ok` (usable replicates).
#'   * `"loadings"`: one row per item x factor x level, with columns
#'     `level`, `factor`, `item`, `loading`, `se`, `ci_lower`, `ci_upper`.
#'     The columns `se`, `ci_lower`, and `ci_upper` are populated only for
#'     `engine = "esem"`. That engine produces loading SEs that account for
#'     the rotation, which is the step that turns the raw solution into a more
#'     interpretable one. They are `NA` for PCA (principal component
#'     analysis) and EFA (exploratory factor analysis). The confidence level
#'     is controlled by `conf_level`.
#'   * `"variance"`: one row per factor x level, with columns
#'     `level`, `factor`, `proportion`, `cumulative`. Both are proportions of
#'     total item variance on a 0-1 scale (multiply by 100 for a percentage).
#'   * `"fit"`: one row per fit statistic x level, with columns `level`,
#'     `statistic`, `value`. For PCA objects the statistics are eigenvalues.
#'     For EFA objects they are `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, and
#'     `BIC`, where `chi` is the likelihood-ratio chi-square ([psych::fa()]'s
#'     `STATISTIC`). So `chi`, `dof`, `p_value`, `RMSEA`, and `TLI` all share
#'     one statistical framing. (psych's residual-based *empirical* chi-square
#'     is a different statistic and is not reported.) For ESEM (exploratory
#'     structural equation modeling) they are `chi`, `dof`, `p_value`, `CFI`,
#'     `TLI`, `RMSEA`, `SRMR`, `BIC`. Three estimators run a scaled test:
#'     `"WLSMV"` and `"ULSMV"` for ordinal items (a few ordered categories,
#'     such as a 1 to 5 rating), and `"MLR"` for continuous ones. Under any
#'     of them the whole row reports lavaan's
#'     mean-and-variance-adjusted ("scaled") variant, so every quantity shares
#'     one scaling. That covers `chi`, `dof`, and `p_value` **and** `CFI`,
#'     `TLI`, and `RMSEA`. This matters most for WLSMV and ULSMV. The naive
#'     chi-square has no valid reference distribution (lavaan's own
#'     `summary()` labels its p-value "Unknown"), and the naive `CFI` and
#'     `TLI` are badly optimistic for ordinal data (Xia & Yang, 2019). The
#'     estimator `"ML"` has no scaled variant, so it reports the naive values
#'     (the correct ones for ML). The statistic `SRMR` has no scaled variant
#'     and is reported as-is. The statistic `BIC` is
#'     `NA` under WLSMV and ULSMV, because a limited-information estimator has
#'     no proper log-likelihood, and it is populated under ML and MLR. Use
#'     `format = "wide"` for one row per **non-anchor** level and one column
#'     per statistic. Non-anchor means k >= 2, because the saturated 1-factor
#'     anchor is dropped, matching `summary()` and `autoplot(what = "fit")`.
#'     Conventional fit cutoffs (Hu & Bentler 1999) are shown as
#'     reference lines in `autoplot(what = "fit")` and inline in `summary()`,
#'     but are not returned as a pass/fail column here. They are contested
#'     thresholds, report-only, and never gate anything (see those functions'
#'     docs). The `format` argument is oriented to the EFA and ESEM model-fit
#'     statistics. For PCA the "statistics" are per-component eigenvalues.
#'   * `"nodes"`: Forbes-extension pruning annotations (requires `prune != "none"`
#'     when the object was created). One row per factor across all levels,
#'     with columns
#'     `id`, `level`, `pruned`, `prune_reason`. Returns an empty data frame with
#'     the same columns when no pruning was applied.
#'   * `"scores"`: long-format per-observation factor scores, where a factor
#'     score is each person's estimated standing on a factor (requires
#'     `keep_scores = TRUE` at fit time or use [augment.ackwards()] for on-the-fly
#'     computation). Columns: `obs` (row index), `level`, `factor`, `score`.
#' @param primary_only For `what = "edges"` only. When `TRUE`, returns just each
#'   factor's primary-parent edge (`is_primary == TRUE`). That is the lineage
#'   tree that the diagram draws as solid arrows. Default `FALSE` (all edges).
#'   Errors for any other value of `what`.
#' @param sort For `what = "edges"` only. One of `"none"` (default, natural order)
#'   or `"strength"` (descending `|r|`). Ignored for all other values of `what`.
#' @param format For `what = "fit"` only. One of `"long"` (default, one row per
#'   statistic x level) or `"wide"` (one row per level, one column per
#'   statistic). Errors for all other values of `what`.
#' @param conf_level For `what = "loadings"` only. Confidence level for the
#'   loading intervals. Default `0.95`. The intervals are computed as
#'   `loading ± qnorm((1 + conf_level) / 2) * se` and are `NA` for engines
#'   that carry no SEs (PCA, EFA). Errors for all other values of `what`.
#' @param ... Ignored.
#'
#' @return A data frame (class `data.frame`).
#'
#' @section Factor labels:
#' If [factor labels][set_factor_labels] have been attached to the object, the
#' output gains display-only label columns: `factor_label` for `what =
#' "loadings"`, `"variance"`, or `"scores"`, and `from_label`/`to_label` for
#' `what = "edges"`. Each carries the label for a labeled factor and `NA`
#' otherwise. These columns are **absent** when no labels are set, so an
#' unlabeled object's output is unchanged. The ID columns (`factor`, `from`,
#' `to`) are never altered.
#'
#' @seealso [glance.ackwards()], [print.ackwards()], [set_factor_labels()]
#'
#' @examples
#' x <- ackwards(sim16, k_max = 5)
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
  what = c("edges", "loadings", "variance", "fit", "nodes", "scores", "factor_cor"),
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
    scores   = .tidy_scores(x),
    factor_cor = .tidy_factor_cor(x)
  )
  # Factor labels (M51): display-only columns, present only when labels have
  # been set (unlabeled objects keep their exact pre-M51 schema; Invariant 5 --
  # the ID columns are never mutated). `nodes`/`fit` carry no factor-ID column.
  labels <- x$meta$factor_labels
  if (!is.null(labels)) {
    if (what %in% c("loadings", "variance", "scores")) {
      out$factor_label <- unname(labels[out$factor])
    } else if (what == "edges") {
      out$from_label <- unname(labels[out$from])
      out$to_label <- unname(labels[out$to])
    } else if (what == "factor_cor") {
      out$factor_a_label <- unname(labels[out$factor_a])
      out$factor_b_label <- unname(labels[out$factor_b])
    }
  }
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
  out <- x$edges$tidy
  # M47: when boot_edges() has run, expose its SE + percentile-CI columns on
  # the edge table. Joined on the directed (from, to) key -- boot rows are in
  # the same order, but match by key so the merge is robust to reordering.
  if (!is.null(x$boot)) {
    be <- x$boot$edges
    key_out <- paste(out$from, out$to, sep = "\r")
    key_be <- paste(be$from, be$to, sep = "\r")
    m <- match(key_out, key_be)
    out$se <- be$se[m]
    out$lo <- be$lo[m]
    out$hi <- be$hi[m]
    out$n_boot_ok <- be$n_boot_ok[m]
  }
  out
}

# One row per unordered factor pair within each level, read from the stored
# `factor_cor` (already permuted and sign-aligned to the stored loadings).
# Pairs are listed in upper-triangle order: (1,2), (1,3), (2,3), ... A level
# with one factor has no pair and contributes no row, so an object whose only
# level is k = 1 returns zero rows with the same four columns.
.tidy_factor_cor <- function(x) {
  empty <- data.frame(
    level = integer(0L),
    factor_a = character(0L),
    factor_b = character(0L),
    cor = numeric(0L),
    stringsAsFactors = FALSE
  )
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    k <- as.integer(ki)
    if (k < 2L) {
      return(NULL)
    }
    Phi <- lev$factor_cor
    idx <- which(upper.tri(Phi), arr.ind = TRUE)
    idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
    data.frame(
      level = k,
      factor_a = lev$labels[idx[, 1L]],
      factor_b = lev$labels[idx[, 2L]],
      cor = Phi[idx],
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(empty)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
#' Returns a one-row data frame of top-level model metadata. For EFA
#' (exploratory factor analysis) and ESEM (exploratory structural equation
#' modeling) objects, fit indices at the deepest converged level are included.
#' The same five columns (`CFI`, `TLI`, `RMSEA`, `SRMR`, `BIC`) are present
#' across all engines. Columns unavailable for a given engine or estimator are
#' `NA`. For example, `CFI` and `SRMR` are `NA` for EFA, and all five are `NA`
#' for PCA (principal component analysis). For ESEM, `BIC` is `NA` under
#' `estimator = "WLSMV"` or `"ULSMV"`, because these limited-information
#' estimators have no proper log-likelihood, and it is populated under `"ML"`
#' or `"MLR"`. Under a scaled-test estimator
#' (`"WLSMV"`, `"ULSMV"`, or `"MLR"`) the `CFI`, `TLI`, and `RMSEA` reported
#' here are the scaled variants (see [tidy.ackwards()] for the rationale).
#'
#' @param x An `ackwards` object.
#' @param ... Ignored.
#'
#' @return A one-row `data.frame`.
#'
#' @seealso [tidy.ackwards()], [print.ackwards()]
#'
#' @examples
#' x <- ackwards(sim16, k_max = 5)
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
