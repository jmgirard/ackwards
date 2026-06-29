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
#'   * `"edges"` *(default)* — one row per directed between-level edge:
#'     `from`, `to`, `level_from`, `level_to`, `r`, `is_primary`, `above_cut`.
#'   * `"loadings"` — one row per item × factor × level:
#'     `level`, `factor`, `item`, `loading`.
#'   * `"loadings_se"` — one row per item × factor × level with the rotation-aware
#'     standard error of each loading: `level`, `factor`, `item`, `se`. Only the
#'     ESEM engine (`engine = "esem"`) produces these; errors informatively for
#'     PCA/EFA objects, which carry no loading standard errors.
#'   * `"variance"` — one row per factor × level:
#'     `level`, `factor`, `variance_pct`, `cumulative_pct`.
#'   * `"fit"` — one row per fit index × level: `level`, `index`, `value`.
#'     For PCA objects the indices are eigenvalues; for EFA objects they are
#'     `chi`, `dof`, `p_value`, `RMSEA`, `TLI`, `BIC`.
#'   * `"nodes"` — Forbes-extension pruning annotations (requires `prune != "none"`
#'     when the object was created). One row per factor across all levels:
#'     `id`, `level`, `pruned`, `prune_reason`. Returns an empty data frame with
#'     the same columns when no pruning was applied.
#'   * `"scores"` — long-format per-observation factor scores (requires
#'     `keep_scores = TRUE` at fit time or use [augment.ackwards()] for on-the-fly
#'     computation). Columns: `obs` (row index), `level`, `factor`, `score`.
#' @param primary_only For `what = "edges"` only. When `TRUE`, returns just each
#'   factor's primary-parent edge (`is_primary == TRUE`) — the lineage tree that
#'   the diagram draws as solid arrows. Default `FALSE` (all edges). Errors for
#'   any other value of `what`.
#' @param sort For `what = "edges"` only. One of `"none"` (default, natural order)
#'   or `"strength"` (descending `|r|`). Ignored for all other values of `what`.
#' @param ... Ignored.
#'
#' @return A data frame (class `data.frame`).
#'
#' @seealso [glance.ackwards()], [print.ackwards()]
#'
#' @examples
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   x <- ackwards(psych::bfi[, 1:25], k_max = 5)
#'   tidy(x) # edges in natural order
#'   tidy(x, sort = "strength") # strongest edges first
#'   tidy(x, primary_only = TRUE) # just the primary-parent lineage
#'   tidy(x, what = "loadings")
#'   tidy(x, what = "variance")
#' }
#'
#' @export
tidy.ackwards <- function(
  x,
  what = c("edges", "loadings", "loadings_se", "variance", "fit", "nodes", "scores"),
  primary_only = FALSE,
  sort = c("none", "strength"),
  ...
) {
  what <- match.arg(what)
  sort <- match.arg(sort)
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
  out <- switch(what,
    edges       = .tidy_edges(x),
    loadings    = .tidy_loadings(x),
    loadings_se = .tidy_loadings_se(x),
    variance    = .tidy_variance(x),
    fit         = .tidy_fit(x),
    nodes       = .tidy_nodes(x),
    scores      = .tidy_scores(x)
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

.tidy_loadings <- function(x) {
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    L <- lev$loadings
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(L)), function(j) {
      data.frame(
        level = k,
        factor = colnames(L)[j],
        item = rownames(L),
        loading = L[, j],
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.tidy_loadings_se <- function(x) {
  has_se <- vapply(x$levels, function(lev) !is.null(lev$loadings_se), logical(1L))
  if (!any(has_se)) {
    cli::cli_abort(c(
      "!" = "Loading standard errors are not available in this {.cls ackwards} \\
             object.",
      "i" = "Rotation-aware loading SEs are produced only by \\
             {.code engine = \"esem\"}; PCA and EFA carry none."
    ))
  }
  rows <- lapply(names(x$levels), function(ki) {
    SE <- x$levels[[ki]]$loadings_se
    if (is.null(SE)) {
      return(NULL)
    }
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(SE)), function(j) {
      data.frame(
        level = k,
        factor = colnames(SE)[j],
        item = rownames(SE),
        se = SE[, j],
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
#' Returns a one-row data frame of top-level model metadata.
#'
#' @param x An `ackwards` object.
#' @param ... Ignored.
#'
#' @return A one-row `data.frame`.
#'
#' @seealso [tidy.ackwards()], [print.ackwards()]
#'
#' @export
glance.ackwards <- function(x, ...) {
  data.frame(
    engine            = x$engine,
    rotation          = x$rotation,
    cor               = x$cor,
    k_max             = x$k_max,
    n_obs             = x$n_obs,
    deepest_converged = x$meta$deepest_converged,
    n_edges           = nrow(x$edges$tidy),
    stringsAsFactors  = FALSE
  )
}
