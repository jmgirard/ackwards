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
#'   * `"variance"` — one row per factor × level:
#'     `level`, `factor`, `variance_pct`, `cumulative_pct`.
#' @param ... Ignored.
#'
#' @return A data frame (class `data.frame`).
#'
#' @seealso [glance.ackwards()], [print.ackwards()]
#'
#' @importFrom generics tidy
#' @export
tidy.ackwards <- function(x, what = c("edges", "loadings", "variance"), ...) {
  what <- match.arg(what)
  switch(what,
    edges    = .tidy_edges(x),
    loadings = .tidy_loadings(x),
    variance = .tidy_variance(x)
  )
}

.tidy_edges <- function(x) {
  x$edges$tidy
}

.tidy_loadings <- function(x) {
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    L   <- lev$loadings
    k   <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(L)), function(j) {
      data.frame(
        level   = k,
        factor  = colnames(L)[j],
        item    = rownames(L),
        loading = L[, j],
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.tidy_variance <- function(x) {
  rows <- lapply(names(x$levels), function(ki) {
    lev <- x$levels[[ki]]
    k   <- as.integer(ki)
    fac_labels <- lev$labels
    var_vals   <- lev$variance[fac_labels]
    cum_val    <- lev$variance["cumulative"]
    data.frame(
      level          = k,
      factor         = fac_labels,
      variance_pct   = round(var_vals * 100, 2),
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
#' @importFrom generics glance
#' @export
glance.ackwards <- function(x, ...) {
  data.frame(
    method            = x$method,
    rotation          = x$rotation,
    cor_type          = x$cor_type,
    k_max             = x$k_max,
    n_obs             = x$n_obs,
    deepest_converged = x$meta$deepest_converged,
    n_edges           = nrow(x$edges$tidy),
    stringsAsFactors  = FALSE
  )
}
