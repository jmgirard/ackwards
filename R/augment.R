#' @importFrom generics augment
#' @export
generics::augment

#' Augment data with factor scores from an ackwards object
#'
#' Appends per-observation factor scores for every level to a data frame.
#' Score columns are named `.m{k}f{j}` (e.g., `.m1f1`, `.m2f1`, `.m2f2`),
#' matching the factor labels used throughout the object.
#'
#' @details
#' **Score computation.** Scores are `S = Z W / sqrt(score_var)`, where
#' `Z = .standardize(data)` (item z-scores), `W` is the per-level weight matrix
#' stored in the object, and `sqrt(score_var)` standardizes by the real
#' score standard deviations (Invariant 1: never assume unit variance).
#' For PCA the method is `"components"`; for EFA/ESEM it is `"tenBerge"`.
#'
#' **Missing data.** Score projection applies weights row-wise and propagates
#' NAs listwise: any observation with at least one missing item variable will
#' produce `NA` scores at every level. This differs from fitting, which uses
#' pairwise-complete correlations. A warning is issued if NA rows are detected.
#' Use `na.omit(data)` before scoring if NA rows are unwanted.
#'
#' **Data source.** If `data` is supplied, scores are always recomputed from
#' it using the stored weights -- this is how to score new observations. If
#' `data` is `NULL` and `keep_scores = TRUE` was set at fit time, the stored
#' scores are returned. If neither is available an informative error is raised.
#'
#' **Scores-only output (`append = FALSE`).** By default (`append = TRUE`) the
#' scores are appended to the supplied `data`, following the broom convention.
#' Set `append = FALSE` to return *only* the score columns -- convenient for
#' feeding scores straight into [cor()], [lm()], or a clustering call without
#' dragging the item columns along. Because `augment()` always preserves row
#' order and row count, `cbind(data, augment(x, data, append = FALSE))`
#' reproduces the appended output exactly. A row that scores `NA` (because it is
#' missing an item -- see *Missing data* above) is still returned in place, so the
#' positional alignment holds even with `NA` scores. For a rejoin that survives
#' *filtering* the scores afterwards, name identifier columns with `id_cols` so
#' they travel with the scores.
#'
#' @param x An `ackwards` object.
#' @param data A data frame or numeric matrix with the same variables (columns)
#'   used to fit `x`. When `NULL` (default), uses pre-stored scores if
#'   available (requires `keep_scores = TRUE` at fit time).
#' @param append Logical. When `TRUE` (default) the score columns are appended
#'   to `data` (or to a `.obs` index when `data` is `NULL`). When `FALSE` only
#'   the score columns are returned (plus any `id_cols`).
#' @param id_cols Optional character vector naming columns of `data` to carry
#'   through alongside the scores when `append = FALSE` (for example a subject
#'   identifier, so scores can be rejoined after filtering). Ignored -- and an
#'   error -- when `append = TRUE` (all columns are already kept) or when `data`
#'   is `NULL` (there are no source columns to carry). `NULL` (default) returns
#'   the bare score columns.
#' @param ... Ignored.
#'
#' @return A data frame. With `append = TRUE`: the supplied `data` (or a `.obs`
#'   index when `data` is `NULL`) with score columns appended. With
#'   `append = FALSE`: only the score columns, optionally prefixed by the
#'   `id_cols`. Row order and count always match the input.
#'
#' @seealso [ackwards()], [tidy.ackwards()]
#'
#' @examples
#' # Score the training data on the fly (no keep_scores = TRUE needed)
#' x <- ackwards(bfi25, k_max = 5)
#' scores_df <- augment(x, data = bfi25)
#' head(scores_df[, startsWith(names(scores_df), ".m")])
#'
#' # Scores-only: just the .m{k}f{j} columns, ready for cor()/lm()
#' scores_only <- augment(x, data = bfi25, append = FALSE)
#' round(cor(scores_only[, c(".m5f1", ".m5f2")]), 2)
#'
#' # Carry an identifier through for a safe post-filter rejoin
#' df <- data.frame(id = seq_len(nrow(bfi25)), bfi25)
#' scored <- augment(x, data = df, append = FALSE, id_cols = "id")
#' head(scored[, c("id", ".m5f1")])
#'
#' # Store at fit time and augment without re-supplying data
#' x2 <- ackwards(bfi25, k_max = 5, keep_scores = TRUE)
#' scores_df2 <- augment(x2)
#'
#' @export
augment.ackwards <- function(x, data = NULL, append = TRUE, id_cols = NULL, ...) {
  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    cli::cli_abort("{.arg append} must be a single {.code TRUE} or {.code FALSE}.")
  }
  if (!is.null(id_cols)) {
    if (!is.character(id_cols)) {
      cli::cli_abort("{.arg id_cols} must be a character vector of column names or {.code NULL}.")
    }
    if (isTRUE(append)) {
      cli::cli_abort(c(
        "!" = "{.arg id_cols} is only used with {.code append = FALSE}.",
        "i" = "With {.code append = TRUE} every column of {.arg data} is already kept."
      ))
    }
    if (is.null(data)) {
      cli::cli_abort(c(
        "!" = "{.arg id_cols} needs {.arg data} to carry columns from.",
        "i" = "Supply {.arg data}, or drop {.arg id_cols} to return bare score columns."
      ))
    }
    missing_id <- setdiff(id_cols, colnames(as.data.frame(data)))
    if (length(missing_id) > 0L) {
      cli::cli_abort(c(
        "!" = "{.arg id_cols} column{?s} not found in {.arg data}: {.val {missing_id}}."
      ))
    }
  }
  scores_list <- if (!is.null(data)) {
    data_mat <- as.matrix(data)
    if (!is.numeric(data_mat)) {
      non_num <- names(which(!vapply(
        as.data.frame(data), is.numeric, logical(1L)
      )))
      cli::cli_abort(c(
        "!" = "{.arg data} must contain only numeric columns.",
        "x" = "Non-numeric column{?s}: {.val {non_num}}"
      ))
    }
    W_ref <- x$levels[[1L]]$scoring$weights
    p_expected <- nrow(W_ref)
    vars_expected <- rownames(W_ref)
    if (!is.null(vars_expected) && !is.null(colnames(data_mat))) {
      missing_vars <- setdiff(vars_expected, colnames(data_mat))
      if (length(missing_vars) > 0L) {
        cli::cli_abort(c(
          "!" = "{.arg data} is missing {length(missing_vars)} variable{?s} \\
                 that the model was fit on.",
          "x" = "Missing: {.val {missing_vars}}"
        ))
      }
      data_mat <- data_mat[, vars_expected, drop = FALSE]
    } else if (ncol(data_mat) != p_expected) {
      cli::cli_abort(c(
        "!" = "{.arg data} has {ncol(data_mat)} column{?s} but the model \\
               was fit on {p_expected}.",
        "i" = "Supply a data frame with the same {p_expected} variables \\
               used at fit time (or with matching column names)."
      ))
    }
    .compute_scores(x$levels, data_mat)
  } else if (!is.null(x$scores)) {
    x$scores
  } else {
    cli::cli_abort(c(
      "!" = "Factor scores are not stored in this {.cls ackwards} object \\
             and no {.arg data} was supplied.",
      "i" = "Either refit with {.code keep_scores = TRUE} or pass the original \\
             data: {.code augment(x, data = your_data)}."
    ))
  }

  n <- nrow(scores_list[[1L]])

  # Base carrier: the passthrough columns the scores get attached to.
  # append = TRUE  -> full data (or a .obs index when data is NULL).
  # append = FALSE -> only the requested id_cols (empty frame when none), so the
  #                   result is scores-only. Row order/count always match input.
  base <- if (isTRUE(append)) {
    if (!is.null(data)) as.data.frame(data) else data.frame(.obs = seq_len(n))
  } else if (!is.null(id_cols)) {
    as.data.frame(data)[, id_cols, drop = FALSE]
  } else {
    data.frame(row.names = seq_len(n))[, integer(0L), drop = FALSE]
  }

  out <- base
  for (ki in names(scores_list)) {
    S <- scores_list[[ki]]
    for (j in seq_len(ncol(S))) {
      out[[paste0(".", colnames(S)[j])]] <- S[, j]
    }
  }
  rownames(out) <- NULL

  out
}
