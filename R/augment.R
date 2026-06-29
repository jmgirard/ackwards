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
#' @param x An `ackwards` object.
#' @param data A data frame or numeric matrix with the same variables (columns)
#'   used to fit `x`. When `NULL` (default), uses pre-stored scores if
#'   available (requires `keep_scores = TRUE` at fit time).
#' @param ... Ignored.
#'
#' @return A data frame. If `data` is supplied it is returned with score
#'   columns appended. If `data` is `NULL` the return is a minimal data frame
#'   with a `.obs` index column followed by score columns.
#'
#' @seealso [ackwards()], [tidy.ackwards()]
#'
#' @examples
#' # Score the training data on the fly (no keep_scores = TRUE needed)
#' x <- ackwards(bfi25, k_max = 5)
#' scores_df <- augment(x, data = bfi25)
#' head(scores_df[, startsWith(names(scores_df), ".m")])
#'
#' # Store at fit time and augment without re-supplying data
#' x2 <- ackwards(bfi25, k_max = 5, keep_scores = TRUE)
#' scores_df2 <- augment(x2)
#'
#' @export
augment.ackwards <- function(x, data = NULL, ...) {
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
  out <- if (!is.null(data)) as.data.frame(data) else data.frame(.obs = seq_len(n))

  for (ki in names(scores_list)) {
    S <- scores_list[[ki]]
    for (j in seq_len(ncol(S))) {
      out[[paste0(".", colnames(S)[j])]] <- S[, j]
    }
  }

  out
}
