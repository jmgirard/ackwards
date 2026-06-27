#' Augment data with factor scores from an ackwards object
#'
#' Appends per-observation factor scores for every level to a data frame.
#' Score columns are named `.m{k}f{j}` (e.g., `.m1f1`, `.m2f1`, `.m2f2`),
#' matching the factor labels used throughout the object.
#'
#' @details
#' **Score computation.** Scores are `S = Z W / sqrt(score_var)`, where
#' `Z = scale(data)` (item z-scores), `W` is the per-level weight matrix
#' stored in the object, and `sqrt(score_var)` standardizes by the real
#' score standard deviations (Invariant 1: never assume unit variance).
#' For PCA the method is `"components"`; for EFA/ESEM it is `"tenBerge"`.
#'
#' **Data source.** If `data` is supplied, scores are always recomputed from
#' it using the stored weights — this is how to score new observations. If
#' `data` is `NULL` and `scores = TRUE` was set at fit time, the stored
#' scores are returned. If neither is available an informative error is raised.
#'
#' @param x An `ackwards` object.
#' @param data A data frame or numeric matrix with the same variables (columns)
#'   used to fit `x`. When `NULL` (default), uses pre-stored scores if
#'   available (requires `scores = TRUE` at fit time).
#' @param ... Ignored.
#'
#' @return A data frame. If `data` is supplied it is returned with score
#'   columns appended. If `data` is `NULL` the return is a minimal data frame
#'   with a `.obs` index column followed by score columns.
#'
#' @seealso [ackwards()], [tidy.ackwards()]
#'
#' @examples
#' \dontrun{
#' # Score the training data on the fly (no scores=TRUE needed)
#' x <- ackwards(psych::bfi[, 1:25], k = 5)
#' scores_df <- augment(x, data = psych::bfi[, 1:25])
#' head(scores_df[, grep("^\\.m", names(scores_df))])
#'
#' # Or store at fit time and augment without re-supplying data
#' x2 <- ackwards(psych::bfi[, 1:25], k = 5, scores = TRUE)
#' scores_df2 <- augment(x2)
#' }
#'
#' @importFrom generics augment
#' @export
augment.ackwards <- function(x, data = NULL, ...) {
  scores_list <- if (!is.null(data)) {
    .compute_scores(x$levels, as.matrix(data))
  } else if (!is.null(x$scores)) {
    x$scores
  } else {
    cli::cli_abort(c(
      "!" = "Factor scores are not stored in this {.cls ackwards} object \\
             and no {.arg data} was supplied.",
      "i" = "Either refit with {.code scores = TRUE} or pass the original \\
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
