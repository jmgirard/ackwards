#' Score new observations with a fitted ackwards model
#'
#' Applies a fitted bass-ackwards model to new data — for example the held-out
#' test split of a cross-validation design — producing factor scores for every
#' level **without retraining**. This is the idiomatic `predict()` front door
#' to the same machinery as [augment.ackwards()]: the call
#' `predict(object, newdata)` returns exactly
#' `augment(object, data = newdata, append = FALSE)`, a data frame holding
#' only the `.m{k}f{j}` score columns, one row per row of `newdata`.
#'
#' @details
#' Under the default `scaling = "fit"`, `newdata` is standardized by the
#' **fit-time** item means/SDs stored in the object before the stored weight
#' matrices are applied, so the new scores land on the same metric as the
#' training solution: an observation's score does not depend on which other
#' observations share its split, and train and test scores are directly
#' comparable. See [augment.ackwards()] (section *Scoring new observations*)
#' for the full semantics, the `scaling = "sample"` alternative, and the
#' non-Pearson-basis caveat.
#'
#' `newdata` must contain the variables the model was fit on (matched by
#' column name, with extra columns ignored; a bare unnamed matrix is matched
#' positionally). Rows with missing items produce `NA` scores (scoring does
#' not impute).
#'
#' @param object An `ackwards` object.
#' @param newdata A data frame or numeric matrix with the same variables
#'   (columns) used to fit `object`. Required — to retrieve scores stored at
#'   fit time, use `augment(object)` instead.
#' @param scaling Which item means/SDs standardize `newdata`: `"fit"`
#'   (default, the training moments) or `"sample"` (`newdata`'s own moments).
#'   See [augment.ackwards()].
#' @param ... Ignored.
#'
#' @return A data frame of factor scores (columns `.m{k}f{j}`, one per factor
#'   per level), with one row per row of `newdata` in the original order.
#'
#' @seealso [augment.ackwards()] for appending scores to the data (and the
#'   full scoring documentation), [ackwards()].
#'
#' @examples
#' # Cross-validation: fit on a training split, score the test split
#' train <- bfi25[1:500, ]
#' test <- bfi25[501:1000, ]
#' x <- ackwards(train, k_max = 5)
#' test_scores <- predict(x, test)
#' head(test_scores)
#'
#' # Identical to the augment() spelling:
#' identical(test_scores, augment(x, data = test, append = FALSE))
#'
#' @export
predict.ackwards <- function(object, newdata, scaling = c("fit", "sample"), ...) {
  scaling <- rlang::arg_match(scaling)
  if (missing(newdata) || is.null(newdata)) {
    cli::cli_abort(c(
      "!" = "{.arg newdata} is required.",
      "i" = "{.fn predict} scores new observations; to retrieve scores stored \\
             at fit time, use {.code augment(object)}."
    ))
  }
  augment.ackwards(object, data = newdata, append = FALSE, scaling = scaling)
}
