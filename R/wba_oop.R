
# S3 constructors ---------------------------------------------------------

new_bar <- function(correlations = list(),
                    loadings = list(),
                    details = list()) {

  stopifnot(is.list(correlations))
  stopifnot(is.list(loadings))
  stopifnot(is.list(details))

  structure(
    list(
      correlations = correlations,
      loadings = loadings,
      details = details
    ),
    class = "ackwards_bar"
  )

}


# S3 methods --------------------------------------------------------------

#' Print Bass-Ackwards Results
#'
#' @param x An object of class "ackwards_bar" from `wba()` or `xba()`.
#' @param digits A
#' @param cut A double indicating which magnitude of loadings to display or
#'   `NULL` to display all loadings (default = `.3`).
#' @param ... Ignored for now.
#'
#' @return Prints and returns x
#' @export print.ackwards_bar
#' @export
#'
#' @examples
#' wba(cor(skiers), 3)
print.ackwards_bar <- function(x, digits = 3, cut = .3, ...) {

  # Assertions
  stopifnot(is.numeric(digits))
  stopifnot(length(digits) == 1)
  stopifnot(floor(digits) == ceiling(digits))
  stopifnot(is.null(cut) || is.numeric(cut))
  stopifnot(is.null(cut) || is.finite(cut))
  stopifnot(is.null(cut) || length(cut) == 1)
  stopifnot(is.null(cut) || cut >= 0)

  # Print header
  cat(
    "Bass-Ackwards Analysis\n",
    paste0("  Method = ", x$details$method),
    paste0("  Engine = ", x$details$engine),
    sep = "\n"
  )

  # Print cross-level correlations
  cat("\nCross-Level Correlations\n\n")

  for (i in 2:length(x$correlations)) {
    print.default(
      round(x$correlations[[i]], digits),
      print.gap = 3L
    )
    cat("\n")
  }

  # Print factor loadings
  cat("\nFactor Loadings\n\n")

  for (i in 1:length(x$loadings)) {
    load_i <- round(x$loadings[[i]], digits)
    if (!is.null(cut)) {
      load_i[abs(load_i) < cut] <- NA
    }
    print.default(
      load_i,
      print.gap = 3L,
      na.print = "."
    )
    cat("\n")
  }
  if (!is.null(cut)) {
    cat(paste0(". = Loading magnitude less than ", cut, "\n\n"))
  }

  invisible(x)
}
