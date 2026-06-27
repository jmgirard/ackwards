#' Suggest a maximum number of factors for bass-ackwards analysis
#'
#' Runs two complementary criteria and reports their recommendations. Neither
#' alone is definitive; the goal is a consensus range to inform your choice of
#' `k` in [ackwards()].
#'
#' **Criteria computed:**
#' * **Parallel analysis** (Horn 1965) -- compares observed eigenvalues to those
#'   from random correlation matrices of the same size. Suggests the number of
#'   components whose eigenvalues exceed the 95th-percentile chance level.
#' * **MAP** (Velicer 1976) -- the Minimum Average Partial criterion. Finds the
#'   number of components that minimises the average squared partial correlation
#'   after the components are partialled out.
#'
#' Both use the same Pearson or Spearman correlation matrix as [ackwards()].
#' Parallel analysis is implemented via [psych::fa.parallel()]; MAP via
#' [psych::vss()].
#'
#' @section Interpreting the output:
#' `k` in [ackwards()] is a **maximum depth**, not a claim that exactly k
#' factors exist. Users commonly set k one or two levels above the consensus to
#' watch higher-level factors fragment -- this is a feature of the method, not
#' overextraction.
#'
#' @param data A data frame or numeric matrix (items in columns, observations in
#'   rows).
#' @param k_max Maximum number of components to test. Defaults to
#'   `min(ncol(data) - 1, 8)`. Increase if you expect a deeper hierarchy.
#' @param cor Correlation basis: `"pearson"` (default) or `"spearman"`. Should
#'   match the `cor` argument you plan to use in [ackwards()].
#' @param n_iter Number of Monte Carlo iterations for parallel analysis. Default
#'   `20`. Reduce to `5` for fast/exploratory runs; increase to `100+` for
#'   publication.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"suggest_k"`. Print it for a formatted summary.
#'   The list contains:
#'   \item{k_parallel}{Recommended k from parallel analysis.}
#'   \item{k_map}{Recommended k from MAP.}
#'   \item{criteria}{Data frame with one row per k: `k`, `map` value,
#'     `pa_suggested` (logical; `TRUE` if k is within the parallel-analysis threshold).}
#'   \item{k_max, n_obs, n_vars, cor}{Metadata.}
#'
#' @section A note on overextraction:
#' Parallel analysis in particular tends to recommend more factors than replicate
#' across independent samples, especially with correlated items (Forbes, 2023).
#' Treat these criteria as a starting range for exploration, not a definitive
#' stopping rule. Setting `k` in [ackwards()] one or two levels above the
#' consensus is intentional -- watching factors fragment is part of the method.
#'
#' @seealso [ackwards()]
#'
#' @references
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000578}
#'
#' Horn, J. L. (1965). A rationale and test for the number of factors in factor
#'   analysis. *Psychometrika*, 30, 179--185.
#'
#' Velicer, W. F. (1976). Determining the number of components from the matrix
#'   of partial correlations. *Psychometrika*, 41, 321--327.
#'
#' @examples
#' \dontrun{
#' suggest_k(psych::bfi[, 1:25])
#' suggest_k(psych::bfi[, 1:25], k_max = 6, n_iter = 5)
#' }
#'
#' @export
suggest_k <- function(data, k_max = NULL, cor = "pearson", n_iter = 20L, ...) {
  rlang::check_installed("psych", reason = "for suggest_k()")

  cor <- rlang::arg_match(cor, c("pearson", "spearman"))

  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
  }
  data_mat <- as.matrix(data)
  if (!is.numeric(data_mat)) {
    cli::cli_abort("{.arg data} must contain only numeric columns.")
  }

  p <- ncol(data_mat)
  n <- nrow(data_mat)

  if (is.null(k_max)) k_max <- min(p - 1L, 8L)
  k_max    <- as.integer(k_max)
  n_iter   <- as.integer(n_iter)

  if (k_max < 1L || k_max >= p) {
    cli::cli_abort(
      "{.arg k_max} must be between 1 and {p - 1L} (number of variables - 1)."
    )
  }

  R <- stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")

  # --- Parallel analysis (Horn) -----------------------------------------------
  cli::cli_progress_step("Running parallel analysis ({n_iter} iterations)...")
  pa <- psych::fa.parallel(
    R,
    n.obs  = n,
    fa     = "pc",
    n.iter = n_iter,
    plot   = FALSE,
    quant  = 0.95
  )
  # Cap at k_max in case pa reports more components than we tested
  k_parallel <- min(pa$ncomp, k_max)

  # --- MAP (Velicer) ----------------------------------------------------------
  cli::cli_progress_step("Running MAP (Velicer)...")
  vss_out <- psych::vss(
    R,
    n      = k_max,
    n.obs  = n,
    rotate = "varimax",
    fm     = "pc",
    plot   = FALSE
  )
  map_vals <- vss_out$map[seq_len(k_max)]
  k_map    <- which.min(map_vals)

  cli::cli_progress_done()

  # --- Build criteria table ---------------------------------------------------
  # pa_suggested: TRUE if k is within the parallel-analysis threshold
  criteria <- data.frame(
    k            = seq_len(k_max),
    map          = map_vals,
    pa_suggested = seq_len(k_max) <= k_parallel,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      k_parallel = k_parallel,
      k_map      = k_map,
      criteria   = criteria,
      k_max      = k_max,
      n_obs      = n,
      n_vars     = p,
      cor        = cor
    ),
    class = "suggest_k"
  )
}

#' Print a suggest_k object
#'
#' @param x A `suggest_k` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.suggest_k <- function(x, ...) {
  cli::cli_h1("Factor / Component Count Suggestion ({.pkg ackwards})")

  cli::cli_dl(c(
    "Variables" = as.character(x$n_vars),
    "n"         = format(x$n_obs, big.mark = ","),
    "Basis"     = x$cor,
    "Tested k"  = paste0("1\u2013", x$k_max)
  ))

  cli::cli_h2("Criteria (k = 1\u2013{x$k_max})")

  cr      <- x$criteria
  pa_sym  <- ifelse(cr$pa_suggested, cli::col_green(cli::symbol$tick), cli::col_grey("-"))
  map_fmt <- formatC(cr$map, digits = 4, format = "f")

  for (i in seq_len(nrow(cr))) {
    cli::cli_text(
      "  {pa_sym[i]} k = {cr$k[i]}:  MAP = {map_fmt[i]}  |  \\
       PA {ifelse(cr$pa_suggested[i], 'suggested', 'not suggested')}"
    )
  }

  cli::cli_h2("Recommendations")

  cli::cli_bullets(c(
    "*" = "Parallel analysis: k {cli::symbol$leq} {x$k_parallel}",
    "*" = "MAP (Velicer):     k = {x$k_map}"
  ))

  lo <- min(x$k_map, x$k_parallel)
  hi <- max(x$k_map, x$k_parallel)
  if (lo == hi) {
    cli::cli_text("{.strong Consensus: k = {lo}}")
  } else {
    cli::cli_text("{.strong Consensus range: k = {lo}\u2013{hi}}")
  }

  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Note: k in ackwards() is a maximum depth. Consider setting k one \\
       or two levels above the consensus to observe factor fragmentation."
    )
  )
  cli::cli_text(
    cli::col_grey(
      "Caution: parallel analysis tends to overextract; many suggested \\
       structures do not replicate (Forbes, 2023). Treat this as a range."
    )
  )

  invisible(x)
}
