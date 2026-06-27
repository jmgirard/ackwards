#' Print an ackwards object
#'
#' Displays a compact summary of the bass-ackwards result using cli formatting.
#' No matrix dumps -- use [tidy.ackwards()] to access values programmatically.
#'
#' @param x An `ackwards` object.
#' @param ... Ignored.
#'
#' @return `x` invisibly.
#'
#' @seealso [tidy.ackwards()], [glance.ackwards()]
#'
#' @importFrom rlang `%||%`
#' @export
print.ackwards <- function(x, ...) {
  cut_show <- x$meta$cut_show %||% 0.3
  cli::cli_h1("Bass-Ackwards Analysis ({.pkg ackwards})")

  # --- Call / settings --------------------------------------------------------
  cli::cli_dl(c(
    "Engine"    = cli::style_bold(x$method),
    "Rotation"  = x$rotation,
    "Basis"     = x$cor_type,
    "n"         = format(x$n_obs, big.mark = ","),
    "k (max)"   = as.character(x$k_max)
  ))

  # --- Per-level table --------------------------------------------------------
  cli::cli_h2("Levels")

  K <- x$k_max
  for (ki in seq_len(K)) {
    lev <- x$levels[[as.character(ki)]]
    cum_pct <- round(lev$variance["cumulative"] * 100, 1)
    conv_sym <- if (isTRUE(lev$converged)) cli::col_green(cli::symbol$tick) else cli::col_red(cli::symbol$cross)
    cli::cli_text(
      "  {conv_sym} {.strong k = {ki}}: {ki} factor{?s}, {cum_pct}% variance"
    )
  }

  # --- Edge summary -----------------------------------------------------------
  tidy <- x$edges$tidy
  n_edges_total <- nrow(tidy)
  n_above <- sum(tidy$above_cut, na.rm = TRUE)

  cli::cli_h2("Edges")
  cli::cli_text(
    "{n_above} of {n_edges_total} edges have |r| {cli::symbol$geq} {cut_show}"
  )

  # --- Caveat (DESIGN.md section 2) ------------------------------------------
  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Note: This is a series of linked solutions, not a fitted hierarchical \\
       model. Cross-level edges are descriptive score correlations."
    )
  )

  invisible(x)
}
