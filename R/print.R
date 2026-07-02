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
  cor_label <- if (is.na(x$cor)) "(user-supplied matrix)" else x$cor
  n_label <- if (is.na(x$n_obs)) "NA" else format(x$n_obs, big.mark = ",")
  cli::cli_dl(c(
    "Engine"    = cli::style_bold(x$engine),
    "Rotation"  = x$rotation,
    "Basis"     = cor_label,
    "n"         = n_label,
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
  # --- Bootstrap edge CIs (M47) -----------------------------------------------
  if (!is.null(x$boot)) {
    pct <- round(100 * x$boot$conf)
    cli::cli_text(
      cli::col_grey(
        "  {pct}% bootstrap CIs on all edges ({x$boot$n_boot} replicates); \\
         see {.code tidy(x)} or {.code x$boot$edges}."
      )
    )
  }

  # --- Pruning summary (Forbes extension; DESIGN.md s14.18) -------------------
  if (!is.null(x$prune)) {
    cli::cli_h2("Pruning")
    rules <- x$prune$rules
    if ("redundant" %in% rules) {
      n_flagged <- sum(x$prune$nodes$pruned & x$prune$nodes$prune_reason == "redundant")
      r_thr <- x$prune$redundancy_r
      phi_note <- if (!is.null(x$prune$redundancy_phi)) {
        paste0(", phi > ", x$prune$redundancy_phi)
      } else {
        ""
      }
      cli::cli_text(
        "  Redundancy (|r| {cli::symbol$geq} {r_thr}{phi_note}): \\
         {n_flagged} node{?s} flagged"
      )
    }
    if ("artifact" %in% rules) {
      n_phi <- if (!is.null(x$prune$phi)) nrow(x$prune$phi) else 0L
      cli::cli_text(
        "  Artifact: Tucker's phi computed for {n_phi} cross-level factor pair{?s}"
      )
      if (!is.null(x$prune$structural)) {
        n_struct <- sum(
          x$prune$structural$few_items | x$prune$structural$orphan |
            x$prune$structural$split_merge,
          na.rm = TRUE
        )
        cli::cli_text(
          "  Structural signals: {n_struct} factor{?s} flagged \\
           (inspect {.code x$prune$structural})"
        )
      }
    }
    if (!is.null(x$prune$manual) && length(x$prune$manual) > 0L) {
      cli::cli_text(
        "  Manual: {length(x$prune$manual)} node{?s} explicitly flagged \\
         ({paste(x$prune$manual, collapse = ', ')})"
      )
    }
    cli::cli_rule()
    cli::cli_text(
      cli::col_grey(
        "Note: Pruning is interpretive relabeling, not re-estimation. \\
         Flagged nodes remain in the object; all edges are preserved. \\
         Inspect with {.code x$prune$nodes} and {.code tidy(x, what = \"nodes\")}."
      )
    )
  }

  # --- Caveat (DESIGN.md section 2) ------------------------------------------
  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Note: This is a series of linked solutions, not a fitted hierarchical \\
       model. Cross-level edges are descriptive score correlations. \\
       Per-level fit indices (EFA/ESEM) describe how well a k-factor model \\
       fits the items at that level -- they do not validate the edges or \\
       the hierarchy itself."
    )
  )

  invisible(x)
}
