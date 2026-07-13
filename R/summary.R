#' Summarise an ackwards object
#'
#' Returns a structured `summary_ackwards` object that, when printed, shows
#' per-level variance and fit indices, a readable lineage list, and (when
#' present) pruning annotations. More verbose than [print.ackwards()]; designed
#' for inspection and reporting.
#'
#' @param object An `ackwards` object.
#' @param ... Ignored.
#'
#' @return An object of class `"summary_ackwards"`, printed via
#'   [print.summary_ackwards()].
#'
#' @seealso [print.ackwards()], [tidy.ackwards()], [glance.ackwards()]
#'
#' @examples
#' x <- ackwards(sim16, k_max = 5)
#' summary(x)
#'
#' @export
summary.ackwards <- function(object, ...) {
  structure(
    list(
      call = object$call,
      engine = object$engine,
      rotation = object$rotation,
      cor = object$cor,
      n_obs = object$n_obs,
      k_max = object$k_max,
      estimator = object$meta$estimator,
      near_singular = isTRUE(object$meta$near_singular),
      min_eigenvalue = object$meta$min_eigenvalue,
      variance = .tidy_variance(object),
      fit = .tidy_fit(object),
      lineage = .summary_lineage(object),
      prune = .prune_digest(object),
      boot = object$boot,
      factor_labels = object$meta$factor_labels # M51
    ),
    class = "summary_ackwards"
  )
}

#' Print a summary_ackwards object
#'
#' @param x A `summary_ackwards` object (produced by `summary()`).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.summary_ackwards <- function(x, ...) {
  cli::cli_h1("Summary: Bass-Ackwards Analysis ({.pkg ackwards})")

  .print_ba_header(x$engine, x$rotation, x$cor, x$n_obs, x$k_max)

  # --- Per-level variance + fit -----------------------------------------------
  cli::cli_h2("Levels")

  is_esem <- x$engine == "esem"
  is_efa <- x$engine == "efa"
  is_pca <- x$engine == "pca"

  fit_tbl <- x$fit # pre-pulled for readability

  for (ki in seq_len(x$k_max)) {
    if (ki > 1L) cli::cli_text("") # blank line between level blocks
    var_rows <- x$variance[x$variance$level == ki, , drop = FALSE]
    # cumulative is a running proportion within the level; max gives the total.
    # Fixed-precision with trailing zeros so per-level figures do not drift
    # (e.g. "20.9%" not "20.91%"/"13.6%").
    cum_total <- sprintf("%.1f", max(var_rows$cumulative) * 100)
    cli::cli_text(
      "{.strong k = {ki}}: {ki} factor{?s}  ({cum_total}% cumulative variance)"
    )
    fit_rows <- fit_tbl[fit_tbl$level == ki, , drop = FALSE]

    for (i in seq_len(nrow(var_rows))) {
      fac <- var_rows$factor[i]
      vpct <- sprintf("%.1f", var_rows$proportion[i] * 100)

      # PCA: eigenvalue statistics are "eigenvalue.<label>" (see engine_pca.R:65)
      suffix <- if (is_pca && nrow(fit_rows) > 0L) {
        eig_row <- fit_rows[
          fit_rows$statistic == paste0("eigenvalue.", fac), ,
          drop = FALSE
        ]
        if (nrow(eig_row) > 0L) {
          paste0("  eigenvalue ", sprintf("%.2f", eig_row$value[1L]))
        } else { # nocov start
          ""
        } # nocov end
      } else {
        ""
      }
      # M51: display "label (id)" when set; `fac` stays the ID for the
      # eigenvalue lookup above.
      fac_disp <- .label_id(fac, x$factor_labels)
      cli::cli_text("  {fac_disp}  {vpct}%{suffix}")
    }

    # EFA / ESEM: append a fit-statistic line per level (not per-factor)
    if ((is_efa || is_esem) && ki > 1L && nrow(fit_rows) > 0L) {
      idx_show <- if (is_esem) {
        c("CFI", "TLI", "RMSEA", "SRMR")
      } else {
        c("RMSEA", "TLI", "chi", "dof")
      }
      fit_sub <- fit_rows[fit_rows$statistic %in% idx_show, , drop = FALSE]
      if (nrow(fit_sub) > 0L) {
        cuts <- .fit_cutoffs()
        parts <- vapply(
          seq_len(nrow(fit_sub)),
          function(j) {
            val <- fit_sub$value[j]
            idx <- fit_sub$statistic[j]
            formatted <- if (idx %in% c("chi", "dof")) {
              paste0(idx, " = ", round(val, 1L))
            } else {
              paste0(idx, " = ", round(val, 3L))
            }
            cut <- cuts[[idx]]
            if (!is.null(cut) && !is.na(val)) {
              ok <- if (cut$direction == "hi") val >= cut$threshold else val <= cut$threshold
              formatted <- paste0(formatted, if (ok) " \u2714" else " \u2718")
            }
            formatted
          },
          character(1L)
        )
        cli::cli_text(
          "  {cli::col_grey(paste(parts, collapse = '  '))}"
        )
      }
    }
  }

  # ESEM under a scaled-test estimator: name the scaled reporting once (M32;
  # see tidy(what = "fit") docs for the full rationale). Not shown for ML,
  # which has no scaled variant and reports naive values.
  if (is_esem && isTRUE(x$estimator %in% c("WLSMV", "ULSMV", "MLR"))) {
    cli::cli_text(
      cli::col_grey(
        "  Fit indices above report lavaan's mean-and-variance-adjusted \\
         (\"scaled\") variant, per the {.val {x$estimator}} estimator."
      )
    )
  }

  # --- Lineage ----------------------------------------------------------------
  cli::cli_h2("Lineage (primary parents)")

  lin <- x$lineage
  if (nrow(lin) == 0L) {
    cli::cli_text(cli::col_grey("  (only one level - no edges to display)"))
  } else {
    for (i in seq_len(nrow(lin))) {
      cli::cli_text("  {lin$parent[i]} {cli::symbol$arrow_right} {lin$children[i]}")
    }
  }

  # --- Durable near-singularity caution (DESIGN.md s6) ------------------------
  if (isTRUE(x$near_singular)) {
    cli::cli_text("") # summary sets the caution off with a blank line
    .print_near_singular(
      x$min_eigenvalue,
      "per-level fit indices and factor scores may be unreliable -- the solution rests on a rank-deficient matrix."
    )
  }

  # --- Bootstrap edge CIs (M47) -----------------------------------------------
  if (!is.null(x$boot)) {
    cli::cli_h2("Bootstrap edge CIs")
    pct <- round(100 * x$boot$conf)
    n_ok_min <- min(x$boot$edges$n_boot_ok)
    cli::cli_text(
      "  {pct}% percentile CIs from {x$boot$n_boot} replicate{?s} \\
       ({n_ok_min}+ usable per edge)"
    )
    cli::cli_text(
      cli::col_grey(
        "  Per-edge se / lo / hi in {.code x$boot$edges} or {.code tidy(x)}."
      )
    )
  }

  # --- Pruning ----------------------------------------------------------------
  if (!is.null(x$prune)) {
    cli::cli_h2("Pruning")
    p <- x$prune
    if ("redundant" %in% p$rules) {
      nodes_str <- if (length(p$redundant) > 0L) {
        paste(p$redundant, collapse = ", ")
      } else {
        "(none)"
      }
      r_thr <- p$redundancy_r
      phi_note <- if (!is.null(p$redundancy_phi)) {
        paste0(", phi > ", p$redundancy_phi)
      } else {
        ""
      }
      cli::cli_text(
        "  Redundant ({p$redundancy_criterion}, |r| {cli::symbol$geq} {r_thr}{phi_note}): \\
         {length(p$redundant)} node{?s} flagged"
      )
      if (length(p$redundant) > 0L) {
        cli::cli_text("  Flagged: {nodes_str}")
      }
    }
    if ("artifact" %in% p$rules) {
      cli::cli_text(
        "  Artifact: Tucker's phi computed for {p$artifact_n} cross-level pair{?s}"
      )
      if (!is.null(p$structural_n)) {
        cli::cli_text(
          "  Structural signals: {p$structural_n} factor{?s} flagged \\
           (inspect {.code x$prune$structural})"
        )
      }
    }
    if (!is.null(p$manual) && length(p$manual) > 0L) {
      cli::cli_text(
        "  Manual: {length(p$manual)} node{?s} explicitly flagged \\
         ({paste(p$manual, collapse = ', ')})"
      )
    }
  }

  # --- Footer: one rule, then prune note (if any) + caveat --------------------
  prune_note <- if (!is.null(x$prune)) {
    "Note: Pruning is interpretive relabeling, not re-estimation. \\
     Flagged nodes remain in the object with all edges preserved."
  }
  .print_honesty_footer(prune_note)

  invisible(x)
}

# Build the lineage data frame: one row per parent with all primary children.
# Uses only adjacent-level primary edges (level_to == level_from + 1) so
# skip-level edges from pairs="all" don't appear in the tree listing.
# fill_primary() converts every remaining NA is_primary to FALSE before the
# object is assembled, so is_primary is never NA here; which() is kept as
# cheap defence against edge tables built outside that path.
# Explicit sort by level then label ensures stable top-down ordering.
.summary_lineage <- function(x) {
  e <- x$edges$tidy
  idx <- which(e$is_primary & (e$level_to == e$level_from + 1L))
  adj_primary <- e[idx, , drop = FALSE]

  if (nrow(adj_primary) == 0L) {
    return(data.frame(
      parent   = character(0L),
      children = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  parents <- unique(adj_primary$from)
  # Sort by level (level_from of first occurrence) then by label within level.
  parent_levels <- adj_primary$level_from[match(parents, adj_primary$from)]
  parents <- parents[order(parent_levels, parents)]

  # Factor labels (M51): display each node as "label (id)" when a label is set,
  # bare ID otherwise. Sorting stays keyed on IDs above, so labels never reorder.
  labels <- x$meta$factor_labels

  rows <- lapply(parents, function(p) {
    kids <- adj_primary$to[adj_primary$from == p]
    data.frame(
      parent = .label_id(p, labels),
      children = paste(.label_id(kids, labels), collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
