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
#' x <- ackwards(bfi25, k_max = 5)
#' summary(x)
#'
#' @export
summary.ackwards <- function(object, ...) {
  structure(
    list(
      call       = object$call,
      engine     = object$engine,
      rotation   = object$rotation,
      cor        = object$cor,
      n_obs      = object$n_obs,
      k_max      = object$k_max,
      variance   = .tidy_variance(object),
      fit        = .tidy_fit(object),
      lineage    = .summary_lineage(object),
      prune      = .summary_prune(object)
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

  cor_label <- if (is.na(x$cor)) "(user-supplied matrix)" else x$cor
  n_label <- if (is.na(x$n_obs)) "NA" else format(x$n_obs, big.mark = ",")
  cli::cli_dl(c(
    "Engine"   = cli::style_bold(x$engine),
    "Rotation" = x$rotation,
    "Basis"    = cor_label,
    "n"        = n_label,
    "k (max)"  = as.character(x$k_max)
  ))

  # --- Per-level variance + fit -----------------------------------------------
  cli::cli_h2("Levels")

  is_esem <- x$engine == "esem"
  is_efa <- x$engine == "efa"
  is_pca <- x$engine == "pca"

  fit_tbl <- x$fit # pre-pulled for readability

  for (ki in seq_len(x$k_max)) {
    var_rows <- x$variance[x$variance$level == ki, , drop = FALSE]
    # cumulative_pct is a running sum within the level; max gives the total.
    cum_total <- max(var_rows$cumulative_pct)
    cli::cli_text(
      "{.strong k = {ki}}: {ki} factor{?s}  ({cum_total}% cumulative variance)"
    )
    fit_rows <- fit_tbl[fit_tbl$level == ki, , drop = FALSE]

    for (i in seq_len(nrow(var_rows))) {
      fac <- var_rows$factor[i]
      vpct <- var_rows$variance_pct[i]

      # PCA: eigenvalue indices are "eigenvalue.<label>" (see engine_pca.R:65)
      suffix <- if (is_pca && nrow(fit_rows) > 0L) {
        eig_row <- fit_rows[
          fit_rows$index == paste0("eigenvalue.", fac), ,
          drop = FALSE
        ]
        if (nrow(eig_row) > 0L) {
          paste0("  eigenvalue ", round(eig_row$value[1L], 2))
        } else {
          ""
        }
      } else {
        ""
      }
      cli::cli_text("  {fac}  {vpct}%{suffix}")
    }

    # EFA / ESEM: append a fit-index line per level (not per-factor)
    if ((is_efa || is_esem) && ki > 1L && nrow(fit_rows) > 0L) {
      idx_show <- if (is_esem) {
        c("CFI", "TLI", "RMSEA", "SRMR")
      } else {
        c("RMSEA", "TLI", "chi", "dof")
      }
      fit_sub <- fit_rows[fit_rows$index %in% idx_show, , drop = FALSE]
      if (nrow(fit_sub) > 0L) {
        parts <- vapply(
          seq_len(nrow(fit_sub)),
          function(j) {
            val <- fit_sub$value[j]
            idx <- fit_sub$index[j]
            if (idx %in% c("chi", "dof")) {
              paste0(idx, " = ", round(val, 1L))
            } else {
              paste0(idx, " = ", round(val, 3L))
            }
          },
          character(1L)
        )
        cli::cli_text(
          "  {cli::col_grey(paste(parts, collapse = '  '))}"
        )
      }
    }
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
        "  Redundant (|r| {cli::symbol$geq} {r_thr}{phi_note}): \\
         {length(p$redundant)} node{?s} flagged"
      )
      if (length(p$redundant) > 0L) {
        cli::cli_text("  Flagged: {nodes_str}")
      }
    }
    if ("artefact" %in% p$rules) {
      cli::cli_text(
        "  Artefact: Tucker's phi computed for {p$artefact_n} cross-level pair{?s}"
      )
    }
    cli::cli_rule()
    cli::cli_text(
      cli::col_grey(
        "Note: Pruning is interpretive relabeling, not re-estimation. \\
         Flagged nodes remain in the object with all edges preserved."
      )
    )
  }

  # --- Caveat -----------------------------------------------------------------
  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Note: This is a series of linked solutions, not a fitted hierarchical \\
       model. Cross-level edges are descriptive score correlations."
    )
  )

  invisible(x)
}

# Build the lineage data frame: one row per parent with all primary children.
# Uses only adjacent-level primary edges (level_to == level_from + 1) so
# skip-level edges from pairs="all" don't appear in the tree listing.
# which() is used instead of direct logical indexing for NA-safety: fill_primary()
# leaves is_primary = NA on skip-level edges, and NA & FALSE = NA (row included).
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

  rows <- lapply(parents, function(p) {
    kids <- adj_primary$to[adj_primary$from == p]
    data.frame(
      parent = p,
      children = paste(kids, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Extract pruning info for display: lists of flagged-node IDs and artefact count.
# rules is carried through so print.summary_ackwards can gate on what was requested
# (e.g. prune="artefact" never flags redundant nodes).
.summary_prune <- function(x) {
  if (is.null(x$prune)) {
    return(NULL)
  }
  nodes <- x$prune$nodes
  redundant <- if (!is.null(nodes)) {
    nodes$id[nodes$pruned & nodes$prune_reason == "redundant"]
  } else {
    character(0L)
  }
  artefact_n <- if (!is.null(x$prune$phi)) nrow(x$prune$phi) else NULL
  list(
    rules          = x$prune$rules,
    redundant      = redundant,
    artefact_n     = artefact_n,
    redundancy_r   = x$prune$redundancy_r,
    redundancy_phi = x$prune$redundancy_phi
  )
}
