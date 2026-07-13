# ── Shared print/summary helpers (M59) ────────────────────────────────────────
# print.ackwards() and print.summary_ackwards() render the same fixed blocks;
# these give each block one source of truth so the two surfaces never drift.

# Engine/rotation/basis/n/k definition-list header.
.print_ba_header <- function(engine, rotation, cor, n_obs, k_max) {
  cor_label <- if (is.na(cor)) "(user-supplied matrix)" else cor
  n_label <- if (is.na(n_obs)) "NA" else format(n_obs, big.mark = ",")
  cli::cli_dl(c(
    "Engine"   = cli::style_bold(engine),
    "Rotation" = rotation,
    "Basis"    = cor_label,
    "n"        = n_label,
    "k (max)"  = as.character(k_max)
  ))
}

# Durable near-singularity caution (DESIGN.md s6). The yellow warning wrapper and
# the ?ackwards pointer are shared; `detail` is the surface-specific middle
# clause (print is terse, summary spells out the rank deficiency), passed in so
# each surface keeps its exact wording.
.print_near_singular <- function(min_eigenvalue, detail) {
  cli::cli_text(cli::col_yellow(
    "{cli::symbol$warning} Near-singular correlation matrix (min eigenvalue \\
     {signif(min_eigenvalue, 2)}): {detail} See {.code ?ackwards} \\
     (\"When to trust the result\")."
  ))
}

# Closing footer: one rule, an optional surface-supplied prune note, then the
# load-bearing "linked solutions, not a fitted hierarchy" honesty caveat (D-001).
# The caveat is one verbatim source of truth; the prune note wording differs
# slightly between surfaces, so it is passed in rather than unified.
.print_honesty_footer <- function(prune_note = NULL) {
  cli::cli_rule()
  if (!is.null(prune_note)) {
    cli::cli_text(cli::col_grey(prune_note))
  }
  cli::cli_text(
    cli::col_grey(
      "Note: This is a series of linked solutions, not a fitted hierarchical \\
       model. Cross-level edges are descriptive score correlations. \\
       Per-level fit indices (EFA/ESEM) describe how well a k-factor model \\
       fits the items at that level -- they do not validate the edges or \\
       the hierarchy itself."
    )
  )
}

# Extract pruning info for display: flagged redundant-node IDs, artifact/phi
# count, structural-signal count, and the redundancy thresholds. Read by BOTH
# print.ackwards() (live) and summary.ackwards() (stored on the summary object),
# so the two surfaces report identical counts. `rules` is carried through so each
# surface can gate on what was requested (rules="artifact" never flags redundant
# nodes). NULL when the object carries no pruning annotations.
.prune_digest <- function(x) {
  if (is.null(x$prune)) {
    return(NULL)
  }
  nodes <- x$prune$nodes
  redundant <- if (!is.null(nodes)) {
    nodes$id[nodes$pruned & nodes$prune_reason == "redundant"]
  } else { # nocov start
    character(0L)
  } # nocov end
  artifact_n <- if (!is.null(x$prune$phi)) nrow(x$prune$phi) else NULL
  structural_n <- if (!is.null(x$prune$structural)) {
    sum(
      x$prune$structural$few_items | x$prune$structural$orphan |
        x$prune$structural$split_merge,
      na.rm = TRUE
    )
  } else {
    NULL
  }
  list(
    rules                = x$prune$rules,
    redundant            = redundant,
    artifact_n           = artifact_n,
    structural_n         = structural_n,
    redundancy_r         = x$prune$redundancy_r,
    redundancy_phi       = x$prune$redundancy_phi,
    redundancy_criterion = x$prune$redundancy_criterion,
    manual               = x$prune$manual
  )
}

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
  .print_ba_header(x$engine, x$rotation, x$cor, x$n_obs, x$k_max)

  # --- Factor labels (M51) ----------------------------------------------------
  # Only surfaced when the user has attached labels; unlabeled objects print
  # exactly as before. Reports coverage; the labels themselves show in the
  # per-node surfaces (summary(), autoplot(), top_items()).
  if (!is.null(x$meta$factor_labels)) {
    n_set <- length(x$meta$factor_labels)
    n_tot <- length(unlist(lapply(x$levels, `[[`, "labels"), use.names = FALSE))
    cli::cli_text(cli::col_grey(
      "  Factor labels: {n_set} of {n_tot} factor{?s} labelled \\
       (see {.fn set_factor_labels})."
    ))
  }

  # --- Per-level table --------------------------------------------------------
  cli::cli_h2("Levels")

  K <- x$k_max
  for (ki in seq_len(K)) {
    lev <- x$levels[[as.character(ki)]]
    cum_pct <- .fmt_pct(lev$variance["cumulative"])
    conv_sym <- if (isTRUE(lev$converged)) {
      cli::col_green(.ok_glyph(TRUE))
    } else {
      cli::col_red(.ok_glyph(FALSE))
    }
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

  # --- Durable near-singularity caution (DESIGN.md s6) ------------------------
  if (isTRUE(x$meta$near_singular)) {
    .print_near_singular(
      x$meta$min_eigenvalue,
      "fit indices and factor scores may be unreliable."
    )
  }

  # --- Pruning summary (Forbes extension; DESIGN.md s14.18) -------------------
  # Counts come from the shared .prune_digest() so print and summary can never
  # disagree on the redundant / phi-note / structural-signal figures (M59). The
  # display wording stays per-surface (print is terser).
  if (!is.null(x$prune)) {
    cli::cli_h2("Pruning")
    pd <- .prune_digest(x)
    if ("redundant" %in% pd$rules) {
      n_flagged <- length(pd$redundant)
      r_thr <- pd$redundancy_r
      phi_note <- if (!is.null(pd$redundancy_phi)) {
        paste0(", phi > ", pd$redundancy_phi)
      } else {
        ""
      }
      crit <- pd$redundancy_criterion
      cli::cli_text(
        "  Redundancy ({crit}, |r| {cli::symbol$geq} {r_thr}{phi_note}): \\
         {n_flagged} node{?s} flagged"
      )
    }
    if ("artifact" %in% pd$rules) {
      n_phi <- pd$artifact_n %||% 0L
      cli::cli_text(
        "  Artifact: Tucker's phi computed for {n_phi} cross-level factor pair{?s}"
      )
      if (!is.null(pd$structural_n)) {
        cli::cli_text(
          "  Structural signals: {pd$structural_n} factor{?s} flagged \\
           (inspect {.code x$prune$structural})"
        )
      }
    }
    if (!is.null(pd$manual) && length(pd$manual) > 0L) {
      cli::cli_text(
        "  Manual: {length(pd$manual)} node{?s} explicitly flagged \\
         ({paste(pd$manual, collapse = ', ')})"
      )
    }
  }

  # --- Footer: one rule, then prune note (if any) + caveat -------------------
  prune_note <- if (!is.null(x$prune)) {
    "Note: Pruning is interpretive relabeling, not re-estimation. \\
     Flagged nodes remain in the object; all edges are preserved. \\
     Inspect with {.code x$prune$nodes} and {.code tidy(x, what = \"nodes\")}."
  }
  .print_honesty_footer(prune_note)

  invisible(x)
}
