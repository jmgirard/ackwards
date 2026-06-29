#' Generate a node-label scaffold for autoplot
#'
#' Returns a named character vector covering all factor IDs in the object,
#' ready to pass to `autoplot(x, node_labels = ...)`. The vector is also
#' printed as an editable `c(...)` literal so you can copy it into a script,
#' fill in substantive labels, and use it directly.
#'
#' The factor IDs are returned in the same left-to-right, top-to-bottom order
#' used by [ba_layout()] and [autoplot.ackwards()], so the printed literal maps
#' directly onto the diagram.
#'
#' @section Style options:
#' * `"id"` *(default)* -- every value equals the factor ID (`"m1f1"`, `"m2f1"`,
#'   ...). This is a round-trip no-op: passing the result to `node_labels` without
#'   editing reproduces the default labels exactly. Useful as the starting point
#'   for adding substantive labels.
#' * `"forbes"` -- values follow the Forbes (2023) convention: level-letter +
#'   within-level index (`"A1"`, `"B1"`, `"B2"`, ...). Level 1 -> `A`, level 2 ->
#'   `B`, level 3 -> `C`, and so on. Within-level indices are assigned in
#'   canonical layout order (left to right). Requires `k_max <= 26` (LETTERS has
#'   26 entries); an error is raised for deeper objects.
#' * `"blank"` -- all values are empty strings. Useful as a starting scaffold
#'   when you want to supply every label from scratch with no defaults showing
#'   through.
#'
#' @param x An `ackwards` object.
#' @param style One of `"id"` (default), `"forbes"`, or `"blank"`. See
#'   **Style options** above.
#'
#' @return A named character vector: names are factor IDs (`"m{k}f{j}"`),
#'   values are the label strings for the chosen `style`. The vector is
#'   suitable for direct use as the `node_labels` argument to
#'   [autoplot.ackwards()]. The `c(...)` literal is also printed to the
#'   console for copy-paste.
#'
#' @seealso [autoplot.ackwards()], [top_items()], [ba_layout()]
#'
#' @examples
#' if (requireNamespace("psych", quietly = TRUE)) {
#'   x <- ackwards(psych::bfi[, 1:25], k_max = 5)
#'
#'   # Start from ID defaults, then fill in your own labels:
#'   labs <- label_template(x)
#'   labs["m5f1"] <- "Neuroticism"
#'
#'   # Forbes letter convention:
#'   label_template(x, style = "forbes")
#' }
#'
#' @export
label_template <- function(x, style = c("id", "forbes", "blank")) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }
  style <- match.arg(style)

  # Get factor IDs in canonical layout order (same as autoplot uses)
  layout_nodes <- ba_layout(x)$nodes
  node_ids <- layout_nodes$id

  values <- switch(style,
    id = node_ids,
    blank = rep("", length(node_ids)),
    forbes = {
      levels_vec <- layout_nodes$level
      max_level <- max(levels_vec)
      if (max_level > 26L) {
        cli::cli_abort(c(
          "style = \"forbes\" requires at most 26 levels (LETTERS has 26 entries).",
          "x" = "This object has {max_level} levels.",
          "i" = "Use style = \"id\" or style = \"blank\" instead."
        ))
      }
      level_seq <- lapply(unique(levels_vec), function(k) {
        ids_at_k <- node_ids[levels_vec == k]
        letter <- LETTERS[k]
        paste0(letter, seq_along(ids_at_k))
      })
      unlist(level_seq, use.names = FALSE)
    }
  )

  out <- stats::setNames(values, node_ids)

  # Print an editable literal the user can copy into their script
  pairs <- paste0('  "', names(out), '" = "', out, '"', collapse = ",\n")
  cli::cli_text("{.fn label_template} scaffold ({style} style):")
  cat(paste0("c(\n", pairs, "\n)\n"))

  invisible(out)
}
