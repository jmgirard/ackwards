# R/factor_labels.R -- persistent factor labels (M51)
#
# Design (DESIGN.md s.14 item 45):
#   - Factor labels are the *factor*-level counterpart to M50's *item* labels
#     (meta$item_labels). The two are kept lexically distinct: nothing here
#     touches item_labels, and every user-facing string says "factor label".
#   - Stored in x$meta$factor_labels (named character vector, names = factor IDs
#     "m{k}f{j}"). Living in meta means it rides through prune()/boot_edges()/
#     augment()/predict() unmodified -- those verbs already copy meta.
#   - Purely additive: an unlabeled object (the default, NULL) behaves exactly
#     as it did pre-M51. Lineage still lives in IDs (Invariant 5); labels are
#     display only and never replace an ID column.

# Internal: every factor ID in the object, in level-then-within-level order
# (same enumeration prune() uses). Kept private -- callers that need canonical
# *layout* order use ba_layout(); this is just the membership set for validation.
.all_factor_ids <- function(x) {
  unlist(lapply(x$levels, `[[`, "labels"), use.names = FALSE)
}

#' Attach persistent factor labels to an ackwards object
#'
#' Store substantive names for factors (e.g. `"Neuroticism"` for `"m5f1"`) on
#' the object itself, so that [print()][print.ackwards], [summary()][summary.ackwards],
#' [tidy()][tidy.ackwards], [autoplot()][autoplot.ackwards], and
#' [top_items()] display them without re-supplying the labels each time. This is
#' the *factor*-label counterpart to item / variable labels (see
#' [top_items()] and `?ackwards`); the two are distinct and never interchanged.
#'
#' Labels are display only -- they never change a factor's stable ID
#' (`m{k}f{j}`), and every lineage / edge / score column continues to key on the
#' ID. A factor with no label falls back to its ID everywhere. The stored labels
#' ride along through [prune()], [boot_edges()], [augment()][augment.ackwards],
#' and [predict()][predict.ackwards] unchanged.
#'
#' @section Updating and clearing:
#' `set_factor_labels()` **merges** into any labels already stored, so you can
#' build them up incrementally. Within a single call:
#' * a normal string sets (or overwrites) that factor's label;
#' * an `NA` or `""` value **removes** just that factor's label;
#' * passing `labels = NULL` clears **all** labels at once.
#'
#' The scaffold printed by [label_template()] is a convenient starting point:
#' copy its `c(...)` literal, fill in the substantive names, and pass it here.
#'
#' @param x An `ackwards` object.
#' @param labels A named character vector mapping factor IDs (`"m{k}f{j}"`) to
#'   label strings, or `NULL` to clear all stored labels. Every name must match
#'   a factor ID in `x` (an unknown ID is an error, not a warning -- the object's
#'   IDs are knowable up front; see [factor_labels()] to read the current set).
#'
#' @return The `ackwards` object, with `meta$factor_labels` updated. Pipeable.
#'
#' @seealso [factor_labels()] to read them back, [label_template()] for a
#'   ready-to-edit scaffold, [autoplot.ackwards()] (`node_labels` overrides a
#'   stored label per node).
#'
#' @examples
#' x <- ackwards(sim16, k_max = 4, engine = "pca")
#'
#' # Start from a scaffold, fill in names, attach them:
#' x <- set_factor_labels(x, c(m4f1 = "Alpha", m4f2 = "Beta"))
#' factor_labels(x)
#'
#' # Merge in more; remove one; everything downstream now shows the labels:
#' x <- set_factor_labels(x, c(m2f1 = "Broad", m4f2 = NA))
#' summary(x)
#'
#' @export
set_factor_labels <- function(x, labels) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }

  if (is.null(labels)) {
    x$meta$factor_labels <- NULL
    return(x)
  }

  # A pure removal call reads naturally as `c(m4f2 = NA)`, which R types as
  # logical -- accept it (coerce to character NA, i.e. "remove"). Any other
  # non-character input is a genuine type error.
  if (is.logical(labels) && all(is.na(labels))) {
    labels <- stats::setNames(rep(NA_character_, length(labels)), names(labels))
  }
  if (!is.character(labels)) {
    cli::cli_abort("{.arg labels} must be a named character vector or {.val NULL}.")
  }
  nm <- names(labels)
  if (is.null(nm) || any(nm == "" | is.na(nm))) {
    cli::cli_abort("{.arg labels} must be a {.emph named} character vector (names are factor IDs).")
  }

  valid_ids <- .all_factor_ids(x)
  unknown <- setdiff(nm, valid_ids)
  if (length(unknown) > 0L) {
    cli::cli_abort(c(
      "Some {.arg labels} name{?s} match no factor ID in {.arg x}: {.val {unknown}}.",
      "i" = "Valid IDs run from {.val {valid_ids[[1]]}} to {.val {valid_ids[[length(valid_ids)]]}}."
    ))
  }

  current <- x$meta$factor_labels %||% stats::setNames(character(0), character(0))

  # Split into removals (NA / "") and sets, so a single call can do both.
  drop <- nm[is.na(labels) | labels == ""]
  keep <- labels[!(is.na(labels) | labels == "")]

  current <- current[setdiff(names(current), drop)]
  current[names(keep)] <- keep

  x$meta$factor_labels <- if (length(current) > 0L) current else NULL
  x
}

#' Read the factor labels stored on an ackwards object
#'
#' @param x An `ackwards` object.
#'
#' @return The named character vector of factor labels (names are factor IDs), or
#'   `NULL` if none have been set. See [set_factor_labels()] to attach them.
#'
#' @seealso [set_factor_labels()]
#'
#' @examples
#' x <- ackwards(sim16, k_max = 4, engine = "pca")
#' factor_labels(x) # NULL -- none set yet
#' x <- set_factor_labels(x, c(m4f1 = "Alpha"))
#' factor_labels(x)
#'
#' @export
factor_labels <- function(x) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }
  x$meta$factor_labels
}

# Internal: format a factor ID as "label (id)" when a label is stored, else the
# bare ID. Vectorised over `ids`. Shared by print()/summary()/top_items() so the
# display form (DESIGN s.14 item 45) lives in one place.
.label_id <- function(ids, labels) {
  if (is.null(labels)) {
    return(ids)
  }
  lab <- unname(labels[ids])
  ifelse(is.na(lab), ids, paste0(lab, " (", ids, ")"))
}
