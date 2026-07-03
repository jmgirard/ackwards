#' Display the salient items for each factor
#'
#' Returns, per level and factor, the items whose absolute loading meets or
#' exceeds `cut`, sorted by descending absolute loading. This gives a concise
#' reading of "what each factor is about" without printing a full item-by-factor
#' matrix, which does not scale well to large `k` or many items.
#'
#' Loadings are signed and reflect the object's primary-parent sign alignment
#' (see [ackwards()]). Items near the cut threshold may appear for one sign
#' orientation but not the other; this is expected and informative.
#'
#' If [factor labels][set_factor_labels] have been attached, the factor
#' dimension is shown as `label (id)` wherever it appears -- the group headers
#' under `by = "factor"` and the body entries under `by = "item"`.
#'
#' @param x An `ackwards` object.
#' @param level Integer vector selecting which level(s) to include. `NULL`
#'   (default) returns all levels.
#' @param cut Absolute-loading threshold. Items with `|loading| >= cut` are
#'   shown. Default `0.3`, matching the `cut_show` plotting default.
#' @param n Maximum number of items to show per factor. `NULL` (default) shows
#'   all items meeting the cut. When set, the top-`n` items by `|loading|` are
#'   kept after applying the cut.
#' @param sort Logical. When `TRUE` (default), items within each group are
#'   ordered by descending `|loading|`. Set to `FALSE` to keep the original
#'   order (useful when items have a meaningful sequence).
#' @param by One of `"factor"` (default) or `"item"`. `"factor"` groups the
#'   listing by factor (the salient items *of* each factor -- "what is this
#'   factor about?"). `"item"` inverts the grouping to list, for each item, the
#'   factors it loads on -- which makes cross-loadings legible ("where does this
#'   item go?"). `n` and `sort` apply within whichever unit `by` selects.
#' @param show_labels Logical. When `TRUE` (default) and the data carried a
#'   variable-label attribute at fit time (see [ackwards()]), items are shown as
#'   `id: label`; items without a label fall back to the bare id. Set to
#'   `FALSE` to always show the bare `m{k}f{j}`-style item ids.
#'
#' @return An object of class `"top_items"`. Print it for a grouped cli listing.
#'   The underlying data frame (one row per selected item) is accessible via
#'   `$data` and contains columns `level`, `factor`, `item`, and `loading` (plus
#'   `label` when labels are available). The values equal the corresponding
#'   `tidy(x, what = "loadings")` rows (after filtering and optional sorting).
#'
#' @seealso [tidy.ackwards()], [label_template()], [set_factor_labels()],
#'   [autoplot.ackwards()]
#'
#' @examples
#' # Fit the raw dataset (not na.omit(), which would drop the column
#' # attributes): bfi25's IPIP item labels are then captured and printed as
#' # `code: label`. `missing = "listwise"` handles the NAs cleanly.
#' x <- ackwards(bfi25, k_max = 5, cor = "polychoric", missing = "listwise")
#' top_items(x)
#' top_items(x, level = 5, cut = 0.4, n = 5)
#'
#' # Invert the grouping to read cross-loadings item-by-item
#' top_items(x, level = 3, cut = 0.25, by = "item")
#'
#' @export
top_items <- function(x, level = NULL, cut = 0.3, n = NULL, sort = TRUE,
                      by = c("factor", "item"), show_labels = TRUE) {
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("{.arg x} must be an {.cls ackwards} object.")
  }
  if (!is.numeric(cut) || length(cut) != 1L || cut < 0 || cut > 1) {
    cli::cli_abort("{.arg cut} must be a single number between 0 and 1.")
  }
  if (!is.null(n) && (!is.numeric(n) || length(n) != 1L || n < 1L)) {
    cli::cli_abort("{.arg n} must be a positive integer or NULL.")
  }
  if (!is.logical(sort) || length(sort) != 1L) {
    cli::cli_abort("{.arg sort} must be TRUE or FALSE.")
  }
  if (!is.logical(show_labels) || length(show_labels) != 1L) {
    cli::cli_abort("{.arg show_labels} must be TRUE or FALSE.")
  }
  by <- rlang::arg_match(by)

  # Effective labels: only carried when the user opts in and the object has them.
  eff_labels <- if (isTRUE(show_labels)) x$meta$item_labels else NULL

  available_levels <- as.integer(names(x$levels))

  if (!is.null(level)) {
    level <- as.integer(level)
    bad <- setdiff(level, available_levels)
    if (length(bad) > 0L) {
      cli::cli_abort(c(
        "!" = "Requested level{?s} not found in {.arg x}: {bad}.",
        "i" = "Available levels: {available_levels}."
      ))
    }
    keep_levels <- as.character(level)
  } else {
    keep_levels <- names(x$levels)
  }

  # Build the full filtered long table first (one row per surviving
  # item x factor x level), retaining the column/row indices so that both
  # grouping orientations and the "keep original order" (sort = FALSE) mode can
  # be honoured downstream regardless of `by`.
  rows <- lapply(keep_levels, function(ki) {
    L <- x$levels[[ki]]$loadings
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(L)), function(j) {
      keep <- which(abs(L[, j]) >= cut)
      if (length(keep) == 0L) {
        return(NULL)
      }
      data.frame(
        level = k,
        factor = colnames(L)[j],
        item = rownames(L)[keep],
        loading = L[keep, j],
        .factor_ord = j,
        .item_ord = keep,
        stringsAsFactors = FALSE
      )
    }))
  })

  df <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(df)) {
    df <- data.frame(
      level = integer(0L),
      factor = character(0L),
      item = character(0L),
      loading = numeric(0L),
      .factor_ord = integer(0L),
      .item_ord = integer(0L),
      stringsAsFactors = FALSE
    )
  }

  # Group by the chosen unit (factor or item); order groups by their original
  # index within the level, then order rows within each group (by descending
  # |loading| when sort = TRUE, else by the other dimension's original order).
  grp_ord <- if (by == "factor") df$.factor_ord else df$.item_ord
  sub_ord <- if (by == "factor") df$.item_ord else df$.factor_ord
  within_key <- if (sort) -abs(df$loading) else sub_ord
  df <- df[order(df$level, grp_ord, within_key, sub_ord), , drop = FALSE]

  # Cap each group to the top-n rows (post-sort).
  if (!is.null(n) && nrow(df) > 0L) {
    grp_val <- if (by == "factor") df$factor else df$item
    grp_key <- paste(df$level, grp_val, sep = "\r")
    rank <- stats::ave(seq_len(nrow(df)), grp_key, FUN = seq_along)
    df <- df[rank <= n, , drop = FALSE]
  }

  df$.factor_ord <- NULL
  df$.item_ord <- NULL
  if (!is.null(eff_labels)) {
    df$label <- unname(eff_labels[df$item])
  }
  rownames(df) <- NULL

  structure(
    list(
      data           = df,
      levels_shown   = as.integer(keep_levels),
      cut            = cut,
      n              = n,
      sort           = sort,
      by             = by,
      item_labels    = eff_labels,
      factor_labels  = x$meta$factor_labels, # M51: shown on factor headers
      engine         = x$engine,
      k_max          = x$k_max
    ),
    class = "top_items"
  )
}

# Format an item id for display: "id: label" when a label is available,
# otherwise the bare id. Leading with the id keeps the list aligned with the
# codes used in the loadings/edge tables. `labs` is the (possibly NULL) named
# label vector.
.format_item_label <- function(id, labs) {
  if (is.null(labs)) {
    return(id)
  }
  lab <- unname(labs[id])
  ifelse(is.na(lab), id, paste0(id, ": ", lab))
}

#' Print a top_items object
#'
#' @param x A `top_items` object (produced by [top_items()]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.top_items <- function(x, ...) {
  by_item <- identical(x$by, "item")
  cli::cli_h1(
    "Salient {if (by_item) 'factors by item' else 'items by factor'} ({.pkg ackwards})"
  )
  cli::cli_dl(c(
    "Engine"  = cli::style_bold(x$engine),
    "Cut"     = paste0("|loading| >= ", x$cut),
    "Top-n"   = if (is.null(x$n)) "all" else as.character(x$n)
  ))

  df <- x$data
  labs <- x$item_labels

  if (nrow(df) == 0L) {
    cli::cli_text(cli::col_grey(
      "No items met the |loading| >= {x$cut} threshold."
    ))
    return(invisible(x))
  }

  for (k in x$levels_shown) {
    df_k <- df[df$level == k, , drop = FALSE]
    if (nrow(df_k) == 0L) next
    cli::cli_h2("Level {k} ({k} factor{?s})")

    # Group header = factor (default) or item; body lists the other dimension.
    group_col <- if (by_item) "item" else "factor"
    body_col <- if (by_item) "factor" else "item"
    groups <- unique(df_k[[group_col]])
    for (gi in seq_along(groups)) {
      grp <- groups[[gi]]
      if (gi > 1L) cli::cli_text("") # blank line between factor/item groups
      df_g <- df_k[df_k[[group_col]] == grp, , drop = FALSE]
      # M51: factor headers show "label (id)" when a factor label is set;
      # item headers keep the item "id: label" form.
      header <- if (by_item) {
        .format_item_label(grp, labs)
      } else {
        .label_id(grp, x$factor_labels)
      }
      cli::cli_text("{.strong {header}}")

      for (i in seq_len(nrow(df_g))) {
        loading_str <- formatC(df_g$loading[i], digits = 3L, format = "f")
        # M51: the body's factor entries (under by = "item") also show
        # "label (id)"; the body's item entries (under by = "factor") show the
        # item "id: label" form.
        entry <- if (by_item) {
          .label_id(df_g[[body_col]][i], x$factor_labels)
        } else {
          .format_item_label(df_g[[body_col]][i], labs)
        }
        cli::cli_text("  {entry}  [{loading_str}]")
      }
    }
  }

  cli::cli_rule()
  cli::cli_text(cli::col_grey(
    "Loadings reflect primary-parent sign alignment. \\
     Use tidy(x, what = \"loadings\") for the full matrix."
  ))

  invisible(x)
}
