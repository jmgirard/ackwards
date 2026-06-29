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
#' @param x An `ackwards` object.
#' @param level Integer vector selecting which level(s) to include. `NULL`
#'   (default) returns all levels.
#' @param cut Absolute-loading threshold. Items with `|loading| >= cut` are
#'   shown. Default `0.3`, matching the `cut_show` plotting default.
#' @param n Maximum number of items to show per factor. `NULL` (default) shows
#'   all items meeting the cut. When set, the top-`n` items by `|loading|` are
#'   kept after applying the cut.
#' @param sort Logical. When `TRUE` (default), items within each factor are
#'   ordered by descending `|loading|`. Set to `FALSE` to keep the original
#'   item order (useful when items have a meaningful sequence).
#'
#' @return An object of class `"top_items"`. Print it for a grouped
#'   level-by-factor cli listing. The underlying data frame (one row per
#'   selected item) is accessible via `$data` and contains columns `level`,
#'   `factor`, `item`, and `loading`. The values equal the corresponding
#'   `tidy(x, what = "loadings")` rows (after filtering and optional sorting).
#'
#' @seealso [tidy.ackwards()], [label_template()], [autoplot.ackwards()]
#'
#' @examples
#' x <- ackwards(bfi25, k_max = 5)
#' top_items(x)
#' top_items(x, level = 5, cut = 0.4, n = 5)
#'
#' @export
top_items <- function(x, level = NULL, cut = 0.3, n = NULL, sort = TRUE) {
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

  rows <- lapply(keep_levels, function(ki) {
    lev <- x$levels[[ki]]
    L <- lev$loadings
    k <- as.integer(ki)
    do.call(rbind, lapply(seq_len(ncol(L)), function(j) {
      loadings_j <- L[, j]
      keep <- abs(loadings_j) >= cut
      if (!any(keep)) {
        return(NULL)
      }
      items <- rownames(L)[keep]
      vals <- loadings_j[keep]
      if (sort) {
        ord <- order(abs(vals), decreasing = TRUE)
        items <- items[ord]
        vals <- vals[ord]
      }
      if (!is.null(n)) {
        take <- seq_len(min(n, length(items)))
        items <- items[take]
        vals <- vals[take]
      }
      data.frame(
        level = k,
        factor = colnames(L)[j],
        item = items,
        loading = vals,
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
      stringsAsFactors = FALSE
    )
  }
  rownames(df) <- NULL

  structure(
    list(
      data           = df,
      levels_shown   = as.integer(keep_levels),
      cut            = cut,
      n              = n,
      sort           = sort,
      engine         = x$engine,
      k_max          = x$k_max
    ),
    class = "top_items"
  )
}

#' Print a top_items object
#'
#' @param x A `top_items` object (produced by [top_items()]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.top_items <- function(x, ...) {
  cli::cli_h1("Salient items by factor ({.pkg ackwards})")
  cli::cli_dl(c(
    "Engine"  = x$engine,
    "Cut"     = paste0("|loading| >= ", x$cut),
    "Top-n"   = if (is.null(x$n)) "all" else as.character(x$n)
  ))

  df <- x$data

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

    for (fac in unique(df_k$factor)) {
      df_f <- df_k[df_k$factor == fac, , drop = FALSE]
      cli::cli_text("{.strong {fac}}")

      for (i in seq_len(nrow(df_f))) {
        loading_str <- formatC(df_f$loading[i], digits = 3L, format = "f")
        cli::cli_text("  {df_f$item[i]}  [{loading_str}]")
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
