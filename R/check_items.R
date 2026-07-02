#' Screen items for problems before factor analysis
#'
#' @description
#' Real-world item sets often contain columns that break or *silently* degrade a
#' factor analysis: an item everyone answered the same way (no variance), an
#' item dominated by one response with a couple of stray answers, or an item
#' with heavy missingness. On the polychoric basis these are especially costly
#' -- a near-empty response category can make [psych::polychoric()] fail
#' outright (see the `correct` argument of [ackwards()]), and a near-constant
#' item can produce a plausible-looking but meaningless factor with no warning.
#'
#' `check_items()` reports these problems *before* you fit, one row per item, so
#' you can collapse rare categories, drop degenerate items, or set
#' `cor`/`correct` deliberately rather than debugging a cryptic error. It only
#' reports; it never changes your data. [ackwards()] runs the same screen
#' internally: it **errors** on a constant item and **warns** on a
#' near-degenerate one, naming the offenders.
#'
#' @param data A data frame or numeric matrix (items in columns).
#' @param cor The correlation basis you plan to use: `"polychoric"` (default),
#'   `"pearson"`, or `"spearman"`. Only affects which problems are flagged --
#'   sparse response categories destabilise the polychoric basis but not the
#'   others.
#'
#' @return A data frame (class `check_items`) with one row per item and columns
#'   `item`, `n_valid`, `pct_missing`, `n_distinct`, `min_count` (smallest
#'   observed category count), `top_prop` (proportion of valid responses in the
#'   most common value), and `flag` -- one of `"ok"`, `"constant"`,
#'   `"near-constant"`, `"sparse category"`, or `"high missing"`. Print it for a
#'   grouped summary with guidance; treat it as a plain data frame for the full
#'   per-item table.
#'
#' @seealso [ackwards()] and its `correct` argument (the polychoric failure this
#'   screens for); [suggest_k()] and [comparability()], the other pre-analysis
#'   diagnostics.
#'
#' @examples
#' # sim16 is clean continuous data -> nothing flagged
#' check_items(sim16, cor = "pearson")
#'
#' # An item dominated by one response with a couple of stray answers is flagged
#' d <- sim16
#' d$bad <- c(rep(1L, nrow(d) - 2L), 2L, 3L)
#' check_items(d)
#'
#' @export
check_items <- function(data, cor = c("polychoric", "pearson", "spearman")) {
  cor <- rlang::arg_match(cor)
  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
  }
  df <- as.data.frame(data)
  if (ncol(df) == 0L) {
    cli::cli_abort("{.arg data} has no columns to check.")
  }
  report <- .screen_items(df, cor)
  structure(report, class = c("check_items", "data.frame"), basis = cor)
}

# Per-item screen shared by check_items() and ackwards(). Returns a plain
# data.frame (one row per column of `df`) with the stats and a single worst-case
# `flag`. Thresholds are deliberately conservative so ordinary Likert items
# (which routinely have a small category) are not flagged; the failure mode that
# motivates this is a *dominant* category with near-empty others.
.screen_items <- function(df, cor) {
  dominant_prop <- 0.95 # top category this common -> near-constant
  min_category <- 5L # smallest observed category below this -> sparse (polychoric)
  max_missing <- 0.20 # more missing than this -> high missing

  rows <- lapply(names(df), function(nm) {
    x <- df[[nm]]
    n_total <- length(x)
    xv <- x[!is.na(x)]
    n_valid <- length(xv)
    pct_missing <- if (n_total > 0L) 1 - n_valid / n_total else 1
    tab <- table(xv)
    n_distinct <- length(tab)
    min_count <- if (n_distinct > 0L) as.integer(min(tab)) else 0L
    top_prop <- if (n_valid > 0L) max(tab) / n_valid else NA_real_

    flag <- if (n_distinct <= 1L) {
      "constant"
    } else if (!is.na(top_prop) && top_prop >= dominant_prop) {
      "near-constant"
    } else if (cor == "polychoric" && min_count < min_category) {
      "sparse category"
    } else if (pct_missing > max_missing) {
      "high missing"
    } else {
      "ok"
    }

    data.frame(
      item = nm,
      n_valid = n_valid,
      pct_missing = round(pct_missing, 3),
      n_distinct = n_distinct,
      min_count = min_count,
      top_prop = round(top_prop, 3),
      flag = flag,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Print an item quality check
#'
#' @param x A `check_items` object (produced by [check_items()]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.check_items <- function(x, ...) {
  cli::cli_h1("Item quality check ({.pkg ackwards})")
  cli::cli_dl(c(
    "Basis"   = cli::style_bold(attr(x, "basis")),
    "Items"   = as.character(nrow(x)),
    "Flagged" = as.character(sum(x$flag != "ok"))
  ))

  flagged <- x[x$flag != "ok", , drop = FALSE]
  if (nrow(flagged) == 0L) {
    cli::cli_text(cli::col_green("{cli::symbol$tick} No item problems detected."))
  } else {
    cli::cli_h2("Flagged items")
    # constant is fatal (red cross); the rest are advisory (yellow bang)
    for (fl in c("constant", "near-constant", "sparse category", "high missing")) {
      items_fl <- flagged$item[flagged$flag == fl]
      if (length(items_fl) > 0L) {
        sym <- if (fl == "constant") {
          cli::col_red(cli::symbol$cross)
        } else {
          cli::col_yellow("!")
        }
        cli::cli_text("{sym} {.strong {fl}}: {.val {items_fl}}")
      }
    }
  }

  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Constant items must be dropped (no variance). A near-constant item (one \\
       response dominates) can yield a meaningless factor; a sparse category can \\
       make {.code cor = \"polychoric\"} fail -- collapse rare categories, try \\
       {.code correct = 0}, or drop the item. Full per-item table: treat this \\
       object as a data frame."
    )
  )
  invisible(x)
}
