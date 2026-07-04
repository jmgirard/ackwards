# Split-half factor comparability (Everett 1983; Goldberg 1990) -- the
# replicability gate for hierarchy depth. See DESIGN.md s.14 item 35.
#
# Mechanics per split:
#   1. Randomly divide the rows into two halves; fit levels 1..k in each half.
#   2. Anchor both half-solutions to the full-sample solution: match each
#      full-sample factor to its best half-solution counterpart (greedy
#      max-|r| with removal; a bijection is well-posed here because both
#      solutions have the same k, unlike parent matching in s.7).
#   3. The comparability coefficient for full-sample factor j is the
#      correlation between the two matched half-solution scores, computed on
#      the POOLED correlation matrix through the same W'RW algebra as every
#      other score correlation in the package (Invariant 1): the cross-
#      solution matrix is literally compute_edges() on a two-element levels
#      list. Tucker's phi on the matched loading columns is reported alongside
#      (report-first, flag-second -- nothing is auto-judged).

#' Split-half factor comparability
#'
#' Measures how well each factor at each level of a bass-ackwards hierarchy
#' **replicates** across random split-halves of the sample -- Everett's (1983)
#' factor comparability coefficients, the depth gate Goldberg's own lab applied
#' to its hierarchies (Goldberg, 1990) and the direct instrument for the
#' overextraction caution in [suggest_k()]: non-replicable structure
#' concentrates in the deeper levels of an overextracted hierarchy
#' (Forbes, 2023).
#'
#' For each of `n_splits` random half-splits, solutions at every level
#' `1..k_max` are fit independently in each half. Each half-solution's factors
#' are matched to the **full-sample** solution's factors (so coefficients are
#' reported under the same `m{k}f{j}` labels you get from [ackwards()]), and
#' the comparability coefficient for a factor is the correlation between its
#' two matched half-solution scores, computed on the pooled correlation matrix
#' via the same `W'RW` algebra used for between-level edges -- applying both
#' halves' scoring weights to the full sample, exactly Everett's procedure.
#' Tucker's congruence coefficient (phi) between the matched half-solution
#' loading columns is reported alongside: comparability asks whether the two
#' halves' *scores* agree; phi asks whether their *loading patterns* agree.
#'
#' @section Interpreting the output:
#' Coefficients near 1 mean the factor re-emerges in independent half-samples;
#' a factor whose comparability is low is sample-idiosyncratic and should not
#' anchor substantive interpretation. Goldberg's lab treated roughly .95 as
#' comfortable and .90 as a floor; these are conventions, not tests, so
#' `comparability()` reports every coefficient and flags nothing. The deepest
#' level at which all factors replicate is a natural **hierarchy floor** for
#' [ackwards()]'s `k_max`; see `vignette("ackwards-girard")` for the full
#' workflow.
#'
#' A level that fails to converge in a half-sample yields `NA` coefficients
#' for that split (convergence is data, not an error); the number of usable
#' splits per factor is reported in `summary$n_splits_ok` and a message
#' summarises any shortfall.
#'
#' @param data A data frame or numeric matrix of observed variables (items in
#'   columns, observations in rows). **Raw data only** -- splitting needs rows,
#'   so a correlation matrix is not accepted. Missing values are handled
#'   pairwise throughout (as in [ackwards()]'s default).
#' @param k_max Maximum number of factors/components to evaluate -- normally
#'   the same value (or one or two above it) you intend to pass to
#'   [ackwards()]. Required.
#' @param engine Extraction engine: `"pca"` (default) or `"efa"`. `"esem"` is
#'   not yet supported here (fitting `2 * n_splits` lavaan hierarchies needs
#'   its own performance treatment); for ESEM workflows, run
#'   `comparability()` with `engine = "efa"` as a structural screen.
#' @param cor Correlation basis: `"pearson"` (default) or `"spearman"`. As
#'   with [suggest_k()], `"polychoric"` is not supported (estimating
#'   polychoric matrices in every half-sample is slow and unstable); users
#'   analysing ordinal data should screen replicability on the Pearson basis
#'   and fit the final model with `cor = "polychoric"` in [ackwards()].
#' @param fm Factor extraction method passed to [psych::fa()]; only used when
#'   `engine = "efa"`. One of `"minres"` (default), `"ml"`, or `"pa"`.
#' @param n_splits Number of random split-half replicates. Default `10L`
#'   (repeated random splits, following Goldberg's practice -- a single split
#'   can mislead by luck of the draw). Each replicate fits `2 * k_max`
#'   solutions, so the default costs 20 hierarchy fits; PCA and EFA are fast
#'   enough that this is typically a few seconds.
#' @param seed Integer seed for reproducible splits. `NULL` (default) uses the
#'   current RNG state.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"comparability"`. Print it for a per-level
#'   summary; call [autoplot()] on it for a diagnostic plot. The list contains:
#'   \item{coefficients}{Data frame with one row per split x level x factor:
#'     `split`, `level`, `factor` (full-sample `m{k}f{j}` label), `r` (score
#'     comparability), `phi` (Tucker's congruence of the matched loading
#'     columns). `NA` when the level did not converge in one of the halves.}
#'   \item{summary}{Data frame with one row per level x factor: `level`,
#'     `factor`, `r_median`, `r_min`, `phi_median`, `phi_min` (across splits),
#'     and `n_splits_ok` (splits in which both halves converged).}
#'   \item{k_max}{Deepest level evaluated. Can be lower than the `k_max` you
#'     asked for when the full-sample fit truncated (non-convergence at deep
#'     levels); the original request is kept in `k_requested`.}
#'   \item{k_requested, n_splits, n_half, engine, cor, fm, n_obs, n_vars,
#'     seed}{Metadata.}
#'
#' @seealso [suggest_k()] for the plausible depth *range* (eigenstructure),
#'   [factorability()] for sampling adequacy before you fit, [prune()] for
#'   factors that perpetuate without differentiating (redundancy), and
#'   [ackwards()] for the extraction itself.
#'
#' @references
#' Everett, J. E. (1983). Factor comparability as a means of determining the
#'   number of factors and their rotation. *Multivariate Behavioral Research*,
#'   18(2), 197--218. \doi{10.1207/s15327906mbr1802_5}
#'
#' Goldberg, L. R. (1990). An alternative "description of personality": The
#'   Big-Five factor structure. *Journal of Personality and Social
#'   Psychology*, 59(6), 1216--1229. \doi{10.1037/0022-3514.59.6.1216}
#'
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' @examples
#' \donttest{
#' cmp <- comparability(bfi25, k_max = 5, n_splits = 5, seed = 1)
#' cmp
#' cmp$summary
#' }
#'
#' @export
comparability <- function(data, k_max, engine = "pca", cor = "pearson",
                          fm = "minres", n_splits = 10L, seed = NULL, ...) {
  cl <- match.call()
  .check_unknown_dots(list(...), "comparability")

  # Targeted messages for the two deliberately unsupported values before the
  # generic arg_match errors (both are documented scope decisions, not typos).
  if (identical(engine, "esem")) {
    cli::cli_abort(c(
      "!" = "{.code engine = \"esem\"} is not supported by {.fn comparability}.",
      "i" = "Fitting {.code 2 * n_splits} lavaan hierarchies needs its own \\
             performance treatment; use {.code engine = \"efa\"} as a \\
             structural screen for ESEM workflows."
    ))
  }
  if (identical(cor, "polychoric")) {
    cli::cli_abort(c(
      "!" = "{.code cor = \"polychoric\"} is not supported by \\
             {.fn comparability} (as with {.fn suggest_k}).",
      "i" = "Estimating a polychoric matrix in every half-sample is slow and \\
             unstable. Screen replicability on the Pearson basis, then fit \\
             the final model with {.code cor = \"polychoric\"} in \\
             {.fn ackwards}."
    ))
  }
  engine <- rlang::arg_match(engine, c("pca", "efa"))
  cor <- rlang::arg_match(cor, c("pearson", "spearman"))
  fm <- rlang::arg_match(fm, c("minres", "ml", "pa"))

  if (!is.numeric(n_splits) || length(n_splits) != 1L || is.na(n_splits) ||
    n_splits < 1L || n_splits != as.integer(n_splits)) {
    cli::cli_abort("{.arg n_splits} must be a single positive integer.")
  }
  n_splits <- as.integer(n_splits)

  # Raw data only: splitting needs rows.
  if (is.matrix(data)) .check_maybe_cov_matrix(data)
  if (.is_cor_matrix(data)) {
    cli::cli_abort(c(
      "!" = "{.fn comparability} requires raw item data.",
      "i" = "Split-half replication resamples rows, which a correlation \\
             matrix does not carry."
    ))
  }
  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
  }
  data_mat <- as.matrix(data)
  if (!is.numeric(data_mat)) {
    cli::cli_abort("{.arg data} must contain only numeric columns.")
  }

  n <- nrow(data_mat)
  p <- ncol(data_mat)
  n_half <- n %/% 2L

  if (n_half <= p) {
    cli::cli_abort(c(
      "!" = "Each half-sample must have more rows than variables \\
             (n/2 = {n_half}, p = {p}).",
      "i" = "Half-sample correlation matrices would be singular, so \\
             half-solutions cannot be scored."
    ))
  }

  # --- Full-sample anchor -----------------------------------------------------
  # A genuine ackwards() fit, so the reported factor labels are the same
  # m{k}f{j} the user sees when they fit the hierarchy themselves (its ordinal /
  # missing-data advisories also fire here, once -- Invariant 6). Its stored $r
  # is the pooled correlation matrix every cross-solution correlation uses.
  x_full <- ackwards(data_mat, k_max = k_max, engine = engine, cor = cor, fm = fm)
  k_eff <- x_full$k_max
  r_pooled <- x_full$r

  # --- Split-half replicates ---------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  cli::cli_progress_step(
    "Fitting {n_splits} split-half replicate{?s} ({engine}, k = 1-{k_eff})..."
  )

  rows <- vector("list", n_splits)
  for (s in seq_len(n_splits)) {
    perm <- sample.int(n)
    half_a <- .fit_half(data_mat[perm[seq_len(n_half)], , drop = FALSE],
      k_max = k_eff, engine = engine, cor = cor, fm = fm
    )
    half_b <- .fit_half(data_mat[perm[(n_half + 1L):n], , drop = FALSE],
      k_max = k_eff, engine = engine, cor = cor, fm = fm
    )

    rows[[s]] <- do.call(rbind, lapply(seq_len(k_eff), function(k) {
      lev_f <- x_full$levels[[as.character(k)]]
      lev_a <- half_a[[as.character(k)]]
      lev_b <- half_b[[as.character(k)]]
      out <- .level_comparability(lev_f, lev_a, lev_b, r_pooled)
      cbind(split = s, out)
    }))
  }
  cli::cli_progress_done()

  coefficients <- do.call(rbind, rows)
  rownames(coefficients) <- NULL

  # --- Aggregate across splits -------------------------------------------------
  agg_stat <- function(v, f) if (all(is.na(v))) NA_real_ else f(v, na.rm = TRUE)
  key <- interaction(coefficients$level, coefficients$factor, drop = TRUE)
  summary_df <- do.call(rbind, lapply(split(coefficients, key), function(d) {
    data.frame(
      level = d$level[1L],
      factor = d$factor[1L],
      r_median = agg_stat(d$r, stats::median),
      r_min = agg_stat(d$r, min),
      phi_median = agg_stat(d$phi, stats::median),
      phi_min = agg_stat(d$phi, min),
      n_splits_ok = sum(!is.na(d$r)),
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- summary_df[order(summary_df$level, summary_df$factor), , drop = FALSE]
  rownames(summary_df) <- NULL

  # Summarise convergence shortfalls once instead of per-fit warnings
  # (per-half engine warnings are muffled inside .fit_half).
  short <- summary_df[summary_df$n_splits_ok < n_splits, , drop = FALSE]
  if (nrow(short) > 0L) {
    lv <- sort(unique(short$level))
    n_ok <- vapply(lv, function(k) {
      min(short$n_splits_ok[short$level == k])
    }, integer(1L))
    cli::cli_inform(c(
      "!" = "Some half-solutions did not converge: {length(lv)} level{?s} \\
             ({lv}) yielded usable coefficients in as few as \\
             {min(n_ok)}/{n_splits} splits.",
      "i" = "Convergence failure in a half-sample is itself replicability \\
             evidence; see {.code $summary$n_splits_ok}."
    ))
  }

  structure(
    list(
      coefficients = coefficients,
      summary = summary_df,
      k_max = k_eff,
      k_requested = as.integer(k_max),
      n_splits = n_splits,
      n_half = n_half,
      engine = engine,
      cor = cor,
      fm = fm,
      n_obs = n,
      n_vars = p,
      seed = seed,
      call = cl
    ),
    class = "comparability"
  )
}

# Fit levels 1..k_max in one half-sample, returning the engine's levels list.
# Engine warnings (Heywood, per-level truncation) are muffled here: across
# 2 * n_splits fits they would repeat unusably, and their aggregate signal is
# exactly what the coefficients + n_splits_ok report. A half whose correlation
# matrix cannot be factored at all (e.g. a zero-variance column) yields an
# empty list (all-NA split).
.fit_half <- function(data_half, k_max, engine, cor, fm) {
  out <- tryCatch(
    suppressMessages(suppressWarnings({
      R_h <- stats::cor(data_half, method = cor, use = "pairwise.complete.obs")
      switch(engine,
        pca = pca_levels(R_h, k_max = k_max, cor = cor),
        efa = efa_levels(R_h,
          k_max = k_max, fm = fm, n_obs = nrow(data_half), cor = cor
        )
      )
    })),
    error = function(e) list(levels = list())
  )
  out$levels
}

# Cross-solution score-correlation matrix between two same-k solutions,
# evaluated on the pooled correlation matrix R. Routed through compute_edges()
# on a two-element levels list so the one edge path (Invariant 1) computes it:
# rows = lev_a's factors, columns = lev_b's factors.
.cross_cor <- function(lev_a, lev_b, R) {
  compute_edges(
    levels = stats::setNames(list(lev_a, lev_b), c("1", "2")),
    R = R,
    edge_method = "algebra",
    pairs = "adjacent"
  )$matrices[["1:2"]]
}

# Greedy-with-removal bijection on |E|: repeatedly take the globally largest
# remaining |r|, assign that row to that column, remove both. Well-posed here
# (square, same k); ties resolve deterministically (first in column-major
# order), and near-ties signal exactly the instability the coefficients report.
.match_square <- function(E) {
  A <- abs(E)
  k <- nrow(A)
  assignment <- rep(NA_integer_, k)
  for (i in seq_len(k)) {
    mx <- which(A == max(A), arr.ind = TRUE)[1L, , drop = TRUE]
    assignment[[mx[[1L]]]] <- as.integer(mx[[2L]])
    A[mx[[1L]], ] <- -Inf
    A[, mx[[2L]]] <- -Inf
  }
  assignment
}

# One level's comparability rows: anchor both half-solutions to the
# full-sample solution, then correlate the matched half-solution scores.
# Signs: each matched half factor is oriented so its correlation with the
# anchoring full-sample factor is positive; the coefficient is signed after
# that alignment (a negative value would itself be diagnostic, not hidden).
.level_comparability <- function(lev_f, lev_a, lev_b, R) {
  k <- lev_f$k
  na_out <- data.frame(
    level = k,
    factor = lev_f$labels,
    r = NA_real_,
    phi = NA_real_,
    stringsAsFactors = FALSE
  )
  if (is.null(lev_a) || is.null(lev_b)) {
    return(na_out)
  }

  e_fa <- .cross_cor(lev_f, lev_a, R)
  e_fb <- .cross_cor(lev_f, lev_b, R)
  e_ab <- .cross_cor(lev_a, lev_b, R)

  # Pathological pairwise missingness can leave NA cells in the pooled R and
  # hence here; matching is then undefined, so report NA for the level rather
  # than crash mid-run (estimability is data, not an error -- Invariant 7).
  if (anyNA(e_fa) || anyNA(e_fb) || anyNA(e_ab)) {
    return(na_out)
  }

  ja <- .match_square(e_fa)
  jb <- .match_square(e_fb)

  out <- na_out
  for (j in seq_len(k)) {
    sa <- sign(e_fa[j, ja[[j]]])
    sb <- sign(e_fb[j, jb[[j]]])
    out$r[[j]] <- sa * sb * e_ab[ja[[j]], jb[[j]]]
    out$phi[[j]] <- sa * sb *
      .tucker_phi(lev_a$loadings[, ja[[j]]], lev_b$loadings[, jb[[j]]])
  }
  out
}

#' Print a comparability object
#'
#' @param x A `comparability` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.comparability <- function(x, ...) {
  cli::cli_h1("Split-Half Factor Comparability ({.pkg ackwards})")

  cli::cli_dl(c(
    "Engine" = cli::style_bold(x$engine),
    "Basis" = x$cor,
    "n" = paste0(
      format(x$n_obs, big.mark = ","),
      " (", format(x$n_half, big.mark = ","), " per half)"
    ),
    "Splits" = as.character(x$n_splits),
    "Levels" = paste0(
      "1-", x$k_max,
      # k_requested absent on pre-follow-up objects; show the suffix only when
      # the full-sample anchor genuinely truncated below the request.
      if (!is.null(x$k_requested) && x$k_requested > x$k_max) {
        paste0(" (requested 1-", x$k_requested, "; full-sample fit truncated)")
      } else {
        ""
      }
    )
  ))

  cli::cli_h2("Comparability by level (median across splits)")

  sm <- x$summary
  for (k in sort(unique(sm$level))) {
    d <- sm[sm$level == k, , drop = FALSE]
    if (all(is.na(d$r_median))) {
      cli::cli_text(
        "  k = {k}:  ",
        cli::col_grey("no usable splits (half-solutions did not converge)")
      )
      next
    }
    weakest <- which.min(d$r_median)
    line <- paste0(
      "  k = ", k, ":  ",
      "median r ", .format_r(stats::median(d$r_median, na.rm = TRUE)),
      ", min r ", .format_r(d$r_median[weakest]),
      " (", d$factor[weakest], ")"
    )
    if (any(d$n_splits_ok < x$n_splits)) {
      line <- paste0(
        line,
        cli::col_grey(paste0(
          "  [", min(d$n_splits_ok), "/", x$n_splits, " splits usable]"
        ))
      )
    }
    cli::cli_text(line)
  }

  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Per-factor detail (incl. Tucker's \u03c6) in \\
       {.code $summary}; per-split values in {.code $coefficients}."
    )
  )
  cli::cli_text(
    cli::col_grey(
      "Conventional benchmarks: {cli::symbol$geq} .95 comfortable, \\
       {cli::symbol$geq} .90 floor (Everett, 1983; Goldberg, 1990) -- \\
       conventions, not tests. Interpret levels whose factors all replicate."
    )
  )

  invisible(x)
}

#' Plot a comparability diagnostic
#'
#' Renders a two-panel ggplot2 diagnostic for a [comparability()] object:
#' score comparability (r) and loading congruence (Tucker's phi) for every
#' factor at every level. Grey points are individual splits; black points are
#' the per-factor medians. Dashed and dotted reference lines mark the
#' conventional .90 / .95 benchmarks -- visual guides, not tests.
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @param object A `comparability` object.
#' @param ... Ignored.
#' @return A `ggplot` object.
#'
#' @seealso [comparability()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   cmp <- comparability(bfi25, k_max = 5, n_splits = 5, seed = 1)
#'   autoplot(cmp)
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
autoplot.comparability <- function(object, ...) {
  rlang::check_installed("ggplot2", reason = "for autoplot.comparability()")

  # Plotmath strings (rendered via label_parsed): the phi glyph as UTF-8 text
  # fails on the pdf device R CMD check uses for examples (mbcsToSbcs), while
  # plotmath draws it from the symbol font on every device.
  panel_r <- "Score~comparability~(italic(r))"
  panel_phi <- "Loading~congruence~(phi)"

  # Long format: one row per split x level x factor x statistic, plus the
  # per-factor medians as a second point type. Factors within a level are
  # dodged horizontally so they do not overprint.
  cf <- object$coefficients
  sm <- object$summary

  dodge_x <- function(level, fac) {
    j <- as.integer(sub("^m\\d+f", "", fac))
    k <- level
    level + (j - (k + 1) / 2) * (0.6 / max(sm$level))
  }

  long_split <- rbind(
    data.frame(
      level = cf$level, factor = cf$factor, value = cf$r,
      panel = panel_r, type = "Split", stringsAsFactors = FALSE
    ),
    data.frame(
      level = cf$level, factor = cf$factor, value = cf$phi,
      panel = panel_phi, type = "Split", stringsAsFactors = FALSE
    )
  )
  long_med <- rbind(
    data.frame(
      level = sm$level, factor = sm$factor, value = sm$r_median,
      panel = panel_r, type = "Median", stringsAsFactors = FALSE
    ),
    data.frame(
      level = sm$level, factor = sm$factor, value = sm$phi_median,
      panel = panel_phi, type = "Median", stringsAsFactors = FALSE
    )
  )
  plot_data <- rbind(long_split, long_med)
  plot_data <- plot_data[!is.na(plot_data$value), , drop = FALSE]
  plot_data$x <- dodge_x(plot_data$level, plot_data$factor)
  plot_data$panel <- factor(plot_data$panel, levels = c(panel_r, panel_phi))
  plot_data$type <- factor(plot_data$type, levels = c("Split", "Median"))

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$x, y = .data$value, color = .data$type, size = .data$type)
  ) +
    ggplot2::geom_hline(
      yintercept = 0.90, linetype = "dashed", color = "grey60", linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0.95, linetype = "dotted", color = "grey60", linewidth = 0.5
    ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_manual(
      values = c(Split = "grey55", Median = "black"), name = NULL
    ) +
    ggplot2::scale_size_manual(
      values = c(Split = 1.3, Median = 2.6), name = NULL
    ) +
    ggplot2::scale_x_continuous(
      name = "Number of factors / components (level)",
      breaks = seq_len(object$k_max)
    ) +
    ggplot2::facet_wrap(~panel, ncol = 1L, labeller = ggplot2::label_parsed) +
    ggplot2::labs(
      y = "Coefficient",
      caption = paste0(
        "Dashed / dotted lines: conventional .90 / .95 benchmarks ",
        "(Everett, 1983; Goldberg, 1990) -- visual guides, not tests"
      )
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
