# Bootstrap confidence intervals on between-level edges (DESIGN.md s.14
# item 36) -- the inferential-honesty companion to prune()'s point-estimate
# rules and the Forbes "strongest edge" claim.
#
# Mechanics per replicate:
#   1. Resample n rows with replacement (indices drawn upfront from `seed`, so
#      each replicate is deterministic given its indices and the serial and
#      parallel paths agree exactly).
#   2. Recompute the correlation matrix with the same basis + missing-data
#      routine as the original fit, and refit levels 1..k_max (muffled --
#      .fit_half() pattern; a failed replicate yields NAs, never an abort).
#   3. Anchor every replicate level to the full-sample solution (greedy
#      max-|r| matching + sign orientation, evaluated on the full-sample R --
#      the M46 comparability machinery). Without this step, label switching
#      and sign flipping across replicates would corrupt the pooled
#      distributions.
#   4. Compute the replicate's edges through compute_edges() on the replicate
#      R (Invariant 1), with the same `pairs` as the object.
# Percentile CIs and bootstrap SEs are then read off each edge's replicate
# distribution.

#' Bootstrap confidence intervals for between-level edges
#'
#' Attaches nonparametric bootstrap standard errors and percentile confidence
#' intervals to every between-level correlation (edge) of a fitted
#' bass-ackwards hierarchy. Every edge [ackwards()] reports is a point
#' estimate; `boot_edges()` quantifies its sampling uncertainty, which matters
#' most where a hard threshold consumes the estimate -- [prune()]'s
#' `|r| >= redundancy_r` redundancy rule, and the Forbes (2023) practice of
#' interpreting the strongest all-pairs edge.
#'
#' For each of `n_boot` replicates, `n` rows are resampled with replacement,
#' the correlation matrix is recomputed with the same basis and missing-data
#' routine used at fit time, and the full hierarchy is refit. Each replicate
#' level is then **anchored to the full-sample solution** -- its factors
#' matched (greedy max-|r| with removal) and sign-oriented against the
#' full-sample factors on the full-sample correlation matrix -- before its
#' edges are computed. Without anchoring, factor label switching and sign
#' flipping across replicates would corrupt the pooled edge distributions;
#' this is the same matching machinery [comparability()] uses.
#'
#' All resample indices are drawn upfront from `seed`, so results are
#' reproducible and identical whether replicates run serially or in parallel.
#' Replicate fits are dispatched through \pkg{future.apply} when it is
#' installed and the user has set a [future::plan()] (serial otherwise, as in
#' [ackwards()]'s ESEM engine).
#'
#' @section What the intervals do and do not fix:
#' Per-edge intervals make sampling uncertainty **visible**: an edge whose
#' interval straddles [prune()]'s `redundancy_r` threshold should not be
#' treated as decisively above or below it. They do **not** correct the
#' selection bias of scanning many edges for the strongest one -- the maximum
#' of hundreds of correlations capitalizes on chance even when every
#' individual interval is honest. Treat the intervals as per-edge error bars,
#' not a familywise inference.
#'
#' @section Failed replicates:
#' A replicate whose hierarchy fails to converge (in full or at some levels)
#' contributes `NA` to the affected edges and is dropped from their
#' distributions -- convergence is data, not an error. The usable replicate
#' count is reported per edge in `n_boot_ok`, and a message summarises any
#' shortfall.
#'
#' @param x An `ackwards` object fit with `engine = "pca"` or `"efa"` on a
#'   `"pearson"` or `"spearman"` basis. ESEM objects are not supported
#'   (refitting `n_boot` lavaan hierarchies needs its own performance
#'   treatment), nor are polychoric-basis objects (estimating a polychoric
#'   matrix in every resample is slow and unstable) or objects fit from a
#'   correlation matrix (resampling needs rows).
#' @param data The raw item data the model was fit on. Required -- the
#'   `ackwards` object deliberately does not store raw data (light core), so
#'   it must be re-supplied here. Columns are matched by name against the
#'   fit; a warning is issued if the data do not look like the fit data.
#' @param n_boot Number of bootstrap replicates. Default `1000L`; larger
#'   values (2000+) are advisable for published interval endpoints (Efron &
#'   Tibshirani, 1993). Each replicate refits the full hierarchy, so cost
#'   scales linearly.
#' @param conf Confidence level for the percentile intervals. Default `0.95`.
#' @param seed Integer seed for reproducible resampling. `NULL` (default)
#'   uses the current RNG state.
#' @param ... Reserved for future arguments.
#'
#' @return `x`, invisibly modified: the `$boot` element is populated with
#'   \item{edges}{Data frame with one row per edge (aligned with
#'     `tidy(x, what = "edges")`): `from`, `to`, `level_from`, `level_to`,
#'     `r` (the full-sample point estimate), `se` (bootstrap standard error),
#'     `lo`, `hi` (percentile interval endpoints), and `n_boot_ok` (usable
#'     replicates for that edge).}
#'   \item{n_boot, conf, seed}{The request.}
#'   After calling `boot_edges()`, `tidy(x, what = "edges")` gains `se`,
#'   `lo`, and `hi` columns, and `print(x)`/`summary(x)` note the interval
#'   coverage.
#'
#' @seealso [prune()] for the thresholded rules the intervals contextualise,
#'   [comparability()] for split-half replicability of the factors
#'   themselves, [tidy.ackwards()] for the augmented edge table.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). *An introduction to the bootstrap*.
#'   Chapman & Hall.
#'
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' @examples
#' \donttest{
#' x <- ackwards(sim16, k_max = 3)
#' x <- boot_edges(x, sim16, n_boot = 100, seed = 1)
#' x$boot$edges
#' head(tidy(x)) # now carries se / lo / hi
#' }
#'
#' @export
boot_edges <- function(x, ...) {
  UseMethod("boot_edges")
}

#' @rdname boot_edges
#' @export
boot_edges.ackwards <- function(x, data, n_boot = 1000L, conf = 0.95,
                                seed = NULL, ...) {
  .check_unknown_dots(list(...), "boot_edges")

  # --- Scope guards (documented decisions, not typos) -------------------------
  if (identical(x$meta$input_type, "cor_matrix")) {
    cli::cli_abort(c(
      "!" = "{.fn boot_edges} requires a model fit on raw item data.",
      "i" = "Bootstrap resampling draws rows, which a correlation matrix \\
             does not carry."
    ))
  }
  if (identical(x$engine, "esem")) {
    cli::cli_abort(c(
      "!" = "{.code engine = \"esem\"} objects are not supported by \\
             {.fn boot_edges}.",
      "i" = "Refitting {.code n_boot} lavaan hierarchies needs its own \\
             performance treatment; fit an {.code engine = \"efa\"} \\
             counterpart to bootstrap its edges."
    ))
  }
  if (identical(x$cor, "polychoric")) {
    cli::cli_abort(c(
      "!" = "Polychoric-basis objects are not supported by {.fn boot_edges} \\
             (as with {.fn comparability}).",
      "i" = "Estimating a polychoric matrix in every resample is slow and \\
             unstable. Bootstrap a Pearson-basis fit of the same data \\
             instead."
    ))
  }

  n_boot <- .check_count(n_boot, "n_boot", min = 2L)
  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) ||
    conf <= 0 || conf >= 1) {
    cli::cli_abort("{.arg conf} must be a single number in (0, 1).")
  }

  # --- Data validation ---------------------------------------------------------
  if (missing(data) || is.null(data)) {
    cli::cli_abort(c(
      "!" = "{.arg data} is required.",
      "i" = "The {.cls ackwards} object does not store raw data; re-supply \\
             the data the model was fit on."
    ))
  }
  data_mat <- .as_numeric_matrix(data)

  # Match columns to the fit (same contract as augment()/predict()).
  W_ref <- x$levels[[1L]]$scoring$weights
  vars_expected <- rownames(W_ref)
  if (!is.null(vars_expected) && !is.null(colnames(data_mat))) {
    missing_vars <- setdiff(vars_expected, colnames(data_mat))
    if (length(missing_vars) > 0L) {
      cli::cli_abort(c(
        "!" = "{.arg data} is missing {length(missing_vars)} variable{?s} \\
               that the model was fit on.",
        "x" = "Missing: {.val {missing_vars}}"
      ))
    }
    data_mat <- data_mat[, vars_expected, drop = FALSE]
  } else if (ncol(data_mat) != nrow(W_ref)) {
    cli::cli_abort(c(
      "!" = "{.arg data} has {ncol(data_mat)} column{?s} but the model \\
             was fit on {nrow(W_ref)}.",
      "i" = "Supply the data the model was fit on (or a matrix with \\
             matching column names)."
    ))
  }

  # Reproduce the fit-time missing-data preparation, then check the supplied
  # data actually look like the fit data (row count + item means). The CIs
  # are only meaningful for the data the point estimates came from.
  missing_eff <- x$meta$missing
  if (identical(missing_eff, "listwise")) {
    data_mat <- data_mat[stats::complete.cases(data_mat), , drop = FALSE]
  }
  n <- nrow(data_mat)
  mu_fit <- x$meta$item_means
  looks_wrong <- !identical(n, as.integer(x$n_obs)) ||
    (!is.null(mu_fit) &&
      max(abs(colMeans(data_mat, na.rm = TRUE) - unname(mu_fit))) > 1e-6)
  if (looks_wrong) {
    cli::cli_warn(c(
      "!" = "{.arg data} does not look like the data the model was fit on \\
             (row count or item means differ).",
      "i" = "Bootstrap intervals describe the sampling uncertainty of the \\
             fitted estimates only when computed from the same data."
    ))
  }

  # --- Upfront resample indices --------------------------------------------
  # All randomness happens here: each replicate is then deterministic given
  # its index vector, so serial and parallel dispatch agree exactly.
  if (!is.null(seed)) set.seed(seed)
  idx_list <- lapply(seq_len(n_boot), function(b) {
    sample.int(n, n, replace = TRUE)
  })

  if (identical(missing_eff, "fiml")) {
    cli::cli_inform(c(
      "i" = "{.code missing = \"fiml\"} fit: each replicate re-estimates the \\
             correlation matrix via {.fn psych::corFiml}, which is slow -- \\
             expect {n_boot} FIML estimations."
    ))
  }

  # Canonical edge order = the object's tidy edge rows (row-major within each
  # matrices key -- the order compute_edges() builds them in).
  keys <- names(x$edges$matrices)
  dims <- lapply(x$edges$matrices, dim)

  cli::cli_progress_step(
    "Fitting {n_boot} bootstrap replicate{?s} ({x$engine}, \\
     k = 1-{x$k_max})..."
  )
  rep_rows <- .boot_lapply(idx_list, function(idx) {
    .boot_replicate(
      idx, data_mat,
      x_levels = x$levels, R_full = x$r, keys = keys, dims = dims,
      engine = x$engine, cor = x$cor, fm = x$meta$fm %||% "minres",
      missing_eff = missing_eff, k_max = x$k_max,
      pairs = x$meta$pairs
    )
  })
  cli::cli_progress_done()

  boot_mat <- do.call(rbind, rep_rows) # n_boot x n_edges

  # --- Percentile CIs + SEs ----------------------------------------------------
  alpha <- (1 - conf) / 2
  n_ok <- colSums(!is.na(boot_mat))
  stat <- function(v, f) if (sum(!is.na(v)) < 2L) NA_real_ else f(v)
  se <- apply(boot_mat, 2L, stat, f = function(v) stats::sd(v, na.rm = TRUE))
  lo <- apply(boot_mat, 2L, stat, f = function(v) {
    stats::quantile(v, alpha, na.rm = TRUE, names = FALSE)
  })
  hi <- apply(boot_mat, 2L, stat, f = function(v) {
    stats::quantile(v, 1 - alpha, na.rm = TRUE, names = FALSE)
  })

  tidy_edges <- x$edges$tidy
  boot_df <- data.frame(
    from = tidy_edges$from,
    to = tidy_edges$to,
    level_from = tidy_edges$level_from,
    level_to = tidy_edges$level_to,
    r = tidy_edges$r,
    se = se,
    lo = lo,
    hi = hi,
    n_boot_ok = as.integer(n_ok),
    stringsAsFactors = FALSE
  )

  if (any(boot_df$n_boot_ok < n_boot)) {
    cli::cli_inform(c(
      "!" = "Some replicates did not converge: usable replicates per edge \\
             range from {min(boot_df$n_boot_ok)} to \\
             {max(boot_df$n_boot_ok)} of {n_boot}.",
      "i" = "Convergence failure in a resample is itself stability \\
             evidence; see {.code $boot$edges$n_boot_ok}."
    ))
  }

  x$boot <- list(
    edges = boot_df,
    n_boot = n_boot,
    conf = conf,
    seed = seed
  )
  x
}

# Parallel dispatch: future.apply when installed (runs under the user's
# future::plan(); sequential by default), serial lapply otherwise -- the M26
# pattern. future.seed = TRUE silences future's RNG advisory; the replicates
# themselves are deterministic given the precomputed indices.
.boot_lapply <- function(X, FUN) {
  if (rlang::is_installed("future.apply")) {
    future.apply::future_lapply(X, FUN, future.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

# One bootstrap replicate: resample -> recompute R -> refit -> anchor to the
# full-sample solution -> edges on the replicate R. Returns a numeric vector
# aligned with the object's tidy edge rows (row-major within each matrices
# key); edges touching an unusable level are NA.
.boot_replicate <- function(idx, data_mat, x_levels, R_full, keys, dims,
                            engine, cor, fm, missing_eff, k_max, pairs) {
  na_out <- rep(NA_real_, sum(vapply(dims, prod, numeric(1L))))

  levels_rep <- tryCatch(
    suppressMessages(suppressWarnings({
      d_b <- data_mat[idx, , drop = FALSE]
      R_b <- if (identical(missing_eff, "fiml")) {
        .corfiml_R(d_b)
      } else {
        stats::cor(d_b, method = cor, use = "pairwise.complete.obs")
      }
      out <- switch(engine,
        pca = pca_levels(R_b, k_max = k_max, cor = cor),
        efa = efa_levels(R_b,
          k_max = k_max, fm = fm, n_obs = nrow(d_b), cor = cor
        )
      )
      list(levels = out$levels, R_b = R_b)
    })),
    error = function(e) NULL
  )
  if (is.null(levels_rep) || length(levels_rep$levels) < 2L) {
    return(na_out)
  }

  anchored <- .anchor_levels(x_levels, levels_rep$levels, R_full)
  usable <- !vapply(anchored, is.null, logical(1L))
  if (sum(usable) < 2L) {
    return(na_out)
  }

  E_rep <- suppressWarnings(compute_edges(
    levels = anchored[usable],
    R = levels_rep$R_b,
    edge_method = "algebra",
    pairs = pairs,
    build_tidy = FALSE
  )$matrices)

  # Flatten in the canonical key order, row-major within each matrix (the
  # order compute_edges() builds tidy rows in); NA-fill keys the replicate
  # could not produce (failed / unmatched levels).
  unlist(lapply(seq_along(keys), function(i) {
    E <- E_rep[[keys[[i]]]]
    if (is.null(E)) rep(NA_real_, prod(dims[[i]])) else as.vector(t(E))
  }), use.names = FALSE)
}

# Anchor replicate levels to the full-sample solution: match each full-sample
# factor to its best replicate counterpart (greedy max-|r| with removal on
# the cross-solution correlations, evaluated on the full-sample R -- the M46
# comparability machinery), then permute and sign-orient the replicate's
# weights/loadings into the full-sample factor order. Returns a list indexed
# like x_levels; NULL entries mark levels the replicate could not anchor
# (missing, or NA cross-correlations under pathological missingness).
.anchor_levels <- function(x_levels, levels_rep, R_full) {
  out <- stats::setNames(
    vector("list", length(x_levels)), names(x_levels)
  )
  for (ki in names(x_levels)) {
    lev_f <- x_levels[[ki]]
    lev_b <- levels_rep[[ki]]
    if (is.null(lev_b)) next
    e_fb <- .cross_cor(lev_f, lev_b, R_full)
    if (anyNA(e_fb)) next
    jb <- .match_square(e_fb)
    s <- sign(e_fb[cbind(seq_along(jb), jb)])
    s[s == 0] <- 1 # an exactly-zero anchor correlation carries no orientation

    anch <- lev_b
    anch$labels <- lev_f$labels
    anch$loadings <- sweep(
      lev_b$loadings[, jb, drop = FALSE], 2L, s, "*"
    )
    colnames(anch$loadings) <- lev_f$labels
    anch$scoring$weights <- sweep(
      lev_b$scoring$weights[, jb, drop = FALSE], 2L, s, "*"
    )
    colnames(anch$scoring$weights) <- lev_f$labels
    anch$scoring$score_var <- lev_b$scoring$score_var[jb]
    out[[ki]] <- anch
  }
  out
}
