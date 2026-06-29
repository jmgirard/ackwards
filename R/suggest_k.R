#' Suggest a maximum number of factors for bass-ackwards analysis
#'
#' Runs five complementary selection criteria and reports their recommendations.
#' No single criterion is definitive; the goal is a consensus range to inform
#' your choice of `k` in [ackwards()].
#'
#' **Criteria computed:**
#' * **PA-PC** (Horn 1965, PC basis) -- parallel analysis on principal-component
#'   eigenvalues. Compares observed eigenvalues to those from random correlation
#'   matrices; suggests retaining components whose eigenvalues exceed the 95th
#'   percentile of chance. Tends to overextract; treat as an upper bound.
#' * **PA-FA** (Horn 1965, FA basis) -- parallel analysis using common-factor
#'   eigenvalues. More conservative than PA-PC and the better match for the EFA
#'   and ESEM engines in [ackwards()].
#' * **MAP** (Velicer 1976) -- Minimum Average Partial criterion. Finds the k
#'   that minimises the average squared partial correlation remaining after
#'   extracting k components. Usually conservative.
#' * **VSS-1 / VSS-2** (Revelle & Rocklin 1979) -- Very Simple Structure fit
#'   at complexities 1 and 2. Finds the k maximising the fit of a very simple
#'   loading structure.
#' * **CD** (Ruscio & Roche 2012; optional) -- Comparison Data. Resamples
#'   from the observed item distributions to generate comparison eigenvalue
#'   profiles; retains factors until adding one no longer improves RMSE beyond
#'   chance. Requires the \pkg{EFAtools} package (install separately).
#'
#' PA (both bases), MAP, and VSS share the same correlation matrix as
#' [ackwards()]. CD operates on the raw data matrix directly (required for
#' resampling) and is skipped gracefully when \pkg{EFAtools} is not installed.
#'
#' @section Interpreting the output:
#' `k` in [ackwards()] is a **maximum depth**, not a claim that exactly k
#' factors exist. Users commonly set k one or two levels above the consensus to
#' watch higher-level factors fragment -- this is a feature of the method, not
#' overextraction.
#'
#' @param data A data frame or numeric matrix (items in columns, observations in
#'   rows). Alternatively, a pre-computed **correlation matrix** may be supplied
#'   (a square, symmetric, numeric matrix with unit diagonal). When a
#'   correlation matrix is supplied, `n_obs` is required (PA and VSS need N),
#'   the `cor` argument is ignored, and the Comparison Data (CD) criterion is
#'   skipped (CD requires raw item distributions for resampling).
#' @param k_max Maximum number of components to test. Defaults to
#'   `min(ncol(data) - 1, 8)`. Increase if you expect a deeper hierarchy.
#' @param cor Correlation basis: `"pearson"` (default) or `"spearman"`. Should
#'   match the `cor` argument you plan to use in [ackwards()]. Ignored when
#'   `data` is a correlation matrix (the basis is already fixed).
#' @param n_obs Number of observations. Required when `data` is a pre-computed
#'   correlation matrix (PA and VSS need N). Ignored when raw data are supplied
#'   (N is determined from `nrow(data)`).
#' @param n_iter Number of Monte Carlo iterations for parallel analysis. Default
#'   `20`. Reduce to `5` for fast/exploratory runs; increase to `100+` for
#'   publication.
#' @param seed Integer seed passed to [set.seed()] before the Comparison Data
#'   (CD) step. `NULL` (default) uses the current RNG state. **Note:** the
#'   parallel-analysis step uses `psych::fa.parallel()`, which does not respond
#'   reliably to `set.seed()` -- PA simulation results will vary across calls
#'   regardless of `seed`.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"suggest_k"`. Print it for a formatted summary;
#'   call [autoplot()] on it for a diagnostic scree/criteria plot. The list
#'   contains:
#'   \item{k_parallel_pc}{Recommended k from PC-based parallel analysis.}
#'   \item{k_parallel_fa}{Recommended k from FA-based parallel analysis
#'     (`NA_integer_` if no FA factor exceeded the random threshold).}
#'   \item{k_map}{Recommended k from MAP.}
#'   \item{k_vss1}{Recommended k from VSS complexity-1.}
#'   \item{k_vss2}{Recommended k from VSS complexity-2.}
#'   \item{k_cd}{Recommended k from Comparison Data (`NA_integer_` when
#'     \pkg{EFAtools} is not installed or CD fails).}
#'   \item{cd_available}{Logical; `TRUE` when \pkg{EFAtools} was found and CD
#'     ran successfully.}
#'   \item{criteria}{Data frame with one row per k: `k`, `ev_obs` (observed PC
#'     eigenvalue), `ev_obs_fa` (observed FA eigenvalue),
#'     `pa_pc_quant` / `pa_fa_quant` (95th-pct simulated eigenvalue for each
#'     basis), `pa_pc_suggested` / `pa_fa_suggested` (logical retention),
#'     `map`, `vss1`, `vss2`.}
#'   \item{k_max, n_obs, n_vars, cor}{Metadata.}
#'
#' @section A note on overextraction:
#' PA-PC in particular tends to recommend more factors than replicate across
#' independent samples, especially with correlated items (Forbes, 2023).
#' PA-FA and CD are more conservative. Treat the full set of criteria as a
#' range: the true k is likely somewhere in the middle.
#'
#' @seealso [ackwards()]
#'
#' @references
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' Horn, J. L. (1965). A rationale and test for the number of factors in factor
#'   analysis. *Psychometrika*, 30, 179--185.
#'
#' Revelle, W., & Rocklin, T. (1979). Very simple structure: An alternative
#'   procedure for estimating the optimal number of interpretable factors.
#'   *Multivariate Behavioral Research*, 14(4), 403--414.
#'
#' Ruscio, J., & Roche, B. (2012). Determining the number of factors to retain
#'   in an exploratory factor analysis using comparison data of a known factorial
#'   structure. *Psychological Assessment*, 24(2), 282--292.
#'
#' Velicer, W. F. (1976). Determining the number of components from the matrix
#'   of partial correlations. *Psychometrika*, 41, 321--327.
#'
#' @examples
#' \donttest{
#' sk <- suggest_k(bfi25)
#' sk
#' autoplot(sk)
#'
#' # Faster exploratory run
#' suggest_k(bfi25, k_max = 6, n_iter = 5)
#'
#' # Correlation-matrix input (CD is skipped; n_obs required)
#' R <- cor(bfi25, use = "pairwise.complete.obs")
#' suggest_k(R, n_obs = 875L)
#' }
#'
#' @export
suggest_k <- function(data, k_max = NULL, cor = "pearson", n_obs = NULL,
                      n_iter = 20L, seed = NULL, ...) {
  # --- Detect R-matrix vs. raw-data input -------------------------------------
  if (is.matrix(data)) .check_maybe_cov_matrix(data)
  input_type <- if (.is_cor_matrix(data)) "cor_matrix" else "data"

  n_iter <- as.integer(n_iter)

  if (input_type == "cor_matrix") {
    # Validate and normalise; synthesise V1..Vp dimnames if absent
    R <- .validate_cor_matrix(as.matrix(data))
    p <- nrow(R)
    n_vars <- p

    # Warn if cor set explicitly (ignored for R input)
    cor_was_set <- !missing(cor) && cor != "pearson"
    if (cor_was_set) {
      cli::cli_warn(
        c(
          "!" = "{.arg cor} is ignored when a correlation matrix is supplied.",
          "i" = "Remove {.arg cor} from your call to suppress this warning."
        ),
        .frequency = "once",
        .frequency_id = "suggest_k_cor_ignored"
      )
    }
    # Warn if n_obs not supplied (also caught below for required check)
    if (is.null(n_obs)) {
      cli::cli_abort(c(
        "!" = "{.arg n_obs} is required when a correlation matrix is supplied.",
        "i" = "Parallel analysis and VSS need the number of observations.",
        "i" = "Supply the number of observations: {.code n_obs = <N>}."
      ))
    }
    if (!is.numeric(n_obs) || length(n_obs) != 1L ||
      n_obs < 1L || n_obs != as.integer(n_obs)) {
      cli::cli_abort("{.arg n_obs} must be a positive integer.")
    }
    n <- as.integer(n_obs)

    if (is.null(k_max)) k_max <- min(p - 1L, 8L)
    k_max <- as.integer(k_max)

    if (k_max < 1L || k_max >= p) {
      cli::cli_abort(
        "{.arg k_max} must be between 1 and {p - 1L} (number of variables - 1)."
      )
    }

    # CD requires raw data for resampling -- gate it off with an info note
    cd_available <- FALSE
    cli::cli_inform(c(
      "i" = "Comparison Data (CD) is skipped when a correlation matrix is \\
             supplied (CD requires raw item distributions for resampling)."
    ))
    cor_stored <- NA_character_
    data_mat <- NULL # no raw data
  } else {
    # --- Raw data path --------------------------------------------------------
    cor <- rlang::arg_match(cor, c("pearson", "spearman"))

    if (!is.null(n_obs)) {
      cli::cli_warn(
        c(
          "!" = "{.arg n_obs} is ignored when raw data are supplied \\
               (N is determined from {.code nrow(data)}).",
          "i" = "Remove {.arg n_obs} from your call to suppress this warning."
        ),
        .frequency = "once",
        .frequency_id = "suggest_k_nobs_ignored"
      )
    }

    if (!is.data.frame(data) && !is.matrix(data)) {
      cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
    }
    data_mat <- as.matrix(data)
    if (!is.numeric(data_mat)) {
      cli::cli_abort("{.arg data} must contain only numeric columns.")
    }

    p <- ncol(data_mat)
    n_vars <- p
    n <- nrow(data_mat)
    cor_stored <- cor

    if (is.null(k_max)) k_max <- min(p - 1L, 8L)
    k_max <- as.integer(k_max)

    if (k_max < 1L || k_max >= p) {
      cli::cli_abort(
        "{.arg k_max} must be between 1 and {p - 1L} (number of variables - 1)."
      )
    }

    R <- stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
    cd_available <- rlang::is_installed("EFAtools") # will be checked later
  }

  # --- Parallel analysis: PC + FA (Horn) -----------------------------------------
  cli::cli_progress_step(
    "Running parallel analysis ({n_iter} iterations, PC + FA)..."
  )
  invisible(utils::capture.output(
    pa <- psych::fa.parallel(
      R,
      n.obs  = n,
      fa     = "both",
      n.iter = n_iter,
      plot   = FALSE,
      quant  = 0.95
    )
  ))
  k_parallel_pc <- min(pa$ncomp, k_max)

  # PA-FA may legitimately return 0 when no observed FA eigenvalue exceeds the
  # random-data threshold. Use NA_integer_ so downstream code can detect and
  # handle the "undetermined" state cleanly rather than emitting "k <= 0".
  nfact_raw <- pa$nfact
  k_parallel_fa <- if (is.null(nfact_raw) || nfact_raw == 0L) {
    NA_integer_
  } else {
    min(as.integer(nfact_raw), k_max)
  }

  # Observed PC and FA eigenvalues; simulated thresholds at the 95th pct.
  # psych uses $pc.values / $pc.sim in newer versions with fa = "both";
  # fall back to $values / $sim (backwards-compat aliases) if absent.
  ev_obs <- (pa$pc.values %||% pa$values)[seq_len(k_max)]
  ev_obs_fa <- pa$fa.values[seq_len(k_max)]
  pa_pc_quant <- (pa$pc.sim %||% pa$sim)[seq_len(k_max)]
  pa_fa_quant <- pa$fa.sim[seq_len(k_max)]

  # --- MAP + VSS (Velicer; Revelle & Rocklin) ---------------------------------
  cli::cli_progress_step("Running MAP and VSS...")
  invisible(utils::capture.output(
    vss_out <- psych::vss(
      R,
      n      = k_max,
      n.obs  = n,
      rotate = "varimax",
      fm     = "pc",
      plot   = FALSE
    )
  ))
  map_vals <- vss_out$map[seq_len(k_max)]
  k_map <- which.min(map_vals)
  vss1_vals <- vss_out$vss.stats$cfit.1[seq_len(k_max)]
  vss2_vals <- vss_out$vss.stats$cfit.2[seq_len(k_max)]
  k_vss1 <- which.max(vss1_vals)
  k_vss2 <- which.max(vss2_vals)

  # --- Comparison Data (Ruscio & Roche 2012) -- optional ---------------------
  # cd_available is FALSE when R-matrix input (gated above); otherwise check
  # whether EFAtools is installed.
  if (input_type == "data") {
    cd_available <- rlang::is_installed("EFAtools")
  }
  k_cd <- NA_integer_
  cd_rmse <- NULL
  if (cd_available) {
    # Warn early when the correlation basis won't carry over to CD.
    if (cor == "spearman") {
      cli::cli_warn(
        c(
          "!" = "CD ({.pkg EFAtools}) always uses Pearson correlations.",
          "i" = "The {.code cor = \"spearman\"} basis applies to MAP/VSS/PA \\
                 but not to CD. CD and the other criteria may diverge."
        )
      )
    }

    # CD needs complete cases for resampling; apply listwise deletion.
    data_complete <- data_mat[stats::complete.cases(data_mat), , drop = FALSE]
    n_dropped <- nrow(data_mat) - nrow(data_complete)
    if (n_dropped > 0L) {
      cli::cli_inform(
        "CD: {n_dropped} row{?s} with missing values removed \\
         ({nrow(data_complete)} complete cases used)."
      )
    }

    cli::cli_progress_step("Running Comparison Data (CD)...")
    if (!is.null(seed)) set.seed(seed)
    cd_out <- tryCatch(
      EFAtools::CD(data_complete, n_factors_max = k_max),
      error = function(e) {
        cli::cli_warn(
          c(
            "!" = "Comparison Data (CD) failed and will be omitted.",
            "i" = "{conditionMessage(e)}"
          )
        )
        NULL
      }
    )
    if (!is.null(cd_out)) {
      k_cd <- min(as.integer(cd_out$n_factors), k_max)
      cd_rmse <- colMeans(cd_out$RMSE_eigenvalues)[seq_len(k_max)]
    } else {
      cd_available <- FALSE
    }
  }

  cli::cli_progress_done()

  # --- Build criteria table ---------------------------------------------------
  criteria <- data.frame(
    k = seq_len(k_max),
    ev_obs = ev_obs,
    ev_obs_fa = ev_obs_fa,
    pa_pc_quant = pa_pc_quant,
    pa_pc_suggested = seq_len(k_max) <= k_parallel_pc,
    pa_fa_quant = pa_fa_quant,
    pa_fa_suggested = if (is.na(k_parallel_fa)) {
      rep(FALSE, k_max)
    } else {
      seq_len(k_max) <= k_parallel_fa
    },
    map = map_vals,
    vss1 = vss1_vals,
    vss2 = vss2_vals,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      k_parallel_pc = k_parallel_pc,
      k_parallel_fa = k_parallel_fa,
      k_map         = k_map,
      k_vss1        = k_vss1,
      k_vss2        = k_vss2,
      k_cd          = k_cd,
      cd_available  = cd_available,
      cd_rmse       = cd_rmse,
      criteria      = criteria,
      k_max         = k_max,
      n_obs         = n,
      n_vars        = n_vars,
      cor           = cor_stored,
      input_type    = input_type
    ),
    class = "suggest_k"
  )
}

#' Print a suggest_k object
#'
#' @param x A `suggest_k` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.suggest_k <- function(x, ...) {
  cli::cli_h1("Factor / Component Count Suggestion ({.pkg ackwards})")

  cor_label <- if (is.na(x$cor)) "(user-supplied matrix)" else x$cor
  cli::cli_dl(c(
    "Variables" = as.character(x$n_vars),
    "n"         = format(x$n_obs, big.mark = ","),
    "Basis"     = cor_label,
    "Tested k"  = paste0("1-", x$k_max)
  ))

  cli::cli_h2("Criteria (k = 1-{x$k_max})")

  cr <- x$criteria
  tick <- cli::col_green(cli::symbol$tick)
  dash <- cli::col_grey("-")
  star <- cli::col_cyan("*")

  map_fmt <- formatC(cr$map, digits = 4, format = "f")
  vss1_fmt <- formatC(cr$vss1, digits = 4, format = "f")
  vss2_fmt <- formatC(cr$vss2, digits = 4, format = "f")

  map_fmt[x$k_map] <- paste0(map_fmt[x$k_map], star)
  vss1_fmt[x$k_vss1] <- paste0(vss1_fmt[x$k_vss1], star)
  vss2_fmt[x$k_vss2] <- paste0(vss2_fmt[x$k_vss2], star)

  for (i in seq_len(nrow(cr))) {
    pc_sym <- if (cr$pa_pc_suggested[i]) tick else dash
    fa_sym <- if (cr$pa_fa_suggested[i]) tick else dash

    if (!x$cd_available) {
      cd_part <- ""
    } else if (!is.na(x$k_cd) && i == x$k_cd) {
      cd_part <- paste0("  CD ", tick, star)
    } else if (!is.na(x$k_cd) && i < x$k_cd) {
      cd_part <- paste0("  CD ", tick)
    } else {
      cd_part <- paste0("  CD ", dash)
    }

    cli::cli_text(
      "  k = {cr$k[i]}:  PA-PC {pc_sym}  PA-FA {fa_sym}  \\
       MAP {map_fmt[i]}  VSS-1 {vss1_fmt[i]}  VSS-2 {vss2_fmt[i]}{cd_part}"
    )
  }

  if (!x$cd_available) {
    if (!is.null(x$input_type) && identical(x$input_type, "cor_matrix")) {
      cli::cli_text(
        cli::col_grey(
          "  + CD skipped (requires raw data; not available for matrix input)."
        )
      )
    } else {
      cli::cli_text(
        cli::col_grey(
          "  + CD requires {.pkg EFAtools} (install to enable)."
        )
      )
    }
  }

  cli::cli_h2("Recommendations")

  # Build recommendations; PA-FA is omitted when undetermined (NA).
  recs <- c("PA-PC" = paste0("k <= ", x$k_parallel_pc))
  if (!is.na(x$k_parallel_fa)) {
    recs["PA-FA"] <- paste0("k <= ", x$k_parallel_fa)
  } else {
    recs["PA-FA"] <- "undetermined (no FA factor exceeded random threshold)"
  }
  recs["MAP"] <- paste0("k = ", x$k_map)
  recs["VSS-1"] <- paste0("k = ", x$k_vss1)
  recs["VSS-2"] <- paste0("k = ", x$k_vss2)
  if (x$cd_available && !is.na(x$k_cd)) {
    recs["CD"] <- paste0("k = ", x$k_cd)
  }
  for (nm in names(recs)) {
    cli::cli_bullets(c("*" = paste0(nm, ": ", recs[[nm]])))
  }

  # Consensus uses only the numeric k values; NA PA-FA and absent CD excluded.
  all_k <- stats::na.omit(
    c(
      x$k_parallel_pc, x$k_parallel_fa, x$k_map, x$k_vss1, x$k_vss2,
      if (x$cd_available) x$k_cd
    )
  )
  lo <- min(all_k)
  hi <- max(all_k)

  if (lo == hi) {
    cli::cli_text("{.strong Consensus: k = {lo}}")
  } else {
    cli::cli_text("{.strong Consensus range: k = {lo}-{hi}}")
  }

  cli::cli_rule()
  cli::cli_text(
    cli::col_grey(
      "Note: k_max in ackwards() is a maximum depth. Setting k_max one or two \\
       levels above the consensus to observe factor fragmentation is intentional."
    )
  )
  cli::cli_text(
    cli::col_grey(
      "Caution: PA-PC tends to overextract; structures may not replicate \\
       (Forbes, 2023). PA-FA and CD are more conservative. Use the range."
    )
  )

  invisible(x)
}

#' Plot a suggest_k diagnostic
#'
#' Renders a ggplot2 diagnostic for a `suggest_k` object. When
#' \pkg{EFAtools} is installed and CD was computed, the plot is a 2x2 grid:
#' Scree/PA (top-left), MAP (top-right), VSS (bottom-left), and CD RMSE
#' (bottom-right). When CD is unavailable, the plot is a single-column
#' three-panel layout (Scree/PA, MAP, VSS). The recommended k for each
#' criterion is marked with a star-shaped point.
#'
#' The scree panel shows both PC and FA observed eigenvalues alongside their
#' respective random-data thresholds. PA-PC compares the blue "Observed (PC)"
#' line to the dashed PA-PC threshold; PA-FA compares the teal "Observed (FA)"
#' line to the dotted PA-FA threshold. Reading the lines on the same panel is
#' informative but the two comparisons are independent.
#'
#' The CD panel plots the mean RMSE between observed and comparison-data
#' eigenvalues at each k; the retention threshold (star) is where this curve
#' first crosses below the comparison-data average.
#'
#' Requires the \pkg{ggplot2} package.
#'
#' @param object A `suggest_k` object.
#' @param ... Ignored.
#' @return A `ggplot` object.
#'
#' @seealso [suggest_k()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sk <- suggest_k(bfi25, n_iter = 5)
#'   autoplot(sk)
#' }
#' }
#'
#' @importFrom rlang .data
#' @export
autoplot.suggest_k <- function(object, ...) {
  rlang::check_installed("ggplot2", reason = "for autoplot.suggest_k()")

  cr <- object$criteria
  k_max <- object$k_max
  show_cd <- isTRUE(object$cd_available) && !is.null(object$cd_rmse)

  # --- Long-format data -------------------------------------------------------
  # Scree panel: 4 series so each PA comparison is shown correctly.
  # PA-PC: Observed (PC) vs PA-PC (95th pct).
  # PA-FA: Observed (FA) vs PA-FA (95th pct).
  scree_data <- data.frame(
    k = rep(cr$k, 4L),
    value = c(cr$ev_obs, cr$pa_pc_quant, cr$ev_obs_fa, cr$pa_fa_quant),
    series = rep(
      c("Observed (PC)", "PA-PC (95th pct)", "Observed (FA)", "PA-FA (95th pct)"),
      each = k_max
    ),
    panel = "Scree / Parallel Analysis",
    stringsAsFactors = FALSE
  )

  map_data <- data.frame(
    k = cr$k,
    value = cr$map,
    series = "MAP",
    panel = "MAP (minimize)",
    stringsAsFactors = FALSE
  )

  vss_data <- data.frame(
    k = rep(cr$k, 2L),
    value = c(cr$vss1, cr$vss2),
    series = rep(c("VSS-1", "VSS-2"), each = k_max),
    panel = "VSS (maximize)",
    stringsAsFactors = FALSE
  )

  plot_data <- rbind(scree_data, map_data, vss_data)
  panel_levels <- c("Scree / Parallel Analysis", "MAP (minimize)", "VSS (maximize)")

  all_series <- c(
    "Observed (PC)", "PA-PC (95th pct)",
    "Observed (FA)", "PA-FA (95th pct)",
    "MAP", "VSS-1", "VSS-2"
  )

  # Add CD panel when available; it gets its own row in the 2x2 grid.
  if (show_cd) {
    cd_data <- data.frame(
      k = seq_len(k_max),
      value = object$cd_rmse,
      series = "CD (RMSE)",
      panel = "CD (RMSE, minimize)",
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, cd_data)
    panel_levels <- c(panel_levels, "CD (RMSE, minimize)")
    all_series <- c(all_series, "CD (RMSE)")
  }

  plot_data$panel <- factor(plot_data$panel, levels = panel_levels)
  plot_data$series <- factor(plot_data$series, levels = all_series)

  # --- Mark optimal k for each criterion with a star point --------------------
  plot_data$is_opt <- FALSE

  # Scree: PA-PC star on PC observed at k_parallel_pc;
  #        PA-FA star on FA observed at k_parallel_fa (omit if NA).
  plot_data$is_opt <- plot_data$is_opt |
    (plot_data$panel == "Scree / Parallel Analysis" &
      plot_data$series == "Observed (PC)" &
      plot_data$k == object$k_parallel_pc)

  if (!is.na(object$k_parallel_fa)) {
    plot_data$is_opt <- plot_data$is_opt |
      (plot_data$panel == "Scree / Parallel Analysis" &
        plot_data$series == "Observed (FA)" &
        plot_data$k == object$k_parallel_fa)
  }

  # MAP: mark the minimum
  plot_data$is_opt <- plot_data$is_opt |
    (plot_data$panel == "MAP (minimize)" & plot_data$k == object$k_map)

  # VSS: mark each criterion's maximum
  plot_data$is_opt <- plot_data$is_opt |
    (plot_data$series == "VSS-1" & plot_data$k == object$k_vss1) |
    (plot_data$series == "VSS-2" & plot_data$k == object$k_vss2)

  # CD: mark its retention threshold on the CD panel (not MAP)
  if (show_cd && !is.na(object$k_cd)) {
    plot_data$is_opt <- plot_data$is_opt |
      (plot_data$series == "CD (RMSE)" & plot_data$k == object$k_cd)
  }

  # --- Scales -----------------------------------------------------------------
  series_color <- c(
    "Observed (PC)" = "#2166AC",
    "PA-PC (95th pct)" = "#999999",
    "Observed (FA)" = "#4DAC26",
    "PA-FA (95th pct)" = "#BBBBBB",
    "MAP" = "#D6604D",
    "VSS-1" = "#1A9850",
    "VSS-2" = "#762A83",
    "CD (RMSE)" = "#B35806"
  )
  series_lt <- c(
    "Observed (PC)" = "solid",
    "PA-PC (95th pct)" = "dashed",
    "Observed (FA)" = "solid",
    "PA-FA (95th pct)" = "dotted",
    "MAP" = "solid",
    "VSS-1" = "solid",
    "VSS-2" = "dashed",
    "CD (RMSE)" = "solid"
  )
  series_shape <- c(
    "Observed (PC)" = 16L,
    "PA-PC (95th pct)" = 17L,
    "Observed (FA)" = 15L,
    "PA-FA (95th pct)" = 4L,
    "MAP" = 16L,
    "VSS-1" = 16L,
    "VSS-2" = 17L,
    "CD (RMSE)" = 16L
  )

  # ncol = 2 when CD panel present (2x2 grid); single column otherwise.
  n_col <- if (show_cd) 2L else 1L

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x        = .data$k,
      y        = .data$value,
      color    = .data$series,
      linetype = .data$series,
      shape    = .data$series
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_point(
      data        = plot_data[plot_data$is_opt, , drop = FALSE],
      size        = 3.5,
      shape       = 8L,
      stroke      = 1.2,
      show.legend = FALSE,
      inherit.aes = TRUE
    ) +
    ggplot2::scale_color_manual(values = series_color, name = NULL) +
    ggplot2::scale_linetype_manual(values = series_lt, name = NULL) +
    ggplot2::scale_shape_manual(values = series_shape, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(k_max)) +
    ggplot2::facet_wrap(~panel, ncol = n_col, scales = "free_y") +
    ggplot2::labs(
      x = "Number of factors / components (k)",
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text       = ggplot2::element_text(face = "bold")
    )

  p
}
