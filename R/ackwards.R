#' Bass-ackwards hierarchical structural analysis
#'
#' Extracts factor/component solutions at levels 1 through `k`, then
#' characterises the hierarchy by computing between-level factor-score
#' correlations. The "hierarchy" is descriptive: edges are score correlations,
#' not a fitted higher-order SEM.
#'
#' @section Defaults and why:
#' * **`engine = "pca"`** — the original Goldberg (2006) method; fastest; never
#'   fails to converge; the Waller (2007) algebra is exact for components.
#' * **`rotation = "varimax"`** — the `T'=T^-1` property of orthogonal
#'   rotation enables the closed-form `W'RW` edge algebra and keeps
#'   within-level factors uncorrelated so cross-level edges reflect only the
#'   hierarchical signal. Matches Goldberg (2006), Kim & Eaton (2015), and
#'   Forbush et al. (2024). Varimax is the only supported rotation; oblique
#'   rotation would confound the between-level signal that is the method's
#'   core output.
#' * **`cor = "pearson"`** — no silent basis switching. If your items look
#'   ordinal (≤ 7 distinct integer values), a cli warning will suggest
#'   `cor = "polychoric"`, which is available for all three engines.
#' * **`align_signs = TRUE`** — unaligned signs make the output unreadable.
#'   Anchor: m1f1 is oriented toward the positive manifold; each subsequent
#'   factor is flipped so its edge to its primary parent is positive.
#' * **`keep_scores = FALSE` / `keep_fits = FALSE`** — memory and privacy.
#'   Scores are O(n × Σk) and often sensitive; raw engine fits can be large.
#'   Both are recomputable from the stored `r` matrix.
#'
#' @param data A data frame or numeric matrix of observed variables (items in
#'   columns, observations in rows). Missing values are handled via pairwise
#'   deletion when computing `R`. Note: `n_obs` passed to the EFA engine is
#'   always `nrow(data)`; under missingness the effective per-correlation N may
#'   be smaller, making chi-square / RMSEA / p-value slightly anti-conservative.
#' @param k_max Maximum number of factors/components to extract. Required; use
#'   [suggest_k()] if uncertain. Sets the depth of the hierarchy: levels
#'   1 through `k_max` are all extracted.
#' @param engine Extraction engine: `"pca"` (default), `"efa"`, or `"esem"`.
#'   `"esem"` uses [lavaan::efa()] with rotation-aware SEs and per-level fit
#'   indices; recommended for the clinical/HiTOP workflow (Kim & Eaton, 2015;
#'   Forbush et al., 2024). Requires lavaan >= 0.6-13.
#' @param fm Factor extraction method passed to [psych::fa()]; only used when
#'   `engine = "efa"`. One of `"minres"` (default, robust OLS), `"ml"`
#'   (maximum likelihood, gives chi-square fit but converges less reliably at
#'   deep levels), or `"pa"` (principal axis). Ignored for `engine = "pca"`.
#' @param cor Correlation basis: `"pearson"` (default), `"spearman"`, or
#'   `"polychoric"`. For PCA/EFA, `"polychoric"` computes a polychoric
#'   correlation matrix via `psych::polychoric()` (requires psych). For ESEM,
#'   it triggers WLSMV estimation via lavaan. A cli warning is emitted when
#'   ordinal-looking columns are detected and `cor` is not `"polychoric"`.
#' @param estimator Estimation method for the ESEM engine. `NULL` (default)
#'   auto-selects: `"WLSMV"` when `cor = "polychoric"`, `"ML"` otherwise.
#'   Pass explicitly to override: `"ULSMV"` (unweighted WLS), `"MLR"`
#'   (robust ML). Ignored for PCA and EFA engines.
#' @param align_signs Logical; sign-align factors to primary-parent lineage?
#'   Default `TRUE`.
#' @param keep_scores Logical; store factor scores in the result? Default
#'   `FALSE` (recomputable via [augment.ackwards()]). When `TRUE`,
#'   per-observation scores are stored in `x$scores` as a named list of
#'   `n × k_j` matrices, one per level, standardized by real score SDs
#'   (see [augment.ackwards()]).
#' @param keep_fits Logical; store raw engine fit objects? Default `FALSE`.
#'   When `TRUE`, the per-level fit objects (psych or lavaan) are stored in
#'   `x$fits` as a named list indexed by level.
#' @param seed Integer seed for stochastic engines (not used by PCA but
#'   captured for reproducibility metadata). Default `NULL`.
#' @param pairs Which level pairs to compute edges for: `"adjacent"` (default,
#'   classic Goldberg — only consecutive levels) or `"all"` (Forbes extension —
#'   every pair of levels). `"all"` reveals associations that span multiple
#'   levels and is required for redundancy pruning. Setting `prune` to anything
#'   other than `"none"` automatically upgrades this to `"all"` with a message.
#' @param prune Character vector controlling Forbes-extension pruning. Default
#'   `"none"` (no pruning). Options:
#'   * `"redundant"` — identify chains of factors connected by primary-parent
#'     links with `|r| >= redundancy_r` (and optionally `phi > redundancy_phi`).
#'     Applies Forbes's (2023) retention rule: keep the bottom node when the
#'     chain reaches level `k` (most specific); keep the top node otherwise.
#'     Pruning is *flag-only*: flagged nodes stay in the object with
#'     `pruned = TRUE` and `prune_reason = "redundant"` in `x$prune$nodes`.
#'   * `"artefact"` — compute Tucker's congruence coefficient (φ) for all
#'     cross-level factor pairs and store in `x$prune$phi` for researcher
#'     inspection. No factors are auto-flagged; artefact identification requires
#'     judgment (Forbes, 2023; Wicherts et al., 2016).
#' @param redundancy_r Scalar in `(0, 1]`. Adjacent primary-parent `|r|`
#'   threshold for redundancy chains. Default `0.9` (Forbes, 2023).
#' @param redundancy_phi Scalar in `(0, 1]` or `NULL` (default). If non-`NULL`,
#'   Tucker's φ must *also* exceed this threshold for a link to be included in
#'   a redundancy chain (conjunctive with `redundancy_r`). `NULL` means only
#'   `redundancy_r` is used. Recommended: `0.95` (Lorenzo-Seva & ten Berge, 2006).
#' @param cut_show Edges with `|r| >= cut_show` are flagged `above_cut` in
#'   `tidy()` output. Default `0.3`.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"ackwards"`. See [print.ackwards()],
#'   [tidy.ackwards()], [glance.ackwards()], and [augment.ackwards()] for
#'   output methods.
#'
#' @seealso [print.ackwards()], [tidy.ackwards()], [glance.ackwards()]
#'
#' @references
#' Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research in
#'   Personality*, 40(4), 347–358. \doi{10.1016/j.jrp.2006.01.001}
#'
#' Waller, N. G. (2007). A general method for computing hierarchical component
#'   structures by Goldberg's bass-ackwards method. *Journal of Research in
#'   Personality*, 41(4), 745–752. \doi{10.1016/j.jrp.2006.08.005}
#'
#' @examples
#' \dontrun{
#' x <- ackwards(psych::bfi[, 1:25], k_max = 5)
#' print(x)
#' tidy(x)
#' glance(x)
#' }
#'
#' @export
ackwards <- function(
  data,
  k_max,
  engine = "pca",
  cor = "pearson",
  fm = "minres",
  estimator = NULL,
  align_signs = TRUE,
  keep_scores = FALSE,
  keep_fits = FALSE,
  seed = NULL,
  pairs = "adjacent",
  prune = "none",
  redundancy_r = 0.9,
  redundancy_phi = NULL,
  cut_show = 0.3,
  ...
) {
  cl <- match.call()

  # --- Input validation -------------------------------------------------------
  engine <- rlang::arg_match(engine, c("pca", "efa", "esem"))
  fm <- rlang::arg_match(fm, c("minres", "ml", "pa"))
  cor <- rlang::arg_match(cor, c("pearson", "spearman", "polychoric"))
  if (!is.null(estimator)) {
    estimator <- rlang::arg_match(estimator, c("ML", "MLR", "WLSMV", "ULSMV"))
  }
  pairs <- rlang::arg_match(pairs, c("adjacent", "all"))
  prune <- rlang::arg_match(prune, c("none", "redundant", "artefact"),
    multiple = TRUE
  )

  if (!is.numeric(redundancy_r) || length(redundancy_r) != 1L ||
    redundancy_r <= 0 || redundancy_r > 1) {
    cli::cli_abort("{.arg redundancy_r} must be a single number in (0, 1].")
  }
  if (!is.null(redundancy_phi) &&
    (!is.numeric(redundancy_phi) || length(redundancy_phi) != 1L ||
      redundancy_phi <= 0 || redundancy_phi > 1)) {
    cli::cli_abort("{.arg redundancy_phi} must be a single number in (0, 1] or {.code NULL}.")
  }

  # Auto-upgrade pairs to "all" when pruning is requested (chains need all-levels
  # edges for the endpoint-r enrichment; artefact phi is most informative across
  # all pairs). Invariant 6: announce the upgrade rather than switching silently.
  if (!identical(prune, "none") && pairs == "adjacent") {
    pairs <- "all"
    cli::cli_inform(c(
      "i" = "{.arg pairs} upgraded to {.val all}: pruning requires all-levels \\
             edges to assess redundancy chains and artefact relationships."
    ))
  }

  # cor = "spearman" + engine = "esem" is semantically inconsistent: lavaan fits
  # Pearson ML on raw data while compute_edges() uses Spearman R for scoring
  # (DESIGN.md §14 known limitations). Warn loudly rather than silently mixing bases.
  if (engine == "esem" && cor == "spearman") {
    cli::cli_warn(
      c(
        "!" = "{.code cor = \"spearman\"} with {.code engine = \"esem\"} \\
               uses inconsistent bases: lavaan fits a Pearson-ML model on raw \\
               data while edges are computed from the Spearman correlation matrix.",
        "i" = "Consider {.code cor = \"polychoric\"} for ordinal data or \\
               {.code cor = \"pearson\"} for a consistent continuous path."
      ),
      .frequency = "once",
      .frequency_id = "ackwards_esem_spearman"
    )
  }

  if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 2L || k_max != as.integer(k_max)) {
    cli::cli_abort("{.arg k_max} must be an integer >= 2 (need at least two levels for a hierarchy).")
  }
  k_max <- as.integer(k_max)

  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
  }
  data_mat <- as.matrix(data)
  if (!is.numeric(data_mat)) {
    cli::cli_abort("{.arg data} must contain only numeric columns.")
  }

  p <- ncol(data_mat)
  n_obs <- nrow(data_mat)

  if (k_max > p) {
    cli::cli_abort("{.arg k_max} ({k_max}) cannot exceed the number of variables ({p}).")
  }

  # --- Ordinal-detection warning (DESIGN.md §9, Invariant 6) -----------------
  # Compute once; reused in meta$ordinal_warned below.
  # Only warn when the user has NOT already opted into the polychoric basis.
  is_ordinal <- detect_ordinal(as.data.frame(data_mat))
  if (is_ordinal && cor != "polychoric") {
    cli::cli_warn(
      c(
        "!" = "One or more columns look like ordinal/Likert items \\
               ({cli::symbol$ellipsis} {.val <= 7} distinct integer values).",
        "i" = "Results use a {.val {cor}} basis. Consider \\
               {.code cor = \"polychoric\"} for ordinal data."
      ),
      .frequency = "once",
      .frequency_id = "ackwards_ordinal_warning"
    )
  }

  # --- Seed capture -----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  # --- Compute correlation matrix + extract levels ----------------------------
  # ESEM: lavaan owns R computation for polychoric (WLSMV uses its own polychoric
  # matrix internally); for continuous paths we pre-compute R and pass it through.
  # PCA/EFA: always compute R here and pass to the engine.
  if (engine == "esem") {
    R_ext <- if (cor != "polychoric") {
      stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
    } else {
      NULL
    }
    estimator_eff <- if (!is.null(estimator)) {
      estimator
    } else if (cor == "polychoric") {
      "WLSMV"
    } else {
      "ML"
    }
    esem_out <- esem_levels(data_mat,
      k_max = k_max, estimator = estimator_eff, cor_type = cor,
      n_obs = n_obs, R_external = R_ext, keep_fits = keep_fits
    )
    levels_list <- esem_out$levels
    fits_stored <- esem_out$fits
    R <- esem_out$r_lv
    if (is.null(R)) R <- R_ext
    if (is.null(R)) {
      R <- stats::cor(data_mat, method = "pearson", use = "pairwise.complete.obs")
    }
  } else {
    # PCA / EFA: compute R then dispatch
    if (cor == "polychoric") {
      rlang::check_installed("psych", reason = "for polychoric correlations")
      poly_out <- tryCatch(
        psych::polychoric(data_mat),
        error = function(e) {
          cli::cli_abort(
            c(
              "!" = "{.fn psych::polychoric} failed: {conditionMessage(e)}",
              "i" = "Check that your data contains integer-like columns with few \\
                   distinct values, or use {.code cor = \"pearson\"}."
            )
          )
        }
      )
      R <- poly_out$rho
      # Guard against non-positive-definite polychoric matrices (DESIGN.md §14 remaining)
      min_eig <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
      if (min_eig <= 0) {
        cli::cli_warn(
          c(
            "!" = "Polychoric correlation matrix is not positive definite \\
                   (min eigenvalue = {round(min_eig, 4)}).",
            "i" = "Applying smoothing via {.fn psych::cor.smooth}."
          )
        )
        R <- psych::cor.smooth(R)
      }
    } else {
      R <- stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
    }
    engine_out <- switch(engine,
      pca = pca_levels(R,
        k_max = k_max, cor_type = cor,
        keep_fits = keep_fits
      ),
      efa = efa_levels(R,
        k_max = k_max, fm = fm, n_obs = n_obs,
        cor_type = cor, keep_fits = keep_fits
      )
    )
    levels_list <- engine_out$levels
    fits_stored <- engine_out$fits
  }

  # --- Handle convergence truncation (PCA never truncates; EFA/ESEM may) -----
  k_eff <- length(levels_list)
  if (k_eff < 2L) {
    cli::cli_abort(
      c(
        "!" = "{toupper(engine)} failed to build at least 2 converged levels \\
               (k_eff = {k_eff}; at least 2 are required for a hierarchy).",
        "i" = "Try fewer factors, or check your data for \\
               perfect multicollinearity or near-singular correlation."
      )
    )
  }
  if (k_eff < k_max) {
    cli::cli_warn(
      c(
        "!" = "Hierarchy truncated: {k_max - k_eff} level{?s} did not converge \\
               (requested k_max = {k_max}, built k = {k_eff}).",
        "i" = "Set {.arg k_max = {k_eff}} to suppress this message."
      )
    )
  }

  # --- Lineage matching -------------------------------------------------------
  # Always adjacent: lineage and sign alignment are defined by parent-child
  # relationships between consecutive levels.
  raw_edges_adj <- compute_edges(
    levels      = levels_list,
    R           = R,
    edge_method = "auto",
    pairs       = "adjacent",
    align       = FALSE,
    cut_show    = cut_show
  )

  lineage <- vector("list", k_eff)
  names(lineage) <- as.character(seq_len(k_eff))
  for (ki in seq_len(k_eff - 1L) + 1L) {
    key <- paste0(ki - 1L, ":", ki)
    lineage[[as.character(ki)]] <- match_parents(raw_edges_adj$matrices[[key]])
  }

  # --- Sign alignment (DESIGN.md §7) -----------------------------------------
  if (align_signs && k_eff >= 1L) {
    loadings_list <- lapply(levels_list, `[[`, "loadings")
    aligned <- align_signs(loadings_list, raw_edges_adj$matrices, lineage)
    for (ki in seq_len(k_eff)) {
      levels_list[[as.character(ki)]]$loadings <- aligned$loadings[[ki]]
      # Flip weights consistently with loadings
      levels_list[[as.character(ki)]]$scoring$weights <-
        flip_weights(
          levels_list[[as.character(ki)]]$scoring$weights,
          aligned$signs[[ki]]
        )
    }
  }

  # --- Final edge set (adjacent or all-levels Forbes extension) ---------------
  # Recompute using aligned levels_list so skip-level edges inherit correct signs.
  # Adjacent-pair matrices match aligned$edges exactly; recomputing is safe.
  final_edges <- compute_edges(
    levels      = levels_list,
    R           = R,
    edge_method = "auto",
    pairs       = pairs,
    align       = FALSE,
    cut_show    = cut_show
  )

  # --- Assemble aligned edge object -------------------------------------------
  edges_obj <- list(
    matrices = final_edges$matrices,
    tidy     = fill_primary(final_edges$tidy, lineage, levels_list)
  )

  # --- Meta -------------------------------------------------------------------
  conv <- vapply(levels_list, `[[`, logical(1L), "converged")
  meta <- list(
    k_requested       = k_max, # what the user asked for (may exceed k_eff)
    converged_levels  = conv,
    deepest_converged = max(which(conv)),
    pairs             = pairs,
    prune             = prune,
    redundancy_r      = redundancy_r,
    redundancy_phi    = redundancy_phi,
    cut_show          = cut_show,
    ordinal_warned    = is_ordinal
  )

  # --- Assemble result --------------------------------------------------------
  x <- new_ackwards(
    call        = cl,
    engine      = engine,
    cor         = cor,
    n_obs       = n_obs,
    k_max       = k_eff, # effective depth (may be < k_max if truncated)
    seed        = seed,
    pkg_version = utils::packageVersion("ackwards"),
    levels      = levels_list,
    edges       = edges_obj,
    lineage     = lineage,
    scores      = NULL,
    fits        = NULL,
    r           = R,
    data        = NULL,
    meta        = meta
  )

  # --- Opt-in storage (scores and raw fits) -----------------------------------
  # Scores use the sign-aligned levels_list so column orientations match loadings.
  if (isTRUE(keep_scores)) {
    x$scores <- .compute_scores(levels_list, data_mat)
  }
  if (isTRUE(keep_fits)) {
    x$fits <- fits_stored
  }

  # --- Forbes pruning (flag-only, never removes levels) -----------------------
  # Pruning is applied post-assembly so .apply_pruning() can read the final
  # aligned edges directly from the object.
  if (!identical(prune, "none")) {
    x$prune <- .apply_pruning(x, prune, redundancy_r, redundancy_phi)
  }

  x
}

# S3 constructor — validates structure and attaches class
new_ackwards <- function(
  call, engine, cor, n_obs, k_max, seed, pkg_version,
  levels, edges, lineage, scores, fits, r, data, meta
) {
  structure(
    list(
      call        = call,
      engine      = engine,
      rotation    = "varimax",
      cor         = cor,
      n_obs       = n_obs,
      k_max       = k_max,
      seed        = seed,
      pkg_version = pkg_version,
      levels      = levels,
      edges       = edges,
      lineage     = lineage,
      scores      = scores,
      fits        = fits,
      r           = r,
      data        = data,
      meta        = meta,
      prune       = NULL # populated by .apply_pruning() when prune != "none"
    ),
    class = "ackwards"
  )
}
