#' Bass-ackwards hierarchical structural analysis
#'
#' Extracts factor/component solutions at levels 1 through `k`, then
#' characterises the hierarchy by computing between-level factor-score
#' correlations. The "hierarchy" is descriptive: edges are score correlations,
#' not a fitted higher-order SEM.
#'
#' @section Defaults and why:
#' * **`method = "pca"`** — the original Goldberg (2006) method; fastest; never
#'   fails to converge; the Waller (2007) algebra is exact for components.
#' * **`rotation = "cfT"`** — orthogonal Crawford-Ferguson (≈ varimax) keeps
#'   within-level factors uncorrelated, so cross-level edges reflect only the
#'   hierarchical signal, not within-level factor covariance. Matches Goldberg.
#' * **`cor = "pearson"`** — no silent basis switching. If your items look
#'   ordinal (≤ 7 distinct integer values), a cli warning will suggest
#'   `cor = "polychoric"` (available in a later release).
#' * **`align = TRUE`** — unaligned signs make the output unreadable. Anchor:
#'   m1f1 is oriented toward the positive manifold; each subsequent factor is
#'   flipped so its edge to its primary parent is positive.
#' * **`scores = FALSE` / `keep_fits = FALSE`** — memory and privacy. Scores
#'   are O(n × Σk) and often sensitive; raw engine fits can be large. Both are
#'   recomputable from the stored `r` matrix.
#'
#' @param data A data frame or numeric matrix of observed variables (items in
#'   columns, observations in rows). Missing values are handled via pairwise
#'   deletion when computing `R`.
#' @param k Maximum number of factors/components. Required; use
#'   `suggest_k()` (future release) if uncertain. Sets the depth of the
#'   hierarchy: levels 1 through `k` are all extracted.
#' @param method Extraction engine: `"pca"` (default). `"efa"` and `"esem"`
#'   are planned for later milestones.
#' @param rotation Rotation family: `"cfT"` (orthogonal CF, default) or
#'   `"cfQ"` (oblique CF). Oblique is theoretically appropriate when
#'   within-level construct correlations matter, but complicates cross-level
#'   edge interpretation. See DESIGN.md §9.
#' @param kappa CF rotation kappa parameter. `NULL` (default) uses
#'   `1 / (2 * p)` where `p` is the number of variables — the value that
#'   reproduces varimax for orthogonal rotation.
#' @param cor Correlation basis: `"pearson"` (default) or `"spearman"`. A
#'   `"polychoric"` option is planned.
#' @param align Logical; sign-align factors to primary-parent lineage?
#'   Default `TRUE`.
#' @param scores Logical; store factor scores in the result? Default `FALSE`
#'   (recomputable from `r` and the weight matrices).
#' @param keep_fits Logical; store raw engine fit objects? Default `FALSE`.
#' @param seed Integer seed for stochastic engines (not used by PCA but
#'   captured for reproducibility metadata). Default `NULL`.
#' @param cut_show Edges with `|r| >= cut_show` are flagged `above_cut` in
#'   `tidy()` output. Default `0.3`.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"ackwards"`. See [print.ackwards()],
#'   [tidy.ackwards()], and [glance.ackwards()] for output methods.
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
#' x <- ackwards(psych::bfi[, 1:25], k = 5)
#' print(x)
#' tidy(x)
#' glance(x)
#' }
#'
#' @export
ackwards <- function(
    data,
    k,
    method    = "pca",
    rotation  = "cfT",
    kappa     = NULL,
    cor       = "pearson",
    align     = TRUE,
    scores    = FALSE,
    keep_fits = FALSE,
    seed      = NULL,
    cut_show  = 0.3,
    ...) {
  cl <- match.call()

  # --- Input validation -------------------------------------------------------
  method   <- rlang::arg_match(method,   c("pca"))
  rotation <- rlang::arg_match(rotation, c("cfT", "cfQ"))
  cor      <- rlang::arg_match(cor,      c("pearson", "spearman"))

  if (!is.numeric(k) || length(k) != 1L || k < 1L || k != as.integer(k)) {
    cli::cli_abort("{.arg k} must be a positive integer.")
  }
  k <- as.integer(k)

  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
  }
  data_mat <- as.matrix(data)
  if (!is.numeric(data_mat)) {
    cli::cli_abort("{.arg data} must contain only numeric columns.")
  }

  p     <- ncol(data_mat)
  n_obs <- nrow(data_mat)

  if (k > p) {
    cli::cli_abort("{.arg k} ({k}) cannot exceed the number of variables ({p}).")
  }

  # --- Ordinal-detection warning (DESIGN.md §9, Invariant 6) -----------------
  if (detect_ordinal(as.data.frame(data_mat))) {
    cli::cli_warn(
      c(
        "!" = "One or more columns look like ordinal/Likert items \\
               ({cli::symbol$ellipsis} {.val <= 7} distinct integer values).",
        "i" = "Results use a Pearson basis. Consider {.code cor = \"polychoric\"} \\
               for ordinal data (available in a future release)."
      ),
      .frequency = "once",
      .frequency_id = "ackwards_ordinal_warning"
    )
  }

  # --- kappa (stored in meta for future cfT/cfQ path; not used in M1 PCA) ----
  if (is.null(kappa)) kappa <- 1 / (2 * p)

  # --- Seed capture -----------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  # --- Compute correlation matrix ---------------------------------------------
  R <- stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")

  # --- Extract levels ---------------------------------------------------------
  levels_list <- switch(method,
    pca = pca_levels(R, k_max = k, rotation = rotation)
  )

  # --- Lineage matching -------------------------------------------------------
  # Build lineage before sign alignment (matching is on |r|, sign-invariant)
  raw_edges <- compute_edges(
    levels  = levels_list,
    R       = R,
    method  = "auto",
    pairs   = "adjacent",
    align   = FALSE,
    cut_show = cut_show
  )

  lineage <- vector("list", k)
  names(lineage) <- as.character(seq_len(k))
  for (ki in seq_len(k - 1L) + 1L) {
    key <- paste0(ki - 1L, ":", ki)
    lineage[[as.character(ki)]] <- match_parents(raw_edges$matrices[[key]])
  }

  # --- Sign alignment (DESIGN.md §7) -----------------------------------------
  if (align && k >= 1L) {
    loadings_list <- lapply(levels_list, `[[`, "loadings")
    aligned <- align_signs(loadings_list, raw_edges$matrices, lineage)
    for (ki in seq_len(k)) {
      levels_list[[as.character(ki)]]$loadings <- aligned$loadings[[ki]]
      # Flip weights consistently with loadings
      levels_list[[as.character(ki)]]$scoring$weights <-
        flip_weights(levels_list[[as.character(ki)]]$scoring$weights,
                     aligned$signs[[ki]])
    }
    edge_matrices <- aligned$edges
  } else {
    edge_matrices <- raw_edges$matrices
  }

  # --- Rebuild tidy edge tibble with aligned matrices -------------------------
  edges_obj <- compute_edges(
    levels   = levels_list,
    R        = R,
    method   = "auto",
    pairs    = "adjacent",
    align    = FALSE,    # already aligned above
    cut_show = cut_show
  )
  # Overwrite matrices with the aligned versions
  edges_obj$matrices <- edge_matrices
  # Rebuild tidy from aligned matrices
  edges_obj$tidy <- .rebuild_tidy(edge_matrices, cut_show)
  edges_obj$tidy <- fill_primary(edges_obj$tidy, lineage, levels_list)

  # --- Meta -------------------------------------------------------------------
  meta <- list(
    converged_levels = vapply(levels_list, `[[`, logical(1L), "converged"),
    deepest_converged = max(which(vapply(levels_list, `[[`, logical(1L), "converged"))),
    kappa             = kappa,
    ordinal_warned    = detect_ordinal(as.data.frame(data_mat))
  )

  # --- Assemble result --------------------------------------------------------
  new_ackwards(
    call        = cl,
    method      = method,
    rotation    = rotation,
    cor_type    = cor,
    n_obs       = n_obs,
    k_max       = k,
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
}

# S3 constructor — validates structure and attaches class
new_ackwards <- function(
    call, method, rotation, cor_type, n_obs, k_max, seed, pkg_version,
    levels, edges, lineage, scores, fits, r, data, meta) {
  structure(
    list(
      call        = call,
      method      = method,
      rotation    = rotation,
      cor_type    = cor_type,
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
      meta        = meta
    ),
    class = "ackwards"
  )
}

# Rebuild a tidy edge data frame from a named list of aligned edge matrices.
.rebuild_tidy <- function(matrices, cut_show = 0.3) {
  rows <- lapply(names(matrices), function(key) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
    ka <- as.integer(parts[1L])
    kb <- as.integer(parts[2L])
    E  <- matrices[[key]]
    from_labs <- rownames(E)
    to_labs   <- colnames(E)
    do.call(rbind, lapply(seq_along(from_labs), function(i) {
      do.call(rbind, lapply(seq_along(to_labs), function(j) {
        data.frame(
          from       = from_labs[i],
          to         = to_labs[j],
          level_from = ka,
          level_to   = kb,
          r          = E[i, j],
          is_primary = NA,
          above_cut  = abs(E[i, j]) >= cut_show,
          stringsAsFactors = FALSE
        )
      }))
    }))
  })
  tidy <- do.call(rbind, rows)
  rownames(tidy) <- NULL
  tidy
}
