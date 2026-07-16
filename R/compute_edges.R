#' Compute between-level factor-score correlations
#'
#' The centerpiece of the bass-ackwards algebra. For any engine whose scoring
#' is a **linear** map `S = Z W`, the cross-level correlation matrix is:
#'
#' ```
#' E(a,b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}
#' ```
#'
#' where `R` is the input correlation matrix and `D_x = diag(W_x' R W_x)` are
#' the **actual** score variances (not assumed to be 1). This avoids
#' materialising scores while remaining exact for PCA, EFA (regression /
#' Bartlett / tenBerge) -- all of which produce linear score maps.
#'
#' When the algebra cannot be used (nonlinear scoring, missing `R`, or the user
#' forces `edge_method = "scores"`), scores are materialised from `data` instead.
#'
#' @param levels Named list (indexed by k) of per-level objects produced by an
#'   engine. Each must contain a `scoring` sub-list with fields `linear`,
#'   `weights`, and `score_var`.
#' @param R Square correlation matrix (p x p). Required for the algebra path.
#' @param edge_method One of `"auto"` (algebra when possible, scores otherwise),
#'   `"algebra"` (force; errors if conditions not met), or `"scores"` (always
#'   materialise).
#' @param pairs `"adjacent"` (classic Goldberg) or `"all"` (Forbes extension).
#' @param data Optional data frame / matrix of raw observations. Required only
#'   when `edge_method = "scores"` or the scores path is triggered.
#' @param use Passed to [stats::cor()] when materialising scores.
#' @param cut_show Edges with `|r| >= cut_show` are flagged `above_cut` in the
#'   tidy tibble.
#' @param build_tidy Build the tidy edge data frame? `FALSE` returns
#'   `tidy = NULL` for matrices-only callers (the pre-alignment lineage pass in
#'   `ackwards()`, `comparability()`'s cross-solution blocks, `boot_edges()`
#'   replicates), which would otherwise build and discard it (M60).
#'
#' @return A list with:
#'   \item{matrices}{Named list of `(k_a x k_b)` edge matrices, keyed
#'     `"k_a:k_b"`.}
#'   \item{tidy}{A data frame with one row per directed edge: `from`, `to`,
#'     `level_from`, `level_to`, `r`, `is_primary`, `above_cut` -- or `NULL`
#'     when `build_tidy = FALSE`.}
#'
#' @keywords internal
compute_edges <- function(
  levels,
  R,
  edge_method = c("auto", "algebra", "scores"),
  pairs = c("adjacent", "all"),
  data = NULL,
  use = "pairwise.complete.obs",
  cut_show = 0.3,
  build_tidy = TRUE
) {
  edge_method <- match.arg(edge_method)
  pairs <- match.arg(pairs)

  K <- length(levels)
  ks <- as.integer(names(levels))

  # Build list of index pairs to compute
  if (pairs == "adjacent") {
    pair_list <- lapply(seq_len(K - 1L), function(i) c(ks[i], ks[i + 1L]))
  } else {
    pair_list <- do.call(c, lapply(seq_len(K - 1L), function(i) {
      lapply(seq(i + 1L, K), function(j) c(ks[i], ks[j]))
    }))
  }

  matrices <- list()

  for (pair in pair_list) {
    ka <- pair[1L]
    kb <- pair[2L]
    key <- paste0(ka, ":", kb)

    la <- levels[[as.character(ka)]]
    lb <- levels[[as.character(kb)]]

    algebra_ok <- isTRUE(la$scoring$linear) &&
      isTRUE(lb$scoring$linear) &&
      !is.null(R)

    use_algebra <- switch(edge_method,
      auto = algebra_ok,
      algebra = {
        if (!algebra_ok) {
          cli::cli_abort(
            "edge_method = 'algebra' requested but conditions not met for pair {ka}:{kb}."
          )
        }
        TRUE
      },
      scores = FALSE
    )

    if (use_algebra) {
      Wa <- la$scoring$weights
      Wb <- lb$scoring$weights
      C <- crossprod(Wa, R %*% Wb)
      sa <- sqrt(.score_var(Wa, R))
      sb <- sqrt(.score_var(Wb, R))
      E <- sweep(sweep(C, 1L, sa, "/"), 2L, sb, "/")
    } else {
      if (is.null(data)) {
        cli::cli_abort(
          "Scores path requires {.arg data}; none supplied for pair {ka}:{kb}."
        )
      }
      Z <- .standardize(as.matrix(data))
      Sa <- Z %*% la$scoring$weights
      Sb <- Z %*% lb$scoring$weights
      E <- stats::cor(Sa, Sb, use = use)
    }

    rownames(E) <- la$labels
    colnames(E) <- lb$labels
    matrices[[key]] <- E
  }

  if (!build_tidy) {
    return(list(matrices = matrices, tidy = NULL))
  }

  # Build tidy edge tibble (before sign alignment, lineage not yet known)
  tidy_rows <- lapply(names(matrices), function(key) {
    parts <- strsplit(key, ":", fixed = TRUE)[[1L]]
    ka <- as.integer(parts[1L])
    kb <- as.integer(parts[2L])
    E <- matrices[[key]]
    from_labs <- rownames(E)
    to_labs <- colnames(E)
    do.call(rbind, lapply(seq_along(from_labs), function(i) {
      do.call(rbind, lapply(seq_along(to_labs), function(j) {
        data.frame(
          from = from_labs[i],
          to = to_labs[j],
          level_from = ka,
          level_to = kb,
          r = E[i, j],
          is_primary = NA, # filled after lineage matching
          above_cut = abs(E[i, j]) >= cut_show,
          stringsAsFactors = FALSE
        )
      }))
    }))
  })

  tidy <- do.call(rbind, tidy_rows)
  rownames(tidy) <- NULL

  list(matrices = matrices, tidy = tidy)
}

# Fill in `is_primary` in the tidy edge data frame once lineage is known.
# lineage: list indexed by k (k >= 2), integer vectors: parent index in k-1
# for each factor in k.
fill_primary <- function(tidy, lineage, levels) {
  ks <- as.integer(names(levels))
  for (k in ks[ks >= 2L]) {
    parents <- lineage[[as.character(k)]]
    la_labs <- levels[[as.character(k - 1L)]]$labels
    lb_labs <- levels[[as.character(k)]]$labels
    for (j in seq_along(lb_labs)) {
      primary_from <- la_labs[parents[j]]
      primary_to <- lb_labs[j]
      mask <- tidy$from == primary_from & tidy$to == primary_to
      tidy$is_primary[mask] <- TRUE
      # Mark non-primary edges in the same column
      other_mask <- tidy$level_from == (k - 1L) &
        tidy$level_to == k &
        tidy$to == primary_to &
        !mask
      tidy$is_primary[other_mask] <- FALSE
    }
  }
  # Any remaining NA is_primary: false (edges not in the primary tree)
  tidy$is_primary[is.na(tidy$is_primary)] <- FALSE
  tidy
}
