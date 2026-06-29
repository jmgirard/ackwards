# Internal utilities — not exported

# Column-wise na.rm-aware standardization: each column is centered and scaled
# using its own non-NA observations. NA values in a column remain NA after
# standardization, so only rows that are missing a particular item produce NA
# in that column of Z — and consequently NA scores only for those rows (via the
# linear map S = Z W). Contrast with base scale(), which returns a fully-NA
# column when any value is missing, making ALL scores NA.
# A zero-variance or all-NA column gets scale factor 1 (avoids Inf/NaN).
.standardize <- function(x) {
  m <- colMeans(x, na.rm = TRUE)
  s <- apply(x, 2, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(x, 2, m, "-"), 2, s, "/")
}

# Format a correlation as an APA-style string: strip the leading zero, pad
# trailing zeros to `digits` decimal places. Does not prepend "-" when the
# magnitude rounds to zero at the requested precision (avoids "-.00").
# Examples (digits = 2): 0.23 -> ".23", -0.3 -> "-.30", 0 -> ".00",
#   1 -> "1.00", -0.003 -> ".00"
.format_r <- function(r, digits = 2L) {
  fmt <- formatC(abs(r), digits = digits, format = "f")
  fmt <- sub("^0\\.", ".", fmt)
  is_zero <- grepl("^[.0]+$", fmt)
  ifelse(r < 0 & !is_zero, paste0("-", fmt), fmt)
}

# Generate standard factor labels for level k with j factors: "m{k}f{j}"
make_labels <- function(k) {
  paste0("m", k, "f", seq_len(k))
}

# Detect whether a data frame likely contains ordinal (Likert-scale) columns.
# Heuristic: a column is flagged if it is integer-like and has <= max_levels
# distinct values. Returns TRUE if any column is flagged.
detect_ordinal <- function(data, max_levels = 7L) {
  is_int_like <- function(x) {
    is.integer(x) || (is.numeric(x) && all(x == floor(x), na.rm = TRUE))
  }
  cols <- vapply(data, function(x) {
    vals <- x[!is.na(x)]
    length(vals) > 0L &&
      is_int_like(x) &&
      length(unique(vals)) <= max_levels
  }, logical(1))
  any(cols)
}

# Primary-parent assignment: for each factor in level b, find the factor in level a
# with the highest |r|. Returns an integer vector of length ncol(E).
# Each child picks its argmax parent independently; multiple children can (and do)
# share a parent — that is normal and expected in the bass-ackwards hierarchy.
# LSAP (bijection) is wrong here: adjacent levels always have n_b = n_a + 1, so
# a bijection would need a padding row that can return an index > nrow(E).
match_parents <- function(E) {
  apply(abs(E), 2, which.max)
}

# Sign-align a list of per-level loadings matrices and a corresponding list of
# edge matrices.
#
# Rules (DESIGN.md §7):
#   1. Anchor m1f1: flip so sum of loadings column is positive.
#      Extend to all factors in level 1: each flipped so sum > 0.
#   2. For each subsequent level: flip each factor so its edge correlation
#      with its primary parent (from match_parents) is positive.
#
# Arguments:
#   loadings_list  — list indexed by k (1..K) of p×k loading matrices
#   edges_list     — list named "k_a:k_b" of (k_a × k_b) edge matrices
#   lineage        — list indexed by k>=2 of integer vectors (parent indices)
#
# Returns: list(loadings = ..., edges = ..., signs = ...) where signs is a list
# of ±1 vectors (one per level) recording the flip applied.
.align_signs <- function(loadings_list, edges_list, lineage) {
  K <- length(loadings_list)
  signs <- vector("list", K)

  # Level 1: anchor so column sum is positive
  L1 <- loadings_list[[1]]
  s1 <- sign(colSums(L1))
  s1[s1 == 0] <- 1L
  loadings_list[[1]] <- sweep(L1, 2, s1, "*")
  signs[[1]] <- s1

  for (k in seq_len(K - 1) + 1L) {
    # Edge matrix between level k-1 and level k (rows = k-1, cols = k)
    key <- paste0(k - 1L, ":", k)
    E <- edges_list[[key]]
    parents <- lineage[[k]] # parent index in level k-1 for each factor in level k
    # Flip factor j in level k if its correlation with its primary parent is negative
    sk <- integer(ncol(E))
    for (j in seq_len(ncol(E))) {
      sk[j] <- if (E[parents[j], j] >= 0) 1L else -1L
    }
    loadings_list[[k]] <- sweep(loadings_list[[k]], 2, sk, "*")
    edges_list[[key]] <- sweep(sweep(E, 2, sk, "*"), 1, signs[[k - 1L]], "*")
    signs[[k]] <- sk
  }

  list(loadings = loadings_list, edges = edges_list, signs = signs)
}

# Apply previously computed sign vectors to a weight matrix W (p × k).
flip_weights <- function(W, sign_vec) {
  sweep(W, 2, sign_vec, "*")
}

# Materialize per-observation factor scores for all levels.
# Z = .standardize(data_mat): standardizes observed items by their sample mean/SD.
# S_k = Z W_k, then each column divided by sqrt(score_var_k) (Invariant 1:
# standardize by real score SDs, never assume unit variance).
# Returns a named list (by level character key) of n × k_j matrices.
#
# Pearson basis: .standardize() gives Pearson z-scores, so W'RW computed on the
# Pearson R matches and empirical score variance = model-implied = 1.
# Non-Pearson basis (polychoric, Spearman): .standardize() still uses Pearson
# standardization, but score_var comes from the non-Pearson R. Empirical score
# SDs will differ from 1.0; a warning is issued to alert the user.
.compute_scores <- function(levels_list, data_mat) {
  bases <- unique(vapply(levels_list, function(lev) lev$scoring$basis, character(1L)))
  if (!all(bases == "pearson")) {
    non_pearson <- bases[bases != "pearson"]
    cli::cli_warn(
      c(
        "!" = "Factor scores are standardized using model-implied SDs from a \\
               {.val {non_pearson}} correlation matrix.",
        "i" = "The raw projection uses {.code .standardize(data)} (Pearson z-scores), \\
               but {.code score_var} comes from the {.val {non_pearson}} R.",
        "i" = "Empirical score SDs will differ from 1.0. \\
               For non-Pearson analyses, between-level edges from \\
               {.fn tidy} are the authoritative associations."
      ),
      .frequency = "once",
      .frequency_id = "ackwards_nonpearson_scores"
    )
  }
  n_na_rows <- sum(!stats::complete.cases(data_mat))
  if (n_na_rows > 0L) {
    cli::cli_warn(c(
      "!" = "{n_na_rows} row{?s} contain missing values and will produce NA scores.",
      "i" = "Score projection applies weights row-wise and propagates NAs \\
             listwise; FIML estimation does not impute missing item responses.",
      "i" = "Use {.code missing = \"listwise\"} when fitting (so the model \\
             and scores share the same complete rows), or call \\
             {.code na.omit(data)} before scoring."
    ))
  }
  Z <- .standardize(data_mat)
  scores <- lapply(names(levels_list), function(ki) {
    lev <- levels_list[[ki]]
    W <- lev$scoring$weights
    score_var <- lev$scoring$score_var
    S <- Z %*% W
    S <- sweep(S, 2, sqrt(score_var), "/")
    colnames(S) <- lev$labels
    rownames(S) <- rownames(data_mat)
    S
  })
  names(scores) <- names(levels_list)
  scores
}

# Validate the missing= argument against engine and estimator.
# Errors clearly when the combination is unsupported (FIML is only valid for
# ESEM with a full-information estimator — ML or MLR). PCA and EFA work on a
# pre-computed correlation matrix; WLSMV/ULSMV are limited-information WLS
# and have no FIML extension.
.resolve_missing <- function(missing, engine, estimator_eff) {
  if (missing != "fiml") {
    return(invisible(NULL))
  }
  if (engine %in% c("pca", "efa")) {
    cli::cli_abort(c(
      "!" = "{.code missing = \"fiml\"} is not supported for \\
             {.code engine = \"{engine}\"}.",
      "i" = "FIML requires a raw-data likelihood; PCA and EFA operate on a \\
             pre-computed correlation matrix.",
      "i" = "Use {.code missing = \"pairwise\"} or \\
             {.code missing = \"listwise\"} instead."
    ))
  }
  if (!is.null(estimator_eff) && estimator_eff %in% c("WLSMV", "ULSMV")) {
    cli::cli_abort(c(
      "!" = "{.code missing = \"fiml\"} is not supported with \\
             {.code estimator = \"{estimator_eff}\"}.",
      "i" = "FIML requires a full-information estimator; \\
             {estimator_eff} is limited-information WLS.",
      "i" = "Use {.code missing = \"fiml\"} with {.code estimator = \"ML\"} \\
             or {.code estimator = \"MLR\"}, or switch to \\
             {.code missing = \"listwise\"}."
    ))
  }
  invisible(NULL)
}

# Detect whether x looks like a correlation matrix (not raw data).
# Heuristic: numeric, square, symmetric, and unit diagonal.
# A false positive on raw data is effectively impossible (requires n x n
# numeric data where n == p and all diagonal elements are exactly 1.0).
.is_cor_matrix <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) return(FALSE)
  if (nrow(x) != ncol(x)) return(FALSE)
  if (!all(abs(diag(x) - 1) < 1e-8)) return(FALSE)
  if (!isSymmetric(unname(x), tol = 1e-8)) return(FALSE)
  TRUE
}

# Validate and normalise a user-supplied correlation matrix. Returns the
# matrix (possibly with synthesised dimnames) or errors with a specific message.
# Non-positive-definite matrices warn and pass through — the engine will error
# naturally if truly degenerate; auto-smoothing would silently alter user input.
.validate_cor_matrix <- function(R) {
  if (!is.matrix(R) || !is.numeric(R)) {
    cli::cli_abort(c(
      "!" = "Supplied correlation matrix must be a numeric matrix.",
      "i" = "Got {.cls {class(R)}}."
    ))
  }
  p <- nrow(R)
  if (p != ncol(R)) {
    cli::cli_abort(c(
      "!" = "Correlation matrix must be square ({p} rows but {ncol(R)} columns)."
    ))
  }
  if (anyNA(R)) {
    cli::cli_abort(c(
      "!" = "Correlation matrix contains {sum(is.na(R))} NA value{?s}.",
      "i" = "Supply a complete (no-NA) correlation matrix."
    ))
  }
  if (!all(is.finite(R))) {
    cli::cli_abort(c(
      "!" = "Correlation matrix contains non-finite values (Inf or NaN)."
    ))
  }
  if (!isSymmetric(unname(R), tol = 1e-8)) {
    cli::cli_abort(c(
      "!" = "Correlation matrix is not symmetric.",
      "i" = "Maximum asymmetry: {max(abs(R - t(R)))}.",
      "i" = "Use {.code R <- (R + t(R)) / 2} to force symmetry."
    ))
  }
  if (!all(abs(diag(R) - 1) < 1e-8)) {
    bad_diag <- which(abs(diag(R) - 1) >= 1e-8)
    cli::cli_abort(c(
      "!" = "Correlation matrix diagonal must be all 1s.",
      "x" = "{length(bad_diag)} diagonal element{?s} differ from 1: \\
             {.val {round(diag(R)[bad_diag], 6)}} (position{?s} {bad_diag}).",
      "i" = "Supply a correlation matrix, not a covariance matrix."
    ))
  }
  off_diag <- R[row(R) != col(R)]
  if (any(abs(off_diag) > 1 + 1e-8)) {
    worst <- max(abs(off_diag))
    cli::cli_abort(c(
      "!" = "Correlation matrix has off-diagonal |r| > 1 (max = {round(worst, 6)}).",
      "i" = "Check that your matrix is a valid correlation matrix."
    ))
  }
  # Synthesise dimnames if absent so loadings/labels work downstream.
  if (is.null(rownames(R))) {
    rn <- paste0("V", seq_len(p))
    rownames(R) <- rn
    colnames(R) <- rn
  }
  # Non-PD: warn and pass through; engine errors naturally if truly degenerate.
  min_eig <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig <= 0) {
    cli::cli_warn(c(
      "!" = "Supplied correlation matrix is not positive definite \\
             (min eigenvalue = {round(min_eig, 4)}).",
      "i" = "PCA/EFA results may be unreliable. Consider regularising \\
             the matrix (e.g., {.fn psych::cor.smooth}) before calling \\
             {.fn ackwards}."
    ))
  }
  R
}

# Validate that x is a well-formed ackwards object (used in tests).
validate_ackwards <- function(x) {
  required <- c(
    "call", "engine", "rotation", "cor", "n_obs", "k_max",
    "seed", "pkg_version", "levels", "edges", "lineage",
    "scores", "fits", "r", "data", "meta", "prune"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    cli::cli_abort("ackwards object missing fields: {.field {missing}}")
  }
  if (!inherits(x, "ackwards")) {
    cli::cli_abort("Object does not have class {.cls ackwards}")
  }
  invisible(x)
}
