# Internal utilities — not exported

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
    is_int_like(x) && length(unique(x[!is.na(x)])) <= max_levels
  }, logical(1))
  any(cols)
}

# Bipartite matching: assign each factor in level b its primary parent in level a.
# Returns an integer vector of length ncol(E): parent index in level a for each
# factor in level b.
# Uses clue::solve_LSAP when available; falls back to greedy which.max.
match_parents <- function(E) {
  abs_E <- abs(E)
  n_a <- nrow(abs_E)
  n_b <- ncol(abs_E)
  if (rlang::is_installed("clue")) {
    # solve_LSAP minimises cost; negate for maximisation
    pad <- max(n_a, n_b)
    cost <- matrix(0, pad, pad)
    cost[seq_len(n_a), seq_len(n_b)] <- -abs_E
    assignment <- as.integer(clue::solve_LSAP(cost))
    assignment[seq_len(n_b)]
  } else {
    apply(abs_E, 2, which.max)
  }
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
align_signs <- function(loadings_list, edges_list, lineage) {
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
    parents <- lineage[[k]]  # parent index in level k-1 for each factor in level k
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

# Validate that x is a well-formed ackwards object (used in tests).
validate_ackwards <- function(x) {
  required <- c(
    "call", "method", "rotation", "cor_type", "n_obs", "k_max",
    "seed", "pkg_version", "levels", "edges", "lineage",
    "scores", "fits", "r", "data", "meta"
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
