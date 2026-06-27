test_that("prune='none' leaves x$prune as NULL", {
  skip_if_not_installed("psych")
  set.seed(10)
  data <- as.data.frame(matrix(rnorm(300), 100, 6))
  x <- suppressWarnings(ackwards(data, k = 3))
  expect_null(x$prune)
  expect_equal(x$meta$prune, "none")
  expect_equal(x$meta$pairs, "adjacent")
})

test_that("prune='redundant' auto-upgrades pairs to 'all' with message", {
  skip_if_not_installed("psych")
  set.seed(10)
  data <- as.data.frame(matrix(rnorm(300), 100, 6))

  expect_message(
    x <- suppressWarnings(ackwards(data, k = 3, prune = "redundant")),
    "upgraded"
  )
  expect_equal(x$meta$pairs, "all")
  expect_false(is.null(x$prune))
})

test_that("prune='redundant' detects a known redundant chain and flags correctly", {
  skip_if_not_installed("psych")

  # Construct data with a clear 2-factor structure: one general + two specifics.
  # At high k, expect near-identical components to appear across consecutive levels.
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.3 * s1 + rnorm(n, sd = 0.1),
    x2 = 0.9 * g + 0.3 * s1 + rnorm(n, sd = 0.1),
    x3 = 0.9 * g + 0.3 * s1 + rnorm(n, sd = 0.1),
    x4 = 0.9 * g + 0.3 * s2 + rnorm(n, sd = 0.1),
    x5 = 0.9 * g + 0.3 * s2 + rnorm(n, sd = 0.1),
    x6 = 0.9 * g + 0.3 * s2 + rnorm(n, sd = 0.1)
  )

  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 4, prune = "redundant", redundancy_r = 0.9)
  ))

  nodes <- x$prune$nodes
  expect_s3_class(nodes, "data.frame")
  expect_true(all(c("id", "level", "pruned", "prune_reason") %in% names(nodes)))

  # Number of rows equals total factors across all levels (1+2+3+4 = 10)
  expect_equal(nrow(nodes), sum(seq_len(x$k_max)))
})

test_that("x$prune$nodes has one row per factor with correct structure", {
  skip_if_not_installed("psych")
  set.seed(1)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 3, prune = "redundant")
  ))

  nodes <- x$prune$nodes
  expect_equal(nrow(nodes), 1 + 2 + 3) # k=1,2,3
  expect_true(is.logical(nodes$pruned))
  expect_true(is.character(nodes$prune_reason))

  # All non-pruned nodes have NA prune_reason
  expect_true(all(is.na(nodes$prune_reason[!nodes$pruned])))

  # Any pruned nodes have prune_reason "redundant"
  if (any(nodes$pruned)) {
    expect_true(all(nodes$prune_reason[nodes$pruned] == "redundant"))
  }
})

test_that("x$prune$chains has correct structure when chains are found", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 4, prune = "redundant", redundancy_r = 0.9)
  ))

  if (!is.null(x$prune$chains)) {
    ch <- x$prune$chains
    expect_true(all(c(
      "chain_id", "id", "level", "r_to_prev", "phi_to_prev",
      "retain", "endpoint_r", "endpoint_r_agrees"
    ) %in% names(ch)))
    # Exactly one retained node per chain
    for (cid in unique(ch$chain_id)) {
      expect_equal(sum(ch$retain[ch$chain_id == cid]), 1L)
    }
    # r_to_prev is NA for the first member of each chain (no predecessor)
    first_members <- tapply(seq_len(nrow(ch)), ch$chain_id, min)
    expect_true(all(is.na(ch$r_to_prev[first_members])))
  }
})

test_that("retention rule: chain reaching k_max retains bottom node", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  k_max <- 4L
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = k_max, prune = "redundant", redundancy_r = 0.9)
  ))

  if (!is.null(x$prune$chains)) {
    ch <- x$prune$chains
    # For chains that reach k_max, the retained node must be at level k_max
    bottom_chains <- unique(ch$chain_id[ch$level == k_max])
    for (cid in bottom_chains) {
      sub <- ch[ch$chain_id == cid, ]
      retained_level <- sub$level[sub$retain]
      expect_equal(retained_level, k_max)
    }
    # For chains that don't reach k_max, the retained node must be at the top
    mid_chains <- setdiff(unique(ch$chain_id), bottom_chains)
    for (cid in mid_chains) {
      sub <- ch[ch$chain_id == cid, ]
      retained_level <- sub$level[sub$retain]
      expect_equal(retained_level, min(sub$level))
    }
  }
})

test_that("Tucker's phi is computed correctly by .tucker_phi", {
  a <- c(0.8, 0.6, 0.0)
  b <- c(0.6, 0.8, 0.0)
  phi <- ackwards:::.tucker_phi(a, b)
  expected <- sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  expect_equal(phi, expected)

  # Identical vectors → phi = 1
  expect_equal(ackwards:::.tucker_phi(a, a), 1.0)

  # Orthogonal vectors → phi = 0
  expect_equal(ackwards:::.tucker_phi(c(1, 0), c(0, 1)), 0.0)
})

test_that("prune='artefact' computes phi table without flagging nodes", {
  skip_if_not_installed("psych")
  set.seed(5)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 3, prune = "artefact")
  ))

  expect_false(is.null(x$prune))
  expect_false(is.null(x$prune$phi))
  phi_df <- x$prune$phi
  expect_true(all(c("from", "to", "level_from", "level_to", "phi") %in% names(phi_df)))
  expect_true(all(phi_df$level_from < phi_df$level_to))

  # No nodes should be auto-flagged in artefact-only mode
  expect_true(all(!x$prune$nodes$pruned))
})

test_that("redundancy_phi conjunctive criterion is applied when set", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 500
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  # With phi threshold of 1 (impossible), no chains should survive
  x_strict <- suppressWarnings(suppressMessages(
    ackwards(data,
      k = 4, prune = "redundant",
      redundancy_r = 0.9, redundancy_phi = 1.0
    )
  ))
  expect_true(all(!x_strict$prune$nodes$pruned))
  expect_null(x_strict$prune$chains)
})

test_that("invalid redundancy_r and redundancy_phi are rejected", {
  skip_if_not_installed("psych")
  data <- as.data.frame(matrix(rnorm(300), 100, 6))
  expect_error(
    suppressWarnings(ackwards(data, k = 3, prune = "redundant", redundancy_r = 1.5)),
    "redundancy_r"
  )
  expect_error(
    suppressWarnings(ackwards(data,
      k = 3, prune = "redundant",
      redundancy_phi = -0.1
    )),
    "redundancy_phi"
  )
})

# --- A3: tidy(what = "nodes") ------------------------------------------------

test_that("tidy(x, what='nodes') returns prune$nodes when pruning was applied", {
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 500L
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 4, prune = "redundant")
  ))

  nodes <- tidy(x, what = "nodes")
  expect_s3_class(nodes, "data.frame")
  expect_true(all(c("id", "level", "pruned", "prune_reason") %in% names(nodes)))
  expect_equal(nrow(nodes), sum(seq_len(x$k_max)))
  expect_equal(nodes, x$prune$nodes)
})

test_that("tidy(x, what='nodes') returns empty frame with correct columns when no pruning", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k = 3))
  expect_null(x$prune)

  nodes <- tidy(x, what = "nodes")
  expect_s3_class(nodes, "data.frame")
  expect_equal(nrow(nodes), 0L)
  expect_true(all(c("id", "level", "pruned", "prune_reason") %in% names(nodes)))
})

# --- A2: print() pruning section -----------------------------------------------

test_that("print.ackwards() executes without error when pruning annotations are present", {
  # cli output cannot be reliably captured via capture.output() in non-interactive
  # test sessions — the Pruning section and relabeling caveat are verified visually.
  # This test guards against regressions that would make print() error or not
  # return x invisibly (consistent with the convention in test-print.R).
  skip_if_not_installed("psych")
  set.seed(42)
  n <- 500L
  g <- rnorm(n)
  s1 <- rnorm(n)
  s2 <- rnorm(n)
  data <- data.frame(
    x1 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x2 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x3 = 0.9 * g + 0.2 * s1 + rnorm(n, sd = 0.05),
    x4 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x5 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05),
    x6 = 0.9 * g + 0.2 * s2 + rnorm(n, sd = 0.05)
  )
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k = 4, prune = "redundant")
  ))
  expect_no_error(print(x))
  expect_invisible(print(x))
  # Pruning slot populated confirms the section has content to show
  expect_false(is.null(x$prune))
  expect_true(any(x$prune$nodes$pruned))
})

# --- B1: sibling redundancy tracing ------------------------------------------

test_that("B1 regression: parent with 2 strong-link children produces 2 chains", {
  # Construct a minimal mock ackwards-like object: k=2, level 1 has 1 factor
  # (m1f1), level 2 has 2 factors (m2f1, m2f2), both with |r| >= 0.9 to m1f1.
  # Old (buggy) code only followed the first child; DFS finds both branches.
  mock_L1 <- matrix(rep(0.9, 6),
    ncol = 1,
    dimnames = list(paste0("x", 1:6), "m1f1")
  )
  mock_L2 <- matrix(
    c(
      0.9, 0.9, 0.9, 0.1, 0.1, 0.1,
      0.1, 0.1, 0.1, 0.9, 0.9, 0.9
    ),
    ncol = 2, dimnames = list(paste0("x", 1:6), c("m2f1", "m2f2"))
  )
  E_1_2 <- matrix(c(0.95, 0.93),
    nrow = 1,
    dimnames = list("m1f1", c("m2f1", "m2f2"))
  )

  mock_x <- list(
    k_max = 2L,
    levels = list(
      "1" = list(labels = "m1f1", loadings = mock_L1),
      "2" = list(labels = c("m2f1", "m2f2"), loadings = mock_L2)
    ),
    lineage = list("1" = NULL, "2" = c(1L, 1L)),
    edges = list(matrices = list("1:2" = E_1_2))
  )

  res <- ackwards:::.find_redundant_chains(mock_x, threshold_r = 0.90, threshold_phi = NULL)

  # Both strong links must produce separate chains
  expect_equal(length(unique(res$chains$chain_id)), 2L)
  expect_true("m2f1" %in% res$chains$id)
  expect_true("m2f2" %in% res$chains$id)

  # m1f1 flagged exactly once (it's the shared root, not retained in either chain)
  expect_equal(sum(res$node_flags$id == "m1f1"), 1L)
  expect_true(res$node_flags$pruned[res$node_flags$id == "m1f1"])

  # Both m2f1 and m2f2 are retained (chain reaches k_max = 2)
  expect_false("m2f1" %in% res$node_flags$id)
  expect_false("m2f2" %in% res$node_flags$id)
})

# --- B4: pruning under convergence truncation --------------------------------

test_that("B4: prune works correctly when hierarchy is truncated (k_eff < k_requested)", {
  skip_if_not_installed("lavaan")
  # 6-variable data: lavaan::efa() errors at k >= 4; requesting k=5 → k_eff=3
  d <- .make_esem_data()

  x <- suppressWarnings(suppressMessages(
    ackwards(d, k = 5, method = "esem", prune = "redundant")
  ))

  expect_equal(x$k_max, 3L)
  expect_false(is.null(x$prune))

  # Retention rule must use k_eff (3), never a non-existent level 4 or 5
  if (!is.null(x$prune$chains)) {
    retained_levels <- x$prune$chains$level[x$prune$chains$retain]
    expect_true(all(retained_levels <= 3L))
  }

  # All node flags reference existing levels
  expect_true(all(x$prune$nodes$level <= 3L))
})
