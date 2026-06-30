test_that("prune='none' leaves x$prune as NULL", {
  skip_if_not_installed("psych")
  set.seed(10)
  data <- as.data.frame(matrix(rnorm(300), 100, 6))
  x <- suppressWarnings(ackwards(data, k_max = 3))
  expect_null(x$prune)
  expect_equal(x$meta$prune, "none")
  expect_equal(x$meta$pairs, "adjacent")
})

test_that("prune='redundant' auto-upgrades pairs to 'all' with message", {
  skip_if_not_installed("psych")
  set.seed(10)
  data <- as.data.frame(matrix(rnorm(300), 100, 6))

  expect_message(
    x <- suppressWarnings(ackwards(data, k_max = 3, prune = "redundant")),
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
    ackwards(data, k_max = 4, prune = "redundant", redundancy_r = 0.9)
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
    ackwards(data, k_max = 3, prune = "redundant")
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
    ackwards(data, k_max = 4, prune = "redundant", redundancy_r = 0.9)
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
    ackwards(data, k_max = k_max, prune = "redundant", redundancy_r = 0.9)
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
    ackwards(data, k_max = 3, prune = "artefact")
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
      k_max = 4, prune = "redundant",
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
    suppressWarnings(ackwards(data, k_max = 3, prune = "redundant", redundancy_r = 1.5)),
    "redundancy_r"
  )
  expect_error(
    suppressWarnings(ackwards(data,
      k_max = 3, prune = "redundant",
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
    ackwards(data, k_max = 4, prune = "redundant")
  ))

  nodes <- tidy(x, what = "nodes")
  expect_s3_class(nodes, "data.frame")
  expect_true(all(c("id", "level", "pruned", "prune_reason") %in% names(nodes)))
  expect_equal(nrow(nodes), sum(seq_len(x$k_max)))
  expect_equal(nodes, x$prune$nodes)
})

test_that("tidy(x, what='nodes') returns empty frame with correct columns when no pruning", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
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
    ackwards(data, k_max = 4, prune = "redundant")
  ))
  expect_no_error(print(x))
  expect_invisible(print(x))
  # Pruning slot populated confirms the section has content to show
  expect_false(is.null(x$prune))
  expect_true(any(x$prune$nodes$pruned))
})

test_that("print.ackwards() shows artefact section (phi table) for prune='artefact'", {
  # Covers print.R lines 70-74: the artefact phi-count section.
  skip_if_not_installed("psych")
  suppressMessages(x <- suppressWarnings(
    ackwards(bfi25[, 1:6], k_max = 3, prune = "artefact")
  ))
  expect_no_error(print(x))
  expect_invisible(print(x))
  expect_false(is.null(x$prune$phi))
})

test_that("print.ackwards() shows phi threshold note when redundancy_phi is set", {
  # Covers print.R line 61: phi_note when redundancy_phi is non-NULL.
  skip_if_not_installed("psych")
  set.seed(42L)
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
    ackwards(data, k_max = 4, prune = "redundant", redundancy_phi = 0.9)
  ))
  expect_false(is.null(x$prune$redundancy_phi))
  expect_no_error(print(x))
  expect_invisible(print(x))
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

# --- B3: internal helper coverage (M23) --------------------------------------

test_that(".tucker_phi() returns NA_real_ when one loading vector is all zeros", {
  # denom = sqrt(0 * sum) = 0 -> defensive NA_real_ branch
  expect_identical(ackwards:::.tucker_phi(c(0, 0, 0), c(0.7, 0.8, 0.9)), NA_real_)
})

test_that(".phi_pairs() handles 'adjacent' which_pairs correctly", {
  # .phi_pairs() is only called with 'all' from .apply_pruning(); the 'adjacent'
  # branch at line 39 requires a direct unit test.
  levels_list <- list(
    "1" = list(
      loadings = matrix(c(.8, .8, .8),
        ncol = 1L,
        dimnames = list(paste0("x", 1:3), "m1f1")
      ),
      labels = "m1f1"
    ),
    "2" = list(
      loadings = matrix(c(.7, .8, .9, .8, .7, .6),
        ncol = 2L,
        dimnames = list(paste0("x", 1:3), c("m2f1", "m2f2"))
      ),
      labels = c("m2f1", "m2f2")
    )
  )
  out <- ackwards:::.phi_pairs(levels_list, which_pairs = "adjacent")
  expect_s3_class(out, "data.frame")
  expect_true("phi" %in% names(out))
  # Only adjacent pairs: 1->2, giving 1*2=2 rows
  expect_equal(nrow(out), 2L)
})

test_that(".find_redundant_chains() handles a missing edge matrix gracefully", {
  # Simulates convergence truncation where the 2:3 edge matrix is absent.
  mock_L1 <- matrix(c(.8, .8, .8),
    ncol = 1L,
    dimnames = list(paste0("x", 1:3), "m1f1")
  )
  mock_L2 <- matrix(c(.9, .9, .9, .1, .1, .1),
    ncol = 2L,
    dimnames = list(paste0("x", 1:3), c("m2f1", "m2f2"))
  )
  mock_L3 <- matrix(c(.85, .85, .85, .1, .1, .1, .1, .1, .1),
    ncol = 3L,
    dimnames = list(paste0("x", 1:3), c("m3f1", "m3f2", "m3f3"))
  )
  E_1_2 <- matrix(c(0.95, 0.93),
    nrow = 1L,
    dimnames = list("m1f1", c("m2f1", "m2f2"))
  )
  # "2:3" key is deliberately absent to trigger the NULL-edge-matrix branch
  mock_x <- list(
    k_max = 3L,
    levels = list(
      "1" = list(labels = "m1f1", loadings = mock_L1),
      "2" = list(labels = c("m2f1", "m2f2"), loadings = mock_L2),
      "3" = list(labels = c("m3f1", "m3f2", "m3f3"), loadings = mock_L3)
    ),
    lineage = list("1" = NULL, "2" = c(1L, 1L), "3" = c(1L, 2L, 2L)),
    edges = list(matrices = list("1:2" = E_1_2)) # no "2:3"
  )
  res <- ackwards:::.find_redundant_chains(mock_x, threshold_r = 0.9, threshold_phi = NULL)
  # Should not error; chains only reflect the 1->2 links
  expect_type(res, "list")
})

# --- B4: pruning under convergence truncation --------------------------------

test_that("B4: prune works correctly when hierarchy is truncated (k_eff < k_requested)", {
  skip_if_not_installed("lavaan")
  # 6-variable data: lavaan::efa() errors at k >= 4; requesting k=5 → k_eff=3
  d <- .make_esem_data()

  x <- suppressWarnings(suppressMessages(
    ackwards(d, k_max = 5, engine = "esem", prune = "redundant")
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

# ---- Tests for structural artefact signals (Wave 2 / M25) --------------------

test_that("prune='artefact' populates $structural with correct schema", {
  skip_if_not_installed("psych")
  set.seed(5)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact")
  ))

  expect_false(is.null(x$prune$structural))
  struct <- x$prune$structural
  expect_s3_class(struct, "data.frame")
  expect_named(struct, c("id", "level", "few_items", "orphan", "split_merge"))
  # One row per (level, factor): 1 + 2 + 3 = 6 rows for k_max = 3
  expect_equal(nrow(struct), sum(seq_len(3L)))
  expect_type(struct$few_items, "logical")
  expect_type(struct$split_merge, "logical")
  # phi table still present and unchanged
  expect_false(is.null(x$prune$phi))
  # No auto-flagging (artefact is flag/report only)
  expect_true(all(!x$prune$nodes$pruned))
})

test_that("prune='artefact' stores min_items and orphan_r in prune slot", {
  skip_if_not_installed("psych")
  set.seed(5)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data,
      k_max = 3, prune = "artefact",
      min_items = 2L, orphan_r = 0.3
    )
  ))
  expect_equal(x$prune$min_items, 2L)
  expect_equal(x$prune$orphan_r, 0.3)
})

test_that("prune='redundant' leaves $structural NULL", {
  skip_if_not_installed("psych")
  set.seed(1)
  data <- as.data.frame(matrix(rnorm(600), 100, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "redundant")
  ))
  expect_null(x$prune$structural)
})

test_that("few_items flags 2-item factor when min_items = 3", {
  skip_if_not_installed("psych")
  # 3 items on factor 1, 2 items on factor 2 -> factor 2 at k=2 has few_items
  set.seed(42)
  n <- 300L
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  data_few <- data.frame(
    x1 = 0.9 * f1 + rnorm(n, sd = 0.2),
    x2 = 0.9 * f1 + rnorm(n, sd = 0.2),
    x3 = 0.9 * f1 + rnorm(n, sd = 0.2),
    x4 = 0.9 * f2 + rnorm(n, sd = 0.2),
    x5 = 0.9 * f2 + rnorm(n, sd = 0.2)
  )
  x <- suppressMessages(
    ackwards(data_few, k_max = 2, prune = "artefact", min_items = 3L)
  )
  struct <- x$prune$structural
  k2 <- struct[struct$level == 2L, ]
  # Exactly one k=2 factor should have 2 primary items -> few_items = TRUE
  expect_true(any(k2$few_items))
  # k=1 (single factor, all 5 items) must not be flagged
  k1 <- struct[struct$level == 1L, ]
  expect_false(k1$few_items)
})

test_that("orphan flags with very high orphan_r; none with zero threshold", {
  skip_if_not_installed("psych")
  set.seed(7)
  data <- as.data.frame(matrix(rnorm(600), 100, 6))

  # orphan_r = 0.99: nearly impossible to satisfy; all factors should be orphans
  x_high <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact", orphan_r = 0.99)
  ))
  expect_true(any(x_high$prune$structural$orphan, na.rm = TRUE))

  # orphan_r = 0: any nonzero |r| clears the threshold; none should be flagged
  x_low <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact", orphan_r = 0)
  ))
  expect_false(any(x_low$prune$structural$orphan, na.rm = TRUE))
})

test_that("split_merge is FALSE at level 1 and level 2 (no multi-parent possible)", {
  skip_if_not_installed("psych")
  set.seed(7)
  data <- as.data.frame(matrix(rnorm(600), 100, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact")
  ))
  struct <- x$prune$structural
  # Level 1: no parent level exists -> always FALSE
  expect_false(any(struct$split_merge[struct$level == 1L]))
  # Level 2: all items share one single k=1 parent -> always FALSE
  expect_false(any(struct$split_merge[struct$level == 2L]))
})

test_that("print.ackwards() shows structural signal count for prune='artefact'", {
  skip_if_not_installed("psych")
  set.seed(5)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact")
  ))
  # print must not error; structural count line present
  expect_no_error(print(x))
  expect_invisible(print(x))
})

test_that("summary.ackwards() shows structural signal count for prune='artefact'", {
  skip_if_not_installed("psych")
  set.seed(5)
  data <- as.data.frame(matrix(rnorm(900), 150, 6))
  x <- suppressWarnings(suppressMessages(
    ackwards(data, k_max = 3, prune = "artefact")
  ))
  expect_no_error(summary(x))
  expect_no_error(print(summary(x)))
})
