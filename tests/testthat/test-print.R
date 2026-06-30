test_that("print.ackwards runs without error and returns x invisibly", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  # cli writes to stderr in non-interactive mode; just verify no error + invisible
  expect_no_error(print(x))
  expect_invisible(print(x))
})

test_that("tidy(x, what = 'edges') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  ed <- generics::tidy(x, what = "edges")
  expect_s3_class(ed, "data.frame")
  expect_true(all(c("from", "to", "level_from", "level_to", "r", "is_primary", "above_cut") %in% names(ed)))
  # 3 levels → 2 adjacent pairs → (1*2) + (2*3) = 8 edges
  expect_equal(nrow(ed), 8L)
  expect_true(all(abs(ed$r) <= 1 + 1e-9))
})

test_that("tidy(x, what = 'loadings') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 2))
  ld <- generics::tidy(x, what = "loadings")
  expect_s3_class(ld, "data.frame")
  expect_true(all(c("level", "factor", "item", "loading") %in% names(ld)))
  # k=1 has 25 rows, k=2 has 50 rows → 75 total
  expect_equal(nrow(ld), 75L)
})

test_that("tidy(x, what = 'variance') returns expected structure", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  vr <- generics::tidy(x, what = "variance")
  expect_s3_class(vr, "data.frame")
  expect_true(all(c("level", "factor", "variance_pct", "cumulative_pct") %in% names(vr)))
  # k=1 has 1 row, k=2 has 2, k=3 has 3 → 6 total
  expect_equal(nrow(vr), 6L)
})

test_that("glance.ackwards returns a one-row data frame with expected columns", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:25], k_max = 3))
  gl <- generics::glance(x)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1L)
  expect_true(all(c(
    "engine", "rotation", "cor", "k_max", "n_obs",
    "deepest_converged", "n_edges"
  ) %in% names(gl)))
  expect_equal(gl$rotation, "varimax")
  expect_equal(gl$k_max, 3L)
  expect_equal(gl$n_obs, 2800L)
  expect_equal(gl$n_edges, 8L)
})

# ── summary.ackwards ──────────────────────────────────────────────────────────

test_that("summary.ackwards returns a summary_ackwards object with expected fields", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:10], k_max = 3))
  s <- summary(x)
  expect_s3_class(s, "summary_ackwards")
  expect_true(all(c(
    "call", "engine", "rotation", "cor", "n_obs", "k_max",
    "variance", "fit", "lineage", "prune"
  ) %in% names(s)))
  expect_equal(s$engine, "pca")
  expect_equal(s$k_max, 3L)
})

test_that("print.summary_ackwards runs without error and returns invisibly", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:10], k_max = 3))
  s <- summary(x)
  expect_no_error(print(s))
  expect_invisible(print(s))
})

test_that("summary fit table carries eigenvalues for PCA", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:10], k_max = 2))
  s <- summary(x)
  # PCA fit indices are named "eigenvalue.<label>"; verify they exist
  fit_k1 <- s$fit[s$fit$level == 1L, ]
  expect_true(any(startsWith(fit_k1$index, "eigenvalue.")))
  expect_true(all(is.finite(fit_k1$value)))
})

test_that("summary lineage: correct structure, ordering, and content", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:10], k_max = 3))
  s <- summary(x)
  lin <- s$lineage
  expect_s3_class(lin, "data.frame")
  expect_true(all(c("parent", "children") %in% names(lin)))
  # k=3 hierarchy has level-1 and level-2 parents: m1f1 + m2f1 + m2f2 = 3 rows
  expect_equal(nrow(lin), 3L)
  # Parents appear in ascending level order (m1fX before m2fX)
  expect_equal(lin$parent[1L], "m1f1")
  # m1f1's children are the two k=2 factors (both primary children of m1f1)
  m1_row <- lin[lin$parent == "m1f1", ]
  m1_children <- trimws(strsplit(m1_row$children, ",")[[1L]])
  expect_setequal(m1_children, c("m2f1", "m2f2"))
})

test_that("summary.ackwards shows pruning info when prune = 'redundant'", {
  skip_if_not_installed("psych")
  suppressWarnings(x <- ackwards(psych::bfi[, 1:10], k_max = 4, prune = "redundant"))
  s <- summary(x)
  expect_false(is.null(s$prune))
  expect_equal(s$prune$rules, "redundant")
  expect_true(!is.null(s$prune$redundancy_r))
})

test_that("summary prune = 'artefact' carries rules and no spurious redundant info", {
  skip_if_not_installed("psych")
  suppressMessages(x <- suppressWarnings(
    ackwards(psych::bfi[, 1:10], k_max = 3, prune = "artefact")
  ))
  s <- summary(x)
  expect_false(is.null(s$prune))
  expect_equal(s$prune$rules, "artefact")
  # artefact-only: no node is flagged as redundant
  expect_equal(length(s$prune$redundant), 0L)
  # phi table exists
  expect_true(!is.null(s$prune$artefact_n) && s$prune$artefact_n > 0L)
  # print does not error (was showing spurious "Redundant: 0 nodes" before fix)
  expect_no_error(suppressMessages(print(s)))
})

test_that("summary prune = c('redundant','artefact') shows both sections", {
  skip_if_not_installed("psych")
  suppressMessages(x <- suppressWarnings(
    ackwards(psych::bfi[, 1:10], k_max = 4, prune = c("redundant", "artefact"))
  ))
  s <- summary(x)
  expect_setequal(s$prune$rules, c("redundant", "artefact"))
  expect_true(!is.null(s$prune$artefact_n))
  expect_no_error(suppressMessages(print(s)))
})

test_that("summary.ackwards handles truncated ESEM hierarchy", {
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  # k=5 requested but p=6 caps lavaan at 3 -> k_eff=3
  x <- suppressWarnings(suppressMessages(ackwards(d, k_max = 5, engine = "esem")))
  s <- summary(x)
  # Summary reflects k_eff, not the requested k
  expect_equal(s$k_max, x$k_max)
  expect_equal(nrow(s$variance), sum(seq_len(x$k_max)))
  expect_no_error(suppressMessages(print(s)))
})

test_that("summary.ackwards works with polychoric basis", {
  skip_if_not_installed("psych")
  d <- .make_ordinal_data()
  suppressMessages(x <- suppressWarnings(
    ackwards(d, k_max = 2, cor = "polychoric")
  ))
  s <- summary(x)
  expect_s3_class(s, "summary_ackwards")
  expect_equal(s$cor, "polychoric")
  expect_no_error(suppressMessages(print(s)))
})

test_that("summary.ackwards works for all three engines", {
  skip_if_not_installed("psych")
  skip_if_not_installed("lavaan")
  d <- .make_esem_data()
  x_pca <- suppressWarnings(ackwards(d, k_max = 2))
  x_efa <- suppressWarnings(ackwards(d, k_max = 2, engine = "efa"))
  x_esem <- suppressWarnings(ackwards(d, k_max = 2, engine = "esem"))
  expect_no_error(summary(x_pca))
  expect_no_error(summary(x_efa))
  expect_no_error(summary(x_esem))
})

# ── Branches added for M23 coverage hardening ──────────────────────────────────

test_that("print.summary_ackwards EFA path shows chi/dof indices", {
  # Covers summary.R lines 102 and 112: EFA idx_show and chi/dof formatting.
  skip_if_not_installed("psych")
  d <- .make_esem_data()
  x_efa <- suppressWarnings(ackwards(d, k_max = 2, engine = "efa"))
  s <- summary(x_efa)
  expect_no_error(suppressMessages(print(s)))
  expect_invisible(suppressMessages(print(s)))
})

test_that("print.summary_ackwards shows '(none)' when redundant pruning finds no flags", {
  # Covers summary.R line 146: redundant pruning enabled but p$redundant is empty.
  # redundancy_r = 1.0 means only |r| >= 1.0 qualifies — impossible → 0 flags.
  skip_if_not_installed("psych")
  set.seed(7L)
  d <- as.data.frame(matrix(rnorm(300L), 50L, 6L))
  x <- suppressWarnings(suppressMessages(
    ackwards(d, k_max = 3, prune = "redundant", redundancy_r = 0.9999)
  ))
  s <- summary(x)
  # No nodes flagged
  expect_equal(length(s$prune$redundant), 0L)
  expect_no_error(suppressMessages(print(s)))
})

test_that("print.summary_ackwards shows phi threshold note when redundancy_phi set", {
  # Covers summary.R line 150: phi_note when redundancy_phi is non-NULL.
  skip_if_not_installed("psych")
  set.seed(42L)
  d <- as.data.frame(matrix(rnorm(300L), 50L, 6L))
  x <- suppressWarnings(suppressMessages(
    ackwards(d, k_max = 3, prune = "redundant", redundancy_phi = 0.9)
  ))
  s <- summary(x)
  expect_false(is.null(s$prune$redundancy_phi))
  expect_no_error(suppressMessages(print(s)))
})

test_that("print.summary_ackwards shows empty lineage message for a no-edge object", {
  # Covers summary.R lines 131, 199-204: both the empty-lineage print line and
  # the .summary_lineage() early return when adj_primary is empty.
  skip_if_not_installed("psych")
  x <- suppressWarnings(ackwards(bfi25[, 1:6], k_max = 3))
  # Patch out all edges so .summary_lineage() returns an empty data frame
  x$edges$tidy <- x$edges$tidy[0L, , drop = FALSE]
  s <- summary(x)
  expect_equal(nrow(s$lineage), 0L)
  expect_no_error(suppressMessages(print(s)))
})
