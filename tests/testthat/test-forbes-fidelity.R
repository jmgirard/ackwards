# Forbes (2023) fidelity guard (M44).
#
# CLAUDE.md's baseline contract: the default output must reproduce Forbes's
# examples exactly. The fixture holds the Spearman correlation matrices of the
# paper's three simulation studies (regenerated from Forbes's public OSF
# script, set.seed(123)) plus expected outputs computed with her reference
# implementation ("ExtendedBassAckwards functions with annotation.R",
# https://osf.io/pcwm8/) -- see attr(fixture, "provenance"). The tests below
# run only ackwards on the fixed inputs and compare against those expected
# values, so no external code or network is needed at test time.
#
# Correspondence conventions (verified during the M44 feasibility study, where
# the 155-variable applied example also matched to 3.9e-14):
#   * Forbes's comp.corr is t(W_a) %*% R %*% W_b with no sign alignment; ours
#     is the same algebra with primary-parent alignment, so |values| must be
#     identical entrywise (W'RW = I for PCA makes her unstandardized products
#     equal our standardized edges).
#   * Her cong is psych::factor.congruence, which rounds to 2 decimals; ours
#     is exact Tucker's phi, so agreement is within 0.005.
#   * Her component labels are letter-level + index in psych::pca column
#     order ("a1"; "b1","b2"; ...); ours are m{k}f{j} in the same order.
#   * Her comp.corr list enumerates pairs as: for c in 2..K, for i in 1..c-1.

.forbes_fixture <- function() {
  readRDS(test_path("fixtures", "forbes2023_sims.rds"))
}

# Map a Forbes component label ("c3") to an ackwards one ("m3f3").
.forbes_to_ackwards <- function(lab) {
  lv <- match(substr(lab, 1L, 1L), letters)
  paste0("m", lv, "f", substr(lab, 2L, nchar(lab)))
}

# Forbes's redundancy "chase" on an ackwards object: from `node`, follow
# primary-parent links upward while |r| >= .9; return the topmost node
# reached (the node itself when its first upward link is below .9).
.chase <- function(x, node) {
  lv <- as.integer(sub("^m(\\d+)f\\d+$", "\\1", node))
  cur <- node
  while (lv > 1L) {
    E <- x$edges$matrices[[paste0(lv - 1L, ":", lv)]]
    j <- match(cur, colnames(E))
    pi <- which.max(abs(E[, j]))
    if (abs(E[pi, j]) < 0.9) break
    cur <- rownames(E)[pi]
    lv <- lv - 1L
  }
  cur
}

test_that("default output reproduces Forbes's simulation examples exactly", {
  skip_if_not_installed("psych")
  sims <- .forbes_fixture()

  for (nm in names(sims)) {
    sim <- sims[[nm]]
    suppressWarnings(suppressMessages(
      x <- cached(ackwards(sim$R, k_max = 4, pairs = "all"))
    ))

    # (1) Between-level correlations: |ours| == |hers| entrywise, all 6 pairs.
    idx <- 0L
    for (c2 in 2:4) {
      for (i in 1:(c2 - 1L)) {
        idx <- idx + 1L
        E_forbes <- sim$comp_corr[[idx]]
        E_ours <- x$edges$matrices[[paste0(i, ":", c2)]]
        expect_equal(
          abs(unname(E_ours)), abs(unname(E_forbes)),
          tolerance = 1e-12,
          label = paste0(nm, " |edges| ", i, ":", c2, " (ackwards)"),
          expected.label = "Forbes reference"
        )
      }
    }

    # (2) Loading congruence: |ours| == |hers| within her 2-dp rounding.
    phi_ours <- ackwards:::.phi_pairs(x$levels, "all")
    idx <- 0L
    for (c2 in 2:4) {
      for (i in 1:(c2 - 1L)) {
        idx <- idx + 1L
        C_forbes <- sim$cong[[idx]]
        sub <- phi_ours[phi_ours$level_from == i & phi_ours$level_to == c2, ]
        C_ours <- matrix(sub$phi, nrow = i, ncol = c2, byrow = TRUE)
        expect_lt(
          max(abs(abs(C_ours) - abs(unname(C_forbes)))),
          0.005 + 1e-12
        )
      }
    }

    # (3) Redundancy chase paths: for every component, following our
    # primary-parent links upward while |r| >= .9 must land on the same
    # component her ChaseCorrPaths reports ("X--null" = no move).
    for (entry in sim$corr_chase) {
      parts <- strsplit(entry, "--", fixed = TRUE)[[1L]]
      from <- .forbes_to_ackwards(parts[1L])
      expected_top <- if (parts[2L] == "null") from else .forbes_to_ackwards(parts[2L])
      expect_identical(
        .chase(x, from), expected_top,
        label = paste0(nm, " chase(", parts[1L], ") (ackwards)"),
        expected.label = paste0("Forbes '", entry, "'")
      )
    }
  }
})

test_that("prune('redundant') flags Forbes's Simulation 1 chains with her retention rule", {
  skip_if_not_installed("psych")
  sims <- .forbes_fixture()
  suppressWarnings(suppressMessages({
    x <- cached(ackwards(sims$sim1$R, k_max = 4, pairs = "all"))
    xp <- prune(x, "redundant")
  }))

  # Her chase found exactly three redundant links: c3--b2, d1--c1, d2--c2.
  # Under her retention rule: the b2-c3 chain stops short of k_max, keeping
  # the top (m2f2); the two chains reaching k_max keep their bottoms
  # (m4f1, m4f2). Flagged = the other chain members.
  flagged <- xp$prune$nodes$id[xp$prune$nodes$pruned]
  expect_setequal(flagged, c("m3f3", "m3f1", "m3f2"))

  ch <- xp$prune$chains
  expect_setequal(ch$id[ch$retain], c("m2f2", "m4f1", "m4f2"))
})

# ---------------------------------------------------------------------------
# AMH applied example (M53).
#
# Forbes's 155-variable "Assessing Mental Health" applied example, k = 10 (OSF
# pcwm8, CC-BY 4.0). The matrix is the exported dataset `forbes2023` (M54);
# fixtures/forbes2023_amh.rds holds only the expected comp_corr/cong computed
# with HER reference implementation from that same md5-pinned matrix (see
# data-raw/forbes2023.R and attr(fixture, "provenance")). As with the simulations,
# only ackwards() runs here; no Forbes code or network at test time.
#
# Redundancy chase: Forbes's ChaseCorrPaths uses the DIRECT (skip-level)
# correlation to a component at each ancestor level -- the criterion our
# prune("redundant") adopts by default since M53 (redundancy_criterion =
# "direct"). On this deep 10-level hierarchy it diverges from an adjacent-hop
# walk on 7 of 54 components (correlation is non-transitive), so this is exactly
# where the direct criterion earns its keep. amh$corr_chase holds her raw
# ChaseCorrPaths output ("X--Y"; "X--null" = no chase) for all 54 components,
# and .direct_chase() below reproduces it from our all-levels edges. See M53.
.amh_fixture <- function() {
  readRDS(test_path("fixtures", "forbes2023_amh.rds"))$amh
}

# Forbes's direct/skip-level chase on an ackwards object: from `node`, at each
# ancestor level take the node with the largest |direct r| and continue while
# |r| >= 0.9 contiguously; return the topmost node reached (the node itself when
# its best direct link to the next level up is < 0.9). This is what
# redundancy_criterion = "direct" traces internally.
.direct_chase <- function(x, node) {
  lv <- as.integer(sub("^m(\\d+)f\\d+$", "\\1", node))
  cur <- node
  for (j in seq_len(lv - 1L)) {
    E <- x$edges$matrices[[paste0(lv - j, ":", lv)]] # rows level lv-j, cols level lv
    col <- E[, node]
    pi <- which.max(abs(col))
    if (abs(col[pi]) < 0.9) break
    cur <- rownames(E)[pi]
  }
  cur
}

test_that("default output reproduces Forbes's AMH applied example (k = 10)", {
  skip_if_not_installed("psych")
  amh <- .amh_fixture()
  K <- amh$k_max
  suppressWarnings(suppressMessages(
    x <- cached(ackwards(forbes2023, k_max = K, pairs = "all"))
  ))

  # (1) Between-level correlations: |ours| == |hers| entrywise, all 45 pairs.
  idx <- 0L
  for (c2 in 2:K) {
    for (i in 1:(c2 - 1L)) {
      idx <- idx + 1L
      E_forbes <- amh$comp_corr[[idx]]
      E_ours <- x$edges$matrices[[paste0(i, ":", c2)]]
      expect_equal(
        abs(unname(E_ours)), abs(unname(E_forbes)),
        tolerance = 1e-12,
        label = paste0("AMH |edges| ", i, ":", c2, " (ackwards)"),
        expected.label = "Forbes reference"
      )
    }
  }

  # (2) Loading congruence: within her factor.congruence 2-dp rounding.
  phi_ours <- ackwards:::.phi_pairs(x$levels, "all")
  idx <- 0L
  for (c2 in 2:K) {
    for (i in 1:(c2 - 1L)) {
      idx <- idx + 1L
      C_forbes <- amh$cong[[idx]]
      sub <- phi_ours[phi_ours$level_from == i & phi_ours$level_to == c2, ]
      C_ours <- matrix(sub$phi, nrow = i, ncol = c2, byrow = TRUE)
      expect_lt(
        max(abs(abs(C_ours) - abs(unname(C_forbes)))),
        0.005 + 1e-12
      )
    }
  }
})

test_that("direct criterion reproduces Forbes's AMH redundancy chase exactly (54/54)", {
  skip_if_not_installed("psych")
  amh <- .amh_fixture()
  suppressWarnings(suppressMessages(
    x <- cached(ackwards(forbes2023, k_max = amh$k_max, pairs = "all"))
  ))

  # Every component's direct chase over our all-levels edges lands on the same
  # node her ChaseCorrPaths reports ("X--null" = the node stays put). This is
  # the criterion prune("redundant") uses by default; here it is the exact
  # published AMH result, including the 7 components where an adjacent walk
  # would diverge (non-transitivity).
  expect_length(amh$corr_chase, 54L)
  for (entry in amh$corr_chase) {
    parts <- strsplit(entry, "--", fixed = TRUE)[[1L]]
    from <- .forbes_to_ackwards(parts[1L])
    expected_top <- if (parts[2L] == "null") from else .forbes_to_ackwards(parts[2L])
    expect_identical(
      .direct_chase(x, from), expected_top,
      label = paste0("AMH direct chase(", parts[1L], ") (ackwards)"),
      expected.label = paste0("Forbes '", entry, "'")
    )
  }
})

test_that("prune('redundant') on AMH: direct default vs adjacent opt-in", {
  skip_if_not_installed("psych")
  amh <- .amh_fixture()
  suppressWarnings(suppressMessages({
    x <- cached(ackwards(forbes2023, k_max = amh$k_max, pairs = "all"))
    xp <- prune(x, "redundant") # default redundancy_criterion = "direct"
    xa <- prune(x, "redundant", redundancy_criterion = "adjacent")
  }))
  ch <- xp$prune$chains

  # The paper's d4 chain: m4f4 -> ... -> m10f4 (Forbes's d4->e4->f5->g5->h5->i4->j4).
  # It reaches k_max, so the retention rule keeps the most-specific bottom node.
  d4 <- ch[ch$chain_id == ch$chain_id[ch$id == "m10f4"], ]
  expect_setequal(
    d4$id,
    c("m4f4", "m5f4", "m6f5", "m7f5", "m8f5", "m9f4", "m10f4")
  )
  expect_identical(d4$id[d4$retain], "m10f4")

  # Direct (Forbes-faithful) decomposition -- regression pin.
  expect_equal(sum(xp$prune$nodes$pruned), 37L)
  expect_setequal(
    unique(ch$id[ch$retain]),
    c(
      "m1f1", "m4f2", "m5f2", "m7f7",
      "m10f1", "m10f2", "m10f3", "m10f4", "m10f5",
      "m10f6", "m10f7", "m10f8", "m10f9"
    )
  )

  # The adjacent opt-in gives a materially different answer on this deep
  # hierarchy (it retains m3f3, which chases further under the direct rule).
  expect_equal(sum(xa$prune$nodes$pruned), 36L)
  expect_true("m3f3" %in% xa$prune$chains$id[xa$prune$chains$retain])
  expect_false("m3f3" %in% ch$id[ch$retain])
})

# ---------------------------------------------------------------------------
# Oblique branch.
#
# Forbes's reference implementation passes `rotate` straight to psych, so it
# fits oblique PCA and EFA too (her fn. 1). fixtures/forbes2023_oblique.rds
# holds its output under rotate = "oblimin" and "promax", for fm = "pca" and
# "minres", on the three simulation matrices above (data-raw/forbes2023-oblique.R
# and attr(fixture, "provenance")). Only ackwards() runs here.
#
# Correspondence conventions:
#   * Her comp.corr is t(W_a) %*% R %*% W_b with her weights, unstandardized.
#     ackwards divides by the real score SDs, so the test compares against
#     D_a^{-1/2} comp.corr D_b^{-1/2}, with D = diag(W'RW) stored per level.
#     D is not assumed: psych::fa's one-factor weights are not unit-variance,
#     so her level-1 EFA rows need it (see `D` in the fixture).
#   * Her signs are psych's; ours are aligned to the primary parent. The test
#     reads each level's flips from the loadings (ours = hers * s) and compares
#     signed values, so a sign error cannot hide behind abs().
#   * fm = "pca" is engine = "pca"; fm = "minres" is engine = "efa" (minres is
#     its default), with n_obs = 5000, which feeds the fit indices only.

test_that("oblique output reproduces Forbes's oblique branch on her simulations", {
  # The fixture was generated with GPArotation's "bb" algorithm, its default
  # since 2026.6-1. The earlier default ("legacy") moves PCA oblimin loadings
  # by 1.3e-6 to 2.8e-6 (sim1, k = 2 to 4, measured 2026-09-30), far above
  # this test's 1e-10.
  skip_if_not_installed("GPArotation", minimum_version = "2026.6-1")
  runs <- readRDS(test_path("fixtures", "forbes2023_oblique.rds"))
  sims <- .forbes_fixture()
  expect_length(runs, 12L)

  # psych::fa() fits an oblimin (GPArotation) rotation from 20 starts, the
  # identity plus 19 random (its n.rotations default), so two unseeded fits
  # agree only to the rotation's convergence tolerance, and her function fits
  # each level twice. Measured 2026-09-30: up to 1.0e-5 on the minres +
  # oblimin edges, and at most 2e-15 on every other run (PCA and promax
  # results do not depend on the RNG). Those
  # runs get an absolute 1e-4; a method error moves edges by 1e-2 or more.
  # Our fits are seeded so the test's outcome repeats.
  for (nm in names(runs)) {
    run <- runs[[nm]]
    engine <- if (run$fm == "pca") "pca" else "efa"
    tol <- if (run$fm == "minres" && run$rotate == "oblimin") 1e-4 else 1e-10
    suppressWarnings(suppressMessages(
      x <- ackwards(sims[[run$sim]]$R,
        k_max = 4, engine = engine, rotation = run$rotate,
        n_obs = 5000, pairs = "all", seed = 1
      )
    ))

    # (1) Loadings and factor correlations, carried by the same flips.
    s <- lapply(1:4, function(k) {
      sign(colSums(unname(x$levels[[k]]$loadings) * run$loadings[[k]]))
    })
    for (k in 1:4) {
      expect_lt(
        max(abs(unname(x$levels[[k]]$loadings) - sweep(run$loadings[[k]], 2, s[[k]], "*"))),
        tol,
        label = paste(nm, "loadings gap, level", k)
      )
      expect_lt(
        max(abs(unname(x$levels[[k]]$factor_cor) - run$Phi[[k]] * tcrossprod(s[[k]]))),
        tol,
        label = paste(nm, "factor_cor gap, level", k)
      )
    }

    # (2) Every level pair: our edges equal her D-standardized comp.corr.
    idx <- 0L
    for (c2 in 2:4) {
      for (i in 1:(c2 - 1L)) {
        idx <- idx + 1L
        E_her <- run$comp_corr[[idx]] /
          sqrt(outer(unname(run$D[[i]]), unname(run$D[[c2]])))
        E_ours <- unname(x$edges$matrices[[paste0(i, ":", c2)]])
        expect_lt(
          max(abs(E_ours - E_her * outer(s[[i]], s[[c2]]))),
          tol,
          label = paste0(nm, " edge gap ", i, ":", c2)
        )
      }
    }

    # (3) Her redundancy chase at .9 (D-036: the oblique chase stays on r),
    # run by her own code on the D-standardized comp.corr. Her raw chase
    # reads unstandardized products, so it can differ wherever D != I; with
    # D = I (every PCA run) the two are the same list.
    # Her chase is the direct (skip-level) one, which .direct_chase() above
    # traces and prune("redundant") uses by default. One edge case departs
    # from her code: for a level-3+ component whose chase is unbroken to
    # level a (run$chase_unbroken), her ChaseCorrPaths() returns "null",
    # because which.min() on a vector with no FALSE counts zero links (see
    # the generator). There the chase does reach a1, so a1 is expected.
    if (run$fm == "pca") expect_identical(run$corr_chase, run$corr_chase_std)
    for (entry in run$corr_chase_std) {
      parts <- strsplit(entry, "--", fixed = TRUE)[[1L]]
      from <- .forbes_to_ackwards(parts[1L])
      expected_top <- if (parts[2L] == "null") from else .forbes_to_ackwards(parts[2L])
      if (parts[1L] %in% run$chase_unbroken) {
        expect_identical(parts[2L], "null")
        expected_top <- "m1f1"
      }
      expect_identical(
        .direct_chase(x, from), expected_top,
        label = paste0(nm, " chase(", parts[1L], ") (ackwards)"),
        expected.label = paste0("Forbes '", entry, "'")
      )
    }
  }

  # The correspondence is load-bearing, not a formality: her level-1 EFA
  # weights are not unit-variance, and on some runs her raw chase stops at a
  # b-level component that the correlations chase on to a1.
  efa <- runs[vapply(runs, function(r) r$fm == "minres", logical(1))]
  expect_true(all(vapply(efa, function(r) abs(r$D[[1]] - 1) > 0.05, logical(1))))
  expect_true(any(vapply(efa, function(r) !identical(r$corr_chase, r$corr_chase_std), logical(1))))
  # The unbroken-chase edge case is exercised (sim3, minres, promax).
  expect_gt(sum(lengths(lapply(runs, `[[`, "chase_unbroken"))), 0L)
})
