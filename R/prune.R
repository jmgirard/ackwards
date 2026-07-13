# R/prune.R -- Forbes (2023) extension: redundancy/artifact pruning
#
# Design (DESIGN.md s.14 items 17--21; M34):
#   - Standalone verb: prune(x, ...) pipes off ackwards(), not an ackwards()
#     argument. Lets a researcher re-prune with new thresholds without
#     re-running the (expensive) extraction.
#   - Flag-only, never remove. Adds pruned/prune_reason annotations; the object
#     retains all levels (preserves Invariant 5).
#   - Redundancy: two criteria (redundancy_criterion), both threshold |r| >=
#     redundancy_r, optionally conjunctive Tucker's phi > redundancy_phi, same
#     retention rule (keep the bottom node if the chain reaches level k_max, the
#     most specific/best-defined; else keep the top node, the broadest
#     manifestation; Forbes 2023, p.3):
#     * "direct" (DEFAULT, M53): faithful to Forbes's ChaseCorrPaths -- chase up
#       via the DIRECT (skip-level) correlation to the component at each ancestor
#       level. This is what the paper does; it reproduces her AMH applied example
#       exactly (test-forbes-fidelity.R).
#     * "adjacent" (pre-M53 behavior, opt-in): trace ADJACENT primary-parent
#       links only. Because correlation is non-transitive, this both over-reaches
#       (a clean adjacent chain whose endpoints are < redundancy_r) and
#       under-reaches (a direct skip-link >= redundancy_r past a weak adjacent
#       hop) versus the direct criterion in deep hierarchies.
#   - Additive enrichments (do not change which nodes are flagged):
#     * Report r AND phi per chain link (report-first/flag-second).
#     * Report direct endpoint r (root-to-leaf, from all-levels edges) and flag
#       whether it clears redundancy_r -- an at-a-glance cross-check on the chain.
#   - Artifact ("artefact" accepted as an alias): compute Tucker's phi for all
#     cross-level factor pairs; never auto-flag. Researcher interprets (Forbes
#     explicitly flags artifact identification as introducing researcher DoF;
#     cf. Wicherts et al. 2016).
#   - Manual pruning: `manual =` unions user-named nodes onto the auto rules
#     (or works standalone with `rules = "none"`); auto reason wins on overlap.
#   - Edges are recomputed fresh on every call via compute_edges(pairs = "all")
#     from the stored levels/R (Invariant 3); x$edges is never mutated
#     (Invariant 1: one edge path).

# Compute Tucker's phi for factor pairs across levels.
# which_pairs: "adjacent" or "all".
# Returns a data frame: from, to, level_from, level_to, phi
.phi_pairs <- function(levels_list, which_pairs = c("all", "adjacent")) {
  which_pairs <- match.arg(which_pairs)
  ks <- as.integer(names(levels_list))
  K <- max(ks)

  if (which_pairs == "adjacent") {
    pair_list <- lapply(seq_len(K - 1L), function(i) c(ks[i], ks[i + 1L]))
  } else {
    pair_list <- do.call(c, lapply(seq_len(K - 1L), function(i) {
      lapply(seq(i + 1L, K), function(j) c(ks[i], ks[j]))
    }))
  }

  rows <- lapply(pair_list, function(pair) {
    ka <- pair[1L]
    kb <- pair[2L]
    La <- levels_list[[as.character(ka)]]$loadings
    Lb <- levels_list[[as.character(kb)]]$loadings
    labs_a <- colnames(La)
    labs_b <- colnames(Lb)
    do.call(rbind, lapply(seq_along(labs_a), function(i) {
      do.call(rbind, lapply(seq_along(labs_b), function(j) {
        phi_val <- .tucker_phi(La[, i], Lb[, j])
        data.frame(
          from = labs_a[i],
          to = labs_b[j],
          level_from = ka,
          level_to = kb,
          phi = phi_val,
          abs_phi = abs(phi_val),
          stringsAsFactors = FALSE
        )
      }))
    }))
  })

  out <- do.call(rbind, rows)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# --- Strong-link builders (one row per redundancy link: shallower -> deeper) ---
# Both return a data frame with columns from_label, to_label, level_from,
# level_to, r_link, phi_link (r_link/phi_link are the ADJACENT-level stats for
# the two consecutive chain nodes, so downstream chain/retention code is
# criterion-agnostic), or NULL when no links meet the threshold.

# ADJACENT criterion: a link is an adjacent primary-parent pair with
# |r| >= threshold_r (and, if set, phi > threshold_phi). This is the pre-M53
# behavior, retained as an opt-in.
.strong_links_adjacent <- function(x, threshold_r, threshold_phi) {
  k_max <- x$k_max
  levels_list <- x$levels
  lineage <- x$lineage

  link_rows <- lapply(seq(2L, k_max), function(k) {
    key <- paste0(k - 1L, ":", k)
    E <- x$edges$matrices[[key]]
    if (is.null(E)) {
      return(NULL)
    }

    lab_from <- rownames(E)
    lab_to <- colnames(E)
    parents_k <- lineage[[as.character(k)]]
    La <- levels_list[[as.character(k - 1L)]]$loadings
    Lb <- levels_list[[as.character(k)]]$loadings

    res <- lapply(seq_along(lab_to), function(j) {
      pi <- parents_k[j]
      r_val <- E[pi, j]
      if (abs(r_val) < threshold_r) {
        return(NULL)
      }

      phi_val <- .tucker_phi(La[, pi], Lb[, j])
      if (!is.null(threshold_phi) && (is.na(phi_val) || phi_val <= threshold_phi)) {
        return(NULL)
      }

      data.frame(
        from_label = lab_from[pi],
        to_label = lab_to[j],
        level_from = k - 1L,
        level_to = k,
        r_link = r_val,
        phi_link = phi_val,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, res)
  })

  sl <- do.call(rbind, link_rows)
  if (is.null(sl) || nrow(sl) == 0L) {
    return(NULL)
  }
  rownames(sl) <- NULL
  sl
}

# DIRECT criterion (default; M53): faithful to Forbes's (2023) ChaseCorrPaths.
# For each component, chase upward using the DIRECT (skip-level) correlation to
# that component at each ancestor level -- at level j take the node with the
# largest |direct r| and continue while |r| >= threshold_r (and, if set, direct
# phi > threshold_phi), contiguously. Because correlation is non-transitive this
# differs from the adjacent walk in deep hierarchies: it reaches an ancestor a
# weak adjacent hop would stop at, and stops where the direct link to the next
# level fails even though adjacent hops continue. Consecutive chase nodes are
# always one level apart, so the emitted links slot straight into the shared
# chain machinery. Verified to reproduce Forbes's AMH applied example exactly
# (test-forbes-fidelity.R).
.strong_links_direct <- function(x, threshold_r, threshold_phi) {
  levels_list <- x$levels
  node_level <- .node_levels(levels_list)
  labs <- names(node_level)

  rows <- list()
  for (node in labs) {
    L <- node_level[[node]]
    if (L < 2L) next
    prev <- node # previous (deeper) node in the chase, at level j + 1
    for (j in seq.int(L - 1L, 1L)) {
      E_dir <- x$edges$matrices[[paste0(j, ":", L)]] # rows level j, cols level L
      if (is.null(E_dir)) break
      dcol <- E_dir[, node]
      p <- which.max(abs(dcol))
      if (abs(dcol[p]) < threshold_r) break
      Pnode <- rownames(E_dir)[p]

      La <- levels_list[[as.character(j)]]$loadings
      Lc <- levels_list[[as.character(L)]]$loadings
      phi_dir <- .tucker_phi(La[, Pnode], Lc[, node])
      if (!is.null(threshold_phi) && (is.na(phi_dir) || phi_dir <= threshold_phi)) break

      # Adjacent-level stats for the consecutive chain nodes (Pnode@j, prev@j+1).
      E_adj <- x$edges$matrices[[paste0(j, ":", j + 1L)]]
      Lb <- levels_list[[as.character(j + 1L)]]$loadings
      rows[[length(rows) + 1L]] <- data.frame(
        from_label = Pnode,
        to_label = prev,
        level_from = j,
        level_to = j + 1L,
        r_link = E_adj[Pnode, prev],
        phi_link = .tucker_phi(La[, Pnode], Lb[, prev]),
        stringsAsFactors = FALSE
      )
      prev <- Pnode
    }
  }

  if (length(rows) == 0L) {
    return(NULL)
  }
  sl <- unique(do.call(rbind, rows))
  rownames(sl) <- NULL
  sl
}

# Find redundant chains from adjacent primary-parent or direct/skip-level links.
# Returns list(node_flags = df | NULL, chains = df | NULL).
.find_redundant_chains <- function(x, threshold_r, threshold_phi,
                                   criterion = c("direct", "adjacent")) {
  criterion <- match.arg(criterion)
  k_max <- x$k_max
  levels_list <- x$levels

  sl <- if (criterion == "direct") {
    .strong_links_direct(x, threshold_r, threshold_phi)
  } else {
    .strong_links_adjacent(x, threshold_r, threshold_phi)
  }
  if (is.null(sl) || nrow(sl) == 0L) {
    return(list(node_flags = NULL, chains = NULL))
  }

  # --- Map label -> level -------------------------------------------------------
  label_to_level <- .node_levels(levels_list)

  # --- Find chain roots (from_label not appearing as any to_label) -------------
  root_labels <- setdiff(unique(sl$from_label), unique(sl$to_label))

  # --- Enumerate all root-to-leaf paths via DFS --------------------------------
  # A root with multiple strong-link children (sibling redundancies) spawns one
  # chain per branch -- each is analysed independently. DFS also handles diamonds
  # (a node reachable from two roots) without special-casing.
  chains_raw <- list()
  for (root in root_labels) {
    stack <- list(root)
    while (length(stack) > 0L) {
      path <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      last <- path[length(path)]
      children <- sl$to_label[sl$from_label == last]
      if (length(children) == 0L) {
        if (length(path) >= 2L) chains_raw <- c(chains_raw, list(path))
      } else {
        for (child in children) {
          stack <- c(stack, list(c(path, child)))
        }
      }
    }
  }
  if (length(chains_raw) == 0L) { # nocov start
    return(list(node_flags = NULL, chains = NULL))
  } # nocov end

  # --- Per-chain metadata: retain label + link stats + endpoint r --------------
  chain_metas <- lapply(seq_along(chains_raw), function(i) {
    ch <- chains_raw[[i]]
    ch_levels <- label_to_level[ch]

    # Retention rule (Forbes 2023, p.3):
    retain_label <- if (ch_levels[length(ch)] == k_max) {
      ch[length(ch)] # chain reaches most-specific level -> keep bottom
    } else {
      ch[1L] # chain ends mid-hierarchy -> keep top (broadest)
    }

    r_vals <- vapply(seq_along(ch), function(idx) {
      if (idx == 1L) {
        return(NA_real_)
      }
      row <- sl[sl$from_label == ch[idx - 1L] & sl$to_label == ch[idx], ]
      if (nrow(row) == 0L) NA_real_ else row$r_link[1L]
    }, numeric(1L))

    phi_vals <- vapply(seq_along(ch), function(idx) {
      if (idx == 1L) {
        return(NA_real_)
      }
      row <- sl[sl$from_label == ch[idx - 1L] & sl$to_label == ch[idx], ]
      if (nrow(row) == 0L) NA_real_ else row$phi_link[1L]
    }, numeric(1L))

    # Endpoint r enrichment: direct r between chain root and leaf from all-levels edges.
    # Disagreement (endpoint_r_agrees = FALSE) means the chain method and the direct
    # correlation give conflicting answers -- flag for researcher attention.
    ep_key <- paste0(ch_levels[1L], ":", ch_levels[length(ch)])
    E_ep <- x$edges$matrices[[ep_key]]
    endpoint_r <- NA_real_
    endpoint_r_agrees <- NA

    if (!is.null(E_ep)) {
      ri <- match(ch[1L], rownames(E_ep))
      ci <- match(ch[length(ch)], colnames(E_ep))
      if (!is.na(ri) && !is.na(ci)) {
        endpoint_r <- E_ep[ri, ci]
        endpoint_r_agrees <- abs(endpoint_r) >= threshold_r
      }
    }

    list(
      ch = ch, ch_levels = ch_levels, retain_label = retain_label,
      r_vals = r_vals, phi_vals = phi_vals,
      endpoint_r = endpoint_r, endpoint_r_agrees = endpoint_r_agrees
    )
  })

  # --- Global retain set -------------------------------------------------------
  # A node retained in ANY chain is never pruned. This prevents contradictions
  # when sibling branches or diamonds produce conflicting per-chain decisions
  # (e.g. a shared root is retained by one branch but would be flagged by another).
  retain_any <- unique(vapply(chain_metas, `[[`, character(1L), "retain_label"))

  # --- Build output data frames ------------------------------------------------
  chain_rows <- vector("list", length(chain_metas))
  for (i in seq_along(chain_metas)) {
    cm <- chain_metas[[i]]
    chain_rows[[i]] <- data.frame(
      chain_id          = i,
      id                = cm$ch,
      level             = cm$ch_levels,
      r_to_prev         = cm$r_vals,
      phi_to_prev       = cm$phi_vals,
      retain            = cm$ch == cm$retain_label,
      endpoint_r        = cm$endpoint_r,
      endpoint_r_agrees = cm$endpoint_r_agrees,
      stringsAsFactors  = FALSE
    )
  }

  chains_df <- do.call(rbind, chain_rows)
  rownames(chains_df) <- NULL

  # Flagged = appears in some chain but never the retained node in any chain
  all_chained <- unique(unlist(lapply(chain_metas, `[[`, "ch")))
  flagged_ids <- setdiff(all_chained, retain_any)

  node_flags_df <- if (length(flagged_ids) > 0L) {
    out <- data.frame(
      id = flagged_ids,
      level = label_to_level[flagged_ids],
      pruned = TRUE,
      prune_reason = "redundant",
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
  } else { # nocov start
    NULL
  } # nocov end

  list(node_flags = node_flags_df, chains = chains_df)
}

# Compute structural artifact signals for each (level, factor).
# Returns a data frame: id, level, few_items, orphan, split_merge.
#
# few_items:   factor is primary parent for fewer than min_items items.
# orphan:      factor's max adjacent-level |r| is below orphan_r (non-replicating).
# split_merge: factor's primary items came from multiple different primary parents
#              at the shallower adjacent level (items merged from separate groups).
.compute_structural_signals <- function(x, min_items, orphan_r) {
  levels_list <- x$levels
  ks <- sort(as.integer(names(levels_list)))

  # Primary factor label for each item at each level.
  primary_factor_of <- lapply(ks, function(k) {
    L <- levels_list[[as.character(k)]]$loadings
    labels_k <- colnames(L)
    idx <- apply(abs(L), 1L, which.max)
    stats::setNames(labels_k[idx], rownames(L))
  })
  names(primary_factor_of) <- as.character(ks)

  rows <- lapply(ks, function(k) {
    lev <- levels_list[[as.character(k)]]
    labels_k <- colnames(lev$loadings)
    pf_k <- primary_factor_of[[as.character(k)]]

    lapply(seq_along(labels_k), function(j_idx) {
      jlab <- labels_k[j_idx]
      primary_items <- names(pf_k)[pf_k == jlab]

      # few_items ---------------------------------------------------------------
      few_items <- length(primary_items) < min_items

      # orphan: max adjacent-level |r| < orphan_r ------------------------------
      adj_r <- numeric(0L)
      if (k > min(ks)) {
        key_up <- paste0(k - 1L, ":", k)
        E_up <- x$edges$matrices[[key_up]]
        if (!is.null(E_up) && jlab %in% colnames(E_up)) {
          adj_r <- c(adj_r, abs(E_up[, jlab]))
        }
      }
      if (k < max(ks)) {
        key_dn <- paste0(k, ":", k + 1L)
        E_dn <- x$edges$matrices[[key_dn]]
        if (!is.null(E_dn) && jlab %in% rownames(E_dn)) {
          adj_r <- c(adj_r, abs(E_dn[jlab, ]))
        }
      }
      # Every factor has >= 1 adjacent level (k_max >= 2, artifact uses
      # pairs = "all"), so the empty-adj_r guard is defensive only.
      orphan <- if (length(adj_r) == 0L) {
        NA # nocov
      } else {
        max(adj_r) < orphan_r
      }

      # split_merge: items came from multiple primary parents at level k-1 ------
      split_merge <- if (k <= min(ks) || length(primary_items) == 0L) {
        FALSE
      } else {
        pf_prev <- primary_factor_of[[as.character(k - 1L)]]
        parents <- unique(pf_prev[primary_items])
        length(parents[!is.na(parents)]) > 1L
      }

      data.frame(
        id = jlab,
        level = k,
        few_items = few_items,
        orphan = orphan,
        split_merge = split_merge,
        stringsAsFactors = FALSE
      )
    })
  })

  out <- do.call(rbind, do.call(c, rows))
  rownames(out) <- NULL
  out
}

# Baseline node table -- all nodes, initially not pruned.
.base_prune_nodes <- function(levels_list) {
  lvl_map <- .node_levels(levels_list)
  out <- data.frame(
    id = names(lvl_map),
    level = unname(lvl_map),
    pruned = FALSE,
    prune_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# Auto-rule pruning engine. `edges_view` is a lightweight stand-in for an
# ackwards object -- list(k_max, levels, lineage, edges = list(matrices = ...))
# -- built by prune.ackwards() from freshly recomputed all-pairs edges, so
# these helpers can stay agnostic to whether x$edges holds adjacent or all
# pairs. `rules` here is the auto-rule subset only (never "none"/"manual").
# Returns the $prune slot (minus $manual, added by the caller).
.apply_pruning <- function(levels_list, edges_view, rules, redundancy_r, redundancy_phi,
                           min_items, orphan_r, redundancy_criterion = "direct") {
  base_nodes <- .base_prune_nodes(levels_list)

  chains_df <- NULL
  phi_df <- NULL

  if ("redundant" %in% rules) {
    chain_result <- .find_redundant_chains(
      edges_view, redundancy_r, redundancy_phi, redundancy_criterion
    )
    chains_df <- chain_result$chains

    if (!is.null(chain_result$node_flags) && nrow(chain_result$node_flags) > 0L) {
      for (ri in seq_len(nrow(chain_result$node_flags))) {
        fl <- chain_result$node_flags[ri, ]
        mask <- base_nodes$id == fl$id & base_nodes$level == fl$level
        base_nodes$pruned[mask] <- TRUE
        base_nodes$prune_reason[mask] <- fl$prune_reason
      }
    }
  }

  structural_df <- NULL
  if ("artifact" %in% rules) {
    # Compute phi for all cross-level pairs for researcher inspection.
    # Pruning is not automated here -- artifact identification requires judgment
    # (Forbes 2023 is explicit: this step introduces researcher DoF).
    phi_df <- .phi_pairs(levels_list, which_pairs = "all")

    # Compute structural artifact signals: flag/report only, never auto-prune.
    structural_df <- .compute_structural_signals(edges_view, min_items, orphan_r)
  }

  n_flagged <- sum(base_nodes$pruned)
  phi_note <- if (!is.null(redundancy_phi)) {
    paste0(" and phi > ", redundancy_phi)
  } else {
    ""
  }

  if ("redundant" %in% rules) {
    if (n_flagged > 0L) {
      cli::cli_inform(c(
        "i" = "Redundancy pruning ({redundancy_criterion} criterion, \\
               |r| {cli::symbol$geq} {redundancy_r}{phi_note}) flagged {n_flagged} node{?s}.",
        "i" = "Nodes are retained in the object; inspect with \\
               {.code x$prune$nodes} and {.code x$prune$chains}."
      ))
    } else {
      cli::cli_inform(c(
        "i" = "Redundancy pruning ({redundancy_criterion} criterion) found no chains \\
               meeting |r| {cli::symbol$geq} {redundancy_r}{phi_note}."
      ))
    }
  }

  if ("artifact" %in% rules) {
    # structural_df is always set above when artifact pruning runs.
    n_struct <- sum(structural_df$few_items | structural_df$orphan |
      structural_df$split_merge, na.rm = TRUE)
    cli::cli_inform(c(
      # cli::symbol has no $phi entry (renders empty); use the escaped glyph.
      "i" = "Artifact mode: Tucker's \u03c6 computed for all \\
             cross-level factor pairs.",
      "i" = "Structural signals computed: {n_struct} factor{?s} flagged \\
             (few_items / orphan / split_merge).",
      "i" = "Inspect {.code x$prune$phi} and {.code x$prune$structural}; \\
             removal is a researcher judgment (Forbes, 2023)."
    ))
  }

  list(
    rules                = rules,
    redundancy_r         = redundancy_r,
    redundancy_phi       = redundancy_phi,
    redundancy_criterion = redundancy_criterion,
    min_items            = min_items,
    orphan_r             = orphan_r,
    nodes                = base_nodes,
    chains               = chains_df,
    phi                  = phi_df,
    structural           = structural_df
  )
}

#' Flag redundant or artifactual factors (Forbes 2023 extension)
#'
#' @description
#' `prune()` never removes anything from an `ackwards` object -- it only
#' annotates factors with `pruned`/`prune_reason` flags in `x$prune$nodes`
#' (flag-only, never removes; the object keeps every level). Because it is a
#' separate, cheap step from extraction, you can re-prune with new thresholds
#' without re-running [ackwards()]:
#' \preformatted{
#'   x <- ackwards(bfi25, k_max = 6, engine = "esem")  # expensive
#'   x |> prune("redundant")                            # cheap, repeatable
#'   x |> prune("redundant", redundancy_r = 0.95)        # no re-extraction
#' }
#' `prune()` is an S3 generic (rather than a plain function) so it coexists
#' with the `prune` generics already defined by recursive-partitioning
#' packages (e.g. `rpart::prune`) regardless of package load order.
#'
#' @param x An `ackwards` object.
#' @param rules Character vector controlling which auto-rules run. Default
#'   `"none"` (no auto rule; combine with `manual` for pure manual pruning, or
#'   call `prune(x)` with no arguments to clear any existing pruning). Options:
#'   * `"redundant"` -- identify chains of factors connected by score
#'     correlations `|r| >= redundancy_r` (and optionally `phi > redundancy_phi`),
#'     using `redundancy_criterion` (default `"direct"`, faithful to Forbes).
#'     Applies Forbes's (2023) retention rule: keep the bottom node when the
#'     chain reaches level `k_max` (most specific); keep the top node otherwise.
#'     Flagged nodes get `pruned = TRUE` and `prune_reason = "redundant"` in
#'     `x$prune$nodes`.
#'   * `"artifact"` (or the alias `"artefact"`, normalized to `"artifact"`) --
#'     compute Tucker's congruence coefficient (phi) for all cross-level factor
#'     pairs and store in `x$prune$phi`, plus structural signals
#'     (`few_items`/`orphan`/`split_merge`) in `x$prune$structural`. No factors
#'     are auto-flagged; artifact identification requires judgment (Forbes,
#'     2023; Wicherts et al., 2016).
#' @param manual Character vector of factor labels (e.g. `c("m4f3", "m4f4")`)
#'   to flag directly, in addition to or instead of an auto rule. Standalone
#'   manual pruning is supported: `prune(x, manual = c("m4f3"))`. Unknown
#'   labels error. A node already flagged by an auto rule keeps that
#'   `prune_reason`; only otherwise-unflagged manual nodes get
#'   `prune_reason = "manual"`.
#' @param redundancy_r Scalar in `(0, 1]`. Score-correlation `|r|` threshold for
#'   redundancy chains. Default `0.9` (Forbes, 2023).
#' @param redundancy_criterion How redundancy chains are traced. One of:
#'   * `"direct"` (default) -- chase upward via the **direct (skip-level)**
#'     correlation between a factor and each ancestor level, continuing while
#'     `|r| >= redundancy_r` contiguously. This is Forbes's (2023) published
#'     `ChaseCorrPaths` rule and the honest operationalization of "the same
#'     construct" (two factor scores share `>= redundancy_r^2` variance
#'     *directly*). It reproduces her AMH applied example exactly.
#'   * `"adjacent"` -- trace **adjacent primary-parent** links only (each
#'     consecutive level `|r| >= redundancy_r`). This was the pre-M53 default;
#'     because correlation is non-transitive it can both over- and under-flag
#'     versus `"direct"` in deep (many-level) hierarchies, so it is retained only
#'     as an opt-in. On shallow/transitive hierarchies the two agree.
#' @param redundancy_phi Scalar in `(0, 1]`, `NULL` (default, auto), or `NA`
#'   (explicit opt-out). When `NULL`:
#'   * `x$engine == "pca"` -- no phi filter. Component scores are
#'     *determinate* (exact linear functions of the data, with no
#'     factor-score indeterminacy), so the score correlation `|r|` is the
#'     true correlation between the components themselves; phi adds nothing
#'     that `|r|` does not already capture.
#'   * `x$engine` is `"efa"` or `"esem"` -- automatically set to `0.95`
#'     (Lorenzo-Seva & ten Berge, 2006). Factor-score indeterminacy off-PCA
#'     means `|r|`-alone is liberal; the conjunctive phi criterion is the
#'     conservative default. A cli message announces the resolved value.
#'   Pass `NA` to disable phi filtering regardless of engine. Pass a numeric
#'   value to override on any engine.
#' @param min_items Minimum number of items for which a factor must be the
#'   primary loader (highest `|loading|`). Factors with fewer than `min_items`
#'   primary items are flagged `few_items = TRUE` in `x$prune$structural`. Only
#'   used when `rules` includes `"artifact"`. Default `3L` -- a factor defined
#'   by one or two items is under-identified and frequently an extraction
#'   artifact rather than a replicable construct (the classic "three-indicator
#'   rule"; Forbes, 2023, Fig. 2).
#' @param orphan_r Threshold for the `orphan` structural signal. A factor whose
#'   maximum **adjacent-level** `|r|` (to the immediately shallower and deeper
#'   levels) falls below `orphan_r` is flagged `orphan = TRUE` in
#'   `x$prune$structural` -- it does not connect to the neighbouring solutions
#'   and so does not replicate across the hierarchy. Only used when `rules`
#'   includes `"artifact"`. Default `0.5` -- a moderate correlation; a factor
#'   that shares less than a quarter of its variance with every neighbour is a
#'   structural outlier worth inspecting.
#' @param ... Reserved for future methods/arguments.
#'
#' @details
#' **Reading `x$prune$chains` under `redundancy_criterion = "direct"`.** The
#' `r_to_prev` and `phi_to_prev` columns report the **adjacent-level**
#' correlation and congruence between consecutive chain members (for continuity
#' of display), but chain *membership* is decided by the **direct (skip-level)**
#' correlation to the chain's deepest factor. A direct chain can therefore
#' legitimately contain a link whose `r_to_prev` is *below* `redundancy_r` -- the
#' stronger direct link is what justified it. The `endpoint_r` column gives the
#' direct root-to-leaf correlation as an at-a-glance cross-check. Under
#' `redundancy_criterion = "adjacent"`, `r_to_prev` *is* the criterion and always
#' meets `redundancy_r`.
#'
#' @return `x`, with `$prune` populated (replacing any prior pruning).
#'
#' @seealso [ackwards()], [tidy.ackwards()] (`what = "nodes"`),
#'   [autoplot.ackwards()] (`drop_pruned`)
#'
#' @references
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
#'   coefficient as a meaningful index of factor similarity.
#'   *Methodology*, 2(2), 57--64. \doi{10.1027/1614-2241.2.2.57}
#'
#' @examples
#' # sim16 has a planted redundant chain + overextraction artifact at k = 5,
#' # so the prune rules always have a finding to show (and no ordinal warning).
#' x <- ackwards(sim16, k_max = 5)
#'
#' xp <- prune(x, "redundant")
#' xp$prune$nodes
#'
#' # Re-prune with a new threshold -- no re-extraction needed
#' prune(x, "redundant", redundancy_r = 0.95)
#'
#' # Manual pruning: standalone, or mixed with an auto rule
#' prune(x, manual = "m4f2")
#' prune(x, "redundant", manual = "m4f2")
#'
#' @export
prune <- function(x, ...) {
  UseMethod("prune")
}

#' @rdname prune
#' @importFrom rlang `%||%`
#' @export
prune.ackwards <- function(x, rules = "none", manual = NULL,
                           redundancy_r = 0.9, redundancy_phi = NULL,
                           redundancy_criterion = c("direct", "adjacent"),
                           min_items = 3L, orphan_r = 0.5, ...) {
  .check_unknown_dots(list(...), "prune")
  levels_list <- x$levels
  all_ids <- unlist(lapply(levels_list, `[[`, "labels"))

  rules <- rlang::arg_match(rules, c("none", "redundant", "artifact", "artefact"),
    multiple = TRUE
  )
  rules[rules == "artefact"] <- "artifact"
  rules <- unique(rules)

  redundancy_criterion <- rlang::arg_match(redundancy_criterion)

  if (!is.null(manual)) {
    if (!is.character(manual)) {
      cli::cli_abort("{.arg manual} must be a character vector of factor labels.")
    }
    unknown <- setdiff(manual, all_ids)
    if (length(unknown) > 0L) {
      cli::cli_abort(c(
        "!" = "{.arg manual} contains unknown factor label{?s}: {.val {unknown}}.",
        "i" = "Valid labels: {.val {all_ids}}."
      ))
    }
    manual <- unique(manual) # flagging is idempotent; keep $manual dup-free
  }

  if (!is.numeric(redundancy_r) || length(redundancy_r) != 1L ||
    redundancy_r <= 0 || redundancy_r > 1) {
    cli::cli_abort("{.arg redundancy_r} must be a single number in (0, 1].")
  }
  # NA is the explicit opt-out ("no phi regardless of engine").
  # NULL is auto (resolved below once rules are known).
  if (!is.null(redundancy_phi) && !isTRUE(is.na(redundancy_phi)) &&
    (!is.numeric(redundancy_phi) || length(redundancy_phi) != 1L ||
      redundancy_phi <= 0 || redundancy_phi > 1)) {
    cli::cli_abort(
      "{.arg redundancy_phi} must be a number in (0, 1], {.code NULL} (auto), or {.code NA} (opt-out)."
    )
  }
  if (!is.numeric(min_items) || length(min_items) != 1L ||
    is.na(min_items) || min_items < 1L || min_items != as.integer(min_items)) {
    cli::cli_abort("{.arg min_items} must be a single positive integer.")
  }
  if (!is.numeric(orphan_r) || length(orphan_r) != 1L ||
    is.na(orphan_r) || orphan_r < 0 || orphan_r > 1) {
    cli::cli_abort("{.arg orphan_r} must be a single number in [0, 1].")
  }
  min_items <- as.integer(min_items)

  auto_rules <- setdiff(rules, "none")

  if (length(auto_rules) == 0L && is.null(manual)) {
    x$prune <- NULL
    return(x)
  }

  # Auto-resolve redundancy_phi = NULL (Invariant 6: announce loud).
  # NA is the explicit opt-out and maps to NULL internally.
  if ("redundant" %in% auto_rules) {
    if (is.null(redundancy_phi)) {
      if (x$engine %in% c("efa", "esem")) {
        redundancy_phi <- 0.95
        cli::cli_inform(c(
          "i" = "{.arg redundancy_phi} auto-set to {.val 0.95} for {.val {x$engine}} engine \\
                 (Lorenzo-Seva & ten Berge, 2006).",
          "i" = "Factor-score indeterminacy off-PCA means {.code |r|}-only \\
                 redundancy is liberal; phi adds a congruence guard.",
          "i" = "To opt out: pass {.code redundancy_phi = NA}."
        ))
      }
      # engine == "pca": keep NULL (|r|-only; component scores are determinate)
    } else if (isTRUE(is.na(redundancy_phi))) {
      redundancy_phi <- NULL # explicit opt-out -> same as PCA default
    }
  }

  if (length(auto_rules) > 0L) {
    # Recompute all-pairs edges fresh from the stored levels/R (Invariant 3;
    # Invariant 1 -- x$edges is never touched). Cheap: W'RW algebra, not
    # re-extraction (DESIGN.md s.3).
    all_edges <- compute_edges(
      levels = levels_list, R = x$r, edge_method = "auto",
      pairs = "all",
      cut_show = x$meta$cut_show %||% 0.3
    )
    edges_view <- list(
      k_max = x$k_max, levels = levels_list, lineage = x$lineage,
      edges = list(matrices = all_edges$matrices)
    )
    result <- .apply_pruning(
      levels_list, edges_view, auto_rules, redundancy_r, redundancy_phi,
      min_items = min_items, orphan_r = orphan_r,
      redundancy_criterion = redundancy_criterion
    )
  } else {
    result <- list(
      rules = "none", redundancy_r = redundancy_r, redundancy_phi = redundancy_phi,
      redundancy_criterion = redundancy_criterion,
      min_items = min_items, orphan_r = orphan_r,
      nodes = .base_prune_nodes(levels_list), chains = NULL, phi = NULL, structural = NULL
    )
  }

  if (!is.null(manual)) {
    idx <- match(manual, result$nodes$id)
    to_add <- idx[!result$nodes$pruned[idx]]
    if (length(to_add) > 0L) {
      result$nodes$pruned[to_add] <- TRUE
      result$nodes$prune_reason[to_add] <- "manual"
    }
  }
  result$manual <- manual

  x$prune <- result
  x
}
