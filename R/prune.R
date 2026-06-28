# R/prune.R — Forbes (2023) extension: redundancy/artefact pruning (internal)
#
# Design (DESIGN.md §14 items 17–21):
#   - Flag-only, never remove. Adds pruned/prune_reason annotations; the object
#     retains all levels (preserves Invariant 5).
#   - Redundancy: trace chains via ADJACENT primary-parent links with
#     |r| >= redundancy_r. Optionally require Tucker's phi > redundancy_phi
#     (conjunctive). Retention rule: keep the bottom node if the chain reaches
#     level k_max (most specific, best-defined); else keep the top node
#     (broadest manifestation, Forbes 2023, p.3).
#   - Additive enrichments (do not change default output):
#     * Report r AND phi per chain link (report-first/flag-second).
#     * Report direct endpoint r from all-levels edges and flag when it
#       disagrees with the chain criterion (correlation is non-transitive:
#       a clean adjacent chain neither implies nor is implied by endpoint identity).
#   - Artefact: compute Tucker's phi for all cross-level factor pairs; never
#     auto-flag. Researcher interprets (Forbes explicitly flags artefact
#     identification as introducing researcher DoF; cf. Wicherts et al. 2016).

# Tucker's congruence coefficient between two loading vectors (Lorenzo-Seva &
# ten Berge, 2006). Formula: phi = sum(a*b) / sqrt(sum(a^2) * sum(b^2))
.tucker_phi <- function(a, b) {
  denom <- sqrt(sum(a^2) * sum(b^2))
  if (denom == 0) {
    return(NA_real_)
  }
  sum(a * b) / denom
}

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

# Find redundant chains using adjacent primary-parent |r| links.
# Returns list(node_flags = df | NULL, chains = df | NULL).
.find_redundant_chains <- function(x, threshold_r, threshold_phi) {
  k_max <- x$k_max
  levels_list <- x$levels
  lineage <- x$lineage

  # --- Build "strong links": adjacent primary-parent pairs meeting threshold ---
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
    return(list(node_flags = NULL, chains = NULL))
  }
  rownames(sl) <- NULL

  # --- Map label → level -------------------------------------------------------
  label_to_level <- stats::setNames(
    unlist(lapply(names(levels_list), function(ki) rep(as.integer(ki), as.integer(ki)))),
    unlist(lapply(levels_list, `[[`, "labels"))
  )

  # --- Find chain roots (from_label not appearing as any to_label) -------------
  root_labels <- setdiff(unique(sl$from_label), unique(sl$to_label))

  # --- Enumerate all root-to-leaf paths via DFS --------------------------------
  # A root with multiple strong-link children (sibling redundancies) spawns one
  # chain per branch — each is analysed independently. DFS also handles diamonds
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
  if (length(chains_raw) == 0L) {
    return(list(node_flags = NULL, chains = NULL))
  }

  # --- Per-chain metadata: retain label + link stats + endpoint r --------------
  chain_metas <- lapply(seq_along(chains_raw), function(i) {
    ch <- chains_raw[[i]]
    ch_levels <- label_to_level[ch]

    # Retention rule (Forbes 2023, p.3):
    retain_label <- if (ch_levels[length(ch)] == k_max) {
      ch[length(ch)] # chain reaches most-specific level → keep bottom
    } else {
      ch[1L] # chain ends mid-hierarchy → keep top (broadest)
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
    # correlation give conflicting answers — flag for researcher attention.
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
  } else {
    NULL
  }

  list(node_flags = node_flags_df, chains = chains_df)
}

# Main pruning dispatcher. Called from ackwards() when prune != "none".
# Returns the $prune slot to attach to the ackwards object.
.apply_pruning <- function(x, prune, redundancy_r, redundancy_phi) {
  levels_list <- x$levels

  # Baseline node table — all nodes, initially not pruned
  base_nodes <- data.frame(
    id = unlist(lapply(levels_list, `[[`, "labels")),
    level = unlist(lapply(
      names(levels_list),
      function(ki) rep(as.integer(ki), as.integer(ki))
    )),
    pruned = FALSE,
    prune_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(base_nodes) <- NULL

  chains_df <- NULL
  phi_df <- NULL

  if ("redundant" %in% prune) {
    chain_result <- .find_redundant_chains(x, redundancy_r, redundancy_phi)
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

  if ("artefact" %in% prune) {
    # Compute phi for all cross-level pairs for researcher inspection.
    # Pruning is not automated here — artefact identification requires judgment
    # (Forbes 2023 is explicit: this step introduces researcher DoF).
    phi_df <- .phi_pairs(levels_list, which_pairs = "all")
  }

  n_flagged <- sum(base_nodes$pruned)
  phi_note <- if (!is.null(redundancy_phi)) {
    paste0(" and phi > ", redundancy_phi)
  } else {
    ""
  }

  if ("redundant" %in% prune) {
    if (n_flagged > 0L) {
      cli::cli_inform(c(
        "i" = "Redundancy pruning (|r| {cli::symbol$geq} {redundancy_r}{phi_note}) \\
               flagged {n_flagged} node{?s}.",
        "i" = "Nodes are retained in the object; inspect with \\
               {.code x$prune$nodes} and {.code x$prune$chains}."
      ))
    } else {
      cli::cli_inform(c(
        "i" = "Redundancy pruning found no chains meeting \\
               |r| {cli::symbol$geq} {redundancy_r}{phi_note}."
      ))
    }
  }

  if ("artefact" %in% prune) {
    cli::cli_inform(c(
      "i" = "Artefact mode: Tucker's {cli::symbol$phi} computed for all \\
             cross-level factor pairs.",
      "i" = "Inspect {.code x$prune$phi} to identify potential artefacts; \\
             removal is a researcher judgment (Forbes, 2023)."
    ))
  }

  list(
    rules          = prune,
    redundancy_r   = redundancy_r,
    redundancy_phi = redundancy_phi,
    nodes          = base_nodes,
    chains         = chains_df,
    phi            = phi_df
  )
}
