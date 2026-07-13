# print/summary snapshot: PCA, converged, no prune

    Code
      snap_print(x)
    Message
      
      -- Bass-Ackwards Analysis (ackwards) -------------------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 4
      
      -- Levels --
      
      v k = 1: 1 factor, 28.2% variance
      v k = 2: 2 factors, 46.5% variance
      v k = 3: 3 factors, 57.5% variance
      v k = 4: 4 factors, 67.7% variance
      
      -- Edges --
      
      9 of 20 edges have |r| >= 0.3
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

---

    Code
      snap_print(summary(x))
    Message
      
      -- Summary: Bass-Ackwards Analysis (ackwards) ----------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 4
      
      -- Levels --
      
      k = 1: 1 factor (28.2% cumulative variance)
      m1f1 28.2% eigenvalue 4.51
      
      k = 2: 2 factors (46.5% cumulative variance)
      m2f1 23.3% eigenvalue 4.51
      m2f2 23.2% eigenvalue 2.93
      
      k = 3: 3 factors (57.5% cumulative variance)
      m3f1 23.0% eigenvalue 4.51
      m3f2 17.5% eigenvalue 2.93
      m3f3 16.9% eigenvalue 1.76
      
      k = 4: 4 factors (67.7% cumulative variance)
      m4f1 17.2% eigenvalue 4.51
      m4f2 16.9% eigenvalue 2.93
      m4f3 16.8% eigenvalue 1.76
      m4f4 16.8% eigenvalue 1.63
      
      -- Lineage (primary parents) --
      
      m1f1 > m2f1, m2f2
      m2f1 > m3f2, m3f3
      m2f2 > m3f1
      m3f1 > m4f3, m4f4
      m3f2 > m4f1
      m3f3 > m4f2
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

# print/summary snapshot: EFA fit-index glyph line

    Code
      snap_print(summary(x))
    Message
      
      -- Summary: Bass-Ackwards Analysis (ackwards) ----------------------------------
      Engine: efa
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 3
      
      -- Levels --
      
      k = 1: 1 factor (25.0% cumulative variance)
      m1f1 25.0%
      
      k = 2: 2 factors (37.8% cumulative variance)
      m2f1 19.1%
      m2f2 18.7%
      chi = 136.4 dof = 26 RMSEA = 0.065 x TLI = 0.912 x
      
      k = 3: 3 factors (41.9% cumulative variance)
      m3f1 19.2%
      m3f2 18.8%
      m3f3 3.9%
      chi = 58.5 dof = 18 RMSEA = 0.047 v TLI = 0.953 v
      
      -- Lineage (primary parents) --
      
      m1f1 > m2f1, m2f2
      m2f1 > m3f1
      m2f2 > m3f2, m3f3
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

# print/summary snapshot: pruning digest (redundant + artifact)

    Code
      snap_print(x)
    Message
      
      -- Bass-Ackwards Analysis (ackwards) -------------------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 4
      
      -- Levels --
      
      v k = 1: 1 factor, 32.3% variance
      v k = 2: 2 factors, 49.6% variance
      v k = 3: 3 factors, 58.6% variance
      v k = 4: 4 factors, 66.2% variance
      
      -- Edges --
      
      9 of 20 edges have |r| >= 0.3
      
      -- Pruning --
      
      Redundancy (direct, |r| >= 0.9): 4 nodes flagged
      Artifact: Tucker's phi computed for 35 cross-level factor pairs
      Structural signals: 3 factors flagged (inspect `x$prune$structural`)
      --------------------------------------------------------------------------------
      Note: Pruning is interpretive relabeling, not re-estimation. Flagged nodes
      remain in the object; all edges are preserved. Inspect with `x$prune$nodes` and
      `tidy(x, what = "nodes")`.
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

---

    Code
      snap_print(summary(x))
    Message
      
      -- Summary: Bass-Ackwards Analysis (ackwards) ----------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 4
      
      -- Levels --
      
      k = 1: 1 factor (32.3% cumulative variance)
      m1f1 32.3% eigenvalue 3.23
      
      k = 2: 2 factors (49.6% cumulative variance)
      m2f1 24.9% eigenvalue 3.23
      m2f2 24.6% eigenvalue 1.73
      
      k = 3: 3 factors (58.6% cumulative variance)
      m3f1 24.3% eigenvalue 3.23
      m3f2 22.6% eigenvalue 1.73
      m3f3 11.6% eigenvalue 0.90
      
      k = 4: 4 factors (66.2% cumulative variance)
      m4f1 22.4% eigenvalue 3.23
      m4f2 17.1% eigenvalue 1.73
      m4f3 15.6% eigenvalue 0.90
      m4f4 11.2% eigenvalue 0.77
      
      -- Lineage (primary parents) --
      
      m1f1 > m2f1, m2f2
      m2f1 > m3f1
      m2f2 > m3f2, m3f3
      m3f1 > m4f2, m4f3
      m3f2 > m4f1
      m3f3 > m4f4
      
      -- Pruning --
      
      Redundant (direct, |r| >= 0.9): 4 nodes flagged
      Flagged: m2f2, m3f1, m3f2, m3f3
      Artifact: Tucker's phi computed for 35 cross-level pairs
      Structural signals: 3 factors flagged (inspect `x$prune$structural`)
      --------------------------------------------------------------------------------
      Note: Pruning is interpretive relabeling, not re-estimation. Flagged nodes
      remain in the object with all edges preserved.
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

# print snapshot: a non-converged level shows the red cross

    Code
      snap_print(x)
    Message
      
      -- Bass-Ackwards Analysis (ackwards) -------------------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 3
      
      -- Levels --
      
      v k = 1: 1 factor, 28.2% variance
      x k = 2: 2 factors, 46.5% variance
      v k = 3: 3 factors, 57.5% variance
      
      -- Edges --
      
      5 of 8 edges have |r| >= 0.3
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

# print/summary snapshot: near-singular caution

    Code
      snap_print(x)
    Message
      
      -- Bass-Ackwards Analysis (ackwards) -------------------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 3
      
      -- Levels --
      
      v k = 1: 1 factor, 28.2% variance
      v k = 2: 2 factors, 46.5% variance
      v k = 3: 3 factors, 57.5% variance
      
      -- Edges --
      
      5 of 8 edges have |r| >= 0.3
      ! Near-singular correlation matrix (min eigenvalue 0.0012): fit indices and
      factor scores may be unreliable. See `?ackwards` ("When to trust the result").
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

---

    Code
      snap_print(summary(x))
    Message
      
      -- Summary: Bass-Ackwards Analysis (ackwards) ----------------------------------
      Engine: pca
      Rotation: varimax
      Basis: pearson
      n: 1,000
      k (max): 3
      
      -- Levels --
      
      k = 1: 1 factor (28.2% cumulative variance)
      m1f1 28.2% eigenvalue 4.51
      
      k = 2: 2 factors (46.5% cumulative variance)
      m2f1 23.3% eigenvalue 4.51
      m2f2 23.2% eigenvalue 2.93
      
      k = 3: 3 factors (57.5% cumulative variance)
      m3f1 23.0% eigenvalue 4.51
      m3f2 17.5% eigenvalue 2.93
      m3f3 16.9% eigenvalue 1.76
      
      -- Lineage (primary parents) --
      
      m1f1 > m2f1, m2f2
      m2f1 > m3f2, m3f3
      m2f2 > m3f1
      
      ! Near-singular correlation matrix (min eigenvalue 0.0012): per-level fit
      indices and factor scores may be unreliable -- the solution rests on a
      rank-deficient matrix. See `?ackwards` ("When to trust the result").
      --------------------------------------------------------------------------------
      Note: This is a series of linked solutions, not a fitted hierarchical model.
      Cross-level edges are descriptive score correlations. Per-level fit indices
      (EFA/ESEM) describe how well a k-factor model fits the items at that level --
      they do not validate the edges or the hierarchy itself.
    Output
      

