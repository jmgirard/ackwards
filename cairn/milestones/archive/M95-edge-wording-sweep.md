# M95: Edge wording sweep across shipped docs and DESIGN §5

**Status:** done (2026-10-01, PR #104 https://github.com/jmgirard/ackwards/pull/104)

**Goal:** Each edge description that the search below finds names the matrix the edge is exact for and credits Waller (2007) only with his own results.

**Outcome:** A paragraph-mode search over the shipped docs found 29 edge descriptions. Each was reworded to say that an edge is exact for, or implied by, the correlation matrix the fit uses. The sites are DESCRIPTION, README, the manuscript, and the help for `ackwards()`, `tidy()`, `prune()`, and `compute_edges()`. Five vignettes changed: engines, intro, ordinal, girard, and forbes. DESCRIPTION no longer offers materialized scores as an edge route (D-038), and it credits Waller with the principal-components result only. The engines PCA paragraph states that result through his transformation matrices and writes `W′RW` at every algebra mention. The `compute_edges()` help gives the scores-branch conditions from the code and says no shipped caller reaches that branch. DESIGN §5.1 cites Eq. 14 and §3 only and names ten Berge with the regression fallback. The §5.2 `NULL if !linear` comments are gone. The §5.3 pseudocode matches `formals(compute_edges)`, with the `"algebra"` abort and no sign-alignment step. A NEWS entry was added. No R code changed.

**Decisions:** none cross-cutting. The plan gate kept the `print()` note's "descriptive score correlations" wording under GP3, and it put DESIGN §5 in this milestone.

**Review:** two passes. Pass 1 failed AC4: two engines algebra mentions did not write `W′RW` (defect return 1, fixed as T8). Pass 2 re-ran all seven criteria and the DoD gate, and all passed. Three lenses gave 18 findings, and none showed a criterion failing. Nine wording slips were fixed before merge (R2-R7, R10-R12). R1 (ESEM matrix ambiguity) joined the cross-check-paths row. R9 (unqualified "exact" in DESIGN outside §5) became a new candidate row. Seven were rejected. Merged on local green under the repo's standing non-release rule.
