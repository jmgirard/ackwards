# M097: Cross-check coverage and the matrix behind each edge

**Status:** done (2026-10-01, PR #106 https://github.com/jmgirard/ackwards/pull/106)

**Goal:** The record of which correlation matrix each edge is exact for, and of which fits the algebra-vs-scores cross-check certifies, matches the code.

**Outcome:** The `ackwards()` and `compute_edges()` help say that each edge is exact for the stored `x$r`. They name the two ESEM settings in which lavaan fits another matrix. One is Spearman with pairwise or listwise deletion. The other is Pearson pairwise on data with missing values. There, ML and MLR fit the complete rows, and WLSMV fits pairwise covariances. The `@param missing` listwise bullet and the varimax bullet were brought into line. Three tests in `tests/testthat/test-compute_edges.R` cover Spearman fits scored from column ranks (plain and tied) and listwise fits scored from complete rows. They run on PCA, EFA, and ESEM over all level pairs at 1e-10. DESIGN's first Known limitations entry lists the four uncovered settings and names the tests. DESIGN names `x$r` at each "exact" edge claim, IP1 included. NEWS entry added. No code or result changed.

**Decisions:** D-039 (IP1 names the matrix each edge is exact for). The plan gate kept README, the vignettes, the manuscript, `?prune`, and `?tidy` on the general wording. AC3 was amended at a mini gate to drop ULSMV.

**Review:** one pass, and all five criteria passed with fresh evidence. A fresh lavaan 0.7.2 sweep of 40 ESEM settings confirmed the help's list. Three lenses gave 16 findings, and none showed a criterion failing. Eight were fixed before merge. They covered two stale `ackwards()` help bullets, four DESIGN "the fit uses" sites, and markers on the M76 and M43 notes. They also covered a "true correlation" overclaim, two Known limitations sentences, and a dangling "instead". One older help bullet went to the ULSMV hotfix row. Six were rejected and one was noted. Two `[high]` hotfix rows came out of T3 (ULSMV pairwise, Spearman with FIML). Merged on local green under the repo's standing non-release rule.
