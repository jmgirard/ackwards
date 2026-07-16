# M60: De-duplicate the setup path (audit bucket 3) — done 2026-07-16

**Goal:** remove the remaining readability-only duplication and dead code in the setup path
(engines + `ackwards()` + comparability/boot_edges) flagged as bucket 3 of the 2026-07-12
audit, zero behavior change. Third and final audit bucket (M58 validation, M59 console/plot).

**Outcome (PR https://github.com/jmgirard/ackwards/pull/60, squash 3954c75):**
- `utils.R` helpers `.variance_explained()` / `.score_var()` are the sole computation sites
  (3 engines + `compute_edges()`'s D^{-1/2}); ESEM ordering key = raw `colSums(L^2)`.
- `.align_signs()` returns loadings+signs only; `compute_edges(build_tidy = FALSE)` for
  matrices-only callers (lineage pass, `.cross_cor()`, `.boot_replicate()`).
- One `eigen()` per input path: `.near_singular_check(min_eig)`; `.corfiml_R()` /
  `.validate_cor_matrix()` return `list(R, min_eig)`.
- Dead code removed: always-true `k_eff >= 1L` guard, reader-less `$meta$converged_levels`;
  `deepest_converged = k_eff` (exact — all engines truncate at first failure, Invariant 7).
- Shared muffled dispatch `.fit_levels_muffled()` for comparability/boot_edges refits.

**Key decisions:** AC6's "net line reduction" clause dropped at a user-approved mini-gate —
documented helper extraction is line-neutral (+9 code lines); the real outcome is
grep-verified single computation sites. Two candidates deliberately not extracted
(positive-manifold flip, shortfall reporter) stay on the ROADMAP candidate list.

**Verification:** suite 2295/0/0/0; dod-gate clean (check 0/0/0, coverage 100%); zero
snapshot diffs; Forbes-fidelity + algebra-vs-scores oracles green; 3-lens independent
review (diff-bug / blame-history / prior-PR) — 0 findings; CI 7/7 green.
