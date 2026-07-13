# M58: Consolidate input-validation helpers & fix two drift bugs

- **Status:** done · **PR:** https://github.com/jmgirard/ackwards/pull/58 (merged 2026-07-13)
- **Priority:** normal · **Depends on:** — · **Principles touched:** —

## Outcome

Removed the input-validation copy-drift the 2026-07-12 audit flagged and fixed **three**
instances of the drift-bug class it produced (one found in review): `suggest_k(n_obs=NA)`,
`ackwards(R, n_obs=NA)`, and `prune(x, <typo>)` now error cleanly instead of swallowing the
arg / dying with base R's "missing value where TRUE/FALSE needed".

New `utils.R` helpers `.as_numeric_matrix()` + `.check_count()` (min-aware, `is.na` guard)
replace 5 coercion triplets + 5 count checks across the verbs. `.tucker_phi` moved to
`utils.R` (shared by prune + comparability); `.node_levels()` unifies prune's 3 label→level
rebuilds; `suggest_k`'s duplicated `k_max` block runs once; `prune.ackwards` `...` now guarded.
Net −40 lines, +5 tests. dod-gate: check 0/0/0, coverage 100%, style/lint/pkgdown clean.

## Key decisions
- `.check_unknown_dots` doc corrected: package-generic methods with reserved dots guard
  (prune/boot_edges); base/tidy generics (print/autoplot) stay permissive.
- Finding 1 (ackwards n_obs=NA twin bug, score 84) fixed on-branch, confined to the numeric
  R-matrix branch. Finding 2 (k_max hoist reorders a "CD skipped" note before an invalid-k_max
  abort; score 66) logged, not actioned — cosmetic, error-path only.
