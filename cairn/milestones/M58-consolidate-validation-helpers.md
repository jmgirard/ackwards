<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section. -->
# M58: Consolidate input-validation helpers & fix two drift bugs

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Principles touched:** —   <!-- invariants not yet numbered as IP/GP; this milestone preserves them (behavior unchanged except the two documented bug fixes) -->
- **Branch/PR:** m58-consolidate-validation-helpers

## Goal

Extract the input-validation boilerplate that the four big verbs each copied (and drifted) into shared `utils.R` helpers, route every call site through them, and fix the two latent bugs that copy-drift produced.

## Scope

**In:**
- New `utils.R` helpers: `.as_numeric_matrix(data, arg)` (the `!is.data.frame && !is.matrix → abort; as.matrix; !is.numeric → abort` triplet) and `.check_count(x, arg, min = 1L)` (the single-positive-integer check, **with an `is.na()` guard**).
- Route every current duplicate through them: matrix coercion at `suggest_k.R`, `comparability.R:167`, `boot_edges.R:170`, `factorability.R:84`, `ackwards.R`; count checks at `suggest_k.R:186` (n_iter), `:217` (n_obs), `comparability.R:152` (n_splits), `boot_edges.R:152` (n_boot), `factorability.R:72` (n_obs).
- **Bug fix 1:** `suggest_k(R, n_obs = NA_real_)` currently reaches `if (... || NA || ...)` and dies with base R's "missing value where TRUE/FALSE needed"; routing through `.check_count`'s `is.na` guard makes it error with the intended "must be a positive integer."
- **Bug fix 2:** `prune.ackwards()` accepts `...` (prune.R:686) but never calls `.check_unknown_dots()` as its siblings do — a mistyped argument is silently swallowed. Add the guard.
- Move `.tucker_phi` from `prune.R:39` to `utils.R` (it is a cross-file utility used by `comparability.R:372`); add `.node_levels(levels_list)` and route `prune.R`'s three divergent label→level rebuilds (`:165`, `:229`, `:443`) through it.
- Hoist `suggest_k.R`'s byte-identical `k_max` default/validate block (`:223-230` ≡ `:293-300`) to run once after the input branch.
- Fold in the stale dangling comment at `ackwards.R:405` ("handled below" — it is not handled there).

**Out:**
- print/summary/autoplot output-layer dedup → M59.
- Readability-only core/engine dedup (per-engine variance/score-var/flip helpers; `ackwards()` doubly-built tidy tibble + discarded `.align_signs()$edges`; `eigen(R)` threading; dead `k_eff>=1L` guard + vestigial convergence meta; `.fit_levels_muffled()` + convergence-shortfall reporter for comparability/boot_edges) → "Remaining dedup pass" candidate.

## Acceptance criteria

- [ ] `.as_numeric_matrix()` and `.check_count()` exist in `utils.R` with direct unit tests covering the edge cases (NA, non-integer, `length != 1`, `< min`); no inline coercion-triplet or hand-rolled positive-integer check remains in the four verbs (grep-verifiable).
- [ ] `suggest_k(R, n_obs = NA_real_)` errors with the "positive integer" message — regression test that **fails before the fix**.
- [ ] `prune(x, <misspelled_arg> = ...)` errors via `.check_unknown_dots()` — regression test that **fails before the fix**.
- [ ] `.tucker_phi` is defined in `utils.R` (not `prune.R`); `.node_levels()` backs prune's three former rebuilds; the duplicated `k_max` block in `suggest_k` runs once — behavior-preserving, with existing `prune`/`comparability`/`suggest_k` tests unchanged and green.
- [ ] `Rscript tools/dod-gate.R` clean (check 0 err/0 warn/0 note, coverage, style, lint, pkgdown).

## Coverage

- AC1 → T1, T2
- AC2 → T1, T2
- AC3 → T3
- AC4 → T4
- AC5 → T5

## Tasks

- [x] T1 — Add `.as_numeric_matrix(data, arg)` and `.check_count(x, arg, min = 1L)` (with `is.na` guard) to `utils.R`; write `test-utils.R` cases for both, including the NA / non-integer / length / min-violation branches.
- [x] T2 — Replace the ~5 matrix-coercion sites and ~5 count-check sites with calls to the new helpers; add the `suggest_k(R, n_obs = NA_real_)` regression test (assert it fails on the current code first, then passes).
- [ ] T3 — Add `.check_unknown_dots(...)` to `prune.ackwards()`; add a regression test that an unknown argument to `prune()` errors (fails first).
- [ ] T4 — Move `.tucker_phi` to `utils.R`; add `.node_levels(levels_list)` and route prune's three rebuilds; hoist the duplicated `suggest_k` `k_max` block; fix the `ackwards.R:405` dangling comment.
- [ ] T5 — `Rscript tools/dod-gate.R` clean; commit with tracking.

## Work log

- 2026-07-12: created by /milestone-plan (bucket 1 + the two latent bugs of the 2026-07-12 codebase audit; scope confirmed at the plan gate).
- 2026-07-12: T1 — added `.as_numeric_matrix()` + `.check_count()` (min-aware; `is.na` guard) to utils.R with direct test-utils.R cases; green.
- 2026-07-12: T2 — routed 5 coercion + 4 count sites (ackwards/suggest_k/comparability/boot_edges/factorability) through the helpers (−40 net lines); added the `suggest_k(n_obs = NA)` regression test. check_items excluded (coerces to data.frame + screens non-numeric, not a numeric-matrix site). Pre-fix failure confirmed by the read buggy code at suggest_k.R:217 (`if(... || NA ...)` → base-R error). 724 pass / 0 fail across affected files.

## Decisions

## Review
