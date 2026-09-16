# M89: Real within-level factor correlations and Φ-partialled edge reporting

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP2, IP3, IP4, IP8, GP2
- **Resolves:** —
- **Surface tier:** user-facing, because it adds `tidy()` columns, a `tidy()` table, and a `summary()` block
- **Branch/PR:** `m089-factor-cor-partialled-edges`

## Goal

Carry each level's real within-level factor correlation through column reordering and sign alignment, show it, and report the Φ-partialled edge quantities (`beta`, `r2`) beside the marginal edges, with default varimax output unchanged.

## Scope

**In:** one shared helper that permutes and sign-flips a within-level factor correlation. PCA and EFA levels read psych's `$Phi` when the fit has one. The ESEM level permutes lavaan's `cor.lv` by its variance sort. `ackwards()` flips `factor_cor` in step with `align_signs`. A `tidy(what = "factor_cor")` table and a `summary()` block. A `beta` column on `tidy(what = "edges")` and an `r2` column on `tidy(what = "variance")` (RR02 Q1, Q2, recommendation 4). A baseline fixture generated on master that proves default output unchanged. Removal of the DESIGN Known-limitations entry on `factor_cor` and `ord`.

**Out:** a user-facing `rotation` argument, oblique scoring weights, oblique variance, the fit-time advisory, and the Forbes oblique oracle go to M90, which is blocked on the design-session decision. An internal oblique switch on the engines also goes to M90, because RR02 recommendation 7 rules out a rotation switch without the scoring and variance kit, even an internal one. Tucker's φ as the documented rotation-consistency diagnostic stays a candidate row gated on the design session. Primary-parent matching, sign anchoring, and the diagram keep the marginal edge. Nothing here changes lineage.

## Acceptance criteria

- [ ] AC1: Default-argument output is unchanged against a committed baseline. `data-raw/baseline-m89.R`, run on master at commit `a645e32` before any code change, stores loadings, weights, `variance`, `factor_cor`, and every edge matrix for default fits of `pca`, `efa`, and `esem` on `sim16` (`k_max = 4`) and of `pca` and `efa` with `cor = "polychoric"` on `bfi25` (`k_max = 3`). The fixture `tests/testthat/fixtures/baseline-m89.rds` carries a `provenance` attribute that names the generator and the commit (IP8) and is registered in `cairn/ORACLES.md`. A test reproduces every stored value within 1e-12, and `tests/testthat/test-forbes-fidelity.R` passes without edit.
- [ ] AC2: One internal helper, `.carry_factor_cor(Phi, ord, signs)`, turns an engine's within-level factor correlation, a column permutation, and a per-factor sign vector into the stored `factor_cor`. A direct test asserts that the result equals `Phi[ord, ord] * tcrossprod(signs)` for k in {2, 4}, for a self-inverse and a non-self-inverse permutation, and for a mixed and an all-negative sign vector. A routing test replaces the helper through `testthat::local_mocked_bindings()` with one that stamps its output, and asserts that a fit through each of the three engines carries the stamp on every level's `factor_cor`.
- [ ] AC3: PCA and EFA levels take Φ from psych's `$Phi` when the fit carries one and the identity otherwise, through a helper `.engine_phi(fit, k)`. A direct test runs it on a fit-shaped list with a non-identity `$Phi`, on one with `$Phi = NULL`, and on a k = 1 fit. The ESEM level's stored `factor_cor` is lavaan's `cor.lv` matrix permuted by `ord` through AC2's helper. `ackwards()` flips the rows and columns of every level's `factor_cor` by the `align_signs` sign vectors, asserted on a copy of a default fit whose level-3 `factor_cor` is non-identity and whose sign vector has one -1.
- [ ] AC4: `tidy(what = "factor_cor")` returns one row per unordered factor pair within each level (columns `level`, `factor_a`, `factor_b`, `cor`) and zero rows with those columns at k = 1. On a default fit every `cor` is 0 within 1e-12. On a copy of a default fit whose level-3 `factor_cor` is replaced by a non-identity matrix with one negative entry, it returns those entries. `summary()` prints a "Within-level factor correlations" block only when some off-diagonal |cor| exceeds 1e-8. Snapshots assert this on the default fit, on a copy with max |cor| = 5e-9 (silent), and on a copy with max |cor| = 2e-8 (prints).
- [ ] AC5: `tidy(what = "edges")` carries `beta`, the Φ-partialled coefficient `B = Φ_s^{-1} E`, where `Φ_s = D^{-1/2} W_a' R W_a D^{-1/2}` comes from the stored weights of the shallower level of each pair (RR02 Q1). This holds for adjacent pairs and, under `pairs = "all"`, for skip-level pairs. `tidy(what = "variance")` carries `r2 = E_j' Φ_s^{-1} E_j` from the adjacent level above for every factor at k >= 2, and `NA` at k = 1. On a default fit `beta` equals `r` and `r2` equals the column sum of `r^2` within 1e-10 for all three engines. On a hand-built case with correlated composites, `beta` matches the standardized coefficients of `lm()` regressing each level-b score on all level-a scores, and `r2` matches that model's R², both within 1e-8. When `Φ_s` cannot be inverted, both columns are `NA` for that pair and a cli warning names the level.
- [ ] AC6: `?tidy.ackwards`, `?summary.ackwards`, and `NEWS.md` describe `beta`, `r2`, and the factor-correlation block, and state that under varimax `beta` equals `r` and `r2` equals the sum of squared `r`. DESIGN.md's Known-limitations entry on `factor_cor` and `ord`, and the matching guard comment in `R/engine_esem.R`, are removed. There is no new export, so `_pkgdown.yml` is unchanged and `pkgdown::check_pkgdown()` passes.
- [ ] AC7: `Rscript tools/dod-gate.R` passes (the profile's verify and check gate: freshness, check, coverage, style, lint, pkgdown).

## Coverage

- AC1 → T1
- AC2 → T2, T3
- AC3 → T2, T3, T4
- AC4 → T5
- AC5 → T6
- AC6 → T7
- AC7 → T7

## Tasks

- [x] T1: Before you touch `R/`, write `data-raw/baseline-m89.R` (default fits per AC1, `provenance` attr with generator path and `git rev-parse HEAD`). Run it on master `a645e32`. Commit the `.rds` under `tests/testthat/fixtures/`. Add the frozen oracle row to `cairn/ORACLES.md`. Write `tests/testthat/test-baseline-m89.R` (1e-12). It must pass on the unchanged code first.
- [x] T2: Add `.carry_factor_cor(Phi, ord, signs)` and `.engine_phi(fit, k)` to `R/utils.R` next to `.variance_explained()` (`R/utils.R:97`), with the AC2 and AC3 direct tests.
- [x] T3: Route the engines. `R/engine_pca.R:58` and `R/engine_efa.R:143` replace `diag(k)` with `.engine_phi(fit, k)`. `R/engine_esem.R:193-200` and `:294-307` permute `cor.lv` by `ord` through the helper and drop the guard comment. Add the `local_mocked_bindings()` routing test.
- [x] T4: In `ackwards()` after `.align_signs()` (`R/ackwards.R:922-929`, where loadings and weights are flipped), flip each level's `factor_cor` through the helper with `ord` set to the identity. Test per AC3.
- [x] T5: Add `"factor_cor"` to `tidy()`'s `what` (`R/tidy.R:117`) with `.tidy_factor_cor()`, and the conditional `summary()` block (`R/summary.R`, `print.summary_ackwards`). Snapshots per AC4.
- [x] T6: Add `.partialled_edges(W_a, R, E)` that returns `beta` and per-column `r2`, with `NA` and one cli warning on a failed `solve()`. Wire `beta` into `.tidy_edges()` (`R/tidy.R:186`) for every stored pair and `r2` into `.tidy_variance()` (`R/tidy.R:323`). Tests per AC5: identity on three engines, `lm()` oracle on a hand-built correlated-composite case, singular path.
- [x] T7: Roxygen for `tidy()` and `summary()`, `NEWS.md`, DESIGN Known-limitations removal, ORACLES rows for the `lm()` live oracle and the varimax identity invariant. Run `devtools::document()` and `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan. Collision sweep: D-034 (supersedes D-002) makes oblique a gated option and leaves this Φ-reporting foundation ungated (RR02 recommendation 4). Absorbs the candidate row "Φ-partialled edge decomposition as reporting". The row stays until completion.
- 2026-09-16: criteria audit ran in full mode (fresh-context [O] reader, 16 findings over M89 and M90, 1 clean). M89 fixes: dropped an internal engine `rotation` switch (it would build the plumbing-only object D-034 rules out). Moved `r2` from the prune-only `nodes` tidier to `variance`. Replaced the PCA-only fidelity suite with a three-engine baseline fixture as AC1's enumerating procedure. Replaced grep-count and test-diff clauses with behavioural tests. Added the `lm()` R² oracle for `r2`. Named the factor-correlation column `cor` because the package already uses `phi` for Tucker congruence. Probed the summary print gate at its threshold. Defined `beta` for skip-level pairs and the singular-Φ path.
- 2026-09-16: plan gate chose columns on existing tidiers (`beta` on edges, `r2` on variance) over a new `what = "partialled"` tidier because the varimax identity is then visible beside `r`. Falsified by users reporting the default edge table as harder to read once `beta` sits beside `r`.
- 2026-09-16: plan gate chose a helper-level unit test plus mocked-binding routing over an internal oblique engine switch for exercising the Φ carry, because the switch produces silently wrong `variance` and weights (RR02 Q5 items 4 and 5) even when unexposed. Falsified by a Φ-carry defect the helper tests cannot reach but an end-to-end oblique fit would (M90's fixtures then catch it).
- 2026-09-16: plan gate chose a master-generated regression fixture over "the fidelity suite passes" as the unchanged-output oracle because the fidelity suite is PCA-only while M89 edits all three engines. Falsified by a default-output change the fixture's five fits do not cover.
- 2026-09-16: /milestone-implement started. Branch cut from master `ee71316` (code tree equal to `a645e32`; the two later commits touched only `cairn/`). No question gate: the plan left no open implementation choice. Column placement chosen in session: `beta` after `r` on the edge table, `r2` after `cumulative` on the variance table.
- 2026-09-16: T1 done. Generator `data-raw/baseline-m89.R`, fixture `baseline-m89.rds` (O13 in ORACLES.md), test `test-baseline-m89.R` passes on the unchanged code (5 fits, 0 failures). Minor amendment to AC1's enumerating procedure: the ESEM fit carries `seed = 1`, because an unseeded ESEM fit differs run to run at about 1e-6 (lavaan rotation random starts, observed twice this session) and cannot serve a 1e-12 oracle. The seed changes no method argument. AC1 wording unchanged.
- 2026-09-16: T2 done. `.carry_factor_cor()` and `.engine_phi()` in `R/utils.R`; direct tests in `test-utils.R` (five permutation/sign cases over k in {2, 4}, psych `$Phi` present/absent/k = 1). psych sets `$Phi` only under oblique rotation (checked on `pca()` and `fa()` with varimax and oblimin this session).
- 2026-09-16: T3 done. PCA/EFA levels build `factor_cor` from `.engine_phi()` through the carry helper; ESEM permutes lavaan's `cor.lv` by `ord` through it and the guard comment is gone. PCA/EFA `factor_cor` stays unnamed so the baseline (O13) compares attribute-for-attribute; tidy() reads level labels. Routing test in `test-factor-cor.R` (mocked helper stamps every level on all three engines). Baseline, efa, esem, pca files: 0 failures.
- 2026-09-16: T4 done. `ackwards()` carries each level's `factor_cor` through the helper with the `align_signs` vectors after the weights flip. Test plants a non-identity Phi at k = 3 (mocked `.engine_phi`) and one -1 in the level-3 sign vector (wrapped `.align_signs`), and asserts the stored matrix equals `Phi * tcrossprod(signs)` and that the loadings really flipped. factor-cor and baseline files: 0 failures.
- 2026-09-16: T5 done. `tidy(what = "factor_cor")` (`.tidy_factor_cor()`, upper-triangle pair order, label columns when factor labels are set) and a `summary()` block gated at max |cor| > 1e-8, one line per pair for each level that has a correlated pair. Tests: default fit all 0, k = 1-only copy zero rows, planted level-3 matrix; snapshots at 5e-9 (silent) and 2e-8 (prints). Existing snapshots unchanged. `k_max = 1` is illegal in `ackwards()`, so the k = 1 case runs on a copy holding only level 1.
- 2026-09-16: T6 done. `.partialled_edges(W_a, R, E, level)` in `R/utils.R` (the `level` argument feeds the warning text); `beta` after `r` on the edge table for every stored pair, `r2` after `cumulative` on the variance table. `test-partialled-edges.R`: varimax identity on pca/efa/esem within 1e-10, skip-level pairs under `pairs = "all"`, lm() oracle on mixed level-2 weights within 1e-8 (beta differs from r there), singular level-2 weights give NA and one warning naming k = 2. Files reading edge/variance tables (print, factor-labels, boot_edges, baseline): 0 failures.
- 2026-09-16: T7 done. Roxygen for `tidy()` (`beta`, `r2`, `"factor_cor"`, label columns, example) and `summary()` (the block and its 1e-8 gate); two NEWS entries; DESIGN Known-limitations entry on `factor_cor`/`ord` removed; ORACLES rows O14 (varimax identity, invariant) and O15 (`lm()` oracle, live). Intro vignette gained two sentences on `beta` and `r2`; all precomputed vignettes regenerated (the baked edge and variance tables now show the new columns). `tools/check-prose.R` clean after shortening four sentences. First `dod-gate.R` run: freshness clean, check 0/0/0, coverage 100%, lint clean, pkgdown complete; failed only on styler restyling uncommitted files, which this commit carries.
- 2026-09-16: second `Rscript tools/dod-gate.R` run on the committed tree: GATE PASSED (freshness, ledger anchors, CI path filters, prose, check 0/0/0, coverage 100%, styler clean, lintr clean, pkgdown index complete).
- 2026-09-16: claim audit: 49 claims read, 0 corrected — NEWS.md, R/tidy.R, R/summary.R, R/utils.R, R/ackwards.R, R/engine_pca.R, R/engine_efa.R, R/engine_esem.R, data-raw/baseline-m89.R, vignettes/ackwards-intro.Rmd.orig, tests. One claim marked unverifiable by the reader (the ~1e-6 unseeded ESEM run-to-run figure in the generator header); it is a dated observation this session made twice before writing it. Reader's incidental note: `.tidy_edges()` would error if a stored edge-matrix cell ever lacked a row in `edges$tidy`; every construction path writes both from one matrix set, so no change made.
- 2026-09-16: all tasks checked, gate clean; status set to review.

## Decisions

## Review
