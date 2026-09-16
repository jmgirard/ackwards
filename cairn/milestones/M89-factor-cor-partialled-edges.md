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

- [x] AC1: Default-argument output is unchanged against a committed baseline. `data-raw/baseline-m89.R`, run on master at commit `a645e32` before any code change, stores loadings, weights, `variance`, `factor_cor`, and every edge matrix for default fits of `pca`, `efa`, and `esem` on `sim16` (`k_max = 4`) and of `pca` and `efa` with `cor = "polychoric"` on `bfi25` (`k_max = 3`). The fixture `tests/testthat/fixtures/baseline-m89.rds` carries a `provenance` attribute that names the generator and the commit (IP8) and is registered in `cairn/ORACLES.md`. A test reproduces every stored value within 1e-12, and `tests/testthat/test-forbes-fidelity.R` passes without edit.
- [x] AC2: One internal helper, `.carry_factor_cor(Phi, ord, signs)`, turns an engine's within-level factor correlation, a column permutation, and a per-factor sign vector into the stored `factor_cor`. A direct test asserts that the result equals `Phi[ord, ord] * tcrossprod(signs)` for k in {2, 4}, for a self-inverse and a non-self-inverse permutation, and for a mixed and an all-negative sign vector. A routing test replaces the helper through `testthat::local_mocked_bindings()` with one that stamps its output, and asserts that a fit through each of the three engines carries the stamp on every level's `factor_cor`.
- [x] AC3: PCA and EFA levels take Φ from psych's `$Phi` when the fit carries one and the identity otherwise, through a helper `.engine_phi(fit, k)`. A direct test runs it on a fit-shaped list with a non-identity `$Phi`, on one with `$Phi = NULL`, and on a k = 1 fit. The ESEM level's stored `factor_cor` is lavaan's `cor.lv` matrix permuted by `ord` through AC2's helper. `ackwards()` flips the rows and columns of every level's `factor_cor` by the `align_signs` sign vectors, asserted on a copy of a default fit whose level-3 `factor_cor` is non-identity and whose sign vector has one -1.
- [x] AC4: `tidy(what = "factor_cor")` returns one row per unordered factor pair within each level (columns `level`, `factor_a`, `factor_b`, `cor`) and zero rows with those columns at k = 1. On a default fit every `cor` is 0 within 1e-12. On a copy of a default fit whose level-3 `factor_cor` is replaced by a non-identity matrix with one negative entry, it returns those entries. `summary()` prints a "Within-level factor correlations" block only when some off-diagonal |cor| exceeds 1e-8. Snapshots assert this on the default fit, on a copy with max |cor| = 5e-9 (silent), and on a copy with max |cor| = 2e-8 (prints).
- [x] AC5: `tidy(what = "edges")` carries `beta`, the Φ-partialled coefficient `B = Φ_s^{-1} E`, where `Φ_s = D^{-1/2} W_a' R W_a D^{-1/2}` comes from the stored weights of the shallower level of each pair (RR02 Q1). This holds for adjacent pairs and, under `pairs = "all"`, for skip-level pairs. `tidy(what = "variance")` carries `r2 = E_j' Φ_s^{-1} E_j` from the adjacent level above for every factor at k >= 2, and `NA` at k = 1. On a default fit `beta` equals `r` and `r2` equals the column sum of `r^2` within 1e-10 for all three engines. On a hand-built case with correlated composites, `beta` matches the standardized coefficients of `lm()` regressing each level-b score on all level-a scores, and `r2` matches that model's R², both within 1e-8. When `Φ_s` cannot be inverted, both columns are `NA` for that pair and a cli warning names the level.
- [x] AC6: `?tidy.ackwards`, `?summary.ackwards`, and `NEWS.md` describe `beta`, `r2`, and the factor-correlation block, and state that under varimax `beta` equals `r` and `r2` equals the sum of squared `r`. DESIGN.md's Known-limitations entry on `factor_cor` and `ord`, and the matching guard comment in `R/engine_esem.R`, are removed. There is no new export, so `_pkgdown.yml` is unchanged and `pkgdown::check_pkgdown()` passes.
- [x] AC7: `Rscript tools/dod-gate.R` passes (the profile's verify and check gate: freshness, check, coverage, style, lint, pkgdown).

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

Review pass 2026-09-16 on branch `m089-factor-cor-partialled-edges` (8 commits ahead of `master` at `ee71316`; `master` had not moved, so no sync merge was needed). Every result below is fresh from this pass.

### Acceptance-criterion evidence

- AC1: `test-baseline-m89.R` reran fresh: 5 tests (pca/efa/esem on sim16, pca/efa polychoric on bfi25), 0 failures, every stored value within 1e-12. `test-forbes-fidelity.R` unchanged by the diff: 5 tests, 0 failures. Fixture `provenance` attribute names generator `data-raw/baseline-m89.R`, tolerance 1e-12, and commit `ee71316`; registered as O13 in `cairn/ORACLES.md`. **Two deviations from the criterion text, recorded verbatim, not reinterpreted:** (a) the text says commit `a645e32`; the fixture was generated at `ee71316`, and `git diff --stat a645e32 ee71316 -- . ':!cairn'` is empty (the one intervening commit touched only `cairn/`), so the code tree is the one the criterion names; (b) the ESEM fit carries `seed = 1`, which the criterion's "default fits" does not name (work-log line of T1 explains why: unseeded lavaan rotation starts differ run to run at ~1e-6). Ticked on the evidence, with both deviations presented at the gate for the maintainer's disposition.
- AC2: `test-utils.R` reran: 42 tests, 0 failures, among them the direct `.carry_factor_cor()` test over k in {2, 4}, self-inverse (`c(2,1)`, `c(4,3,2,1)`) and cycle (`c(2,3,4,1)`) permutations, mixed and all-negative sign vectors, asserting `Phi[ord, ord] * tcrossprod(signs)`. `test-factor-cor.R` reran: 5 tests, 0 failures, among them the `local_mocked_bindings()` routing test that stamps the helper's output and finds the stamp on every level of pca, efa, and esem fits.
- AC3: `test-utils.R` `.engine_phi()` test passes on a fit-shaped list with a non-identity `$Phi`, with `$Phi = NULL`, with no `$Phi` slot, and at k = 1. `R/engine_esem.R` builds `factor_cor` from `lavInspect(fit, "cor.lv")` permuted by `ord` through the helper (diff read). `test-factor-cor.R` sign-flip test passes: planted non-identity level-3 Phi plus one −1 in the level-3 sign vector yields `Phi * tcrossprod(c(1, -1, 1))` and the flipped factor's loadings are negated.
- AC4: `test-factor-cor.R` tidy tests pass: columns `level, factor_a, factor_b, cor`; choose(k, 2) rows per level; zero typed rows on a k = 1-only copy; all |cor| < 1e-12 on a default fit; planted level-3 entries returned in upper-triangle order. `test-print-snapshot.R` run with `NOT_CRAN=true`: 10 tests, 0 failures, 0 skipped, including the default fit (no block), the 5e-9 copy (silent) and the 2e-8 copy (block prints).
- AC5: `test-partialled-edges.R` reran: 4 tests, 0 failures. Varimax identity `beta == r` and `r2 == colSums(r^2)` within 1e-10 on pca, efa, esem; skip-level `beta` present and equal to `r` under `pairs = "all"`; `lm()` oracle (O15) matches standardized coefficients and R² within 1e-8 on the mixed-weights case where `beta` differs from `r`; singular level-2 weights give `NA` and a warning naming `k = 2`, exactly one warning from the helper.
- AC6: diff read: `?tidy.ackwards` documents `beta`, `r2`, `"factor_cor"`, the varimax identities, and the `NA`-with-warning path; `?summary.ackwards` documents the block and its 1e-8 gate; `NEWS.md` carries two entries stating the identities. DESIGN Known-limitations entry on `factor_cor`/`ord` removed; `engine_esem.R` guard comment removed. `NAMESPACE` and `_pkgdown.yml` unchanged in the diff; `pkgdown::check_pkgdown()` passed inside the gate.
- AC7: `Rscript tools/dod-gate.R` (TESTTHAT_CPUS=8) on the committed tree: GATE PASSED — vignette-freshness clean, prose clean, check 0 errors / 0 warnings / 0 notes, coverage 100.00%, styler clean, lintr clean, pkgdown reference index complete. Working tree clean afterwards.

Driving RR: none, so no projection-vs-outcome pairs.

### Consistency gate

- `cairn_validate.py`: all checks passed (exit 0); 17 advisory warnings, all pre-existing (M90 sizing tripwire, M84 wrapped work-log lines).
- No IP/GP changed, so `cairn_impact.py` was skipped.
- `devtools::document()` produced no diff. `README.Rmd`/`README.md` untouched by the diff. `pkgdown::check_pkgdown()` clean. NEWS has the milestone's two entries with no milestone numbers. `data-raw` already `.Rbuildignore`d; check reported 0 notes.

### Independent review (three lenses, fresh context)

- [S] blame-history: zero findings. Verified the removed ESEM guard comment and labelling (relocated to `.label_phi()`), the `diag(k)` replacement, the M47 boot-column append order, D-018's 0–1 scale, D-034/D-035.
- [S] prior-review-record: "no prior-review evidence" — archived Review sections empty on these files, LESSONS silent, PR-comment probe returned an empty array.
- [O] diff-bug: 16 ranked findings; algebra, carry sites, join key, boot-column order, `r2` adjacency rule, label columns, and helper signatures all confirmed correct. Findings and dispositions:
  1. PCA/EFA carry assumes psych sorts `$Phi` with the loadings; unreachable by any test while only varimax runs. — **follow-up → M90** (its oblique fixtures exercise a real non-identity `$Phi`; work-log line added to M90).
  2. `.partialled_edges()` builds Φ_s from weights and R while `E` may come from materialized scores under `edge_method = "scores"`; roxygen does not say so. — **fix now** (one roxygen sentence on `tidy()`; DESIGN Known-limitations algebra-vs-scores entry gains the clause).
  3. Under `pairs = "all"` a singular level warns once per stored pair, not once. — **fix now** (warn once per level in `.tidy_edges()`, test extended).
  4. `summary()` now runs the partialled algebra and can warn on a singular level. — **reject**: `summary()` reports the variance table, and a singular level is worth its warning there too; the one-warning fix in 3 keeps it to one line.
  5. Block gate at 1e-8 vs two-decimal display prints `.00` rows; cosmetic alignment. — **reject**: the 1e-8 gate is AC4's text and two decimals is the package's `r` display convention; a `.00` row says the pair is below display precision, which is true.
  6. `.tidy_edges()` would error on an `NA` match index. — **reject**: latent; every construction path writes `tidy` and `matrices` from one matrix set (claim audit reached the same conclusion).
  7. Column whitelist in `.tidy_edges()` would silently drop a future upstream column. — **fix now** (splice `beta` after `r` instead of a whitelist; behaviour-preserving).
  8. `r2` roxygen says "same 0–1 scale" as the item-variance proportions but the denominator is the factor's score variance. — **fix now** (one clause).
  9. Rotation Φ (`factor_cor`) and score Φ_s (behind `beta`) are undistinguished in docs. — **follow-up → M90** (the distinction only bites once oblique lands; work-log line added to M90).
  10. `stopifnot()` in `.carry_factor_cor()` rather than `cli_abort()`. — **reject**: internal invariant on an engine-supplied matrix, not a user condition; the profile's cli rule is for user-facing conditions.
  11. `ord` not validated as a permutation. — **reject**: internal callers only pass `order()` output or `seq_len(k)`.
  12. PCA/EFA `factor_cor` unnamed while ESEM's is labelled. — **follow-up → M90** (labelling PCA/EFA now would change the O13 attribute comparison; M90 makes `factor_cor` load-bearing and can label all three; work-log line added to M90).
  13. Routing test cannot detect a wrong `ord`. — **follow-up → M90**, same remedy as 1.
  14. Near-singular Φ_s (rcond ~1e-15) may pass `solve()` without a warning. — **follow-up**: candidate row (search-first: no existing row; `meta$near_singular` covers R, not Φ_s).
  15. AC1's commit text (`a645e32`) differs from the fixture's (`ee71316`). — recorded under AC1 above; disposition at the gate.
  16. Under every supported configuration the new columns equal existing ones. — **reject**: the plan-gate intent, with its falsifier already recorded in the work log.
- Return floor: no finding demonstrates a criterion failing inside its procedure's domain; no load-bearing user-facing defect. No status change from the findings.
