# M96: Varimax and oblimin side by side on bfi25 in the engines vignette

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP5
- **Resolves:** —
- **Surface tier:** user-facing — a shipped vignette section
- **Branch/PR:** m096-rotation-comparison

## Goal

If users fit an oblique rotation in place of varimax, show them on one real dataset what changes and what stays the same.

## Scope

**In:** Rewrite the oblique example in the "Orthogonal or oblique rotation" section of `vignettes/ackwards-engines.Rmd.orig`. In this file, "the section" means the lines from `## Orthogonal or oblique rotation` up to, but not including, `## Missing data`. The new example fits EFA with varimax and with oblimin to the vignette's `bfi` object (`na.omit(bfi25)`, setup chunk). It uses the polychoric basis, as the rest of the vignette does. It replaces the `sim16` chunk. Also in scope: a guard test file, the regenerated `.Rmd` and figures, and a NEWS entry.

**Out:**
- A worked comparison in the manuscript becomes a candidate row.
- An exported helper that matches factors across two fits becomes a candidate row. The vignette uses `psych::factor.congruence()` (psych is in Imports).
- PCA, ESEM, promax, and geomin comparisons are not planned. The section's prose already names those options.
- No change goes to `R/`, `man/`, the intro vignette, or the manuscript.

Facts measured at plan time on 2026-10-01 at master `bd55297`, with EFA, polychoric, `na.omit(bfi25)`, and `k_max = 5`:
- Oblimin within-level factor correlations reach .34.
- Each oblimin factor matches one varimax factor, with Tucker congruence of .94 or more. At level 5, the two fits swap the numbers of `m5f2` and `m5f3`.
- After the match, the primary-parent tree is identical.
- Varimax has 0 secondary edges with `above_cut` TRUE (`cut_show = 0.3`). Oblimin has 3: `m3f1→m4f3` (`beta` .016), `m4f1→m5f2` (`beta` .0015), and `m4f3→m5f4` (`beta` .16).

The implement phase re-derives these facts. The prose is written from that run, not from this list.

## Acceptance criteria

- [x] AC1: The section of `vignettes/ackwards-engines.Rmd.orig` fits `ackwards()` twice to `bfi` with `engine = "efa"`, `cor = "polychoric"`, and `k_max = 5`. One fit uses the default rotation, and one uses `rotation = "oblimin"`. A search for `sim16` over the section's lines returns nothing.
- [x] AC2: In the same section, the generated `vignettes/ackwards-engines.Rmd` shows five displays from the two fits:
  - (a) The oblimin fit's within-level factor correlations, from `tidy(x, what = "factor_cor")`.
  - (b) A one-to-one, per-level match of each oblimin factor to a varimax factor by Tucker congruence of loadings, with the congruence values.
  - (c) The primary-parent edges of both fits, with each oblimin ID shown beside its matched varimax ID.
  - (d) For each fit, the secondary edges with `above_cut` TRUE, with `r` and `beta`.
  - (e) One `autoplot()` figure per fit.
- [x] AC3: `tests/testthat/test-vignette-rotation.R` checks four findings on the same two fits, and the section's prose states each one:
  - (1) The largest absolute within-level factor correlation of the oblimin fit.
  - (2) The one-to-one factor match at each level, with its smallest congruence, and the level whose IDs differ between the fits.
  - (3) After the match, the two primary-parent trees are equal.
  - (4) Which secondary edges have `above_cut` TRUE in each fit, with the sign and size of `beta` for each.
  Each check names the factors or edges it asserts. The file passes under `NOT_CRAN=true` with `testthat::test_file()`.
- [x] AC4: `Rscript tools/check-vignette-freshness.R` exits 0. With the md5 stamp line removed, every hunk of `git diff master -- vignettes/ackwards-engines.Rmd` falls inside the section. The diff has no cli `[Nms]` timing change, no U+FE0E check mark, and no `div id=` change. `git diff --stat master -- vignettes/` lists only `ackwards-engines.Rmd`, `ackwards-engines.Rmd.orig`, and `vignettes/assets/ackwards-engines-*` figures from the section's chunks.
- [x] AC5: `git diff --stat master -- R man manuscript` is empty.
- [x] AC6: NEWS.md gains one entry under `# ackwards (development version)` that names the new rotation comparison. `Rscript tools/dod-gate.R` exits 0.

## Coverage

- AC1 → T2
- AC2 → T2, T3
- AC3 → T1, T2
- AC4 → T3
- AC5 → T2, T3
- AC6 → T4

## Tasks

- [x] T1: Fit both models and re-derive the four AC3 findings. Write `tests/testthat/test-vignette-rotation.R` first, before the prose. Use `cached()` fits, and call `skip_if_not_installed("GPArotation")` before any fit. Assert which factor and which edge, never counts alone (follow `test-vignette-m24.R`).
- [x] T2: Rewrite the section in the `.Rmd.orig`. Replace the `sim16` chunk with chunks for AC2 (a) to (e). Guard the chunks on `has_gpa`, and guard the figure chunks on ggplot2 too. Write the prose from the T1 output, then reread the kept "How to decide" paragraph in its new place (LESSONS M92). Run `Rscript tools/check-prose.R` on the file (LESSONS M85).
- [x] T3: Run `Rscript vignettes/precompute.R`. Revert churn in the other vignettes and assets, and revert timing, check-mark, and `div id` noise in this vignette (LESSONS M61, M75, M87). Run the freshness check, then the stamp-less diff check (LESSONS M94).
- [x] T4: Add the NEWS entry. Run `Rscript tools/dod-gate.R`.

## Work log

- 2026-10-01: created by /milestone-plan.
- 2026-10-01: criteria audit (full mode, fresh Opus reader) returned six findings, all fixed before the gate. The correlation basis was unpinned (now polychoric). The section bounds were undefined. AC2 needed rewording. AC3's open "every sentence" domain became four named findings. A skip promise about the test setup moved to T1. AC4 gained a noise check. A redundant prose-check criterion was dropped.
- 2026-10-01: plan gate chose to extend the engines vignette section over a new dedicated vignette. All rotation guidance then stays in one place, with no new pkgdown entry. Falsified by a section too long for one reader to follow, or by a user who cannot find the comparison.
- 2026-10-01: plan chose EFA on the polychoric basis over PCA or the Pearson basis. The vignette uses polychoric everywhere else, and EFA shows larger within-level correlations (.34, against .27 for PCA on Pearson). Falsified by an EFA oblimin fit on bfi25 that does not converge on another platform.

- 2026-10-01: implement started on branch m096-rotation-comparison. Question gate skipped, nothing open.
- 2026-10-01: T1 done. Re-derived all plan facts; oblimin with no seed and seeds 1, 2, 42, 123, 2026 gave the same matches and above-cut edges (differences under 1e-5), so the fit uses `seed = 1` for exact reproducibility. Wrote `test-vignette-rotation.R` (4 tests, 19 expectations). A plant with varimax in place of oblimin fails all four. Suite: 3607 pass, 0 fail.
- 2026-10-01: T2 done. The section now has six chunks for AC2 (a) to (e) and moves "How to decide" to the end. The prose was read against the rendered output. Added a note that `beta` can pass 1 (1.02 for m3f1 to m4f1) and a note that the two diagrams order some nodes differently. Both figures have captions and alt text. `check-prose.R` is clean.
- 2026-10-01: T3 done. Precompute ran clean. Churn in the girard, intro, suggest-k, and visualization vignettes and three assets was reverted. The engines `.Rmd` was rebuilt by a scratch script from master's file with the new stamp line and the freshly knitted section spliced in, because gt tables outside the section changed every element ID. The freshness check passes. The diff has 3 hunks (the stamp and two inside the section) and 0 timing, 0 U+FE0E, and 0 `div id=` lines.
- 2026-10-01: T4 done. Added the NEWS entry. `tools/dod-gate.R` exited 0: check had 0 errors, 0 warnings, and 0 notes. Coverage was 100%. Style, lint, prose, and the pkgdown index were clean.
- 2026-10-01: claim audit: 50 claims read, 4 corrected — vignettes/ackwards-engines.Rmd.orig, tests/testthat/test-vignette-rotation.R, NEWS.md
- 2026-10-01: The audit fixes were the varimax alt text (the top factor has no parent), two test comments, and the NEWS phrase "No code or result changed". It also gave an optional "In oblimin IDs" note. Its two test gaps were closed: test 1 now pins the top-three order, and test 4 pins both small `beta` values. The same reader re-read all eight items once and cleared them.
- 2026-10-01: After the fixes, precompute and the splice ran again. The diff kept 3 hunks and 0 noise lines. The plant still fails all 4 tests. The suite passed 3606 with 0 failures. `tools/dod-gate.R` exited 0 again. Status set to review.
- 2026-10-01: review gate: maintainer chose to fix eight reviewer findings and reject the rest. Fixes landed in `2fb1298`, and the gate passed again.
- step-7 approval: m096-rotation-comparison approved for merge

## Decisions

## Review

Evidence gathered 2026-10-01 on branch head `3a9ccdd`. Master did not move after the branch was cut (`47ac2cf`).

- AC1: The section has 150 lines. It fits `ackwards(bfi, k_max = 5, engine = "efa", cor = "polychoric")` as `x_var`, and the same call with `rotation = "oblimin", seed = 1` as `x_obl`. `grep -c sim16` over the section returns 0. Master's file had 4 hits, and the 3 hits left in the file are outside the section.
- AC2: A read of the rendered section in `vignettes/ackwards-engines.Rmd` shows all five displays. (a) `tidy(x_obl, what = "factor_cor")` prints 20 rows, sorted by absolute size. (b) The `matches` table prints 15 rows with a `congruence` column, and the `stopifnot()` one-to-one guard runs. (c) The merged `tree` table prints 14 primary edges, with `from_obl` and `to_obl` beside the varimax IDs. `anyNA(tree)` prints FALSE. (d) `above_cut_secondary()` prints 0 rows for varimax and 3 rows with `r` and `beta` for oblimin. (e) Two figures, `rotation-plot-varimax-1.png` and `rotation-plot-oblimin-1.png`, are linked with captions and alt text.
- AC3: `testthat::test_file("tests/testthat/test-vignette-rotation.R")` under `NOT_CRAN=true` ran 4 tests and 18 expectations, with 0 failures and 0 skips. Test 1 pins the top three pairs `m5f3-m5f4`, `m3f1-m3f3`, `m4f1-m4f3` at .34, .33, .30. Test 2 pins the one-to-one match, the swap of `m5f2` and `m5f3`, and the smallest congruence .94 at level 5. Test 3 pins equal sorted edge keys after the match. Test 4 pins no varimax edge and three oblimin edges, `m3f1->m4f3`, `m4f1->m5f2`, and `m4f3->m5f4`. It pins each `beta` as positive, at .016, .0015, and .16. The prose states each of the four. A fresh plant put `rotation = "varimax"` in place of `"oblimin"`, and all 4 tests went red. Tests 1 to 3 failed, and test 4 errored on the missing edge.
- AC4: `Rscript tools/check-vignette-freshness.R` exited 0. `git diff -U0 master -- vignettes/ackwards-engines.Rmd` has 5 hunks. One is the stamp line at line 13. The other 4 touch master lines 1501 to 1539 and new lines 1501 to 1708. The section spans master lines 1475 to 1540 and new lines 1475 to 1708, so all 4 hunks are inside it. The diff has 0 `ms]` timing lines, 0 U+FE0E characters, and 0 `div id=` lines. `git diff --name-status master -- vignettes/` lists the `.Rmd`, the `.Rmd.orig`, and two added figures, `ackwards-engines-rotation-plot-varimax-1.png` and `ackwards-engines-rotation-plot-oblimin-1.png`, from the section's two plot chunks.
- AC5: `git diff --stat master -- R man manuscript` printed 0 lines, and all three directories exist.
- AC6: `git diff master -- NEWS.md` adds one bullet under `# ackwards (development version)`, titled "Varimax and oblimin side by side in the engines article". It names the BFI-25 comparison and has no milestone number. `Rscript tools/dod-gate.R` exited 0 and printed GATE PASSED. Check had 0 errors, 0 warnings, and 0 notes. Coverage was 100%. Vignette freshness, prose, style, lint, and the pkgdown index were clean. The gate left the tree unchanged.
- Consistency gate: `cairn_validate.py` exited 0, with 16 work-log format advisories, all in the M84 file. `devtools::document()` produced no diff. The branch touches no DESIGN principle, no README, and no top-level file, so `cairn_impact.py`, `build_readme()`, and the `.Rbuildignore` check do not apply. `pkgdown::check_pkgdown()` and `devtools::check()` passed inside the AC6 gate run.

Independent review, 2026-10-01. Three fresh reviewers ran: diff-bug (Opus), blame-history (Sonnet), and prior review (Sonnet). No finding shows a criterion failing. Each finding is listed with its rank, and dispositions follow at the merge gate.

- R-P (prior review): no prior-review evidence. The GitHub probe found no human review comments, and no archived lesson is broken.
- D1 (diff, rank 1): the prose calls the .943 match "the same factor", but the package's own congruence cutoff is 0.95 (`redundancy_phi`).
- D2 (diff, rank 2): the prose names one `beta` above 1 (m3f1 → m4f1, 1.02). m2f2 → m3f2 is 1.002 and prints as 1.00.
- D3 (diff, rank 3): "Oblimin rotates from random starts" holds only for psych versions with `n.rotations`. DESCRIPTION has no psych floor. The `?ackwards` `seed` doc says the same thing.
- D4 (diff, rank 4): the match shows `abs()` congruence. On other data, a reflected factor then matches with no visible sign. All matched values here are positive.
- D5 (diff, rank 5): when a reader reuses the chunk and two factors pick the same match, the one-to-one `stopifnot()` stops with no next step.
- D6 (diff, rank 6): "set by `cut_show = 0.3`" reads as if an argument was passed. 0.3 is the default.
- D7 (diff, rank 7): the prose explains the two near-zero `beta` values but not why m4f3 → m5f4 keeps .16.
- D8 (diff, rank 8): test 4 pins `signif(beta, 2) == 0.0015` on a true value 2.1e-5 from a rounding boundary. Seeds moved it by at most 4e-6.
- D9 (diff, rank 9): the T1 work-log line says 19 expectations. The file has 18, as the AC3 line records.
- B1 (history, rank 1): the removed sentence "The fit's message says what the edges mean" pointed readers at the oblique fit's advisory. The advisory still prints, but no sentence points to it.
- B2 (history, rank 2): the section shows equal trees and then says "keep varimax" for hierarchy questions, with no example where oblique moves a parent.
- B3 (history, rank 3): "How to decide" now sits under the `### Varimax and oblimin on the same data` heading, after the figures.
- B4 (history, rank 4): the vignette no longer shows an oblique PCA fit or the `sim16` ground truth. The plan chose this.
- B5 (history, rank 5): the `beta` wording is new but agrees with `R/tidy.R`.

Triage at the gate, 2026-10-01. The maintainer chose to fix eight findings and reject the rest.

- Fixed in `2fb1298`. D1: the prose calls the pair a close match, under the .95 `prune()` default. D2: it names m2f2 → m3f2 at 1.002. D4: the table keeps the sign, and test 2 asserts every match is positive. D6: it says the default `cut_show` of 0.3. D7: it explains the .34 to .16 drop through m4f1. D8: `expect_lt` with a 1e-4 tolerance. B1: a paragraph after the fits names the advisory's points. B3: a `### How to decide` heading.
- Rejected: D3, because the wording matches the `?ackwards` `seed` doc and a psych version floor is a separate question. D5, because the stop is meant to fail loudly. B2, because the equal trees are what the data show. B4, because the plan chose EFA on bfi25. B5, because it reports no conflict.
- Noted: D9. The work log is history, and the AC3 line records the count.
- Re-verification after the fixes: the test file ran 4 tests and 21 expectations with 0 failures. Test 3 also pins the two primary `beta` values over 1. The prose check and the freshness check passed. Against `HEAD~1`, the rendered diff changed only prose lines, one code line, and the stamp, and no printed output changed. Against master, the noise counts stayed at 0, and all content hunks fall inside the section (new lines 1475 to 1724). `tools/dod-gate.R` exited 0 again, with check at 0 errors, 0 warnings, and 0 notes, and coverage at 100%. The new positive-congruence assertion was not plant-tested, because no plant of a flipped factor was built.
