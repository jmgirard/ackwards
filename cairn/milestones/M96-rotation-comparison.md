# M96: Varimax and oblimin side by side on bfi25 in the engines vignette

- **Status:** in-progress
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

- [ ] AC1: The section of `vignettes/ackwards-engines.Rmd.orig` fits `ackwards()` twice to `bfi` with `engine = "efa"`, `cor = "polychoric"`, and `k_max = 5`. One fit uses the default rotation, and one uses `rotation = "oblimin"`. A search for `sim16` over the section's lines returns nothing.
- [ ] AC2: In the same section, the generated `vignettes/ackwards-engines.Rmd` shows five displays from the two fits:
  - (a) The oblimin fit's within-level factor correlations, from `tidy(x, what = "factor_cor")`.
  - (b) A one-to-one, per-level match of each oblimin factor to a varimax factor by Tucker congruence of loadings, with the congruence values.
  - (c) The primary-parent edges of both fits, with each oblimin ID shown beside its matched varimax ID.
  - (d) For each fit, the secondary edges with `above_cut` TRUE, with `r` and `beta`.
  - (e) One `autoplot()` figure per fit.
- [ ] AC3: `tests/testthat/test-vignette-rotation.R` checks four findings on the same two fits, and the section's prose states each one:
  - (1) The largest absolute within-level factor correlation of the oblimin fit.
  - (2) The one-to-one factor match at each level, with its smallest congruence, and the level whose IDs differ between the fits.
  - (3) After the match, the two primary-parent trees are equal.
  - (4) Which secondary edges have `above_cut` TRUE in each fit, with the sign and size of `beta` for each.
  Each check names the factors or edges it asserts. The file passes under `NOT_CRAN=true` with `testthat::test_file()`.
- [ ] AC4: `Rscript tools/check-vignette-freshness.R` exits 0. With the md5 stamp line removed, every hunk of `git diff master -- vignettes/ackwards-engines.Rmd` falls inside the section. The diff has no cli `[Nms]` timing change, no U+FE0E check mark, and no `div id=` change. `git diff --stat master -- vignettes/` lists only `ackwards-engines.Rmd`, `ackwards-engines.Rmd.orig`, and `vignettes/assets/ackwards-engines-*` figures from the section's chunks.
- [ ] AC5: `git diff --stat master -- R man manuscript` is empty.
- [ ] AC6: NEWS.md gains one entry under `# ackwards (development version)` that names the new rotation comparison. `Rscript tools/dod-gate.R` exits 0.

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

## Decisions

## Review
