<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". -->
# M55: Address CRAN 0.1.0 feedback and resubmit as 0.1.1

- **Status:** in-progress
- **Priority:** high
- **Depends on:** —
- **Branch/PR:** m55-cran-resubmission

## Goal

Fix the two items from the 2026-07-12 CRAN reviewer feedback (unexplained
acronyms in DESCRIPTION; non-suppressible console output in
`label_template()`) and prepare the package for resubmission from master as
version 0.1.1.

## Scope

**In:** DESCRIPTION acronym expansion; `label_template()` refactor to a
classed visible return + `print()` method (house idiom; gate decision
2026-07-12); test/doc/vignette updates that follow; version bump to 0.1.1,
NEWS 0.1.1 section, `cran-comments.md` rewrite; release-grade check bar
(full CI matrix + fresh win-builder R-devel — CLAUDE.md release exception
applies: merge only on full green matrix). Absorbs the ROADMAP "0.1.0 CRAN
release tail" candidate (its patch-branch-from-tag guidance was superseded
at the plan gate: resubmit from master, including M54).

**Out:** the interactive `devtools::submit_cran()` call, `v0.1.1` tag, and
post-acceptance README update → owner-only post-merge tail (ROADMAP
candidate row). Wiring `forbes2023` into a vignette → existing candidate.
Any wider console-output sweep → not needed (verified 2026-07-12: all other
cli output sits in `print.*` methods, CRAN's allowed exception).

## Acceptance criteria

- [ ] AC1: DESCRIPTION `Description:` expands PCA, EFA, and ESEM on first
      use (per CRAN cookbook, description_issues.html#explaining-acronyms).
- [ ] AC2: `label_template()` returns its named character vector **visibly**
      with class `ackwards_labels`; a top-level call prints the header +
      editable `c(...)` literal via `print.ackwards_labels()`; assignment
      (`labs <- label_template(x)`) emits no console output
      (`expect_silent()` test). No `cat()`/`cli_text()` remains in the
      function body.
- [ ] AC3: The returned object still works downstream: subassignment
      (`labs["m5f1"] <- ...`) and use as `autoplot(x, node_labels = ...)`
      behave as before (tests pass).
- [ ] AC4: Package-wide grep shows no `cat()`/`print()`/cli console writes
      outside `print`/`summary`/plot methods.
- [ ] AC5: DESCRIPTION Version is 0.1.1; NEWS.md has a 0.1.1 section noting
      both CRAN-feedback changes (incl. the `label_template()` visible-return
      behavior change); `cran-comments.md` rewritten as a resubmission
      addressing both reviewer points.
- [ ] AC6: Precomputed vignettes regenerated via `vignettes/precompute.R`;
      `label_template()` sections of `ackwards-visualization` and
      `ackwards-interpret` render correctly.
- [ ] AC7: `Rscript tools/dod-gate.R` passes (check 0 err/0 warn,
      coverage, style, lint, pkgdown); full GitHub Actions CI matrix green;
      fresh win-builder R-devel run with 0 err/0 warn (misspelling NOTE
      for surnames/package name acceptable, already documented).

## Coverage

- AC1 → T1
- AC2 → T2, T3
- AC3 → T3
- AC4 → T2 (grep evidence at review)
- AC5 → T5
- AC6 → T4
- AC7 → T6

## Tasks

- [x] T1: Expand acronyms in DESCRIPTION Description text (PCA, EFA, ESEM;
      wording already exists in cran-comments.md "misspelled words" note).
- [x] T2: Refactor `R/label_template.R`: move header + `c(...)` literal
      rendering (lines 84–86) into `print.ackwards_labels()`; return
      `structure(out, class = c("ackwards_labels", "character"))` visibly;
      roxygen for the method via `@rdname label_template` (no new pkgdown
      topic); `devtools::document()`.
- [x] T3: Update `tests/testthat/test-label_template.R` (31 expectations):
      visible classed return, `expect_silent()` on assignment, print-method
      output test, downstream subassignment + `autoplot` node_labels round
      trip.
- [x] T4: Re-run `Rscript vignettes/precompute.R`; verify the live
      `ackwards-interpret.Rmd` chunk output; commit regenerated `*.Rmd` +
      `vignettes/assets/`.
- [x] T5: Version → 0.1.1; NEWS.md 0.1.1 section (fold the current dev
      section); rewrite `cran-comments.md` for resubmission (respond to both
      reviewer points explicitly; refresh platform table).
- [ ] T6: Run `Rscript tools/dod-gate.R`; open PR; wait for **full** CI
      matrix (release exception — no local-green merge); run
      `devtools::check_win_devel()` and confirm the emailed result with the
      owner before the merge gate.

## Work log

- 2026-07-12: created by /milestone-plan from CRAN reviewer feedback
  (2 items); gate decisions: resubmit from master as 0.1.1; classed
  return + print method; check bar = CI matrix + win-builder devel.
- 2026-07-12: T1 done — PCA/EFA/ESEM expanded in DESCRIPTION Description.
- 2026-07-12: T2+T3 done — label_template() now returns a visible classed
  vector (style kept as attribute); print.ackwards_labels() renders the
  header + c(...) literal; tests updated (42 expectations green, incl.
  expect_silent on assignment and autoplot round trip).
- 2026-07-12: T4 done — precomputed vignettes regenerated; scaffold no
  longer prints when nested in autoplot() (the CRAN complaint), top-level
  rendering unchanged; rest of diff is timing/PNG noise. Live interpret
  vignette needs no edit (top-level call auto-prints via method).
- 2026-07-12: T5 done — Version 0.1.1; NEWS folds dev section under a
  0.1.1 heading led by the two CRAN-feedback changes; cran-comments.md
  rewritten as a point-by-point resubmission response (platform table to
  be confirmed by T6 gate + CI + win-builder before merge).

## Decisions

## Review
