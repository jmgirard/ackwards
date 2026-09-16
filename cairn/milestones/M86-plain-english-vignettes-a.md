# M86: Plain-English pass — vignettes A (intro, suggest-k, engines, visualization)

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M85
- **Driving RR:** —
- **Principles touched:** GP3, IP6
- **Resolves:** —
- **Surface tier:** user-facing — vignettes are the pkgdown articles users read
- **Branch/PR:** `m086-plain-english-vignettes-a`

## Goal

Rewrite the intro, suggest-k, engines, and visualization vignette sources into plain
English that passes the M85 prose checker without changing any code, inline code, or
rendered output.

## Scope

**In:** prose rewrites of `vignettes/ackwards-intro.Rmd.orig`,
`vignettes/ackwards-suggest-k.Rmd.orig`, `vignettes/ackwards-engines.Rmd.orig`, and
`vignettes/ackwards-visualization.Rmd.orig`; re-run of `vignettes/precompute.R` with
timing-only churn in untouched vignettes reverted (LESSONS M61, M75); glosses for the M85
term list at first use in each vignette; the dod-gate prose step widened to these four
files.

**Out:** the other five vignettes → M87; any edit to fenced chunks, chunk options, or inline
`` `r ` `` spans (guarded by AC2); changes to `tools/check-prose.R` beyond adding an
abbreviation or banned phrase the corpus shows is missing (each such change is a work-log
line and keeps M85's tests green); term-list additions → appended to
`tools/prose-terms.txt` in this milestone with a work-log line.

## Acceptance criteria

- [ ] AC1: `Rscript tools/check-prose.R vignettes/ackwards-intro.Rmd.orig
      vignettes/ackwards-suggest-k.Rmd.orig vignettes/ackwards-engines.Rmd.orig
      vignettes/ackwards-visualization.Rmd.orig` exits 0 on the branch head.
- [ ] AC2: `check_code_unchanged()` from `tools/check-prose.R`, run against the merge base
      with `master`, reports no differing fenced-chunk line (options and code) and no
      differing inline `` `r ` `` span in the four batch sources, and each such span still
      sits on one line or is byte-identical to its master form when it crosses a line break.
- [ ] AC3: Every term in `tools/prose-terms.txt` is glossed in plain words at its first use
      in each batch vignette where it appears; the review reads the gloss sites from
      `grep -n` of each term over the four sources.
- [ ] AC4: These claims survive with unchanged meaning at the site named: the intro's framing
      that levels are complementary linked solutions whose edges are score correlations, not
      a fitted hierarchical model (`ackwards-intro`, GP3); the suggest-k vignette's statement
      that `suggest_k()` returns several criteria and a range, never one number
      (`ackwards-suggest-k`, D-013); the engines vignette's statement that a non-converging
      level warns and is skipped rather than aborting the run (`ackwards-engines`, IP7); and
      the ordinal-warning advice that the package never switches basis silently
      (`ackwards-engines`, IP6).
- [ ] AC5: `Rscript vignettes/precompute.R` has been re-run, the regenerated `.Rmd` and
      `vignettes/assets/` committed, untouched vignettes carry no diff against the merge
      base, and `Rscript tools/dod-gate.R` exits 0 with its prose step widened to the four
      batch sources; NEWS.md carries a documentation entry.

## Coverage

- AC1 → T1, T2, T3, T4
- AC2 → T5
- AC3 → T1, T2, T3, T4
- AC4 → T1, T3, T5
- AC5 → T6

## Tasks

- [ ] T1: Rewrite `ackwards-intro.Rmd.orig` (93 sentences, 34 em dashes on master); check.
- [ ] T2: Rewrite `ackwards-suggest-k.Rmd.orig` (148 sentences, 37 em dashes); check.
- [ ] T3: Rewrite `ackwards-engines.Rmd.orig` (115 sentences, 42 em dashes); check.
- [ ] T4: Rewrite `ackwards-visualization.Rmd.orig` (47 sentences, 31 em dashes); check.
- [ ] T5: Run `check_code_unchanged()` against the merge base; re-read the AC4 sites against
      `git show master:<file>`; fix drift.
- [ ] T6: `Rscript vignettes/precompute.R`; revert timing-only churn in the other five
      vignettes; widen the dod-gate prose step; NEWS entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan as the first vignette batch of the plain-English sweep; depends on M85 for the checker and term list.
- 2026-09-16: criteria audit ([O], fresh context, full mode) ran on the shared vignette-batch wording; its generated-`.Rmd` diff clause was dropped as unenforceable (stamp, PNG, and gt-id churn), and an inline-span guard was added.

## Decisions

## Review
