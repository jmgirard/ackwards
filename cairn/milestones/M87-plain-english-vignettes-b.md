# M87: Plain-English pass — vignettes B (girard, forbes, forbes2023, ordinal, interpret)

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M85, M88
- **Driving RR:** —
- **Principles touched:** GP3, IP6, IP9
- **Resolves:** —
- **Surface tier:** user-facing — vignettes are the pkgdown articles users read
- **Branch/PR:** `m087-plain-english-vignettes-b`

## Goal

Rewrite the girard, forbes, forbes2023, ordinal, and interpret vignette sources into plain
English that passes the M85 prose checker without changing any code, inline code, or
rendered output, and widen the dod-gate prose step to the checker's full default domain.

## Scope

**In:** prose rewrites of `vignettes/ackwards-girard.Rmd.orig`,
`vignettes/ackwards-forbes.Rmd.orig`, `vignettes/ackwards-forbes2023.Rmd.orig`,
`vignettes/ackwards-ordinal.Rmd.orig`, and the live `vignettes/ackwards-interpret.Rmd`;
re-run of `vignettes/precompute.R` with timing-only churn in untouched vignettes reverted
(LESSONS M61, M75); glosses for the term list at first use in each vignette; the dod-gate
prose step switched to `check_prose()` with no path argument.

**Out:** the four batch-A vignettes → M86 (either batch may land first; both depend only on
M85); any edit to fenced chunks, chunk options, or inline `` `r ` `` spans (guarded by AC2);
changes to `tools/check-prose.R` beyond adding an abbreviation or banned phrase the corpus
shows is missing (checker hardening → M88) (each a work-log line, M85's tests kept green); term-list additions →
appended to `tools/prose-terms.txt` with a work-log line.

## Acceptance criteria

- [ ] AC1: `Rscript tools/check-prose.R vignettes/ackwards-girard.Rmd.orig
      vignettes/ackwards-forbes.Rmd.orig vignettes/ackwards-forbes2023.Rmd.orig
      vignettes/ackwards-ordinal.Rmd.orig vignettes/ackwards-interpret.Rmd` exits 0 on the
      branch head.
- [ ] AC2: `check_code_unchanged()` from `tools/check-prose.R`, run against the merge base
      with `master`, reports no differing fenced-chunk line (options and code) and no
      differing inline `` `r ` `` span in the five batch sources, and each such span still
      sits on one line or is byte-identical to its master form when it crosses a line break.
- [ ] AC3: Every term in `tools/prose-terms.txt` is glossed in plain words at its first use
      in each batch vignette where it appears; the review reads the gloss sites from
      `grep -n` of each term over the five sources.
- [ ] AC4: These claims survive with unchanged meaning at the site named: the forbes and
      forbes2023 vignettes' statement that the reproducing settings for Forbes (2023) are
      available and which ones they are (`redundancy_criterion = "direct"`, `k_max = 10` on
      the AMH matrix; IP9); the girard vignette's statement that edges are score correlations
      and the hierarchy is descriptive, not a fitted model (GP3); the ordinal vignette's
      statement that the Pearson default is kept and polychoric is opt-in with a warning
      (IP6); and the interpret vignette's IPIP-label lesson, including every item label it
      quotes verbatim from `bfi25`.
- [ ] AC5: `Rscript vignettes/precompute.R` has been re-run, the regenerated `.Rmd` and
      `vignettes/assets/` committed, untouched vignettes carry no diff against the merge
      base, and `Rscript tools/dod-gate.R` exits 0 with its prose step running
      `check_prose()` over the full default domain when M86 has merged, or over M85's paths
      plus the five batch sources when it has not; NEWS.md carries a documentation entry.

## Coverage

- AC1 → T1, T2, T3, T4, T5
- AC2 → T6
- AC3 → T1, T2, T3, T4, T5
- AC4 → T1, T2, T3, T4, T5, T6
- AC5 → T7

## Tasks

- [x] T1: Rewrite `ackwards-girard.Rmd.orig` (89 sentences, 32 em dashes on master); check.
- [x] T2: Rewrite `ackwards-forbes.Rmd.orig` (94 sentences, 43 em dashes); check.
- [x] T3: Rewrite `ackwards-forbes2023.Rmd.orig` (43 sentences, 15 em dashes); check.
- [x] T4: Rewrite `ackwards-ordinal.Rmd.orig` (65 sentences, 20 em dashes); check.
- [x] T5: Rewrite the live `ackwards-interpret.Rmd` (66 sentences, 27 em dashes); it is not
      precomputed, so knit it once locally to confirm it still renders; check.
- [ ] T6: Run `check_code_unchanged()` against the merge base; re-read the AC4 sites against
      `git show master:<file>`; fix drift.
- [ ] T7: `Rscript vignettes/precompute.R`; revert timing-only churn in untouched vignettes;
      set the dod-gate prose step per AC5; NEWS entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan as the second vignette batch of the plain-English sweep; depends on M85 for the checker and term list.
- 2026-09-16: criteria audit ([O], fresh context, full mode) ran on the shared vignette-batch wording; its generated-`.Rmd` diff clause was dropped as unenforceable (stamp, PNG, and gt-id churn), and an inline-span guard was added.
- 2026-09-16: /milestone-plan (M88 gate) added `Depends on: M88` so the vignette batch runs against the hardened checker and the opt-in `DOD_CODE_UNCHANGED=1` gate step.
- 2026-09-16: /milestone-implement started; branch `m087-plain-english-vignettes-b` cut from pushed master. Question gate skipped: nothing open (batch-A conventions reused: page ranges as `--`, Schmid-Leiman hyphenated, gloss in the first prose paragraph that uses a term).
- 2026-09-16: T1 girard rewritten; checker 0 reports, `--code-unchanged master` OK. Glossed: factor, parallel analysis, factor score, split-half, redundancy, PCA, polychoric, ordinal, EFA, loading. `rotation` and `congruence` appear only inside bibliography entry titles (treated as not appearing, like a link title). Step 4 heading's em dash became a comma (anchor changes).
- 2026-09-16: T2 forbes rewritten; checker 0 reports, code guard OK. Glossed: factor, factor score, redundant, artifactual, loading, congruence, rotation, PCA/component, EFA, ESEM, split-half. `polychoric` appears only inside chunks. The em dash before the conditional `r if (top_clean)` span became a space (the span's own text supplies the connective). The direct-criterion paragraph and its "the rule Forbes's own code uses" claim kept (AC4).
- 2026-09-16: T3 forbes2023 rewritten; checker 0 reports, code guard OK. Glossed: factor, redundant, PCA/component, loading. The `k_max = 10` / `n_obs = 3175` fit, the `redundancy_criterion = "direct"` default, and the "reproduces Forbes's published chase exactly, all 54 components" statement kept (AC4).
- 2026-09-16: T4 ordinal rewritten; checker 0 reports, code guard OK. Glossed: ordinal, loading, factor, parallel analysis, polychoric, rotation, ESEM, factor score, EFA. The "Automatic detection" section (Pearson default, warning, `cor = "polychoric"` opt-in) kept (AC4). Recommendation-table cells changed `—` to parentheses. Blockquote `>` markers count as words in the checker, so the score-computation note was split finer.
- 2026-09-16: T5 interpret rewritten; checker 0 reports; knitted once to the scratchpad with `rmarkdown::render()` (renders). `check_code_unchanged()` reads only `*.Rmd.orig`, so the live file's 79 chunk lines and spans were compared to `master:` with `.code_lines_rmd()` directly: identical. Glossed: factor, loading. The `E4: Make friends easily` label quote kept verbatim (AC4). A sentence that opens with a code span never splits in the checker (the span becomes a space, so no capital follows the period): three such sentences were re-opened with a word.

## Decisions

## Review
