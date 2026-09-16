# M87: Plain-English pass — vignettes B (girard, forbes, forbes2023, ordinal, interpret)

- **Status:** review
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

- [x] AC1: `Rscript tools/check-prose.R vignettes/ackwards-girard.Rmd.orig
      vignettes/ackwards-forbes.Rmd.orig vignettes/ackwards-forbes2023.Rmd.orig
      vignettes/ackwards-ordinal.Rmd.orig vignettes/ackwards-interpret.Rmd` exits 0 on the
      branch head.
- [x] AC2: `check_code_unchanged()` from `tools/check-prose.R`, run against the merge base
      with `master`, reports no differing fenced-chunk line (options and code) and no
      differing inline `` `r ` `` span in the five batch sources, and each such span still
      sits on one line or is byte-identical to its master form when it crosses a line break.
- [x] AC3: Every term in `tools/prose-terms.txt` is glossed in plain words at its first use
      in each batch vignette where it appears; the review reads the gloss sites from
      `grep -n` of each term over the five sources.
- [x] AC4: These claims survive with unchanged meaning at the site named: the forbes and
      forbes2023 vignettes' statement that the reproducing settings for Forbes (2023) are
      available and which ones they are (`redundancy_criterion = "direct"`, `k_max = 10` on
      the AMH matrix; IP9); the girard vignette's statement that edges are score correlations
      and the hierarchy is descriptive, not a fitted model (GP3); the ordinal vignette's
      statement that the Pearson default is kept and polychoric is opt-in with a warning
      (IP6); and the interpret vignette's IPIP-label lesson, including every item label it
      quotes verbatim from `bfi25`.
- [x] AC5: `Rscript vignettes/precompute.R` has been re-run, the regenerated `.Rmd` and
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
- [x] T6: Run `check_code_unchanged()` against the merge base; re-read the AC4 sites against
      `git show master:<file>`; fix drift.
- [x] T7: `Rscript vignettes/precompute.R`; revert timing-only churn in untouched vignettes;
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
- 2026-09-16: T6: `Rscript tools/check-prose.R --code-unchanged master` OK on the branch head (four `.Rmd.orig` sources; interpret compared by hand, see T5). The five AC4 sites re-read against `git show master:<file>`: girard Step 5 bullet (score correlations, descriptive, not a fitted model), forbes direct-criterion paragraph, forbes2023 `k_max = 10` fit plus direct default plus "reproduces ... exactly, all 54 components", ordinal "Automatic detection" (Pearson default, warning, polychoric opt-in), interpret `E4: Make friends easily` (matches `attr(bfi25$E4, "label")`). No drift found.
- 2026-09-16: T7: `Rscript vignettes/precompute.R` re-run (exit 0); intro, suggest-k, engines, visualization `.Rmd` and their assets reverted with `git checkout --` (M61). In the four regenerated batch outputs only stamps, one gt div id, cli timing lines, and the check-mark variation selector differ; PNGs re-rendered. `tools/dod-gate.R` prose step now calls `check_prose()` with no path argument (full default domain; M86 has merged). NEWS: third documentation entry added, and the M85 entry's "The vignettes are unchanged" sentence, no longer true for this release, now points at the two vignette entries. `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R`: GATE PASSED (prose clean, code-unchanged clean, check 0 err/0 warn/0 note, coverage 100%, style/lint clean, pkgdown index complete).
- 2026-09-16: claim audit: 41 claims read, 2 corrected — vignettes/ackwards-girard.Rmd.orig, ackwards-forbes.Rmd.orig, ackwards-forbes2023.Rmd.orig, ackwards-ordinal.Rmd.orig, ackwards-interpret.Rmd, NEWS.md, tools/dod-gate.R. Corrected: the girard table cell had called `suggest_k()`'s criteria variance-explained rules (MAP, VSS, CD are not), now "computed from the item correlation matrix"; the new redundancy gloss said the partner sits at a "deeper" level, but a chain that stops short of `k_max` keeps its top node, so all three files now say "another level". Both re-read once by the same reader: correct. The reader also called the ESEM gloss ("adds standard errors") weak but not wrong; kept for consistency with batch A. Precompute re-run for the three sources, untouched vignettes reverted, gate re-run: GATE PASSED (check 0 err/0 warn/0 note, coverage 100%).
- 2026-09-16: all tasks done, gate clean; status → review.

## Decisions

## Review

- 2026-09-16 review start: branch `m087-plain-english-vignettes-b` at bb64b14; `origin/master` unchanged since the branch was cut (no merge needed); no PR exists; tree clean.
- AC1: `Rscript tools/check-prose.R <five sources>` at bb64b14 printed "Prose OK" and exited 0. PASS.
- AC2: `Rscript tools/check-prose.R --code-unchanged master` printed "Code unchanged OK" and exited 0 (four `.Rmd.orig` sources). The live interpret file is outside that guard, so `.code_lines_rmd()` was run on `master:vignettes/ackwards-interpret.Rmd` and the branch file: 79 chunk-and-span lines each, texts identical. Every `` `r `` span in the five sources opens and closes on one line (0 unclosed-span lines per file). PASS.
- AC3: each of the 16 terms was grepped over the five sources and its first prose use read. Glossed at first use: girard (factor L31, parallel analysis L41, factor score L52, split-half L56, redundancy L73, PCA/component L92, polychoric/ordinal L95-97, EFA L229, loading L266); forbes (factor and factor score L21-23, redundant L26, loading L240, congruence L309, rotation L386, PCA/component L392, EFA L441, ESEM L442, split-half L603); forbes2023 (factor and redundancy L21, PCA/component L53, loading L164); ordinal (ordinal L24, loading and factor L27-29, parallel analysis L50, polychoric L56, rotation L72, ESEM L288, factor score L352, EFA L368); interpret (factor L27, loading L30). Terms not present in prose: varimax and FIML (no vignette), rotation and congruence in girard (bibliography titles only), polychoric in forbes (chunks only), and the rest per file. PASS.
- AC4: sites re-read on the branch and against `git show master:`. girard L271-274 keeps "edges are score correlations. It is descriptive, not a fitted hierarchical model"; forbes L161-166 keeps `redundancy_criterion = "direct"` as the default star criterion and "the rule Forbes's own code uses"; forbes2023 L64 keeps `k_max = 10, pairs = "all", n_obs = 3175`, L144 the direct default, L210-211 "reproduces Forbes's published chase exactly, all 54 components"; ordinal L78-93 keeps the Pearson default, the warning, and `cor = "polychoric"` as the opt-in that suppresses it; interpret quotes one item label, `E4: Make friends easily` (L88), identical to master and to `attr(bfi25$E4, "label")`. Meaning unchanged at all sites. PASS.
- AC5: the five regenerated `.Rmd` files and eight `vignettes/assets/` PNGs are committed on the branch (T7 and the claim-audit commits); `git diff origin/master --name-only -- vignettes/` lists only the batch files, so the untouched vignettes carry no diff. `DOD_CODE_UNCHANGED=1 TESTTHAT_CPUS=8 Rscript tools/dod-gate.R` at bb64b14 exited 0: vignette-freshness clean, prose clean, code-unchanged clean, check 0 errors / 0 warnings / 0 notes, coverage 100%, styler and lintr clean, pkgdown reference index complete. The gate's prose step calls `check_prose()` with no path argument (full default domain; M86 is merged). NEWS.md carries the "Plain-English vignettes, second batch" documentation entry. PASS.
- Consistency gate: `cairn_validate.py` all checks passed (16 pre-existing work-log format advisories on M84, not gate failures); no IP/GP changed, so `cairn_impact.py` was skipped; the r-package toolchain checks (document no-diff, check clean, pkgdown index, NEWS entry) are covered by the gate run above. No new top-level files. PASS.
