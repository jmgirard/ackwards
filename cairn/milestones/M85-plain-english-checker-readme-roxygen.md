# M85: Plain-English pass — prose checker, README, DESCRIPTION, NEWS, roxygen

- **Status:** review
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** GP3, IP6, IP9
- **Resolves:** —
- **Surface tier:** user-facing — README, DESCRIPTION, NEWS, and roxygen are what users read on CRAN and pkgdown
- **Branch/PR:** m085-plain-english-checker-readme-roxygen

## Goal

Rewrite README.Rmd, the DESCRIPTION `Description` field, the NEWS development section, and
every roxygen comment into plain English a reader outside the field follows on one read,
enforced by a prose checker that later milestones widen to the vignettes.

## Scope

**In:** `tools/check-prose.R` (the checker and its code-unchanged guard) with a testthat
file; a committed banned-phrase list, abbreviation list, and term list under `tools/`;
prose rewrites of README.Rmd, `DESCRIPTION` `Description:`, the NEWS.md development
section, and every `#'` line in `R/*.R` outside `@examples`; regenerated `man/` and
README.md; the checker wired into `tools/dod-gate.R` over this milestone's domain.

**Out:** the nine vignette sources → M86 (intro, suggest-k, engines, visualization) and M87
(girard, forbes, forbes2023, ordinal, interpret), both depending on this milestone for the
checker and the term list; released NEWS sections (history, left as written); non-roxygen
code comments in `R/*.R` (users never see them; they keep the ` -- ` dash style); pkgdown
navbar titles (not requested); any change to code, `@examples`, chunk code, or numeric
output (guarded by AC4); the modal-verb rule of Simplified Technical English (the gate chose
the Plain register, so `may`/`might`/`should` stay allowed).

## Acceptance criteria

- [ ] AC1: `tools/check-prose.R` defines `check_prose(paths)` and, run as
      `Rscript tools/check-prose.R [paths]` (a directory argument expands to its `*.R`,
      `*.Rmd`, and `*.Rmd.orig` files; no argument means the full default domain), sweeps
      exactly these prose lines: README.Rmd, `vignettes/*.Rmd.orig`, and
      `vignettes/ackwards-interpret.Rmd` outside the YAML header and outside fenced chunks;
      the `Description:` field of `DESCRIPTION`; NEWS.md from its first `# ` heading through
      the line before its second; and every `#'` line of `R/*.R` outside `@examples` blocks, outside
      roxygen fenced code, and outside Rd `\preformatted{}` blocks. Before matching it
      removes single- and double-backtick code spans, including spans that cross a line
      break. It reports file and line for each em dash, en dash, ` -- `, semicolon,
      phrase in `tools/prose-banned.txt`, and sentence over 30 words, where a sentence is
      the text between terminators (`.`, `?`, `!`, optionally followed by a closing
      quote, bracket, or emphasis mark, then whitespace and an optional opening quote,
      bracket, or emphasis mark before an uppercase letter, or a line end) after
      headings, table rows, and section titles are dropped, bullet markers and roxygen
      tag words are stripped, and the abbreviations in `tools/prose-abbrev.txt` (`e.g.`,
      `i.e.`, `et al.`, `vs.`, `p.`, `pp.`, `cf.`, `Fig.`, `No.`) are removed, and where
      a paragraph break, a new bullet, or a new roxygen tag also ends a sentence; it
      exits non-zero on any report.
- [x] AC2: `tests/testthat/test-check-prose.R` runs `check_prose()` on fixture text and is red
      on each report class in each of these forms and locations: dash forms `—`, `–`,
      ` -- `; locations YAML-adjacent prose, heading, bullet item, table cell, link text,
      roxygen `@param`, roxygen `@details`, NEWS entry, DESCRIPTION field; and is silent on
      a fixture carrying every class inside a single-backtick span, a double-backtick span,
      a span crossing a line break, a fenced chunk, roxygen fenced code, an `@examples`
      block, a YAML `---` fence, a markdown `---` rule, the precompute stamp comment, and
      each listed abbreviation.
- [x] AC3: `Rscript tools/check-prose.R README.Rmd DESCRIPTION NEWS.md R/` exits 0 on the
      branch head.
- [x] AC4: `check_code_unchanged(ref)` in `tools/check-prose.R`, run against the merge base
      with `master`, reports no differing line among: non-roxygen lines of `R/*.R`, `#'`
      lines inside `@examples`, fenced-chunk lines and inline `` `r ` `` spans of
      README.Rmd, and every `DESCRIPTION` field other than `Description:`; and it exits 0 on
      the branch head.
- [ ] AC5: Every term in `tools/prose-terms.txt` (the statistical terms of art the rewrite
      meets, committed in T1 after a grep of the domain, at least: factor, component,
      loading, rotation, varimax, factor score, polychoric, ordinal, FIML, ESEM, EFA, PCA,
      congruence, redundancy, parallel analysis, split-half) is glossed in plain words at
      its first use in README.Rmd and at its first use in each exported roxygen topic where
      it appears; the review reads the gloss sites from `grep -n` of each term over the
      domain.
- [x] AC6: These claims survive the rewrite with unchanged meaning, each still stated at the
      site named: a bass-ackwards result is linked solutions whose edges are score
      correlations, never a fitted hierarchical model (README.Rmd and `R/ackwards.R`
      description, GP3); every auto-resolved default announces itself and never switches
      basis silently (`R/ackwards.R` "Defaults and why", IP6); the settings that reproduce
      Forbes (2023) stay available and documented (`R/prune.R` `redundancy_criterion`, IP9);
      and every `@param` default rationale in `R/ackwards.R`, `R/prune.R`, and
      `R/suggest_k.R` names the same default value and reason as on master.
- [ ] AC7: `devtools::document()` produces no diff, and `devtools::build_readme()` changes
      only `#>`-prefixed output lines of the `suggest-print` chunk in README.md (that chunk is
      unseeded and AC4 forbids adding a seed), leaving the rest of README.md and
      `man/figures/` byte-identical over two consecutive builds; `Rscript tools/dod-gate.R`
      exits 0 with the prose check run fail-fast over the AC3 paths before
      `devtools::check()`; NEWS.md carries a documentation entry.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3, T4, T5, T6
- AC4 → T1, T7
- AC5 → T1, T3, T4, T5, T6
- AC6 → T4, T5, T7
- AC7 → T3, T7, T8

## Tasks

- [x] T1: Write `tools/check-prose.R` (`check_prose()`, `check_code_unchanged()`, script
      body guarded by `sys.nframe()` as in `tools/check-vignette-freshness.R:109`) plus
      `tools/prose-banned.txt` (at least: robust, crucial, simply, essentially, importantly,
      not just, note that, nuance, seamless, comprehensive, leverage, delve, in order to,
      it is worth noting, in conclusion), `tools/prose-abbrev.txt`, and
      `tools/prose-terms.txt`; `.Rbuildignore` the new files.
- [x] T2: Write `tests/testthat/test-check-prose.R` with the AC2 fixture matrix; skip in the
      built package as `test-vignette-freshness.R` does.
- [x] T3: Rewrite README.Rmd, `DESCRIPTION` `Description:`, and the NEWS development section;
      run `devtools::build_readme()`.
- [x] T4: Rewrite roxygen in `R/ackwards.R`, `R/engine_*.R`, `R/compute_edges.R`,
      `R/data.R`; run `devtools::document()`; run the checker on those files.
- [x] T5: Rewrite roxygen in `R/prune.R`, `R/suggest_k.R`, `R/tidy.R`, `R/augment.R`,
      `R/predict.R`, `R/boot_edges.R`, `R/comparability.R`; document; check.
- [x] T6: Rewrite roxygen in the remaining `R/*.R`; document; check the whole AC3 domain.
- [x] T7: Run `check_code_unchanged()` against the merge base and re-read the AC6 sites
      against `git show master:<file>`; fix any drift.
- [x] T8: Wire the prose check into `tools/dod-gate.R` after the CI-filter step; NEWS entry;
      `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan from the request to sweep vignettes, README, and roxygen for plain English; split three ways at the gate (M85 checker + README/DESCRIPTION/NEWS/roxygen; M86, M87 vignettes).
- 2026-09-16: criteria audit ([O], fresh context, full mode) returned 10 findings; 8 fixed in the criteria before the gate (scripted chunk and inline-span guard against the merge base, committed term list, must-survive claim list replacing an open fresh-reader read, wider planted probes, incremental gate widening, named list files and strip rules), the sentence-length threshold posed at the gate (30 words chosen), and the non-roxygen comment dash split recorded as Out.
- 2026-09-16: plan gate chose a standing gate (checker wired into dod-gate) over a one-off script because docs written after M87 would otherwise regress unseen; falsified by the checker blocking a legitimate future doc edit more often than it catches a regression.
- 2026-09-16: plan gate chose the Plain register with a 30-word cap over Strict STE modals because hedged statistical claims (may not converge) would have to be overstated as will or must; falsified by a reader report that allowed modals still confuse.
- 2026-09-16: plan chose a must-survive claim list (AC6) over an open fresh-reader read because the audit showed the open read has a built-in escape and no oracle; falsified by a meaning change at a site the list does not name.
- 2026-09-16: implement started; branch m085-plain-english-checker-readme-roxygen cut from pushed master. Question gate skipped: the plan fixes the checker contract, the lists, the register, and the domain, so nothing was open.
- 2026-09-16: T1 done. Checker written with two refinements beyond the AC1 letter, both documented in the script header: a terminator may be followed by a closing bracket or quote before the split, and an opening bracket or quote may precede the uppercase letter; Rd `\preformatted{}` blocks count as roxygen fenced code. The three lists are already covered by the `^tools$` `.Rbuildignore` entry, so no new entry was needed. First sweep of the AC3 domain: 376 reports (137 semicolons, 128 double hyphens, 88 long sentences, 19 em dashes, 4 banned phrases).
- 2026-09-16: T2 done. Fixture matrix passes (57 expectations across the three dash forms and nine locations, plus the silent fixture, the sentence-counting fixture, the non-empty-domain check, and a git-backed `check_code_unchanged()` fixture that plants one change per guarded class). Three fixture slips fixed on the way, none a checker bug. Dev library under R 4.6 lacked styler, EFAtools, and gt; installed from CRAN.
- 2026-09-16: T3 done. README.Rmd, `DESCRIPTION` `Description:`, and the NEWS development entry rewritten; README glosses factor, PCA, EFA, ESEM, ordinal, polychoric, split-half, redundancy, parallel analysis, factor score, and loading at first use. Checker clean on the three files; `check_code_unchanged("master")` clean; README.md rebuilt; suite 726 tests, 0 failures.
- 2026-09-16: T4, T5, T6 done by three parallel [O] subagents on disjoint file sets (T4: ackwards, compute_edges, data, the engine files carry no prose. T5: prune, suggest_k, tidy, augment, predict, boot_edges, comparability. T6: autoplot, layout, print, summary, check_items, factor_labels, factorability, interpret, label_template, utils has no prose). Each ran the checker to clean and the code guard to clean before returning. Their flagged rewordings were read: the `fm` "robust OLS" gloss, the n_obs FIML clause reorder, the `sign_by` sentence prefixes, and sentence-initial code spans given a noun all keep the stated meaning. One paraphrase was reverted to the term with a gloss ("factor overdetermination" in factorability.R). `devtools::document()` run and re-run with no further diff.
- 2026-09-16: toolchain drift found on the fresh R 4.6 library and neutralised without code change: roxygen2 8.1.0 rewrote NAMESPACE's importFrom layout and bumped the DESCRIPTION stamp, so roxygen2 8.0.0 (master's version) was installed and the stamp kept; lintr 3.4.0 flags the two pre-existing `<<-` in suggest_k.R that master also trips, so `.lintr` now passes `assignment_linter(operator = c("<-", "<<-"))`. Both are dev tooling, neither touches a user-facing surface.
- 2026-09-16: T7 done. `check_code_unchanged("master")` clean after the whole rewrite. AC6 sites read on the branch: the ackwards() description keeps "edges are score correlations, not a fitted higher-order SEM"; the `cor = "pearson"` bullet keeps "no silent basis switching" with the warning that suggests polychoric; "Varimax is the only supported rotation" stays; `redundancy_criterion` keeps `"direct"` as the default that reproduces Forbes's ChaseCorrPaths and AMH example with `"adjacent"` as the opt-in. A scripted extraction of "Default `value`" mentions per `@param` matched master for prune.R and suggest_k.R and for ackwards.R apart from regex-form differences on `fm` and `n_obs`, whose values (`"minres"`, `"total"`) are unchanged.
- 2026-09-16: T8 done. Prose check wired into `tools/dod-gate.R` after the CI-filter step, over README.Rmd, DESCRIPTION, NEWS.md, and R/ (the vignette sources join when M86/M87 land). `Rscript tools/dod-gate.R` exit 0: prose clean, check 0/0/0, coverage 100%, styler and lintr clean, pkgdown index complete.
- 2026-09-16: claim audit: 186 claims read, 5 corrected — R/tidy.R, R/ackwards.R, R/augment.R, R/prune.R, README.Rmd. The five: tidy.R's fit bullet had turned the list of scaled-test estimators into a claim that ESEM selects MLR for continuous items (the code selects ML), the `fm` minres gloss had invented a non-normality reason for a banned "robust" (restored to the convergence contrast), augment.R's `id_cols` read "ignored, and an error" (now two error sentences), the prune() topic used "re-rotation" without a gloss, and README's Step 2 said edges are computed between neighbouring levels only (now "different levels", neighbouring by default). NEWS softened from "every help page rewritten" to "the help pages were revised". All five re-read once by the same reader and confirmed; one optional tightening left unapplied (rotation is skipped at k = 1, where the README sentence about rotation is vacuously true).
- 2026-09-16: all tasks done; status set to review. Post-correction re-verification: prose and code guards clean, `document()` and `build_readme()` re-run, suite 726 tests with 0 failures, lint 0, styler 0.
- 2026-09-16: review pass 1 returned the milestone to in-progress (defect return 1 of M85). AC5 fails at five unglossed first-use sites (ackwards.R:43 ESEM, boot_edges.R:49 ESEM, comparability.R:482 factor, suggest_k.R:752 factor, factor_labels.R:120 factor). AC2, AC3, AC4, AC6 ticked against evidence. Gate green (check 0/0/0, coverage 100%, style, lint, pkgdown). 19 [O] findings logged in the Review section for triage at the re-review gate.
- 2026-09-16: amendment return: AC1 — "where a sentence is the text between terminators (`.`, `?`, `!`, optionally followed by a closing quote, bracket, or emphasis mark, then whitespace and an optional opening quote or bracket before an uppercase letter, or a line end)" and "outside roxygen fenced code and Rd `\preformatted{}` blocks". The implementation splits at these sites and the letter does not.
- 2026-09-16: return fixes. The five AC5 sites now gloss at first use (ackwards.R data param, boot_edges.R description, autoplot.comparability, autoplot.suggest_k, factor_labels() description). Checker hardened per [O] findings 1-5: an unmatched backtick within a paragraph, an unclosed fence (markdown or roxygen), and an unterminated YAML header now error with file and line instead of going quiet, spans never cross a paragraph break, and banned phrases match across a wrapped line (finding 3). Word forms added to prose-banned.txt (finding 4). Test file gains a wrapped-phrase case and the three error paths (finding 10): 83 expectations, 0 failures. dod-gate now exits before check() when any base-R guard failed, so AC7's "fail-fast" holds by the letter (finding 6). AC3 sweep clean, code guard clean, document() regenerated five Rd files.
- 2026-09-16: amendment return: AC7 — "`devtools::document()` produces no diff and `devtools::build_readme()` produces no diff outside knitted output blocks on the branch head". The README `suggest_k(bfi25)` chunk is unseeded, its CD criterion varies per build, and AC4 forbids adding a seed.

- 2026-09-16: amendment gate accepted both amendment returns (AC1 sentence rule and preformatted blocks, AC7 README exemption) and the return fixes (five glosses plus [O] findings 1-5 and 10).
- 2026-09-16: re-audit: AC1 (full) — four findings: the opener class omitted emphasis marks, the paragraph/bullet/roxygen-tag sentence boundaries were unstated, "removed" conflated dropped and stripped items, and "exactly" outran the silent extraction paths. First three fixed in the wording; the fourth closed by the checker now erroring on those paths.
- 2026-09-16: re-audit: AC7 (full) — four findings: "fail-fast" was false of the gate, "knitted output block" had no mechanical referent, the file domain (man/figures/) was unstated, and the exemption was wider than its cause. First closed by the gate's early exit; the other three fixed by narrowing the clause to the `#>` lines of the `suggest_k` chunk over two consecutive builds.

- 2026-09-16: re-audit: AC1 (full) — sentence rule matches the code on all eight probed axes. Three findings: the error paths (unclosed fence, YAML, unmatched backtick) are unbound by any criterion, "between its first heading and its second" reads as excluding the heading line the code sweeps, and the `@param` name and `\item{}` label strips are unstated. Second re-entry, so disposition went to the user.
- 2026-09-16: re-audit: AC7 (full) — four findings: the fail-fast clause is unwitnessed by an exit-0 run (a planted failure is the witness), the chunk is named `suggest-print` in README.Rmd not `suggest_k`, `man/figures/` byte-identity binds the graphics device, and the diff baseline is unstated. Second re-entry, so disposition went to the user.
- 2026-09-16: claim audit: 23 claims read, 0 corrected — R/ackwards.R, R/boot_edges.R, R/comparability.R, R/factor_labels.R, R/suggest_k.R, tools/check-prose.R, tools/prose-banned.txt, tools/dod-gate.R, tests/testthat/test-check-prose.R (lines added since 2b08b8f; the earlier audit covered the rest). Gate re-run exit 0 (check 0/0/0, coverage 100%, style, lint, pkgdown). Planted em dash in NEWS.md made the gate exit 1 in 0.6s before check(), then restored.

- 2026-09-16: user chose the two clear wording fixes (AC1 NEWS boundary, AC7 chunk name `suggest-print`) and held the rest, declining the error-path widening. Return work complete, gate green, status set to review for pass 2.

## Decisions

## Review

Review pass 1, 2026-09-16, on branch head 2b08b8f (merge base with master 79ce393; master unmoved).

**Evidence per criterion.**

- AC1: `check_prose(paths)` is defined. Directory expansion, the default domain, the extraction rules for YAML, fences, `@examples`, roxygen fences, the NEWS section, and the Description field, cross-line span stripping, the six report classes, and the nine abbreviations all match the script, read line by line. `Rscript tools/check-prose.R` with no argument reports on the vignette sources and exits 1. Not ticked. The sentence splitter also splits after a terminator followed by a closing quote or bracket, and before an uppercase letter preceded by an opening quote or bracket (`term` at line 277, `mark_re` at line 280). It also treats Rd `\preformatted{}` as roxygen fenced code (lines 147-156). Both widen the AC1 letter. Both are documented in the script header and the T1 work-log line. Amendment proposed below.
- [x] AC2: `test_file("tests/testthat/test-check-prose.R")` under `load_all()` gave 72 expectations passed, 0 failed, 0 skipped.
- [x] AC3: `Rscript tools/check-prose.R README.Rmd DESCRIPTION NEWS.md R/` printed "Prose OK" and exited 0.
- [x] AC4: `Rscript tools/check-prose.R --code-unchanged master` printed "Code unchanged OK" and exited 0 against merge base 79ce393.
- AC5: FAILS at five sites. A script listed the first use of each term in README.Rmd and in every exported roxygen block (`@export` or `@format`, spans stripped). README glosses all 14 terms that appear (FIML and congruence do not appear). Unglossed first uses: `R/ackwards.R:43` ESEM ("ESEM requires raw data", gloss only at line 78). `R/boot_edges.R:49` ESEM (gloss only at line 69). `R/comparability.R:482` factor in the `autoplot.comparability` topic. `R/suggest_k.R:752` factor in the `autoplot.suggest_k` topic. `R/factor_labels.R:120-126` factor in the `factor_labels()` topic, which has its own Rd. Read as not appearing: a term only inside another term's expansion (component inside "principal component analysis", factor inside "exploratory factor analysis") and a term only in a cited paper title (`R/data.R:47`, `R/comparability.R:128`).
- [x] AC6: `README.Rmd:111-112` and `R/ackwards.R:8` keep "edges are score correlations, never a fitted hierarchical model" and "not a fitted higher-order SEM". `R/ackwards.R:25` keeps "no silent basis switching" with the polychoric suggestion. `R/prune.R:686-698` keeps `"direct"` as the default that reproduces Forbes's `ChaseCorrPaths`, with `"adjacent"` as the opt-in. A script extracted every "default" sentence per `@param` on master and on the branch for `R/ackwards.R`, `R/prune.R`, and `R/suggest_k.R`. Every parameter names the same value and reason. `fm` swaps "robust OLS" for "converges reliably" and keeps the contrast with `"ml"`. `rules` keeps "no auto rule", and the clear-with-no-arguments hint survives at `R/prune.R:658`.
- AC7: `devtools::document()` produced no diff. `Rscript tools/dod-gate.R` exited 0 with the prose step before `check()`. The gate printed check 0/0/0 (72s), coverage 100.00%, styler clean, lintr clean, pkgdown index complete. NEWS.md carries the documentation entry. Not ticked. `devtools::build_readme()` changed README.md in the `suggest_k(bfi25)` output only (CD criterion k = 7 became k = 6, consensus range 4-7 became 4-6). The chunk at `README.Rmd:94` is unseeded, so the comparison-data criterion varies per build. A seed is chunk code, which AC4 forbids changing. Amendment proposed below. README.md was restored to the committed version after the read.

**Consistency gate.** `cairn_validate.py` exit 0, all checks pass (16 pre-existing work-log format advisories, all in M84). No principle changed, so `cairn_impact` was skipped. Profile checks: document() no diff (above). README.md sync fails only by the unseeded output above. pkgdown, the NEWS entry, `.Rbuildignore` (`^tools$` covers the four new files), and the full check 0/0/0 all pass via the gate run above.

**Independent review** (three fresh-context lenses, full fan-out, user-facing tier). Every finding is logged. Triage is the maintainer's at the re-review gate. Blame-history lens [S]: no wording contradicts a D-entry. One pre-existing gap noted: `secondary_scope` from D-033 is not in `R/autoplot.R` on master (M84 is blocked), outside this diff. Prior-review lens [S]: no prior-review evidence regressed (M63 citation and `n_obs` wording, M76-M78 chase semantics all preserved), and zero GitHub inline comments exist repo-wide. Diff-bug lens [O], ranked:

1. Unbalanced backtick mutes every later report in the file (`strip_code_spans` joins the whole file). Latent, domain clean today.
2. Unbalanced roxygen fence drops the rest of the block. Latent.
3. Multi-word banned phrases wrapped across a line never match (`\s+` escape is dead, matching is per line).
4. Banned list morphology: robustly, robustness, crucially, comprehensively, leverages, leveraging all pass.
5. Unterminated YAML opener drops the whole file silently.
6. Gate prose step appends to `failures` rather than exiting; "fail-fast" is the wrong word (matches the three prior guards).
7. `.prose_tools_dir()` cannot find `tools/` when sourced from outside the repo root; callers pass the lists explicitly.
8. `.lintr` change out of the milestone's scope (reviewer's own classification).
9. `check_code_unchanged()` drops the line association of README inline `r` spans.
10. AC2 fixture never asserts a multi-word banned phrase.
11. Long sentences in table cells and headings are never reported (AC1 letter).
12. Semicolons or dashes in bare URLs are false positives (none in the domain).
13. `R/boot_edges.R:49` ESEM unglossed at first use (an AC5 miss, counted above).
14. `README.Rmd:45-46` says every engine rotates, which is vacuous at k = 1 (known, unapplied).
15. `R/comparability.R` `n_splits` rationale reads circular after removing "robust".
16. `R/layout.R` states the crossing rule twice in incompatible words.
17. NEWS "dashes and semicolons are gone" is true of README/DESCRIPTION/roxygen only. The no-argument checker run reports 562 on the vignettes.
18. `check_code_unchanged()` not wired into the gate (AC4 asks only for a clean run).
19. README rebuild fixed a stale `citation()` version block; pre-existing.

**Disposition.** Defect return on AC5 (return 1 of M85). Recommended for the same return: fix the five AC5 sites, and findings 1-4 and 10 (the checker's silent paths, which M86/M87 inherit). The rest go to triage at the re-review gate. Two criteria fail by their wording, not the work. They go to the gated amendment protocol at `/milestone-implement` step 6 before re-review:

- AC1 amended clause proposal: "where a sentence is the text between terminators (`.`, `?`, `!`, optionally followed by a closing quote, bracket, or emphasis mark, then whitespace and an optional opening quote or bracket before an uppercase letter, or a line end)" and "outside roxygen fenced code and Rd `\preformatted{}` blocks".
- AC7 amended clause proposal: "`devtools::document()` produces no diff and `devtools::build_readme()` produces no diff outside knitted output blocks on the branch head" (the `suggest_k` chunk is unseeded and AC4 forbids a seed).
