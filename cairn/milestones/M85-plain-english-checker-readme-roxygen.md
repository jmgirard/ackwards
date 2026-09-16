# M85: Plain-English pass — prose checker, README, DESCRIPTION, NEWS, roxygen

- **Status:** in-progress
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
      the `Description:` field of `DESCRIPTION`; NEWS.md between its first `# ` heading and
      its second; and every `#'` line of `R/*.R` outside `@examples` blocks and outside
      roxygen fenced code. Before matching it removes single- and double-backtick code
      spans, including spans that cross a line break. It reports file and line for each
      em dash, en dash, ` -- `, semicolon, phrase in `tools/prose-banned.txt`, and
      sentence over 30 words, where a sentence is the text between terminators (`.`, `?`,
      `!` followed by whitespace and an uppercase letter or a line end) after headings,
      bullet markers, table rows, and the abbreviations in `tools/prose-abbrev.txt`
      (`e.g.`, `i.e.`, `et al.`, `vs.`, `p.`, `pp.`, `cf.`, `Fig.`, `No.`) are removed;
      it exits non-zero on any report.
- [ ] AC2: `tests/testthat/test-check-prose.R` runs `check_prose()` on fixture text and is red
      on each report class in each of these forms and locations: dash forms `—`, `–`,
      ` -- `; locations YAML-adjacent prose, heading, bullet item, table cell, link text,
      roxygen `@param`, roxygen `@details`, NEWS entry, DESCRIPTION field; and is silent on
      a fixture carrying every class inside a single-backtick span, a double-backtick span,
      a span crossing a line break, a fenced chunk, roxygen fenced code, an `@examples`
      block, a YAML `---` fence, a markdown `---` rule, the precompute stamp comment, and
      each listed abbreviation.
- [ ] AC3: `Rscript tools/check-prose.R README.Rmd DESCRIPTION NEWS.md R/` exits 0 on the
      branch head.
- [ ] AC4: `check_code_unchanged(ref)` in `tools/check-prose.R`, run against the merge base
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
- [ ] AC6: These claims survive the rewrite with unchanged meaning, each still stated at the
      site named: a bass-ackwards result is linked solutions whose edges are score
      correlations, never a fitted hierarchical model (README.Rmd and `R/ackwards.R`
      description, GP3); every auto-resolved default announces itself and never switches
      basis silently (`R/ackwards.R` "Defaults and why", IP6); the settings that reproduce
      Forbes (2023) stay available and documented (`R/prune.R` `redundancy_criterion`, IP9);
      and every `@param` default rationale in `R/ackwards.R`, `R/prune.R`, and
      `R/suggest_k.R` names the same default value and reason as on master.
- [ ] AC7: `devtools::document()` and `devtools::build_readme()` produce no diff on the branch
      head; `Rscript tools/dod-gate.R` exits 0 with the prose check run fail-fast over the
      AC3 paths before `devtools::check()`; NEWS.md carries a documentation entry.

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
- [ ] T4: Rewrite roxygen in `R/ackwards.R`, `R/engine_*.R`, `R/compute_edges.R`,
      `R/data.R`; run `devtools::document()`; run the checker on those files.
- [ ] T5: Rewrite roxygen in `R/prune.R`, `R/suggest_k.R`, `R/tidy.R`, `R/augment.R`,
      `R/predict.R`, `R/boot_edges.R`, `R/comparability.R`; document; check.
- [ ] T6: Rewrite roxygen in the remaining `R/*.R`; document; check the whole AC3 domain.
- [ ] T7: Run `check_code_unchanged()` against the merge base and re-read the AC6 sites
      against `git show master:<file>`; fix any drift.
- [ ] T8: Wire the prose check into `tools/dod-gate.R` after the CI-filter step; NEWS entry;
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

## Decisions

## Review
