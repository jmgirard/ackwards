# M88: Prose-checker hardening (`tools/check-prose.R`)

- **Status:** review
- **Priority:** high
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** —
- **Resolves:** —
- **Surface tier:** internal — `tools/check-prose.R` is dev tooling over in-repo documentation sources, and no package user runs it
- **Branch/PR:** `m088-prose-checker-hardening`

## Goal

Fix the six deferred M85 and M86 review findings on the prose checker and its code-unchanged
guard before M87 runs. The checker then resolves its own lists, checks table cells and
headings, ignores URLs, sees inline spans in document order, reports a one-sided vignette
source, and runs from the DoD gate on request.

## Scope

**In:** `tools/check-prose.R` (`.prose_tools_dir()`, `.sentence_reports()`, `.char_reports()`,
`.code_lines_rmd()`, the `.Rmd.orig` loop of `check_code_unchanged()`). An opt-in
`check_code_unchanged()` step in `tools/dod-gate.R` behind `DOD_CODE_UNCHANGED=1`. The
planted-defect tests in `tests/testthat/test-check-prose.R`. Lineage: the `[high]`
prose-checker candidate row (M85 review items 1-4, M86 review items 5-7).

**Out:** any prose rewrite → M87 (vignettes B, which depends on this milestone). Additions to
`tools/prose-*.txt` → M87's work log. The `Rscript` load route of the list defaults (already
resolved through `--file=`, unchanged here). A NEWS entry (internal tooling, no user-visible
change). Banned-phrase or sentence-cap policy changes (none requested).

## Acceptance criteria

- [x] AC1: With the working directory set to a temporary directory outside the repo,
      `read_prose_list("prose-banned.txt")` with no `dir` argument reads
      `tools/prose-banned.txt` from the checker's own directory. If the checker is loaded by
      `source(<path>)` or by `sys.source(<path>, envir)`, this holds.
- [x] AC2: In the markdown fixture of `tests/testthat/test-check-prose.R`, a table cell
      holding 31 words and a heading holding 31 words are each reported as
      `sentence over 30 words` at their own line. A table row whose every cell holds under
      30 words and a heading under 30 words produce no report. Each cell (text between `|`
      separators) and each heading (its text after the `#` marks) counts as one sentence
      under `max_words`. `Rscript tools/check-prose.R README.Rmd DESCRIPTION NEWS.md R
      vignettes/ackwards-intro.Rmd.orig vignettes/ackwards-suggest-k.Rmd.orig
      vignettes/ackwards-engines.Rmd.orig vignettes/ackwards-visualization.Rmd.orig` still
      exits 0 on the branch head.
- [x] AC3: In the same fixture, an em dash, an en dash, and a semicolon placed inside a bare
      `http(s)://` URL, inside an autolink `<http(s)://…>`, and inside a markdown link target
      `](…)` produce no `em dash`, `en dash`, or `semicolon` report. The same three marks in
      the link text `[…]` of that link are each reported.
- [x] AC4: `tools/dod-gate.R` runs `check_code_unchanged("master")` when the environment
      variable `DOD_CODE_UNCHANGED` equals `1`. It prints each returned problem as a note and
      adds one gate failure naming the problem count. If the variable is unset, it prints a
      one-line skip note and runs no guard.
- [x] AC5: In the git fixture of `tests/testthat/test-check-prose.R`, an inline `` `r` ``
      span moved, byte-identical, from before a fenced chunk to after it is reported by
      `check_code_unchanged()`. The reports on the fixture's four planted `.Rmd.orig` edits
      (chunk option, chunk body, inline span text, moved span) each name the working-tree
      line of the first differing item. The items are the fenced-chunk lines and inline `r`
      spans in document order.
- [x] AC6: In the same git fixture, renaming `vignettes/v.Rmd.orig` on the work branch
      yields two `check_code_unchanged()` problems, one naming the old path and one the new.
      Each uses the form of the R-file loop (`<file> exists on only one side of <base>.`).
      The rule covers `vignettes/*.Rmd.orig`, and `README.Rmd` stays hard-coded.
- [x] AC7: `Rscript tools/dod-gate.R` exits 0 on the branch head with `DOD_CODE_UNCHANGED`
      unset.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T5
- AC5 → T4
- AC6 → T4
- AC7 → T6

## Tasks

- [x] T1: `.prose_tools_dir()` (`tools/check-prose.R:34`): resolve the script path at source
      time. Walk the call frames for `source()`'s `ofile` and `sys.source()`'s `file`, then
      fall back to `--file=`. On no hit, stop with a message rather than return `tools`
      (M82 lesson: fail closed). Test the two AC1 routes from a temp working directory.
- [x] T2: `.sentence_reports()` (`:340`): instead of dropping `|` rows and headings, split a
      table row into cells and strip heading marks. Feed each as its own flushed sentence. Keep the separator-only row (`|---|`) dropped. Write the AC2 tests and run
      the AC2 command line.
- [x] T3: `.char_reports()` (`:291`): before the dash and semicolon flags, blank the three
      URL forms of AC3 (bare, autolink, link target) with a same-length non-space filler so
      line numbers and word counts hold. Banned phrases and sentence counts see the blanked text
      too. Write the AC3 tests, positive and silent.
- [x] T4: `.code_lines_rmd()` (`:502`): return a data.frame(line, text) of chunk lines and
      inline spans in document order. Make the `.Rmd.orig` loop (`:571`) report a file
      present on one side only, as the R loop does at `:549`. Report the working-tree line. Extend the git fixture with the moved-span and rename plants.
- [x] T5: `tools/dod-gate.R` (`:68-110`): add the opt-in guard step after the prose step,
      reading `DOD_CODE_UNCHANGED`, with the note and failure form of AC4. Update the header
      comment. Verify by hand once with the variable set on a planted `R/` edit in a scratch
      commit (reverted), and once unset. Record both outcomes in the work log.
- [x] T6: Run the new tests against `git show master:tools/check-prose.R` (M82 lesson) and
      record the failure count in the work log. Then run `styler`, `lintr`, and
      `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan from the `[high]` prose-checker candidate row (M85 review items 1-4, M86 review items 5-7). M87 amended to depend on this milestone.
- 2026-09-16: criteria audit ([O], fresh context, reduced mode) returned seven findings, all fixed before writing. AC1 narrowed to the two sourced routes. AC2 defines per-cell and per-heading counting. AC3 drops ` -- ` (a space ends a URL). AC4's review procedure moved to T5. AC5's universal narrowed to the fixture's four plants with the reported line defined. AC6 states two problems over `vignettes/*.Rmd.orig`. AC7 lost a redundant test clause.
- 2026-09-16: plan gate chose an opt-in `DOD_CODE_UNCHANGED=1` gate step over always-on wiring because every code milestone fails an always-on gate, and over not wiring it because M87 needs a repeatable run. Falsified by a prose-only milestone that ships a code edit with the variable unset, an edit the guard catches when set.
- 2026-09-16: plan gate chose one 30-word cap for cells and headings over a separate shorter cap because a second number adds policy without a corpus case demanding it. Falsified by cells that pass 30 words yet read as too long in review.
- 2026-09-16: plan chose blanking URLs before all reports over skipping only the dash and semicolon flags because a URL's words otherwise inflate a sentence count. Falsified by a URL whose blanking hides a genuine prose defect on the same line.
- 2026-09-16: /milestone-implement started; branch `m088-prose-checker-hardening` cut from pushed master. Question gate skipped: the plan leaves no implementation choice open.
- 2026-09-16: T1 done. The script path is resolved once at load time by walking the frames for `source()`'s `ofile` and `sys.source()`'s `file`, then `--file=`; no hit makes `.prose_tools_dir()` stop. New test loads the checker both ways from a temp working directory; prose test file 9/9 clean, Rscript route from `/tmp` exits 0.
- 2026-09-16: T2 done. `.sentence_units()` expands a heading into one standalone unit and a table row into one unit per cell (separator rows dropped); each flushes as its own sentence. New test: 31-word heading and cell reported at lines 1 and 7, 29-word ones silent. Test file 97 pass, 0 fail. The AC2 sweep command exits 0.
- 2026-09-16: T3 done, with a minor amendment. `.blank_urls()` runs in `check_prose()` right after span stripping, so every report sees URL-free text. A space filler split `[text](target),` into three tokens and surfaced a 31-word artefact in the engines vignette, so the filler is a same-length run of `x` (a URL stays one token, as before); trailing sentence punctuation on a bare URL is kept. Test file 102 pass, 0 fail; the AC2 sweep command exits 0.
- 2026-09-16: T4 done. `.code_lines_rmd()` returns a data.frame(line, text) in document order; a report names the item and the working-tree line (`item N differs from the merge base (line L: text)`). The `.Rmd.orig` loop reports a one-sided file in the R-loop form. Five existing item-number expectations moved to document order. Fixture gains the moved-span plant (item 1, line 5) and the rename plant (two problems). Test file 107 pass, 0 fail; `--code-unchanged master` on the branch exits 0.
- 2026-09-16: T5 done. Opt-in step after the prose step in `tools/dod-gate.R`; header comment updated. Hand run with `DOD_CODE_UNCHANGED=1` on a planted comment line appended to `R/ackwards.R` (working-tree edit, reverted with `git checkout --`): one `code-unchanged:` note naming line 779, gate failed before check() with `code-unchanged guard: 1 problem(s)`, exit 1. The unset run is the T6 full gate.
- 2026-09-16: T6 done. Test file run against `git show master:tools/check-prose.R`: 16 failures, 90 pass (every new assertion fails on the old checker); on the branch 107 pass, 0 fail. Two test lines wrapped for `line_length_linter`. `Rscript tools/dod-gate.R` with the variable unset: skip note printed, check 0/0/0, coverage 100.00%, styler and lintr clean, pkgdown index complete, exit 0.
- 2026-09-16: claim audit: not owed — internal tier.
- 2026-09-16: all tasks checked; status set to review.

## Decisions

## Review

- 2026-09-16 AC1: verified. A scratch script with the working directory set to a temp directory outside the repo (`tools/prose-banned.txt` absent there) loaded `tools/check-prose.R` by `source()` and by `sys.source()`. Both resolved the tools directory to `/Users/jmgirard/github/ackwards/tools`, and `read_prose_list("prose-banned.txt")` returned the trimmed, comment-free list identical to the file. `test-check-prose.R` on the branch head: 107 pass, 0 fail.
- 2026-09-16 AC2: verified. A scratch fixture held a 31-word heading (line 1), a 29-word heading, a 31-word cell (line 7), and a row of two 29-word cells. It produced exactly two reports, `sentence over 30 words (31)` at lines 1 and 7. The eight-path sweep command from the criterion exited 0 with `Prose OK`.
- 2026-09-16 AC3: verified. A scratch fixture placed an em dash, an en dash, and a semicolon inside a bare URL, an autolink, and a link target. Those lines produced no report. The same marks in the link text of line 4 produced exactly three reports at line 4: `em dash`, `en dash`, and `semicolon`.
- 2026-09-16 AC5: verified by the branch-head run of `test-check-prose.R` (107 pass, 0 fail). Its git fixture asserts the moved-span plant is reported as item 1 at line 5. The four planted edits report items 1, 3, 2, and 3 with working-tree lines 4, 6, 5, and 6. Items are chunk lines and inline spans in document order per `.code_lines_rmd()`.
- 2026-09-16 AC6: verified by the same run. The fixture's `git mv` of `vignettes/v.Rmd.orig` to `w.Rmd.orig` yields exactly two problems in the form `<file> exists on only one side of`, one per path. The reverse rename yields none. `README.Rmd` stays hard-coded in `check_code_unchanged()`.
- 2026-09-16 AC4 (set case): verified. A comment line was appended to `R/ackwards.R` in the working tree. `DOD_CODE_UNCHANGED=1 Rscript tools/dod-gate.R` then printed one `code-unchanged:` note naming line 779 and failed before check() with `code-unchanged guard: 1 problem(s)`, exit 1. The plant was reverted. The unset case is recorded with AC7.
- 2026-09-16 AC4 (unset case) and AC7: verified. `Rscript tools/dod-gate.R` on the branch head with the variable unset printed `code-unchanged: skipped (set DOD_CODE_UNCHANGED=1 to run the guard)` and ran no guard. It ended `GATE PASSED` (check 0/0/0, coverage 100.00%, style and lint clean, pkgdown index complete), exit 0. The tree was clean afterwards, so `document()` and `styler` produced no diff.
- 2026-09-16 consistency gate: `cairn_validate.py` exit 0, all checks pass (16 pre-existing work-log format advisories, all in M84). No principle changed, so `cairn_impact` was skipped. NAMESPACE, man/, README, DESCRIPTION, and `.Rbuildignore` are untouched by the branch. No NEWS entry, per the scope (internal tooling).
- 2026-09-16 review lens [S] blame-history: no findings. The removed heading/table drop and the removed `"tools"` fallback both came from M85 with no recorded rationale. The only sourcing caller (`tools/dod-gate.R`, `sys.source`) is on a supported route. No other consumer parses the old report text. No D-entry concerns the checker or the gate.
- 2026-09-16 review lens [S] prior-review record: no findings. Each deferred M85/M86 item is addressed as claimed. The GitHub inline-comment probe returned an empty list, so the PR-thread walk was skipped. Noted: the archives list seven deferred items with two overlaps, against the Goal's "six".
- 2026-09-16 review lens [O] diff-bug: 14 ranked findings. Its verdict: the criteria are met on their own fixtures. Finding 1 is a regression to fix before merge, and findings 2 to 5 are holes the fixtures miss. Triage follows, one line per finding.
- 2026-09-16 finding 1 (fix now): if README.Rmd was absent on both sides, `check_code_unchanged()` still reported `README.Rmd exists on only one side`, because the name is listed unconditionally. Fixed: a both-sides-absent file is skipped. Test: both branches drop README.Rmd and the check returns nothing.
- 2026-09-16 finding 2 (fix now): a link target with nested parentheses or a quoted title was blanked only up to the first `(`. Its remaining marks leaked as prose. Fixed: the target pattern admits one paren level and a title. Test: Wikipedia-style and titled targets holding all three marks are silent.
- 2026-09-16 finding 3 (fix now): the bare-URL class ran past `}`, so `\href{url}{text}` blanked the first word of the text. Fixed: `{` and `}` end a bare URL. Test: a semicolon in the `\href` text is reported.
- 2026-09-16 finding 4 (fix now): cells split on the escaped pipe `\|`, so a long cell can evade the cap. Fixed: cells split on unescaped pipes only. Test: a 41-word cell with one `\|` is reported at 41.
- 2026-09-16 finding 5 (fix now): `base::sys.source()` was not recognised as the loader, so the checker stopped. Fixed: the loader is matched by call name in bare or `base::` form. Test: the `base::sys.source` route resolves the list.
- 2026-09-16 finding 6 (fix now): the script path was kept relative and normalised at call time. A later working-directory change therefore moved the tools directory. Fixed: the path is made absolute at load time. Test: a relative load from the repo root, then a working-directory change, still resolves the list.
- 2026-09-16 finding 7 (fix now): any frame holding a variable `ofile` was accepted as the loader. Fixed with finding 5: the `ofile` branch requires a `source()` call. Covered by the AC1 `source()` route.
- 2026-09-16 finding 9 (fix now): a heading or cell with an internal period was split at it. A 31-word heading `Foo bar. Baz …` therefore went unreported, against AC2's "counts as one sentence". Fixed: a standalone unit is never split. Test: that heading is reported at 31. Classified fix-now rather than a return. AC2's named fixture and sweep both passed, so the failure lay outside the named procedure's domain. The repair is a code fix, not a criterion change.
- 2026-09-16 finding 8 (rejected): the `--file=` fallback returns the running script's path. Under `Rscript tools/check-prose.R` no loader frame exists, so this is the CLI's primary route, not a guess.
- 2026-09-16 finding 10 (rejected): the one-sided message cannot say added versus removed. AC6 fixes the message form to the R-loop's. A vignette source added on a prose-only branch is a code change the guard exists to flag.
- 2026-09-16 finding 11 (rejected): the gate hard-codes `master` and counts a checker error as one problem. `master` is this repo's default branch (CLAUDE.md). A tool error surfacing as a gate failure is the fail-closed behavior wanted. AC4 has no automated test by plan. The Coverage maps it to the hand-verified T5, and `tools/` is outside covr.
- 2026-09-16 finding 12 (rejected): `doi:`, `www.`, `ftp://`, and `mailto:` are not blanked. AC3 names `http(s)` only, and the sweep over DESCRIPTION passes today. A corpus case reopens it.
- 2026-09-16 finding 13 (rejected): `item removed after line 0` when the working tree has no code items. Cosmetic, on a branch AC5 does not constrain.
- 2026-09-16 finding 14 (rejected): `c()` growth in two loops. Irrelevant at the corpus size of a dev tool.
- 2026-09-16 fix-now verification: `test-check-prose.R` 115 pass, 0 fail on the fixed branch. The same test file against the pre-fix checker: 5 fail, 1 error, so every new assertion discriminates. The AC2 sweep command exits 0. styler and lintr clean on both files.
- 2026-09-16 AC4 and AC7 re-verified on the fix-now head (3e4d2c7): `Rscript tools/dod-gate.R` with the variable unset printed the skip note and ended `GATE PASSED`, exit 0. Check 0/0/0, coverage 100.00%, style and lint clean, pkgdown index complete, tree clean. AC1, AC2, AC3, AC5, and AC6 are re-verified by the 115-pass test run above and the sweep exit 0.