# M88: Prose-checker hardening (`tools/check-prose.R`)

- **Status:** in-progress
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

- [ ] AC1: With the working directory set to a temporary directory outside the repo,
      `read_prose_list("prose-banned.txt")` with no `dir` argument reads
      `tools/prose-banned.txt` from the checker's own directory. If the checker is loaded by
      `source(<path>)` or by `sys.source(<path>, envir)`, this holds.
- [ ] AC2: In the markdown fixture of `tests/testthat/test-check-prose.R`, a table cell
      holding 31 words and a heading holding 31 words are each reported as
      `sentence over 30 words` at their own line. A table row whose every cell holds under
      30 words and a heading under 30 words produce no report. Each cell (text between `|`
      separators) and each heading (its text after the `#` marks) counts as one sentence
      under `max_words`. `Rscript tools/check-prose.R README.Rmd DESCRIPTION NEWS.md R
      vignettes/ackwards-intro.Rmd.orig vignettes/ackwards-suggest-k.Rmd.orig
      vignettes/ackwards-engines.Rmd.orig vignettes/ackwards-visualization.Rmd.orig` still
      exits 0 on the branch head.
- [ ] AC3: In the same fixture, an em dash, an en dash, and a semicolon placed inside a bare
      `http(s)://` URL, inside an autolink `<http(s)://…>`, and inside a markdown link target
      `](…)` produce no `em dash`, `en dash`, or `semicolon` report. The same three marks in
      the link text `[…]` of that link are each reported.
- [ ] AC4: `tools/dod-gate.R` runs `check_code_unchanged("master")` when the environment
      variable `DOD_CODE_UNCHANGED` equals `1`. It prints each returned problem as a note and
      adds one gate failure naming the problem count. If the variable is unset, it prints a
      one-line skip note and runs no guard.
- [ ] AC5: In the git fixture of `tests/testthat/test-check-prose.R`, an inline `` `r` ``
      span moved, byte-identical, from before a fenced chunk to after it is reported by
      `check_code_unchanged()`. The reports on the fixture's four planted `.Rmd.orig` edits
      (chunk option, chunk body, inline span text, moved span) each name the working-tree
      line of the first differing item. The items are the fenced-chunk lines and inline `r`
      spans in document order.
- [ ] AC6: In the same git fixture, renaming `vignettes/v.Rmd.orig` on the work branch
      yields two `check_code_unchanged()` problems, one naming the old path and one the new.
      Each uses the form of the R-file loop (`<file> exists on only one side of <base>.`).
      The rule covers `vignettes/*.Rmd.orig`, and `README.Rmd` stays hard-coded.
- [ ] AC7: `Rscript tools/dod-gate.R` exits 0 on the branch head with `DOD_CODE_UNCHANGED`
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
- [ ] T3: `.char_reports()` (`:291`): before the dash and semicolon flags, blank the three
      URL forms of AC3 (bare, autolink, link target) with a same-length space run so line
      numbers and word counts hold. Banned phrases and sentence counts see the blanked text
      too. Write the AC3 tests, positive and silent.
- [ ] T4: `.code_lines_rmd()` (`:502`): return a data.frame(line, text) of chunk lines and
      inline spans in document order. Make the `.Rmd.orig` loop (`:571`) report a file
      present on one side only, as the R loop does at `:549`. Report the working-tree line. Extend the git fixture with the moved-span and rename plants.
- [ ] T5: `tools/dod-gate.R` (`:68-110`): add the opt-in guard step after the prose step,
      reading `DOD_CODE_UNCHANGED`, with the note and failure form of AC4. Update the header
      comment. Verify by hand once with the variable set on a planted `R/` edit in a scratch
      commit (reverted), and once unset. Record both outcomes in the work log.
- [ ] T6: Run the new tests against `git show master:tools/check-prose.R` (M82 lesson) and
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

## Decisions

## Review
