# M86: Plain-English pass — vignettes A (intro, suggest-k, engines, visualization)

- **Status:** review
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
`` `r ` `` spans (guarded by AC2); changes to `tools/check-prose.R` beyond (a) adding an
abbreviation or banned phrase the corpus shows is missing and (b) widening
`check_code_unchanged()` to read `vignettes/*.Rmd.orig` with its existing README
chunk-and-span extraction (each such change is a work-log line, adds a test where it adds
behavior, and keeps M85's tests green); term-list additions → appended to
`tools/prose-terms.txt` in this milestone with a work-log line.

## Acceptance criteria

- [x] AC1: `Rscript tools/check-prose.R vignettes/ackwards-intro.Rmd.orig
      vignettes/ackwards-suggest-k.Rmd.orig vignettes/ackwards-engines.Rmd.orig
      vignettes/ackwards-visualization.Rmd.orig` exits 0 on the branch head.
- [x] AC2: `check_code_unchanged()` from `tools/check-prose.R`, run against the merge base
      with `master`, reports no differing fenced-chunk line (options and code) and no
      differing inline `` `r ` `` span in the four batch sources, and each such span still
      sits on one line or is byte-identical to its master form when it crosses a line break.
- [x] AC3: Every term in `tools/prose-terms.txt` is glossed in plain words at its first use
      in each batch vignette where it appears; the review reads the gloss sites from
      `grep -n` of each term over the four sources.
- [x] AC4: These claims survive with unchanged meaning at the site named: the intro's framing
      that levels are complementary linked solutions whose edges are score correlations, not
      a fitted hierarchical model (`ackwards-intro`, GP3); the suggest-k vignette's statement
      that `suggest_k()` returns several criteria and a range, never one number
      (`ackwards-suggest-k`, D-013); the engines vignette's statement that the between-level
      edges are descriptive correlations between factor scores with no sampling distribution
      of their own, and that per-level fit qualifies one level and never the hierarchy as a
      whole (`ackwards-engines`, GP3); and the engines vignette's missing-data statement that
      `missing = "fiml"` announces the PCA/EFA route with a message, and errors rather than
      switching basis or estimator silently, for WLSMV/ULSMV on the ESEM route and for a
      non-Pearson PCA/EFA basis (`ackwards-engines`, IP6).
- [x] AC5: `Rscript vignettes/precompute.R` has been re-run, the regenerated `.Rmd` and
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

- [x] T1: Rewrite `ackwards-intro.Rmd.orig` (93 sentences, 34 em dashes on master); check.
- [x] T2: Rewrite `ackwards-suggest-k.Rmd.orig` (148 sentences, 37 em dashes); check.
- [x] T3: Rewrite `ackwards-engines.Rmd.orig` (115 sentences, 42 em dashes); check.
- [x] T4: Rewrite `ackwards-visualization.Rmd.orig` (47 sentences, 31 em dashes); check.
- [x] T5: Widen `check_code_unchanged()` to the vignette sources with a planted-edit test;
      run it against the merge base; re-read the AC4 sites against `git show master:<file>`;
      fix drift.
- [x] T6: `Rscript vignettes/precompute.R`; revert timing-only churn in the other five
      vignettes; widen the dod-gate prose step; NEWS entry; `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan as the first vignette batch of the plain-English sweep; depends on M85 for the checker and term list.
- 2026-09-16: criteria audit ([O], fresh context, full mode) ran on the shared vignette-batch wording; its generated-`.Rmd` diff clause was dropped as unenforceable (stamp, PNG, and gt-id churn), and an inline-span guard was added.
- 2026-09-16: implement started on `m086-plain-english-vignettes-a`; question gate skipped (no open choice); the simple-english plugin hook reports pre-existing style hits in the two tracking files, left as they are because cairn history is never rewritten.
- 2026-09-16: substantive amendment (mini gate, user selected): AC4's two engines claims did not exist on master (no warn-and-skip sentence, no never-switches-basis sentence); replaced with the two engines claims the vignette carries (descriptive edges and per-level fit, GP3; the FIML route message and errors, IP6).
- 2026-09-16: re-audit: AC4 (full) — mis-scoped IP6 clause: the WLSMV/ULSMV error belongs to the ESEM estimator path, not the PCA/EFA FIML route; repair adopted.
- 2026-09-16: re-audit: AC4 (full) — residual scope ambiguity in the WLSMV/ULSMV clause; split-clause repair adopted at a user gate (second re-audit line is the stop; no further reader).
- 2026-09-16: substantive amendment (mini gate, user selected): `check_code_unchanged()` read only README.Rmd, R/, and DESCRIPTION, so AC2's instrument could not reach the vignette sources; Scope Out now permits widening it to `vignettes/*.Rmd.orig`, T5 carries the widening plus a planted-edit test, AC2 wording unchanged.
- 2026-09-16: T1 done: intro rewritten; `check-prose.R` clean on the file; all 79 fenced-chunk lines and both inline spans byte-identical to master; 13 terms glossed at first prose use (FIML and congruence absent); the varimax note now says the package "does not offer" oblique rotation instead of "deliberately not offered", matching D-034.
- 2026-09-16: T2 done: suggest-k rewritten; checker clean; 108 chunk lines and 11 inline spans byte-identical to master (the en dash between two spans on the old line 405 became " to ", outside both spans); three criterion headings lost their em dashes; 8 terms glossed (rotation, varimax, factor score, FIML, PCA, congruence, redundancy absent from prose); the D-013 sentence stands at line 69.
- 2026-09-16: T3 done: engines rewritten; checker clean; 188 chunk lines and 1 inline span byte-identical to master (a gt label line retyped with a plain space where master has a no-break space was restored from master by bytes); a gloss paragraph before the comparison table covers PCA, EFA, ESEM, polychoric, ordinal, FIML, and loading, since table cells hit the grep first; the stale in-page anchor `#fiml-for-continuous-pcaefa-via-a-fiml-correlation-matrix` now points at the heading's real anchor `#fiml-for-continuous-pcaefa`; both AC4 engines sites stand (lines 268-271 and 466-469).
- 2026-09-16: T4 done: visualization rewritten; checker clean; 135 chunk lines byte-identical to master, no inline spans; 17 argument headings lost their em dashes; factor and parallel analysis glossed (the only two terms in its prose); the AC1 command over all four sources exits 0.
- 2026-09-16: T5 done: `check_code_unchanged()` now runs its README chunk-and-span extraction over every `vignettes/*.Rmd.orig` on either side of the merge base; `test-check-prose.R` plants a chunk-option edit, a chunk-body edit, and an inline-span edit in a vignette source and sees each reported by file and item (9 expectations, 0 failures); a planted `k_max` edit in the real intro source turned the CLI red (item 23) before the restore; `--code-unchanged master` exits 0 on the branch; styler and lintr clean on both files; the AC4 intro and suggest-k sites re-read against `git show master:` (intro lines 51-54 and 30-32, suggest-k line 69).
- 2026-09-16: T6 done: precompute re-run (45 s); the four generated `.Rmd` regenerated with new stamps; the other four precomputed vignettes and every `vignettes/assets/` PNG reverted to master (re-render churn only; interpret is live and untouched); in the four regenerated files every `#>` output difference is a cli timing line or the check-mark glyph gaining a variation-selector byte from the current cli build, no number changed; dod-gate prose domain widened to the four sources; NEWS documentation entry added and checker-clean; `Rscript tools/dod-gate.R` exit 0 (freshness clean, prose clean, check 0/0/0, coverage 100%, style/lint clean, pkgdown index complete).
- 2026-09-16: claim audit: 163 claims read, 2 corrected — vignettes/ackwards-intro.Rmd.orig, vignettes/ackwards-visualization.Rmd.orig, tools/check-prose.R (the redundancy gloss described the non-default parent walk and now uses `prune()`'s own definition; the guard header now says it reads sources present on both sides of the merge base); both corrections re-read once by the same reader and confirmed; precompute re-run, the two regenerated `.Rmd` recommitted, churn elsewhere reverted; `Rscript tools/dod-gate.R` exit 0 again.
- 2026-09-16: all tasks done; status set to review.

## Decisions

## Review

- 2026-09-16 AC1: `Rscript tools/check-prose.R` over the four batch sources on branch head 44cea84 exits 0 ("Prose OK"). PASS.
- 2026-09-16 AC2: `Rscript tools/check-prose.R --code-unchanged master` exits 0 (merge base d17b646, master unmoved). A grep for an unclosed `` `r `` span on any line finds none in the four sources on either side, and span counts match master (2, 11, 1, 0). PASS.
- 2026-09-16 AC3: `grep -n -i -w` of each of the 16 terms over the four sources. Intro: 13 terms present in prose, each glossed at first prose use (factor 21, parallel analysis 23, factor score 39, ordinal 93, PCA/EFA/ESEM/component/loading 116-121, polychoric 124, varimax/rotation 137-139, redundancy 287, split-half 401). Suggest-k: 9 terms, glossed at 27, 50, 57, 64, 79, 111-113, 181-182, 327 (PCA appears only inside a code span). Engines: 9 terms, glossed at 21-25 and the 32-41 paragraph before the table. Visualization: factor 22 and parallel analysis 389 glossed (loading and polychoric appear only in code). PASS.
- 2026-09-16 AC4: intro lines 28-33 (complementary solutions) and 53-56 (score correlations, not a confirmatory hierarchical model) stand. Suggest-k lines 68-69 ("reports a consensus range, never a single number") stand. Engines lines 269-272 (descriptive score correlations with no sampling distribution of their own) and 357-358 (per-level fit qualifies each level, not the hierarchy as a whole) stand. Engines lines 462-470 state that the FIML PCA/EFA route announces itself with a message and errors for WLSMV/ULSMV and for a non-Pearson PCA/EFA basis. PASS.
- 2026-09-16 AC5: `git diff d17b646..HEAD --name-only` lists only the four batch `.Rmd.orig`, their four regenerated `.Rmd`, the two prose tools, the prose test, NEWS.md, and cairn files (no other vignette and no `vignettes/assets/` file differs). `Rscript tools/dod-gate.R` exit 0 on head 44cea84: vignette-freshness clean, prose clean over the four sources, check 0 errors 0 warnings 0 notes, coverage 100%, styler and lintr clean, pkgdown index complete. NEWS.md carries the "Plain-English vignettes, first batch" documentation entry with no milestone number. PASS.
- 2026-09-16 consistency gate: `cairn_validate.py` exit 0 (all checks pass, 16 advisory work-log warnings, all in M84). No principle changed, impact report skipped. Toolchain slot: `document()` produces no diff (tree clean after the gate's check), README.Rmd untouched, `pkgdown::check_pkgdown()` passed inside the gate, NEWS entry present, no new top-level file, full check 0/0/0. PASS.
