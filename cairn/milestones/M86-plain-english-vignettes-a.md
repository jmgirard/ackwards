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
`` `r ` `` spans (guarded by AC2); changes to `tools/check-prose.R` beyond (a) adding an
abbreviation or banned phrase the corpus shows is missing and (b) widening
`check_code_unchanged()` to read `vignettes/*.Rmd.orig` with its existing README
chunk-and-span extraction (each such change is a work-log line, adds a test where it adds
behavior, and keeps M85's tests green); term-list additions → appended to
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
      (`ackwards-suggest-k`, D-013); the engines vignette's statement that the between-level
      edges are descriptive correlations between factor scores with no sampling distribution
      of their own, and that per-level fit qualifies one level and never the hierarchy as a
      whole (`ackwards-engines`, GP3); and the engines vignette's missing-data statement that
      `missing = "fiml"` announces the PCA/EFA route with a message, and errors rather than
      switching basis or estimator silently, for WLSMV/ULSMV on the ESEM route and for a
      non-Pearson PCA/EFA basis (`ackwards-engines`, IP6).
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

- [x] T1: Rewrite `ackwards-intro.Rmd.orig` (93 sentences, 34 em dashes on master); check.
- [x] T2: Rewrite `ackwards-suggest-k.Rmd.orig` (148 sentences, 37 em dashes); check.
- [x] T3: Rewrite `ackwards-engines.Rmd.orig` (115 sentences, 42 em dashes); check.
- [x] T4: Rewrite `ackwards-visualization.Rmd.orig` (47 sentences, 31 em dashes); check.
- [ ] T5: Widen `check_code_unchanged()` to the vignette sources with a planted-edit test;
      run it against the merge base; re-read the AC4 sites against `git show master:<file>`;
      fix drift.
- [ ] T6: `Rscript vignettes/precompute.R`; revert timing-only churn in the other five
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

## Decisions

## Review
