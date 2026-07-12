# M55: Address CRAN 0.1.0 feedback and resubmit as 0.1.1 — done 2026-07-12

**Goal:** Fix the two CRAN reviewer items on the 0.1.0 submission and ready
a 0.1.1 resubmission from master.

**Outcome:**
- DESCRIPTION Description spells out PCA / EFA / ESEM at first use.
- `label_template()` no longer writes to the console with `cat()`: it
  returns its named vector **visibly** with class `ackwards_labels`, and the
  editable `c(...)` scaffold renders via a new `print.ackwards_labels()`
  method — top-level calls unchanged, assignment/nested use now silent
  (CRAN's suppressibility requirement). Verified sole offending call site.
- Version → 0.1.1; NEWS 0.1.1 section; `cran-comments.md` rewritten as a
  point-by-point resubmission. Precomputed vignettes regenerated.

**Verification:** 7/7 ACs; consistency gate clean; 2-lens review (Opus
diff-bug + Sonnet blame-history) → 0 findings; dod-gate 0/0/0 + 100%
coverage (0 err/0 warn/0 note); full CI matrix + win-builder R-devel green.

**Key decisions:** resubmit from master as 0.1.1 (not a patch branch from
the v0.1.0 tag — superseded the pre-M55 release-tail guidance at the plan
gate); classed-return + print-method design chosen as the house idiom over
`verbose=`/`message()`.

**PR:** https://github.com/jmgirard/ackwards/pull/55
