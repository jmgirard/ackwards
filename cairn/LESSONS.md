# Lessons

Durable, append-only repo lessons (build quirks, testing tricks, gotchas).
One per line: `- YYYY-MM-DD (M<NN>): <lesson>`.

- 2026-07-12 (M54): cairn_validate's ISO-date check flags any slash-separated N/N/N token (e.g. a 0/0/0 clean-check shorthand) as a malformed date — spell it out ("0 err/0 warn/0 note"). Re-proved by M60's archive summary the same day this lesson was pruned for cap space; keep it.
- 2026-07-12 (M54): to bundle a dataset that also backs a fidelity test, drive `data/<name>.rda`
  and the test fixture from one md5-pinned `data-raw/` download so they cannot drift; have the test
  read the matrix from the exported dataset (verifies the shipped object, not a fixture copy).
- 2026-07-12 (M55): CRAN's "cannot suppress console messages" complaint is fixed the R-idiomatic
  way by returning a small **classed** object and moving the rendering into its `print()` method —
  top-level display is unchanged while assignment/nesting go silent. Cheaper and more idiomatic than
  a `verbose=` arg. `cli_text()` writes to stderr (the message stream), so test its output with
  `capture.output(..., type = "message")` and the `cat()` literal with plain `capture.output()`.
- 2026-07-12 (M57): a fixture's *prose* provenance claim ("regenerated from OSF, set.seed(123)") is
  not reproducibility — the M44 sims matrices were provably unrecoverable (no fx/Phi/seed/n/method
  reproduced them; even same R/psych version strings) because no generator was committed. `set.seed`
  is deterministic on a fixed env but the specific draw is lost once the generation path is; pin a
  committed `data-raw/` generator, not a claim. `test-oracle-provenance.R` now enforces this.
- 2026-07-13 (M58): when a milestone fixes a bug caused by copy-drift, grep the ENTIRE codebase for
  the buggy pattern, not just the audit's enumerated sites — the plan scoped `suggest_k`'s
  `n_obs = NA` bug but missed the byte-identical twin in `ackwards()`'s R-matrix `n_obs` check;
  the independent review caught it (`grep -n 'n_obs != as.integer'` would have too).
- 2026-07-13 (M59): line-based coverage can mark an untested branch as covered when `if/else` sits
  on one line — covr counts the line hit if either arm runs. Splitting it during a refactor can
  surface a real, pre-existing gap (here `print`'s non-converged red-cross arm); cover it with a
  test, don't re-collapse the line to game the number.
- 2026-07-16 (M60): don't write "net line reduction" ACs for dedup milestones — extracting *documented* helpers from 3–8-line duplicated blocks is line-neutral (+9 code lines here); make the AC grep-verified single-computation-sites instead.
- 2026-07-16 (M60): cairn_validate's references check wants INDEX.md entries filename-first (`- <name>.md — …`); a citekey-first line reads as a missing entry.
- 2026-07-16 (M62): anchor a grep-verified AC for a wording fix to the *defective context*, not
  the shared phrase — a repo-wide "bass-ackwards method" grep would flag Waller (2007)'s correct
  title; the same phrase can be right in one citation and wrong in another.
- 2026-07-16 (M61): `vignettes/precompute.R` regenerates ALL vignettes — unrelated ones churn
  with pure run noise (cli `[Nms]` timings, gt's random HTML element IDs, PNG re-renders);
  `git checkout --` the untouched vignettes/assets before committing so the PR diff stays scoped.
- 2026-07-16 (M63): when correcting a false absolute claim in user-facing text, enumerate every
  consumer of the value first — the first fix ("n_obs = metadata **only**") traded one false claim
  for a smaller one (n_obs also feeds the factorability N:p screen); the review caught it.
- 2026-07-16 (M63): re-verify a committed reference note against Crossref before propagating it
  repo-wide (supersedes the M56 Waller-title lesson) — goldberg2006.md carried a wrong issue
  number (40(3) vs published 40(4)) while the roxygen it was about to overwrite was right.
- 2026-07-16 (M64): when summarizing a source into a D-entry, any added detail beyond the source is
  a new claim to verify — review caught an invented `lavPredict()` enumeration (extends M63's lesson).
- 2026-07-17 (M65): a file-content hash guard must normalize line endings before hashing — git autocrlf checks out text files as CRLF on Windows, so a raw `tools::md5sum` stamp computed on LF (macOS/Linux) false-fails there; `readLines` + `writeLines(sep="\n")` to a binary connection canonicalizes to LF on both the stamp-write and stamp-verify sides. Windows-only CI caught it.
- 2026-07-17 (M66): to exercise a workflow that only triggers on `schedule`/`workflow_dispatch`, add a temporary branch-scoped `push:` trigger (dispatch registers only once the file is on the default branch) and revert it before review; dispatch-only inputs (`dry_fail`) can't be set on a push, so drive the failure path with a temporary `|| github.event_name == 'push'` on the force-fail step. (Aside: YAML `|` strips the `run:` block indent, so an indented heredoc `EOF` still lands at column 0 at runtime — a recurring reviewer false-positive; verify against the rendered output, not the raw file.)
- 2026-07-19 (M67): a correction is itself a new claim and needs the same verification bar — both independent-review findings on the reference-verification pass were errors *introduced by the correction*, each the very error class being corrected elsewhere in the same diff (a fresh unsourced priority attribution; an unmarked quotation elision). Re-read a correction against the source before committing it, exactly as you would the text it replaces.
- 2026-07-19 (M67): reference pages ingested in one commit can cross-contaminate facts — saucier2005 carried saucier1997's sample N (201), a number absent from its own source, because both pages were written in the same sitting. When verifying a page, treat numbers shared with its commit-siblings as suspect and check each against its own source.
- 2026-07-19 (M67): to check verbatim quotations against a shelf PDF, `pdftotext` + a Python substring search is far cheaper than reading page images, and it settles absence claims outright (e.g. "varimax" occurs zero times). Do NOT grep the flattened text with wide regex windows like `.{130}` — catastrophic backtracking on one long line hangs; use `str.find` plus slicing.
