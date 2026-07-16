# M63: User-facing text accuracy pass — done 2026-07-16

**PR:** https://github.com/jmgirard/ackwards/pull/64 (squash d318067)

**Goal:** correct the two user-facing text defects the M62 review surfaced —
Goldberg (2006) citations missing the published title, and the PCA
matrix-input `n_obs` message promising fit statistics PCA never computes.

**Outcome:** all 12 "Doing it all" full-reference sites (roxygen/man + 8
vignette files) carry the Crossref-verified published title verbatim; the
`n_obs` cli message + `@param`/Details roxygen now say `n_obs` is stored as
metadata and enables the N-based sampling-adequacy checks (a new regression
test in test-cor-input.R locks the wording). Two NEWS bullets; no behavior
changes.

**Key decisions/catches:**
- Crossref re-check at implement time found the *reference note*
  (goldberg2006.md) had the wrong issue number (40(3) → 40(4)); the roxygen
  was right. Fixed the note.
- Independent review (scored 90, empirically confirmed): the first corrected
  wording "metadata **only**" was itself false — `n_obs` also feeds the
  factorability N:p screen. Reworded in all four text sites; dod-gate re-run
  green (check 0 err/0 warn/0 note, coverage 100%).
