# M65: Mechanize the precomputed-vignette staleness guard — done 2026-07-17

**Outcome:** editing a `vignettes/*.Rmd.orig` without re-running `precompute.R` is now a mechanical
failure, not a CLAUDE.md prose warning. `precompute.R` stamps each generated `.Rmd` with an HTML
comment carrying the LF-normalized md5 of its `.Rmd.orig`; `tools/check-vignette-freshness.R`
(base-R, source-able) re-derives it and fails on a stale, missing, or orphaned stamp (both
directions). Enforced in three places: the CI source-checkout step (`R-CMD-check.yaml`, before the
tarball check), `tools/dod-gate.R` (fail-fast), and a skip-aware testthat wrapper
(`test-vignette-freshness.R`, runs in source, skips under the tarball where `.orig` are absent).
CLAUDE.md prose repointed at the guard; stale vignette count fixed (eight of nine precomputed).

**Key decisions / notes:**
- AC2 amended at the question gate: count-agnostic wording (8 `.Rmd.orig`, not the planned "seven").
- Discovered sub-task: added freshness enforcement to `dod-gate.R` (the testthat wrapper skips
  under `check()`'s tarball, so the local gate needed a direct call).
- All 8 committed stamps already equalled the LF-normalized md5, so no vignette regeneration was
  needed for the CRLF fix — checker/precompute.R only.

**Review:** three-lens fan-out clean (no findings); scorer no-op. Windows CI caught a real bug —
`git autocrlf` CRLF checkout made raw-byte md5 false-STALE on all 8 vignettes; fixed by
LF-normalizing before hashing (see LESSONS 2026-07-17). All 6 ACs met; CI fully green incl. Windows.

**PR:** https://github.com/jmgirard/ackwards/pull/68 (squash `443e352`).
