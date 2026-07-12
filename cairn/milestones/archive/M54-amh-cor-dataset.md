# M54: Export `forbes2023` as a bundled dataset — done 2026-07-12

**Goal.** Ship Forbes's (2023) 155-variable "Assessing Mental Health" Spearman correlation
matrix as an exported dataset, extending the M53 test-only fixture.

**Outcome.** `forbes2023` is now bundled + documented; `ackwards(forbes2023, k_max = 10)`
reproduces Forbes's published applied-example hierarchy directly. One md5-pinned OSF generator
(`data-raw/forbes2023.R`) writes **both** `data/forbes2023.rda` and the slimmed fidelity fixture
(`forbes2023_amh.rds`, 117→14 KB — matrix no longer duplicated); the Forbes fidelity test now
reads the matrix from the *exported* dataset (stronger: verifies the shipped object, ≤1e-12 across
all 45 level-pairs). CC-BY 4.0 matrix in an MIT package: Forbes added as data-scoped `cph`,
declared in `LICENSE.note`.

**Key decisions.**
- Single-source the matrix (one generator → dataset + fixture); oracle byte-identical to M53's.
- Named `forbes2023` (author-year) — user override of the plan-gate `amh_cor`.
- `forbes2023` classified as a fidelity/reproduction dataset, outside DESIGN §14.41(b)'s
  declined-third-*teaching*-dataset scope (owner-approved; §14.41(b) + roxygen clarified).

**Verification.** check clean, 0 err/0 warn/0 note (no license NOTE), coverage 100%, style/lint/pkgdown clean, full CI
matrix green. Independent review (2 reviewers + scorer): oracle verified byte-for-byte; 3 findings
(2 fixed, 1 reconciled).

**PR.** https://github.com/jmgirard/ackwards/pull/54
