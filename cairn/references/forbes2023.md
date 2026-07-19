# forbes2023

**Full citation.** Forbes, M. K. (2023). Improving hierarchical models of
individual differences: An extension of Goldberg's bass-ackward method.
*Psychological Methods.* https://doi.org/10.1037/met0000546

**Provenance.** Ingested 2026-07-12 by M57, primarily from the OSF project
below — every oracle value here is code-derived by running Forbes's own
reference implementation, not transcribed from a printed page. The article PDF
is also on the shelf at `cairn/references/sources/forbes2023.pdf` (local only;
gitignored). Pagination: journal pages (advance online publication).
Extraction: verified 2026-07-16 — the numeric claims are re-derived against the
source artifacts by the M44 + M53 fidelity suite, which reproduces the paper's
three simulations and the AMH applied example to 1.3e-14 — observed 2026-07-19.

**OSF project.** `pcwm8` — <https://osf.io/pcwm8/>. Licensed **CC-BY 4.0
International** (<https://creativecommons.org/licenses/by/4.0/>). The project was
briefly CC-BY-NC in early July 2026; Forbes switched it to CC-BY on learning of
the NonCommercial implications, so the AMH matrix can be bundled (legacy M53).

## Why this is the primary source

CLAUDE.md's baseline contract: ackwards' default output must **reproduce
Forbes's examples exactly**. Forbes (2023) footnote 3 names this package
(`github.com/jmgirard/ackwards`) as the reference implementation of the extended
bass-ackward method. Fidelity to her algorithm is therefore test-backed, not
aspirational (M44 + M53).

## Artifacts used as oracles (provenance, not page-transcription)

The oracle **values** are not hand-transcribed from the article; they are
computed by **Forbes's own reference implementation** run on her published
inputs, then frozen as fixtures (so no Forbes code or network runs at test
time). Each artifact below is pinned by OSF guid (+ md5 where an exact byte
match matters).

| Artifact | OSF guid | Role | md5 |
|---|---|---|---|
| `corSpearman_AMH.csv` (155×155 Spearman matrix) | `s9bjz` | The exported `forbes2023` dataset **and** the AMH fidelity oracle input | `c1dd9eca009c2738c268487179d43e87` |
| `ExtendedBassAckwards functions with annotation.R` | `7jfkw` | Her reference implementation — computes the expected `comp.corr` / `cong` / `ChaseCorrPaths` values | `a3e85df897d2a4a4310b9c45dcb068d9` |
| `R script for simulations and applied example_R2.R` | `ztngp` | Regenerates the three simulation Spearman matrices under `set.seed(123)` via `psych::sim.structure` | `8de5ec351992e8fb11c0b8d7015f3329` |

## Which tests / oracles trace to it

- **Oracle O1 — AMH applied example** (k = 10, 45 level-pairs + 54-component
  redundancy chase). Fixture `tests/testthat/fixtures/forbes2023_amh.rds`;
  generator `data-raw/forbes2023.R`; asserted in
  `tests/testthat/test-forbes-fidelity.R`. Matched to 1.3e-14 (M53).
- **Oracle O2 — simulation studies** (three studies, k = 4, 6 level-pairs each).
  Fixture `tests/testthat/fixtures/forbes2023_sims.rds`; generator
  `data-raw/oracle-forbes-sims.R` (M57); asserted in the same test file.

See `cairn/ORACLES.md` for the full registry and the classification of every
oracle in the suite.

## Correspondence conventions (verified M44)

- Her `comp.corr` is `t(W_a) %*% R %*% W_b` with no sign alignment; ours is the
  same algebra with primary-parent alignment, so `|values|` are identical
  entrywise (`W'RW = I` for PCA makes her unstandardized products equal our
  standardized edges).
- Her `cong` is `psych::factor.congruence`, rounded to 2 dp; ours is exact
  Tucker's φ — agreement within 0.005.
- Her `ChaseCorrPaths` uses the **direct/skip-level** correlation to each
  ancestor (non-transitive on deep hierarchies) — adopted as
  `prune("redundant")`'s default `redundancy_criterion = "direct"` (M53).

## Open questions

- None outstanding. Nothing on this page needs anchoring to a printed page or
  table: every oracle value is code-derived from the OSF artifacts above (run by
  her reference implementation), not transcribed — so there is no
  verbatim-critical printed value to re-verify against the article PDF.
