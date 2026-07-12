# Oracle registry

The catalogue of every oracle in the test suite. **The standard (CLAUDE.md
Invariant 8):** every numeric result is verified against ≥2 independent oracle
*types*, and **no unsourced or unreproducible reference value ships** — an oracle
value is either published/closed-form, an independent implementation, a
reproducible seeded artifact, or a cross-route invariant. This file is the map
from each oracle to the test that asserts it and the source/provenance that
backs it.

**Asserted-state is single-sourced to the test file.** The `Status` column names
the test that asserts an oracle — that grep-verifiable fact is the truth, not
this table. An oracle appears here once it is asserted; a not-yet-written one is
planned in its milestone, never here.

## Oracle types

- **frozen** — an external reference value computed once (network/expensive) and
  committed as a fixture with a reproducible `data-raw/` generator + structured
  `provenance` attr. Guarded by `test-oracle-provenance.R`.
- **live** — an *independent implementation* (`psych`, `lavaan`) recomputed at
  test time. This is the **stronger** form and is deliberately **not** frozen:
  a frozen copy is a regression pin, not a cross-check — if the independent
  implementation ever diverges we want the test to fire. Only freeze a live
  oracle if it becomes expensive or network-bound.
- **invariant** — two independent ackwards routes that must agree (the §5.4 /
  D-004 algebra-vs-scores cross-check and kin). No external source; the agreement
  itself is the oracle.
- **closed-form** — a published/definitional formula recomputed in the test with
  deliberately dumb explicit code.

## Registry

| ID | Type | What it pins | Status (asserting test) | Source | Provenance |
|---|---|---|---|---|---|
| **O1** | frozen | AMH applied example, k=10: 45 level-pair edges, 54-component redundancy chase | `test-forbes-fidelity.R` ("…AMH applied example") | Forbes (2023) via her `ExtendedBassAckwards` reference impl — [`references/forbes2023.md`](references/forbes2023.md) | `data-raw/forbes2023.R` (OSF `s9bjz`/`7jfkw`, md5-pinned) → `fixtures/forbes2023_amh.rds` |
| **O2** | frozen | Three simulation studies, k=4: 6 edges + 9-component chase each | `test-forbes-fidelity.R` ("…simulation examples") | Forbes (2023) `sim.structure` DGP + her reference impl — [`references/forbes2023.md`](references/forbes2023.md) | `data-raw/oracle-forbes-sims.R` (OSF `ztngp`/`7jfkw`, md5-pinned; `set.seed(123)`, deterministic) → `fixtures/forbes2023_sims.rds` |
| **O3** | live | PCA between-level edges reproduce the original method | `test-pca.R:30` | `psych::bassAckward` (Revelle) | recomputed at test time on `psych::bfi` |
| **O4** | live | EFA between-level edges | `test-efa.R:223` | `psych::bassAckward(fm="minres")` | recomputed at test time (documents psych's k=1 mis-standardization) |
| **O5** | live | EFA per-level fit statistics | `test-efa.R:148` | `psych::fa` (LR χ² statistic) | recomputed at test time |
| **O6** | live | `factorability()` KMO + Bartlett (Pearson and polychoric) | `test-factorability.R:24, 170` | `psych::KMO`, `psych::cortest.bartlett`, `psych::polychoric` | recomputed at test time on a fixed matrix |
| **O7** | live | ESEM scaled fit (χ²/CFI scaled variant) | `test-esem.R:699` | `lavaan::efa` + `lavaan::fitMeasures` (MLR) | recomputed at test time |
| **O8** | live | FIML-missing edge R (PCA/EFA and ESEM paths) | `test-missing.R:493, 507` (+ saturated `lavInspect` :379) | `psych::corFiml`; lavaan FIML saturated model | recomputed at test time |
| **O9** | invariant | algebra vs materialized-scores edges agree for linear engines (§5.4 / D-004) | `test-compute_edges.R:1, 141`; `test-scores.R:180`; `test-efa.R:246` | none — internal cross-route agreement | two ackwards edge routes, `< 1e-6` |
| **O10** | invariant | `comparability()` pooled-R algebra vs correlating actual scores | `test-comparability.R:92, 109` | none — internal cross-route agreement | two ackwards routes |
| **O11** | closed-form | `boot_edges()` fixed-weights CI = Fisher-z analytic; full-pipeline SE ≥ fixed-weights | `test-boot_edges.R:141` | Fisher-z transform (closed form) | recomputed in-test |
| **O12** | closed-form | `suggest_k()` definitional identities (VSS optima = `which.max`; PA-PC = parallel boundary) | `test-suggest_k.R:117, 138` | definitional | recomputed in-test |

## Policy notes

- **Deliberately live (not frozen): O3–O8.** These are deterministic, sub-second,
  in-package (`psych` is an Import; `lavaan` a Suggest). Freezing them would trade
  a live independent-implementation cross-check for a self-referential pin. Revisit
  only if one becomes expensive or network-bound — the same trigger that justified
  freezing O1/O2 (network/OSF).
- **≥2 independent types per estimator (Invariant 8).** The PCA/EFA edge path is
  pinned by both a live external impl (O3/O4) and the algebra-vs-scores invariant
  (O9); ESEM by O7 + O9; the Forbes contract by O1/O2 (frozen external) + O9
  (invariant). No estimator rests on a single oracle type.
- **New fixtures.** Any new `tests/testthat/fixtures/*.rds` must carry a structured
  top-level `provenance` attr naming its `data-raw/` generator and `source` —
  mechanically enforced by `test-oracle-provenance.R`.
