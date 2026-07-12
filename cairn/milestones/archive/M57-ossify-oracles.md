# M57: Ossify oracles — reproducible, catalogued oracle discipline — done 2026-07-12

**Goal.** Make every oracle catalogued, classified, and reproducible — close the one un-reproducible
fixture, register the rest, codify the standard — adapting `intraclass`'s oracle system to cairn.

**Outcome.** New `cairn/ORACLES.md` registry catalogues all 12 oracles by type (O1/O2 frozen Forbes;
O3–O8 live `psych`/`lavaan` independent-impl; O9/O10 algebra-vs-scores + comparability invariants;
O11/O12 closed-form), each with asserting test:line + provenance. New committed generator
`data-raw/oracle-forbes-sims.R` **replaces** the unreproducible `forbes2023_sims.rds` matrices with
deterministic realizations of Forbes's exact `set.seed(123)`/`sim.structure` recipe (fidelity 3.7e-15).
CLAUDE.md **Invariant #8** + DESIGN §13 cross-ref; Forbes (2023) ingested as
`cairn/references/forbes2023.md`; guard test `test-oracle-provenance.R` blocks any future fixture
lacking a structured `provenance` attr naming its `data-raw/` generator + source. No exported behavior
changed (`data/forbes2023.rda` byte-identical); DoD gate green (0/0/0, cov 100%).

**Key decisions.** Option 1 — catalogue + surgical freeze: live `psych`/`lavaan` oracles stay live (a
frozen copy is a regression pin, not a cross-check). AC2 amended at gate: the pre-M57 sims draw was
provably unrecoverable, so regenerate reproducibly rather than match it. Invariant #8 is the interim
home; fold into a numbered DESIGN IP/GP at the pending `/design-interview` pass.

**Review.** All 6 ACs verified fresh. Two independent reviewers found no substantive defects (full
suite 2254 pass/0 fail; sims replacement confirmed a bug-fix to a false M44 provenance claim, not a
reversal). One minor finding — stale md5 cells in `references/forbes2023.md` — fixed.

PR: https://github.com/jmgirard/ackwards/pull/57 · merged on local-green (non-release).
