<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M57: Ossify oracles — reproducible, catalogued oracle discipline

- **Status:** in-progress   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Branch/PR:** m57-ossify-oracles   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal

Make every oracle in the suite catalogued, classified, and reproducible — close the one
un-reproducible fixture, register the rest, and codify the standard — adapting `intraclass`'s
oracle system to ackwards/cairn.

## Scope

**In:**

- **A `cairn/ORACLES.md` registry** that lists *every* oracle in the test suite, each classified by
  type — **live independent-impl** (`psych::fa`, `psych::KMO`, `psych::cortest.bartlett`,
  `psych::bassAckward`, `lavaan::efa`), **frozen fixture** (with the reason it is frozen),
  **closed-form/published**, or **cross-check invariant** (the §5.4 algebra-vs-scores oracle, D-004)
  — with, per entry: **Status** (the asserting test file, grep-verifiable), **Source** (citation or
  registry ref), **Provenance** (generator script or live-computation note). Models intraclass's
  `project/REFERENCES.md` "Oracle registry".
- **A committed generator `data-raw/oracle-forbes-sims.R`** that reproduces
  `tests/testthat/fixtures/forbes2023_sims.rds` end-to-end: the three simulation Spearman matrices
  from Forbes's OSF script under `set.seed(123)`, expected between-level values from her reference
  implementation, md5-pin + `stopifnot` integrity guards, and a **structured `provenance` attr**
  (replacing the current prose-string provenance) — mirroring `data-raw/forbes2023.R` and the M54
  one-download-drives-both pattern. Regenerate the fixture from it.
- **CLAUDE.md Invariant #8** (oracle-backed numerics: ≥2 independent oracle types; no unsourced or
  unreproducible reference value ships), cross-referenced from DESIGN §13, with a note flagging it
  for formal IP numbering when the design-interview candidate runs.
- **A Forbes (2023) reference summary** `cairn/references/forbes2023.md` (+ `INDEX.md` row):
  citation, OSF/DOI, md5, and which oracles/tests trace to it.
- **A guard test** `tests/testthat/test-oracle-provenance.R`: every `tests/testthat/fixtures/*.rds`
  carries a `provenance` attr and a named `data-raw/` generator (or explicit registry exemption);
  fails on an un-sourced fixture.

**Out:**

- **Full freeze of the live-computed oracles** (`psych::fa`/`KMO`/`bartlett`/`bassAckward`,
  `lavaan::efa`) into fixtures → deliberately *not done*; they stay **live** independent-implementation
  cross-checks (the stronger oracle type — a frozen copy is a regression pin, not an oracle). Rationale
  is recorded in the registry. Revisit only if one becomes expensive/network-bound.
- **Formal IP/GP numbering** of the oracle principle → the pending `/design-interview` candidate;
  Invariant #8 is the interim home and carries a forward-note.
- **New oracles for currently-unoracled behavior / new numeric methods** → out; this milestone
  catalogues and makes reproducible what exists, it does not extend numeric coverage.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [ ] **AC1 — Registry complete & classified.** `cairn/ORACLES.md` exists; every oracle asserted in
      the suite maps to exactly one registry entry and every entry names a real asserting test (no
      orphans either way, verified by grepping the oracle tests against the registry); each entry
      carries type, Status, Source, Provenance.
- [ ] **AC2 — Sims fixture reproducible.** `data-raw/oracle-forbes-sims.R` **deterministically**
      regenerates all three simulation matrices from Forbes's exact `fx`/`Phi`/`set.seed(123)`/
      `sim.structure` recipe (re-running reproduces them to `< 1e-12`) and computes expected
      `comp_corr`/`cong`/`corr_chase` from them via her md5-pinned sourced reference implementation;
      the fixture carries a **structured top-level `provenance` attr**; `test-forbes-fidelity.R`
      passes against the regenerated fixture (ackwards matches Forbes's expected values to `< 1e-12`).
      *(Amended 2026-07-12, gate: the committed sims matrices are **replaced** by the reproducible
      realizations — the prior random draw was unrecoverable; the redundancy topology is preserved so
      the test's hardcoded prune literals still hold.)*
- [ ] **AC3 — Principle codified.** CLAUDE.md Invariants gains #8 (oracle-backed numerics, ≥2
      independent oracle types, no unsourced/unreproducible reference value); DESIGN §13
      cross-references it; a forward-note flags it for the IP/GP design-interview pass.
- [ ] **AC4 — Forbes ingested.** `cairn/references/forbes2023.md` exists (full citation, DOI/OSF,
      md5, oracles that trace to it) with an `INDEX.md` row.
- [ ] **AC5 — Provenance guard test.** `test-oracle-provenance.R` asserts every
      `fixtures/*.rds` has a `provenance` attr and a named generator/exemption; it is demonstrated to
      fail on a deliberately un-sourced fixture and pass on the real ones.
- [ ] **AC6 — Gate green.** `Rscript tools/dod-gate.R` clean (check 0 err/0 warn/0 note, tests,
      style, lint, pkgdown).

## Coverage
<!-- owner: plan · create/amend-via-gate; each AC → task(s) by positional number -->

- AC1 → T3
- AC2 → T2
- AC3 → T5
- AC4 → T1
- AC5 → T4
- AC6 → T6

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] **T1 — Ingest Forbes (2023).** Write `cairn/references/forbes2023.md` (citation, DOI
      `10.1037/met0000546`, OSF `pcwm8`, matrix md5 `c1dd9eca…`, the reference-impl file, which
      tests/oracles trace to it) + one `INDEX.md` row. Primary-source-first (tracking-rules
      "Primary sources rule").
- [x] **T2 — Committed sims generator.** Author `data-raw/oracle-forbes-sims.R` reproducing the
      three sim Spearman matrices from Forbes's OSF simulation script under `set.seed(123)`, computing
      expected `comp_corr`/`cong`/`corr_chase` with her reference implementation, with md5/`stopifnot`
      guards and a structured `provenance` attr; regenerate `fixtures/forbes2023_sims.rds`; confirm
      `devtools::test(filter = "forbes-fidelity")` passes unchanged. Add `.Rbuildignore` coverage if
      needed. Model: `data-raw/forbes2023.R`.
- [x] **T3 — Write `cairn/ORACLES.md`.** Sweep all `tests/testthat/test-*.R`, enumerate every
      oracle, classify by type, and record Status/Source/Provenance per entry; include the header
      statement of the standard and the deliberate live-oracle policy (the Scope "Out" rationale).
- [ ] **T4 — Provenance guard test.** Add `tests/testthat/test-oracle-provenance.R` iterating
      `fixtures/*.rds`, asserting a `provenance` attr + a named generator (map fixture → generator via
      a small in-test table or a registry read); prove it fails on a scratch un-sourced fixture, then
      remove the scratch.
- [ ] **T5 — Codify Invariant #8.** Add invariant #8 to CLAUDE.md's Invariants list; add a one-line
      cross-reference in DESIGN §13; add the forward-note for the IP/GP design-interview candidate.
- [ ] **T6 — DoD gate.** `Rscript tools/dod-gate.R` green; no NEWS/pkgdown reference change (no
      exported-behavior change — internal test infra + tracking/dev docs only); confirm the new
      `data-raw/` + `cairn/` files are `.Rbuildignore`d as appropriate.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan (promoted from the 2026-07-12 "ossify oracles" candidate;
  scope confirmed Option 1 — catalogue + surgical freeze, live oracles stay live).
- 2026-07-12: started; branch m57-ossify-oracles cut from master (CRAN 0.1.1 marker committed to
  master first). Confirmed OSF reachable + sims fixture provenance names the exact regen recipe
  ('R script for simulations and applied example_R2.R', set.seed(123), psych::sim.structure) — T2
  feasible as planned. Decided the guard contract (T4): each fixture's `provenance` must be a list
  with `generator` (chr) + `source` (chr); checks structure not file existence (data-raw is
  .Rbuildignore'd, absent at R CMD check).
- 2026-07-12: T1 done — cairn/references/forbes2023.md + INDEX row; oracle values are code-derived
  from OSF artifacts (not page-transcribed), so no PDF vendored.
- 2026-07-12: T2 done — data-raw/oracle-forbes-sims.R (md5-pins her functions guid 7jfkw + sim
  script guid ztngp; transcribes the 3 fx/Phi DGPs; set.seed(123)/sim.structure/spearman; expected
  values via her ExtendedBassAckwards + FindRedundantComp). Regenerated fixture (deterministic: two
  runs maxdiff 0), top-level structured provenance. forbes-fidelity suite green; AC2 amendment
  verified (prune literals unchanged). data-raw already .Rbuildignore'd.
- 2026-07-12: T3 done — cairn/ORACLES.md registry. Full sweep found 12 oracles across 4 types:
  O1/O2 frozen (Forbes AMH+sims), O3–O8 live independent-impl (psych::bassAckward/fa/KMO/bartlett,
  lavaan::efa/fitMeasures, psych::corFiml FIML path), O9/O10 cross-route invariants (§5.4/D-004
  algebra-vs-scores; comparability), O11/O12 closed-form (Fisher-z; suggest_k identities). Each row:
  type/Status(test:line)/Source/Provenance + the deliberate-live-oracle policy.

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

- 2026-07-12 (T2, gate-approved): the committed `forbes2023_sims.rds` matrices were **not
  reproducible** — no fx/Phi/seed/n/method recipe reproduces them (even n=1e5 spearman stays ~0.036
  off); the exact historical draw is lost. Chose to **regenerate the three matrices reproducibly**
  from Forbes's exact recipe (deterministic: two runs maxdiff 0) and recompute expected values via
  her impl, replacing the murky artifact. Fidelity holds (ackwards vs her impl 3.7e-15); redundancy
  topology preserved (prune literals in test-forbes-fidelity.R unchanged). Provenance placed as a
  structured **top-level** attr (uniform target for the T4 guard). Milestone-local (fixture-specific).

## Review
<!-- owner: review · exclusive -->
