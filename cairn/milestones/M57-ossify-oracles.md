<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M57: Ossify oracles — reproducible, catalogued oracle discipline

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate; M<xx>, M<yy> or — -->
- **Branch/PR:** m57-ossify-oracles · https://github.com/jmgirard/ackwards/pull/57   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal

Make every oracle in the suite catalogued, classified, and reproducible — close the one
un-reproducible fixture, register the rest, and codify the standard — adapting `intraclass`'s
oracle system to ackwards/cairn.

## Scope

**In:**

- **`cairn/ORACLES.md` registry** — lists *every* oracle in the suite, classified (live / frozen /
  closed-form / cross-check invariant), each with Status (asserting test), Source, Provenance.
  Models intraclass's `project/REFERENCES.md`.
- **`data-raw/oracle-forbes-sims.R`** — committed generator that reproduces
  `fixtures/forbes2023_sims.rds` end-to-end (Forbes's OSF sim recipe, `set.seed(123)`; expected
  values via her reference impl; md5-pin + `stopifnot` guards; structured `provenance`).
- **CLAUDE.md Invariant #8** (oracle-backed numerics), cross-ref'd from DESIGN §13, with a
  forward-note for formal IP numbering at the design-interview.
- **Forbes (2023) reference summary** `cairn/references/forbes2023.md` (+ `INDEX.md` row).
- **Guard test** `tests/testthat/test-oracle-provenance.R` — every `fixtures/*.rds` carries a
  `provenance` attr naming a `data-raw/` generator + source; fails on an un-sourced fixture.

**Out:**

- **Full freeze of live oracles** (`psych`/`lavaan`) → deliberately *not done*; they stay live
  independent-impl cross-checks (a frozen copy is a regression pin, not an oracle). Revisit only if
  one becomes expensive/network-bound.
- **Formal IP/GP numbering** of the principle → the pending `/design-interview` candidate (Invariant
  #8 is the interim home).
- **New oracles / new numeric methods** → out; this catalogues and makes reproducible what exists.

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [x] **AC1 — Registry complete & classified.** `cairn/ORACLES.md` exists; every oracle asserted in
      the suite maps to exactly one registry entry and every entry names a real asserting test (no
      orphans either way, verified by grepping the oracle tests against the registry); each entry
      carries type, Status, Source, Provenance.
- [x] **AC2 — Sims fixture reproducible.** `data-raw/oracle-forbes-sims.R` **deterministically**
      regenerates all three simulation matrices from Forbes's exact `fx`/`Phi`/`set.seed(123)`/
      `sim.structure` recipe (re-running reproduces them to `< 1e-12`) and computes expected
      `comp_corr`/`cong`/`corr_chase` from them via her md5-pinned sourced reference implementation;
      the fixture carries a **structured top-level `provenance` attr**; `test-forbes-fidelity.R`
      passes against the regenerated fixture (ackwards matches Forbes's expected values to `< 1e-12`).
      *(Amended 2026-07-12, gate: the committed sims matrices are **replaced** by the reproducible
      realizations — the prior random draw was unrecoverable; the redundancy topology is preserved so
      the test's hardcoded prune literals still hold.)*
- [x] **AC3 — Principle codified.** CLAUDE.md Invariants gains #8 (oracle-backed numerics, ≥2
      independent oracle types, no unsourced/unreproducible reference value); DESIGN §13
      cross-references it; a forward-note flags it for the IP/GP design-interview pass.
- [x] **AC4 — Forbes ingested.** `cairn/references/forbes2023.md` exists (full citation, DOI/OSF,
      md5, oracles that trace to it) with an `INDEX.md` row.
- [x] **AC5 — Provenance guard test.** `test-oracle-provenance.R` asserts every
      `fixtures/*.rds` has a `provenance` attr and a named generator/exemption; it is demonstrated to
      fail on a deliberately un-sourced fixture and pass on the real ones.
- [x] **AC6 — Gate green.** `Rscript tools/dod-gate.R` clean (check 0 err/0 warn/0 note, tests,
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
- [x] **T4 — Provenance guard test.** Add `tests/testthat/test-oracle-provenance.R` iterating
      `fixtures/*.rds`, asserting a `provenance` attr + a named generator (map fixture → generator via
      a small in-test table or a registry read); prove it fails on a scratch un-sourced fixture, then
      remove the scratch.
- [x] **T5 — Codify Invariant #8.** Add invariant #8 to CLAUDE.md; cross-ref DESIGN §13; forward-note for the IP/GP design-interview candidate.
- [x] **T6 — DoD gate.** `Rscript tools/dod-gate.R` green; no NEWS/pkgdown change (no exported-
      behavior change — internal test infra + docs only); new `data-raw/` files `.Rbuildignore`d.

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-12: created by /milestone-plan (promoted from the "ossify oracles" candidate; scope Option 1 — catalogue + surgical freeze, live oracles stay live).
- 2026-07-12: started; branch cut from master (CRAN 0.1.1 marker committed first); confirmed OSF reachable + T2 feasible; fixed the T4 guard contract (fixture `provenance` = list with `generator`+`source`, structure-checked since data-raw is .Rbuildignore'd).
- 2026-07-12: T1 done — cairn/references/forbes2023.md + INDEX row (oracle values code-derived from OSF artifacts, no PDF vendored).
- 2026-07-12: T2 done — data-raw/oracle-forbes-sims.R (md5-pins guids 7jfkw+ztngp; 3 fx/Phi DGPs; expected values via her impl); regenerated fixture deterministic (maxdiff 0), top-level provenance; forbes-fidelity green, prune literals unchanged.
- 2026-07-12: T3 done — cairn/ORACLES.md: 12 oracles, 4 types (O1/O2 frozen, O3–O8 live, O9/O10 invariant, O11/O12 closed-form), each type/Status(test:line)/Source/Provenance + live-oracle policy.
- 2026-07-12: T4 done — test-oracle-provenance.R guard + normalized AMH fixture to top-level provenance via data-raw/forbes2023.R (data/forbes2023.rda byte-identical); teeth demonstrated.
- 2026-07-12: T5 done — CLAUDE.md Invariant 8 + DESIGN §13 cross-ref + IP/GP forward-note.
- 2026-07-12: T6 done — DoD gate green (check 0/0/0, coverage 100%, style+lint clean, pkgdown index); status → review.

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

_Reviewed 2026-07-12 (same session). PR #57._

**AC evidence (fresh):**
- AC1 — 12 registry entries; all cited test:line anchors resolve to real oracle assertions; reverse sweep found no unregistered oracle-bearing test (2 grep hits are lavaan-k-limit comments).
- AC2 — generator run twice → max|R diff| 0 (deterministic); `test-forbes-fidelity.R` green; fixture carries top-level structured provenance.
- AC3 — CLAUDE.md Invariant 8, DESIGN §13 cross-ref, IP/GP forward-note all present (grep-confirmed).
- AC4 — `cairn/references/forbes2023.md` + INDEX row present.
- AC5 — `test-oracle-provenance.R` green (both fixtures pass; teeth verified on scratch fixtures).
- AC6 — DoD gate green (check 0/0/0, coverage 100%, styler/lintr clean, pkgdown index); fresh `document()` no-diff + `pkgdown::check_pkgdown()` pass.

**Consistency gate:** `cairn_validate` all-pass; Coverage map AC1–AC6 → T3/T2/T5/T1/T4/T6 (all exist); no DESIGN IP/GP changed (skip impact); no user-facing change (no NEWS); no stray-file NOTE.

**Independent review (two lenses):** no substantive defects. [O] diff-bug verified the DGP
transcription against the live OSF script (both md5 pins match bit-for-bit), guard-test teeth, and
every registry anchor. [S] blame-history re-ran the full suite (2254 pass/0 fail), byte-confirmed
`data/forbes2023.rda` unchanged + AMH payload identical, and verified the sims replacement is a
bug-fix to a *false* M44 provenance claim (not a reversal), with D-004/D-018 cited correctly.
Triage: 1 minor finding (score ~85) — stale md5 cells in `references/forbes2023.md` — **fixed**
(filled `7jfkw`/`ztngp` md5s). No findings scored <80 to log.
