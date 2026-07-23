# M72: Correct the stale Forbes contract line + author the Goldberg/Forbes departures ledger

- **Status:** planned
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP9, GP1
- **Branch/PR:** —

## Goal

Correct the pre-D-031 contract wording in `forbes2023.md` and author a consolidated `cairn/references/` synthesis note (the departures ledger) cataloguing every ackwards departure from Goldberg (2006) / Forbes (2023) with its documented rationale and empirical/mathematical support, so the "documented reason + ideally empirical support for departing from the canonical sources" bar is auditable in one place.

## Scope

**In:**
- Correct `cairn/references/forbes2023.md:23–24` — replace the superseded "default output must reproduce Forbes's examples exactly" / "CLAUDE.md's baseline contract" framing with the IP9/D-031 capability framing (reproduction is a permanent *capability*, not a default lock-in), marked in place as a current-knowledge correction (`(corrected M72)`); git holds the original.
- New synthesis note `cairn/references/source-departures.md` from `templates/synthesis-note.md`: Provenance, Scope + tracking disclaimer, Evidence snapshot, neutral characterization, an ID'd ledger table, Disposition, Open questions.
- The ledger covers the **default-level departures** (5 identified) + the 2 source-*matching* defaults, each with source behavior, ackwards behavior, rationale location, and empirical/mathematical support (or explicit "philosophy-only").
- Additive capabilities beyond the sources (EFA/ESEM engines, FIML, out-of-sample scoring, `boot_edges`, `comparability`) get a single "governed by GP1, not departures" note — **not** itemized.
- Any departure lacking empirical support → a search-first ROADMAP `candidate` row; documented honestly in the ledger regardless.
- `INDEX.md` line for the synthesis note.

**Out:**
- Itemizing the additive extensions as ledger entries → the GP1 note instead.
- Strengthening IP9/GP1 into a formal principle change (user declined) → a `candidate` if wanted later.
- The 5 application notes → M71.
- Any code change → docs-only, `cairn/`-only.

## Acceptance criteria

- [ ] AC1: `forbes2023.md` no longer asserts "default output must reproduce Forbes's examples exactly" as a live contract; the claim is replaced with the IP9/D-031 capability framing and marked `(corrected M72)`; `grep` for the old wording returns nothing live (only git history).
- [ ] AC2: a synthesis note exists at `cairn/references/source-departures.md`, authored from `templates/synthesis-note.md`, carrying every required section (Provenance with `Ingested`/`Extraction:` lines, Scope + tracking disclaimer, Evidence snapshot, characterization, an ID'd ledger table, Disposition, Open questions).
- [ ] AC3: the ledger enumerates each of the 5 default-level departures (required `k_max` vs auto-stop; W′RW algebra vs score-then-correlate; tenBerge vs components; primary-parent sign alignment vs unaligned `comp.corr`; exact Tucker φ vs rounded congruence) + the 2 matches (`cut_show = 0.3` = Goldberg .30; `redundancy_criterion = "direct"` = Forbes `ChaseCorrPaths`); each row carries a stable ID and states source behavior, ackwards behavior, rationale location, empirical support.
- [ ] AC4: every departure claim cites a verifiable location (DESIGN §9 row / `D-0NN` / citekey) and, where it claims empirical support, names the oracle/test/source (tong2025; Waller 2007 + IP2 test; fidelity suite 1.3e-14; φ within 0.005); additive extensions noted as GP1-governed, not itemized.
- [ ] AC5: any departure lacking empirical support has a search-first ROADMAP `candidate` row; if none lack it, the Disposition states that explicitly.
- [ ] AC6: `INDEX.md` gains a filename-first line for the note; `cairn_validate` exits 0 (`references index<->disk` PASS); diff is docs-only, `cairn/`-only.

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T2, T3
- AC4 → T3
- AC5 → T4
- AC6 → T5

## Tasks

- [ ] T1: Correct `forbes2023.md` — replace the pre-D-031 "must reproduce exactly" contract framing with IP9/D-031 capability wording, marked `(corrected M72)`; verify the old wording greps clean.
- [ ] T2: Author `source-departures.md` from the template — Provenance/Scope/Evidence snapshot/characterization + the ID'd ledger table skeleton with all 5 departures + 2 matches.
- [ ] T3: Fill each ledger row's rationale location + empirical support, verifying each citation/anchor against its source (DESIGN §9 lines, D-entries, tong2025/Waller/fidelity suite); add the GP1-extensions note; set each row's tag.
- [ ] T4: For any departure lacking empirical support, add a search-first ROADMAP `candidate` row; else state "all documented departures carry support" in the Disposition.
- [ ] T5: Add the `INDEX.md` line; run `cairn_validate` (exit 0); confirm the diff is docs-only `cairn/`-only.

## Work log

- 2026-07-23: created by /milestone-plan.

## Decisions

## Review
