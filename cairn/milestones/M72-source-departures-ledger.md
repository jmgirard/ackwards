# M72: Correct the stale Forbes contract line + author the Goldberg/Forbes departures ledger

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP9, GP1
- **Branch/PR:** m72-source-departures-ledger · https://github.com/jmgirard/ackwards/pull/75

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
- **Maintenance hook** (amendment 2026-07-23, user-approved) keeping the ledger current: (a) a Maintenance clause in `source-departures.md` tying any new/changed departure — which per IP9 always carries a D-entry — to a ledger row in the same change; (b) a one-line pointer in the DESIGN §9 defaults table to the ledger; (c) an anchor-integrity check (`tools/check-ledger-anchors.R` + a testthat wrapper) asserting every `D-0NN` / `IPn`/`GPn` / `[[citekey]]` / cited R file the ledger names resolves, skipping in the built package where `cairn/` is absent.

**Out:**
- Itemizing the additive extensions as ledger entries → the GP1 note instead.
- Strengthening IP9/GP1 into a formal principle change (user declined) → a `candidate` if wanted later.
- The 5 application notes → M71.
- A fully mechanical "is the ledger complete" guard → infeasible (detecting an unrecorded departure needs semantic judgment); the anchor check verifies only that cited anchors resolve, the process clause covers completeness.
- Beyond the anchor test + its wrapper, no other code change — the ledger/correction/pointer stay `cairn/`-only.

## Acceptance criteria

- [ ] AC1: `forbes2023.md` no longer asserts "default output must reproduce Forbes's examples exactly" as a live contract; the claim is replaced with the IP9/D-031 capability framing and marked `(corrected M72)`; a grep of `cairn/references/` for the old contract assertion returns nothing (git history + this milestone file's task description aside).
- [ ] AC2: a synthesis note exists at `cairn/references/source-departures.md`, authored from `templates/synthesis-note.md`, carrying every required section (Provenance with `Ingested`/`Extraction:` lines, Scope + tracking disclaimer, Evidence snapshot, characterization, an ID'd ledger table, Disposition, Open questions).
- [ ] AC3: the ledger enumerates each of the 5 default-level departures (required `k_max` vs auto-stop; W′RW algebra vs score-then-correlate; tenBerge vs components; primary-parent sign alignment vs unaligned `comp.corr`; exact Tucker φ vs rounded congruence) + the 2 matches (`cut_show = 0.3` = Goldberg .30; `redundancy_criterion = "direct"` = Forbes `ChaseCorrPaths`); each row carries a stable ID and states source behavior, ackwards behavior, rationale location, empirical support.
- [ ] AC4: every departure claim cites a verifiable location (DESIGN §9 row / `D-0NN` / citekey) and, where it claims empirical support, names the oracle/test/source (tong2025; Waller 2007 + IP2 test; fidelity suite 1.3e-14; φ within 0.005); additive extensions noted as GP1-governed, not itemized.
- [ ] AC5: any departure lacking empirical support has a search-first ROADMAP `candidate` row; if none lack it, the Disposition states that explicitly.
- [ ] AC6: `INDEX.md` gains a filename-first line for the note; `cairn_validate` exits 0 (`references index<->disk` PASS); the `cairn/` documentation changes carry no package-code edits beyond the anchor test (T7); the full DoD gate (`Rscript tools/dod-gate.R`) passes (check 0/0/0, coverage maintained, style/lint/pkgdown clean).
- [ ] AC7: the maintenance hook exists — a Maintenance clause in `source-departures.md` + a one-line pointer in the DESIGN §9 defaults table; `tools/check-ledger-anchors.R` passes on the current ledger **and** fails when a cited anchor is broken (mutation-verified); `tests/testthat/test-ledger-anchors.R` runs the check and skips when `cairn/` is absent (built package).

## Coverage

- AC1 → T1
- AC2 → T2
- AC3 → T2, T3, T6
- AC4 → T3, T6
- AC5 → T4
- AC6 → T5, T8
- AC7 → T7

## Tasks

- [x] T1: Correct `forbes2023.md` — replace the pre-D-031 "must reproduce exactly" contract framing with IP9/D-031 capability wording, marked `(corrected M72)`; verify the old wording greps clean.
- [x] T2: Author `source-departures.md` from the template — Provenance/Scope/Evidence snapshot/characterization + the ID'd ledger table skeleton with all 5 departures + 2 matches.
- [x] T3: Fill each ledger row's rationale location + empirical support, verifying each citation/anchor against its source (DESIGN §9 lines, D-entries, tong2025/Waller/fidelity suite); add the GP1-extensions note; set each row's tag.
- [x] T4: For any departure lacking empirical support, add a search-first ROADMAP `candidate` row; else state "all documented departures carry support" in the Disposition.
- [x] T5: Add the `INDEX.md` line; run `cairn_validate` (exit 0); confirm the diff is docs-only `cairn/`-only.
- [x] T6: (review fixes) Reword E3's "see M1" misdirect (finding 1); drop the dangling "DESIGN §9" anchor from M1's rationale (finding 2); record the 58%/71% recovery figure in `tong2025.md` with its anchor and update that note's Open-questions so E1 is fully traceable.
- [ ] T7: Build the maintenance hook — Maintenance clause in `source-departures.md`; one-line pointer in the DESIGN §9 defaults table; `tools/check-ledger-anchors.R` (parse the ledger's `D-0NN`/`IPn`/`GPn`/`[[citekey]]`/R-file anchors, assert each resolves) + `tests/testthat/test-ledger-anchors.R` (skips when `cairn/` absent); mutation-verify the check fails on a broken anchor.
- [ ] T8: Full DoD gate (`Rscript tools/dod-gate.R`) + `cairn_validate`; fix any fallout.

## Work log

- 2026-07-23: created by /milestone-plan.
- 2026-07-23: set in-progress; branch m72-source-departures-ledger cut from master.
- 2026-07-23: T6 — fixed findings 1 (E3 "see M1" reworded; the PCA-matches fact moved into E3's ackwards column) and 2 (dropped dangling "DESIGN §9" from M1). Traceability: rendered tong2025 p.14 (Results) — the 58%/71% figure is on **p. 14**, not p. 11 as the diff-bug reviewer stated; recorded it verbatim in tong2025.md Extracted values + anchored E1 to p. 14; updated tong2025.md Open-questions (headline rates now recorded, per-condition detail not).
- 2026-07-23: T5 — added INDEX.md "Design-provenance synthesis" section + source-departures.md line; cairn_validate exit 0 (references index<->disk PASS); diff cairn/-only. All tasks done → review.
- 2026-07-23: T4 — no departure is a depart-gap (E1–E5 all supported), so no candidate row spawned; Disposition states this explicitly (satisfies AC5's else-branch).
- 2026-07-23: T2+T3 — authored source-departures.md: 5 departures (E1–E5) + 2 matches (M1–M2), each with rationale location + empirical/math support and a tag (all E1–E5 `depart-supported`). Verified every anchor exists (R/ackwards.R:34/:292; D-004/007/010/017/031 + extension D-003/020/021/022/023; IP1/2/4/9 + GP1; 6 linked notes) and every empirical figure traces to source (tong2025 58%/71%, fidelity 1.3e-14, φ 0.005, chase 54/54). GP1-extensions note added.
- 2026-07-23: T1 — corrected forbes2023.md "Why this is the primary source" to the IP9/D-031 capability framing (was the superseded "default output must reproduce Forbes's examples exactly"), marked (corrected M72). Refined AC1 to scope its grep to cairn/references/ (the milestone file legitimately quotes the old wording). Old assertion greps clean in references/.

## Decisions

## Review

**Fresh evidence per acceptance criterion** (2026-07-23, PR #75):

- AC1 ✓ — `grep` of `cairn/references/` for "reproduce … examples exactly" returns nothing; forbes2023.md now carries the IP9/D-031 capability framing + `(corrected M72)` marker (2 hits).
- AC2 ✓ — `source-departures.md` exists with every required template section: Provenance (`Ingested`/`Extraction:`), Scope + "reference, not an authority" disclaimer, Evidence snapshot, characterization, an ID'd ledger table, Disposition, Open questions (all grep-confirmed present).
- AC3 ✓ — ledger holds 5 departure rows (E1–E5) + 2 match rows (M1–M2), each with a stable ID, source behavior, ackwards behavior, rationale location, and support column.
- AC4 ✓ — every departure cites a verifiable location (DESIGN §9 / D-0NN / citekey) and names its empirical/mathematical support (tong2025 58%/71%; Waller proof + IP2; fidelity 1.3e-14; φ 0.005); all anchors verified to exist at T3 (D-entries, IP1/2/4/9+GP1, `R/ackwards.R:34`/`:292`, 6 linked notes); GP1-extensions note present (line 35).
- AC5 ✓ — no departure is a `depart-gap` (E1–E5 all `depart-supported`), so no candidate row spawned; Disposition states this explicitly; ROADMAP diff adds no candidate row.
- AC6 ✓ — INDEX.md "Design-provenance synthesis" line present; `cairn_validate` exit 0 (`references index<->disk` PASS); `git diff --name-only master..HEAD` is entirely under `cairn/`.

**Consistency gate (r-package `consistency-gate` slot + universal cairn checks):**
- `cairn_validate` exit 0 — every check PASS (incl. `coverage complete`, `references index<->disk`, `weight caps`).
- `cairn_impact` skipped — Principles touched IP9/GP1 are *worked under*, not changed; the diff edits no DESIGN.md principle text.
- `devtools::check()` — **not re-run, deliberately:** the diff is entirely under `cairn/`, which is `.Rbuildignore`d (`^cairn$`), so it is build-excluded and cannot affect `check()`; last green at M70 (dac7f2b, same day). No `R/`, `man/`, DESCRIPTION, NAMESPACE, README, or vignette touched → `document()` no-diff trivially, pkgdown unaffected.
- NEWS.md — justified skip: `cairn/` is internal tracking, not user-facing; no behavior/API/doc-page change.

**Independent fresh-context review (3 lenses + scorer), pass 1:**
- [S] blame-history (Sonnet): **0 findings** — forbes2023.md correction is the D-031-mandated alignment (not an undo); every ledger row consistent with its cited D-entry.
- [S] prior-review (Sonnet): **0 findings** — no regression of the M67/M69 correction-verification, M67 cross-contamination, or M60/M68 format lessons; figures trace to source; `gh` probe found no PR threads.
- [O] diff-bug (Opus): independently confirmed all empirical figures (tong2025 58%/71% at p. 11 + Table 1; fidelity 1.3e-14; φ 0.005; chase 54/54), every anchor, and the forbes2023 correction. **3 findings**, scored by a fresh [S] scorer:
  - Finding 1 (score 92, **actioned** T6): E3's "see M1" cross-reference misdirects — M1 is the unrelated `cut_show` row.
  - Finding 2 (score 88, **actioned** T6): M1's rationale cites "DESIGN §9" for `cut_show = 0.3`, but DESIGN.md doesn't document the show-cut — dangling anchor (other two anchors correct).
  - Finding 3 (score 35, **logged, not actioned**): E3's classification as a Goldberg departure is borderline (tenBerge is factor-engine-only); scorer + reviewer agree it is intentional, transparently-scoped structure, not a defect — left as-is.
- Also actioned in T6 (traceability, surfaced by blame-history): `tong2025.md`'s Open-questions said the recovery rates "were not transcribed", but E1 now cites 58%/71% — record that figure in tong2025.md with its p. 11 anchor.
- Milestone returned to in-progress to action findings 1–2 + the traceability fix + the user-approved maintenance hook (T6–T8); ACs re-verified in pass 2 below.
