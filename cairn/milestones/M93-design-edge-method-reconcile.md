# M93: Reconcile DESIGN's edge_method text with the code

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** —
- **Driving RR:** —
- **Principles touched:** IP1, IP2
- **Resolves:** —
- **Surface tier:** internal — DESIGN.md and DECISIONS.md are in-repo design records, and no exported behavior or user doc changes
- **Branch/PR:** m93-design-edge-method-reconcile

## Goal

`cairn/DESIGN.md` describes the edge scores route as it is: an internal seam of `compute_edges()` that the algebra-vs-scores tests use, not a route a user can choose.

## Scope

**In:**

- Reword IP1 and IP2, with a D-entry that narrows D-004.
- Rewrite the scores-route text in §5: the intro, the §5.2 descriptor, the §5.3 pseudocode comment and "Two situations" list, and §5.4.
- Rewrite the §9 `edge_method` row as internal, the §9 scores-method row's EAP clause to match D-007, and the Documentation standard line on `edge_method = "scores"`.
- Shorten the first Known limitations entry to the limitation plus a pointer to IP2.

**Out:**

- A user-facing option for sample-realized edges → the `[low]` candidate row this plan adds.
- The manuscript's fallback-to-scores sentence → the manuscript accuracy candidate row.
- `R/`, `man/`, vignettes, and NEWS. The sweep found no user-facing promise of a scores route there, and the `compute_edges()` roxygen is accurate for an internal function.
- D-004 and D-007 text, which is history (IP4).

## Acceptance criteria

- [ ] AC1: In `cairn/DESIGN.md`, IP1 states that every call of `compute_edges()` in `R/` passes `edge_method = "auto"` or `"algebra"` with no data, so every shipped edge comes from the `W'RW` algebra. IP1 no longer says "or when the user asks", and its claim matches what `git grep -n "compute_edges(" -- R` lists on the branch head. IP2 states that the scores route stays inside the internal `compute_edges()` as the second route of the algebra-vs-scores agreement tests. IP2 keeps its requirement that a standing test asserts agreement for every linear engine, and it no longer presents `edge_method = "scores"` as available to users. Both bullets keep their numbers.
- [ ] AC2: `cairn/DECISIONS.md` gains one entry that records the IP1 and IP2 rewording with its rationale and the evidence that would reverse it. Its heading names D-004 as the entry it narrows, and its body cites D-007 and D-031.
- [ ] AC3: Two parts of `cairn/DESIGN.md` are checked. The first part is the "Design principles" section, §5, §9, and the "Known limitations" section, each read in full. The second part is every other line that `grep -n -i -E 'edge_method|scores.? route|scores.?path|materiali|user asks|user wants|user would prefer' cairn/DESIGN.md` or `grep -n -E '\bEAP\b' cairn/DESIGN.md` returns. No text in either part presents the scores route, `edge_method`, or EAP scoring as a setting of an exported function.
- [ ] AC4: The §9 `edge_method` row says that `edge_method` belongs to the internal `compute_edges()` and is not an argument of any exported function. The §9 scores-method row states EAP as out of scope (D-007), not as an opt-in.
- [ ] AC5: The first "Known limitations" entry states the limitation: the algebra-vs-scores cross-check does not cover the polychoric or the `missing = "fiml"` PCA/EFA paths, and why. It points to IP2 for how edges are built, keeps its "Corrected M91" mark, and is shorter than the same entry on `master` (`wc -c` of each extracted entry).
- [ ] AC6: The milestone's diff against `master` (`git diff --stat master...m93-design-edge-method-reconcile`) changes files under `cairn/` only. `python3 /Users/jmgirard/.claude/skills/cairn/scripts/cairn_validate.py` passes on the branch head.

## Coverage

- AC1 → T2, T6
- AC2 → T1
- AC3 → T3, T4, T5, T6
- AC4 → T4
- AC5 → T5
- AC6 → T6

## Tasks

- [x] T1: Append D-038 to `cairn/DECISIONS.md`, with a heading that names D-004 as narrowed. Its subject: the scores route is an internal cross-check seam, not a user option. Context: no exported function offers the route, and EAP is out of scope (D-007). Decision: the IP1 and IP2 rewording, with the §9 row kept and marked internal. Consequences: D-004's `_Source` pointer to the §9 row still resolves, and the IP change follows D-031's procedure. A user need for sample-realized edges (the new candidate row) reverses it. (RB tripwire: ip-touching)
- [x] T2: Reword IP1 and IP2 (`cairn/DESIGN.md:116-123`). Check IP1's caller claim against `git grep -n "compute_edges(" -- R` before you write it. (RB tripwire: ip-touching)
- [ ] T3: Rewrite the §5 scores-route text: the intro (`:219`), the §5.2 descriptor's `"EAP"` method value (`:247`), the §5.3 pseudocode comment (`:277`) and "Two situations" list (`:289-295`), and §5.4 (`:299-302`). Keep the algebra derivation unchanged.
- [ ] T4: Rewrite the §9 `edge_method` row (`:429`) as internal, the scores-method row's EAP clause (`:428`) to match D-007, and the Documentation standard line (`:449-450`).
- [ ] T5: Shorten the first Known limitations entry (`:588-599`) to the limitation and a pointer to IP2. Keep the "Corrected M91" mark.
- [ ] T6: Read the four AC3 sections in full and run AC3's two searches on the branch head. Fix any straggler, run `cairn_validate`, and record the sweep result in one work-log line.

## Work log

- 2026-10-01: created by /milestone-plan, promoted from the candidate row "DESIGN still describes `edge_method = "scores"` as a user choice" (added 2026-09-30, M91 review F5).
- 2026-10-01: criteria audit (full mode, ip-touching tag) returned 8 findings. Fixed before the gate: case-insensitive `EAP` matched "cheap", AC3 counted grep hits instead of named sections, AC3 bound a recording act, IP1 named only `ackwards()`, the false clauses were not required gone, and D-031 was called narrowed. The gate settled D-004's §9 source pointer.
- 2026-10-01: plan gate chose fixing the DESIGN text over exposing a user scores route on `ackwards()` because the route is new public API that no user asked for; falsified by a user request for sample-realized edges.
- 2026-10-01: plan gate chose keeping the §9 `edge_method` row marked internal over removing it because D-004's `_Source` cites that row; falsified by a reader taking the marked row as a user default.
- 2026-10-01: plan gate chose a `[low]` candidate row for sample-realized edges over a rejection in the D-entry, which keeps deferral a ROADMAP fact; falsified by a principled reason the package never offers it.
- 2026-10-01: implement started on branch m93-design-edge-method-reconcile. Code read: the six `compute_edges()` calls in `R/` pass "auto" or "algebra" with no data, and the PCA, EFA, and ESEM engines all set `linear = TRUE`.
- 2026-10-01: implement gate approved the IP1, IP2, and D-038 wording as shown, with no escalation. It left the two promoted candidate rows (M93, M94) for post-merge hygiene.
- 2026-10-01: T1 done. D-038 appended to `cairn/DECISIONS.md` as approved at the gate. Its heading names D-004 as narrowed, and its body cites D-007 and D-031.
- 2026-10-01: T2 done. IP1 and IP2 reworded in `cairn/DESIGN.md` to the gate-approved text. IP1's caller claim matches `git grep -n "compute_edges(" -- R` on the branch.

## Decisions

## Review
