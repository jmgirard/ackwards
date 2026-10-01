# M93: Reconcile DESIGN's edge_method text with the code

- **Status:** review
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

- [x] AC1: In `cairn/DESIGN.md`, IP1 states that every call of `compute_edges()` in `R/` passes `edge_method = "auto"` or `"algebra"` with no data, so every shipped edge comes from the `W'RW` algebra. IP1 no longer says "or when the user asks", and its claim matches what `git grep -n "compute_edges(" -- R` lists on the branch head. IP2 states that the scores route stays inside the internal `compute_edges()` as the second route of the algebra-vs-scores agreement tests. IP2 keeps its requirement that a standing test asserts agreement for every linear engine, and it no longer presents `edge_method = "scores"` as available to users. Both bullets keep their numbers.
- [x] AC2: `cairn/DECISIONS.md` gains one entry that records the IP1 and IP2 rewording with its rationale and the evidence that would reverse it. Its heading names D-004 as the entry it narrows, and its body cites D-007 and D-031.
- [x] AC3: Two parts of `cairn/DESIGN.md` are checked. The first part is the "Design principles" section, §5, §9, and the "Known limitations" section, each read in full. The second part is every other line that `grep -n -i -E 'edge_method|scores.? route|scores.?path|materiali|user asks|user wants|user would prefer' cairn/DESIGN.md` or `grep -n -E '\bEAP\b' cairn/DESIGN.md` returns. No text in either part presents the scores route, `edge_method`, or EAP scoring as a setting of an exported function.
- [x] AC4: The §9 `edge_method` row says that `edge_method` belongs to the internal `compute_edges()` and is not an argument of any exported function. The §9 scores-method row states EAP as out of scope (D-007), not as an opt-in.
- [x] AC5: The first "Known limitations" entry states the limitation: the algebra-vs-scores cross-check does not cover the polychoric or the `missing = "fiml"` PCA/EFA paths, and why. It points to IP2 for how edges are built, keeps its "Corrected M91" mark, and is shorter than the same entry on `master` (`wc -c` of each extracted entry).
- [x] AC6: The milestone's diff against `master` (`git diff --stat master...m93-design-edge-method-reconcile`) changes files under `cairn/` only. `python3 /Users/jmgirard/.claude/skills/cairn/scripts/cairn_validate.py` passes on the branch head.

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
- [x] T3: Rewrite the §5 scores-route text: the intro (`:219`), the §5.2 descriptor's `"EAP"` method value (`:247`), the §5.3 pseudocode comment (`:277`) and "Two situations" list (`:289-295`), and §5.4 (`:299-302`). Keep the algebra derivation unchanged.
- [x] T4: Rewrite the §9 `edge_method` row (`:429`) as internal, the scores-method row's EAP clause (`:428`) to match D-007, and the Documentation standard line (`:449-450`).
- [x] T5: Shorten the first Known limitations entry (`:588-599`) to the limitation and a pointer to IP2. Keep the "Corrected M91" mark.
- [x] T6: Read the four AC3 sections in full and run AC3's two searches on the branch head. Fix any straggler, run `cairn_validate`, and record the sweep result in one work-log line.

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
- 2026-10-01: T3 done. §5 intro, §5.2 method values (read from the engines: components, tenBerge, regression), §5.3 comment, the "Two situations" list (now "Where the `scores` route runs"), and §5.4 rewritten. The ESEM polychoric-weights claim was checked against `R/engine_esem.R`. The algebra derivation is unchanged.
- 2026-10-01: T4 done. The §9 `edge_method` row is marked internal and says it is not an argument of any exported function. The scores-method row states EAP as out of scope (D-007). The Documentation standard line no longer names `edge_method = "scores"`.
- 2026-10-01: T5 done. The first Known limitations entry states the uncovered polychoric and FIML paths and why. It points to IP1 and IP2 and keeps the "Corrected M91" mark. It went from 1102 to 686 bytes. Both figures are `wc -c` of the extracted entry, on `master` and on the branch.
- 2026-10-01: T6 done. The principles section, §5, §9, and Known limitations were read in full. Every hit of AC3's two searches falls inside those four sections, and no straggler needed a fix. `cairn_validate` passes, with 16 advisory warnings, all in M84's work log. The diff against `master` touches `cairn/` only.
- 2026-10-01: claim audit: not owed — internal tier
- 2026-10-01: all tasks done. Status set to review.
- 2026-10-01: review in progress. AC1 to AC6 verified and ticked against Review evidence. The full package check and the diff reviewer are still running.

## Decisions

## Review

Review run 2026-10-01 on branch head 3899d8c. The branch already contains `origin/master` (008a0b2), so no sync merge was needed. No PR exists yet.

- AC1 evidence: `cairn/DESIGN.md:116-125`. IP1 states that every `compute_edges()` call in `R/` passes `"auto"` or `"algebra"` with no data and that every shipped edge comes from `W'RW`. The words "or when the user asks" are gone. `git grep -n "compute_edges(" -- R` lists six calls (`ackwards.R:948`, `:1031`, `boot_edges.R:374`, `comparability.R:336`, `layout.R:162`, `prune.R:942`). Each was read: four pass `"auto"`, two pass `"algebra"`, and none passes data. The PCA, EFA, and ESEM engines set `linear = TRUE`. IP2 places the scores route inside the internal `compute_edges()` as the second route of the agreement tests. It says that no exported function offers the route, and it keeps the standing-test requirement. Both bullets keep the numbers IP1 and IP2.
- AC2 evidence: `git diff master...HEAD -- cairn/DECISIONS.md` adds one `### D-` heading, D-038 at `cairn/DECISIONS.md:274`. Its heading reads "narrows D-004". The Decision paragraph records the IP1 and IP2 rewording, and the Context paragraph gives the rationale. The Consequences paragraph names the reversing evidence, a user need for sample-realized edges. The body cites D-007 in Context, D-031 in Decision, and both again in `_Source`.
- AC3 evidence: the four sections were read in full on the branch head. They are "Design principles" (`cairn/DESIGN.md:107-179`), §5 (`:218-309`), §9 with its Documentation standard (`:423-460`), and "Known limitations" (`:592-625`). The first search returns 16 lines and the `\bEAP\b` search returns 2. All 18 hits fall inside those four sections, so no line outside them needed a separate read. No text presents the scores route, `edge_method`, or EAP as a setting of an exported function. The §5.3 pseudocode marks the scores branch as "agreement tests only". The §9 heading reads "users will not override these", but that clause predates M93 and names no scores route.
- AC4 evidence: the §9 `edge_method` row (`cairn/DESIGN.md:435`) is labeled "internal: `compute_edges()` only". Its rationale opens "Not an argument of any exported function" and says that `edge_method` belongs to the internal `compute_edges()`. The scores-method row (`:434`) reads "EAP out of scope (D-007)" in its Default cell. Its rationale says "EAP is out of scope, not an option (D-007, declined M28)".
- AC5 evidence: the first "Known limitations" entry (`cairn/DESIGN.md:594-601`) names the uncovered `cor = "polychoric"` and `missing = "fiml"` PCA/EFA paths. It gives the reason: the algebra uses the polychoric or `psych::corFiml()` matrix, but the scores route standardizes the raw data. It points to IP1 and IP2 for how edges are built, and it keeps the "Corrected M91" mark. An awk script extracted the entry, from its `- ` line to the next `- ` line. The entry measures 686 bytes on the branch and 1102 bytes on `master` (`wc -c`).
- AC6 evidence: `git diff --stat master...m93-design-edge-method-reconcile` lists four files, all under `cairn/`: DECISIONS.md, DESIGN.md, ROADMAP.md, and this milestone file. A filter for paths outside `cairn/` returns none. `cairn_validate.py` exits 0 with all checks passed and 16 advisory warnings, all in M84's work log.

Consistency gate:

- `cairn_validate.py`: exit 0, as recorded under AC6.
- `cairn_impact.py --changed` lists 28 IP1 and 26 IP2 references. The DESIGN.md and D-038 references are M93's own text. The other references are the D-031 entry, the M79 and M72 archives, and the M84 and M94 files. They agree with the new wording, because IP1's single edge path through `compute_edges()` is unchanged. `ROADMAP.md:25` is the candidate row that M93 was promoted from, which the post-merge hygiene pass removes.
- `devtools::document()` produced no diff. `pkgdown::check_pkgdown()` found no problems.
- The diff touches no NEWS.md, README, `.Rbuildignore`, `R/`, `man/`, NAMESPACE, or vignette file, so no NEWS entry is owed.
- `devtools::check()` with `TESTTHAT_CPUS=8`: Status OK, with 0 errors, 0 warnings, and 0 notes.

Independent review: internal tier with a diff under `cairn/` only, so one fresh-context Opus diff reviewer ran. It found no false claim about the code and no text that presents the scores route, `edge_method`, or EAP as a user setting. It reported 11 findings, ranked. The disposition after each one is the proposal put to the maintainer at the merge gate.

- F1: the rewrite dropped a true fact from `DESIGN.md:594-601`. The old entry said that `beta` and `r2` build Φ_s from the stored weights and the fit's R. So they do not share the algebra-vs-scores split. Proposed: fix now, one sentence restored.
- F2: the "Corrected M91" note corrects a claim about `beta` that the entry no longer mentions. Proposed: fix now, together with F1.
- F3: the uncovered-path list omits `cor = "spearman"`, ESEM FIML, and pairwise-missing Pearson data. The scores branch at `R/compute_edges.R` correlates standardized raw data with Pearson `cor()`, so a Spearman R differs in basis (read on the branch). The gap predates M93. Proposed: follow-up candidate row.
- F4: §5.3 says "Under missing data the two differ", but they also differ on complete data under a polychoric or Spearman R (`DESIGN.md:297-298`). M93 wrote this sentence. Proposed: fix now.
- F5: the §5.3 pseudocode signature is stale. It shows `align` and `use = "pairwise"`, omits `cut_show` and `build_tidy`, aligns signs inside the function, and does not show the `"algebra"` abort. This predates M93. Proposed: follow-up candidate row, grouped with F8 and F9.
- F6: IP1 is written as a census of current callers rather than as a rule. Proposed: reject, because AC1 requires that form and the implement gate approved the wording.
- F7: D-038's Context names EAP as the reason the scores route was kept, but the old §5.3 also kept it for a user-requested route. Proposed: fix now, one clause added.
- F8: §5.1 calls the EFA and ESEM paths "regression-scored", but both are tenBerge-scored with regression as the fallback. This predates M93. Proposed: follow-up, grouped with F5.
- F9: the §5.2 comments "NULL if !linear" describe a case that no engine produces. Proposed: follow-up, grouped with F5.
- F10: the ROADMAP hygiene stamp calls the two rows "promoted" while they stay as candidates, and row 26 points to row 25. Proposed: reject, because the post-merge hygiene pass removes both promoted rows and replaces the stamp. The M94 file does not cite row 25.
- F11: the internal `compute_edges()` roxygen still describes the scores triggers. Proposed: reject, because the plan put `R/` and `man/` out of scope and the roxygen describes what the function itself accepts.
