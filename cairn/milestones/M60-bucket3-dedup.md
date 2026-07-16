<!-- Section ownership + write-modes: see tracking-rules.md "Milestone-file
     section ownership". A phase skill never rewrites another phase's section.
     Per-section owners are tagged below. -->
# M60: De-duplicate the setup path (audit bucket 3)

- **Status:** review   <!-- owner: transitioning skill · mirror-update; cairn/ROADMAP.md is the authority -->
- **Priority:** normal   <!-- owner: plan · create/amend-via-gate; high | normal | low -->
- **Depends on:** —   <!-- owner: plan · create/amend-via-gate -->
- **Principles touched:** —   <!-- no formal DESIGN IP/GP yet; operates under CLAUDE.md Invariants 1/2/7 (preserves, does not change them) -->
- **Branch/PR:** m60-bucket3-dedup · https://github.com/jmgirard/ackwards/pull/60   <!-- owner: implement (branch) / review (PR URL) · create -->

## Goal
<!-- owner: plan · create -->

Remove the remaining readability-only duplication and dead code in the setup path
(engines + `ackwards()` + comparability/boot_edges) that the 2026-07-12 audit flagged
as bucket 3, with zero behavior change.

## Scope
<!-- owner: plan · create/amend-via-gate -->

Third and final bucket of the 2026-07-12 codebase audit (M58 = validation helpers,
M59 = console/plot builders). Readability-only: no correctness or performance goal
(DESIGN §3 "measure before optimizing" — the setup-path recomputes are not hotspots),
so the eigen change is framed as removing a redundant computation for clarity, not as
an optimization. Behavior-preserving throughout; existing tests + the algebra-vs-scores
edge oracle are the safety net (Invariants 1/2). No exported behavior, no numeric result,
and no documented object contract changes — internal helpers and dead-code removal only.

**In:**
- **Shared per-engine helpers** in `utils.R`: `.variance_explained(L, p, labels)`
  (replaces the near-identical variance vectors at `engine_pca.R:44-49`,
  `engine_efa.R:114-119`, `engine_esem.R:223-228`, plus ESEM's second compute at
  `engine_esem.R:180-181`) and `.score_var(W, R)` (replaces `engine_pca.R:42`,
  `engine_efa.R:112`, `engine_esem.R:221`).
- **Drop dead edge builds:** `.align_signs()` stops returning (and mutating) its
  discarded `$edges` component (`utils.R:244,248`); `ackwards()`'s first
  `compute_edges()` call (`ackwards.R:832-838`) stops building the tidy tibble it
  throws away (matrices-only path / `build_tidy = FALSE` on `compute_edges()`).
  The second `compute_edges()` (`final_edges`, `:865`) stays — it covers `pairs="all"`
  skip-level edges `.align_signs` never produces.
- **Thread one smallest-eigenvalue:** each input path computes R's smallest eigenvalue
  at most once. Give `.near_singular_check` (`utils.R:527`) an optional precomputed
  `min_eig`; supply it from the polychoric (`ackwards.R:740/748`), FIML
  (`utils.R:372`), and cor-matrix (`utils.R:489`) checks that already computed it.
- **Remove dead code:** the always-true `k_eff >= 1L` guard (`ackwards.R:848` →
  `if (align_signs)`); the vestigial all-TRUE `$meta$converged_levels` vector
  (`ackwards.R:883`, zero readers); simplify `deepest_converged` (`:884`) from
  `max(which(conv))` to `k_eff` (identical value). `deepest_converged` itself stays —
  it is the "convergence summary" DESIGN §6 promises and is read by `glance()`.
- **Shared muffled fit dispatch:** extract `.fit_levels_muffled(R, engine, k_max, cor, fm)`
  — the `suppressMessages(suppressWarnings(switch(engine, pca=…, efa=…)))` core shared
  by comparability's `.fit_half` (`comparability.R:278-292`) and boot_edges'
  `.boot_replicate` (`boot_edges.R:314-331`). Callers keep their own R construction
  (pairwise vs FIML), resample/split step, and return/sentinel shape.

**Out:**
- Positive-manifold flip helper → **candidate** (the k=1 test is identical but each
  engine's follow-up diverges — PCA negates `fit$weights`, EFA reuses a `flip` flag in
  its fallback, ESEM must leave `L_se` alone — and it overlaps `.align_signs()`; a
  shared helper saves ~2 lines/engine while hiding the coupling).
- Shortfall-reporter helper for comparability/boot_edges → **candidate** (same two-bullet
  cli template, but divergent wording, nouns, and computed stats; would need every string
  parameterized for little gain).
- Any behavior/perf/output change; any new exported object; NEWS entry (internal-only).

## Acceptance criteria
<!-- owner: plan · create/amend-via-gate; review reads, never reinterprets -->

- [x] AC1: `.variance_explained()` and `.score_var()` exist in `utils.R` and are the
      sole computation sites — a repo-wide grep finds no inline `colSums(L^2)/p`-style
      variance vector or `diag(crossprod(W, … %*% W))` score-var outside the two helpers
      (LESSONS M58: grep the whole codebase, not just the three engines). All three
      engines call them; existing PCA/EFA/ESEM tests green with no value change.
- [x] AC2: `.align_signs()` no longer returns or mutates `$edges`; `ackwards()`'s first
      `compute_edges()` call builds no discarded tidy tibble. Existing edge/sign tests
      green **and** the algebra-vs-scores edge cross-check oracle still passes (Invariant 2).
- [x] AC3: Each input path (raw polychoric, FIML, cor-matrix, pearson/spearman) computes
      R's smallest eigenvalue at most once, and the near-singular / non-PD conditions fire
      identically to before — existing warning/error tests for those paths green, zero
      snapshot diff.
- [x] AC4: The `k_eff >= 1L` guard and `$meta$converged_levels` are gone;
      `$meta$deepest_converged` retains its exact prior value and `glance()` output is
      unchanged (existing glance tests in `test-efa.R`/`test-pca.R`/`test-esem.R` green).
- [x] AC5: `.fit_levels_muffled()` is the shared muffled-dispatch core called by both
      comparability's half-fit and boot_edges' replicate-fit; `comparability()` and
      `boot_edges()` outputs are unchanged (their existing tests + oracles green).
- [x] AC6: `Rscript tools/dod-gate.R` clean (check 0 err/0 warn/0 note, coverage no
      regression, style/lint/pkgdown), zero snapshot diffs, and no duplicated computation
      sites remain (grep-verified via AC1/AC2/AC5). No new exported object; no NEWS/doc
      change. *(Amended 2026-07-16 at the mini-gate: net-line-reduction clause dropped — see Decisions.)*

## Coverage
<!-- owner: plan · create/amend-via-gate -->

- AC1 → T1
- AC2 → T2
- AC3 → T3
- AC4 → T4
- AC5 → T5
- AC6 → T6

## Tasks
<!-- owner: plan (create) / implement (check-off, minor edits) -->

- [x] T1: Add `.variance_explained(L, p, labels)` and `.score_var(W, R)` to `utils.R`;
      rewire `engine_pca.R`, `engine_efa.R`, `engine_esem.R` (incl. ESEM's second variance
      compute at `:180-181`) to call them. Grep-verify no inline copies remain anywhere.
- [x] T2: Drop `$edges` from `.align_signs()` (return list + the `edges_list` mutation,
      `utils.R:244,248`); add a matrices-only path to `compute_edges()` and use it for the
      first call in `ackwards.R:832-838` so the discarded tidy tibble is never built.
- [x] T3: Add an optional precomputed `min_eig` param to `.near_singular_check`
      (`utils.R:527`); thread the smallest eigenvalue already computed on each input path
      (polychoric `ackwards.R:740/748`, FIML `utils.R:372`, cor-matrix `utils.R:489`) so it
      is computed once per path. Confirm the near-singular/non-PD messages are byte-identical.
- [x] T4: Remove the `k_eff >= 1L` guard (`ackwards.R:848`) and the `converged_levels`
      vector (`:883`); set `deepest_converged = k_eff` (`:884`). Verify `glance()` output
      and the `deepest_converged < 2` glance path (`test-pca.R:218-233`) are unchanged.
- [x] T5: Extract `.fit_levels_muffled(R, engine, k_max, cor, fm)` (the muffled
      `switch(engine, pca=…, efa=…)` dispatch); rewire `.fit_half` (`comparability.R`) and
      `.boot_replicate` (`boot_edges.R`) to call it, each keeping its own R construction,
      resample/split, and return/sentinel handling.
- [x] T6: Run `Rscript tools/dod-gate.R`; confirm 0 err/0 warn/0 note, no coverage
      regression, zero snapshot diffs, algebra-vs-scores oracle green, and no duplicated
      computation sites remain (grep). No NEWS/NAMESPACE change expected — confirm
      (internal-Rd regeneration for the `build_tidy` @param is expected and fine).

## Work log
<!-- owner: any skill · append-only; one line per entry; absolute dates -->

- 2026-07-13: created by /milestone-plan (promotes the bucket-3 dedup candidate; two weak extractions dropped to candidates per plan gate).
- 2026-07-16: implement started on branch m60-bucket3-dedup; T1 done — both helpers in utils.R, engines + compute_edges D^{-1/2} rewired (repo-wide grep clean), ESEM ordering key = raw colSums(L^2) so variance is computed once post-sort; 450 tests green.
- 2026-07-16: T2 done — .align_signs() returns loadings+signs only; compute_edges(build_tidy = FALSE) used by the lineage pass + .cross_cor()/.boot_replicate() (minor amendment, same dead-build pattern); M35 unit test reworked to assert $signs; 250 tests green.
- 2026-07-16: T3 done — .near_singular_check(min_eig); .corfiml_R()/.validate_cor_matrix() return list(R, min_eig) (callers take $R; FIML nocov branch recomputes post-smooth); polychoric threads its smoothed value; messages untouched; 688 tests green incl. snapshots.
- 2026-07-16: T4 done — dead guard, converged_levels, and conv vapply removed; deepest_converged = k_eff (Invariant 7; truncate-at-first-failure verified in all three engines); 476 tests green.
- 2026-07-16: T5 done — .fit_levels_muffled() extracted (signature gains n_obs, minor task edit); callers keep their own suppressed R construction so cor()/corFiml chatter stays muffled as before; 156 tests green.
- 2026-07-16: AC6 amended at a user-approved mini-gate (line-count clause dropped; see Decisions); comment-only density trim + compute_edges.Rd regen.
- 2026-07-16: T6 done — dod-gate PASSED (check 0 err/0 warn/0 note, coverage 100%, style/lint clean, pkgdown complete, zero snapshot diffs, oracles green); NAMESPACE/NEWS untouched; status → review.
- 2026-07-16: review hygiene — same-day work-log entries compressed to the mandated one-line form (cap FAIL) and references/INDEX.md forbes2023 line reformatted to the validator's filename-first format (pre-existing drift).

## Decisions
<!-- owner: implement / review · append-only; milestone-local -->

- 2026-07-16 (mini-gate, user-approved): AC6's "net line reduction vs pre-M60 `R/`" clause
  dropped — measured +39 total / +9 non-comment lines vs master; extracting *documented*
  helpers from 3–8-line duplicated blocks is line-neutral, and forcing a reduction would
  strip helper docs below repo comment density. Real outcome = single computation sites
  (grep-verified, AC1/AC2/AC5). Options offered: drop / ±15 band / strip comments / pause.

## Review
<!-- owner: review · exclusive -->

### Acceptance-criteria evidence (2026-07-16, fresh, by command)

Full parallel suite on the branch: **2295 pass / 0 fail / 0 warn / 0 skip** (devtools::test,
TESTTHAT_CPUS=8). Code identical to dod-gate commit 73d36e9 (`git diff 73d36e9..HEAD --
R/ man/ tests/ DESCRIPTION NAMESPACE` = 0 lines), so the gate run's results remain valid.

- **AC1** ✓ — `grep -rn 'colSums(.*\^2) */ *p\|diag(crossprod(' R/ | grep -v utils.R` → no
  hits; both helpers live in utils.R; all three engines + compute_edges call them; engine
  tests green in the suite run above with no value change (no snapshot/oracle diffs).
- **AC2** ✓ — grep finds no `edges_list[[key]] <-` mutation or `edges =` in `.align_signs()`'s
  return; `build_tidy = FALSE` on ackwards.R:843 (first compute_edges call); edge/sign tests
  and the algebra-vs-scores cross-check oracle (test-compute_edges/scores/esem) green.
- **AC3** ✓ — six `eigen(only.values)` sites, each ≤1 per input path (polychoric pre/post-
  smooth in ackwards.R:744/752; FIML utils.R:387/395; cor-matrix utils.R:521; pearson/
  spearman fallback utils.R:562 inside the check); warning/error tests green; zero snapshot
  diffs (dod-gate + suite).
- **AC4** ✓ — `grep 'k_eff >= 1L\|converged_levels\|max(which(conv))' R/` → no hits;
  `deepest_converged = k_eff`; glance tests in test-pca/efa/esem green (incl. the
  `deepest_converged < 2` path).
- **AC5** ✓ — `.fit_levels_muffled()` defined in utils.R:406, called by comparability.R:284
  and boot_edges.R:324 only; comparability + boot_edges tests/oracles green in the suite.
- **AC6** ✓ (as amended) — dod-gate PASSED: check 0 err/0 warn/0 note, coverage 100%,
  style/lint clean, pkgdown index complete; zero snapshot diffs; dedup greps clean (AC1/2/5
  above); no new export (NAMESPACE/NEWS/README 0-line diff vs master; only internal
  compute_edges.Rd regenerated).

### Consistency gate (2026-07-16)

- `cairn_validate.py` exit 0 after two hygiene fixes on the branch: M60 work-log entries
  compressed to one-line form (plan-owned body 172 → <150) and references/INDEX.md
  forbes2023 line reformatted filename-first (pre-existing drift; validator regex requires
  it). 83 dangling-ID advisories are pre-existing legacy references — non-gating.
- Coverage completeness: PASS (validator; AC1–6 → T1–6, all present).
- cairn_impact: skipped — no DESIGN principle changed.
- Toolchain (r-package slot): `document()` no diff; pkgdown "No problems found"; README/NEWS
  correctly untouched (internal-only); no new top-level files; full `check()` 0/0/0 via
  dod-gate on identical code.
