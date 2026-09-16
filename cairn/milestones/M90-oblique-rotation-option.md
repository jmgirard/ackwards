# M90: Oblique rotation as a documented non-default option

- **Status:** blocked
- **Priority:** normal
- **Depends on:** M89
- **Driving RR:** RR02
- **Principles touched:** IP1, IP2, IP4, IP6, IP8, IP9, GP1, GP4, GP5
- **Resolves:** —
- **Surface tier:** user-facing, because it adds a `rotation` argument on `ackwards()` and changes its outputs
- **Branch/PR:** —

## Goal

Offer oblique rotation as a documented, non-default `rotation` argument on every engine, with real within-level factor correlations, correlation-preserving oblique scoring, oblique-valid variance, a loud fit-time advisory, and a Forbes oracle, while default varimax output stays unchanged.

## Scope

**In:** the `rotation` argument and its per-engine validation. GPArotation as a guarded Suggests (D-035). Oblique ten Berge weights and the oblique regression fallback. The oblique variance formula and the ESEM sort key. Applying the design-session decision to primary-parent matching, sign anchoring, and the diagram. The cli advisory and the `prune()` stance. The Forbes `ExtendedBassAckwards` oblique fidelity fixture. Docs, vignette, DESIGN §9 row, NEWS.

**Out:** the Φ carry helper, `tidy(what = "factor_cor")`, and the `beta` and `r2` columns belong to M89, and this milestone consumes them. Tucker's φ as the documented rotation-consistency diagnostic stays a candidate row. The D-032 gap-tolerant chase question and the Figure 6B retention rule stay their own candidate rows, tied to the same design session. Oblique `boot_edges()` and `comparability()` semantics beyond inheriting the fixed engines become a candidate row (RR02 Q5 item 12).

**Blocker:** BC5 requires a design-session D-entry that fixes which quantity (marginal `E` or Φ-partialled `B`) drives primary-parent matching, sign anchoring, and the diagram, and BC6 requires the redundancy stance. No such entry exists as of 2026-09-16. `/milestone-implement` must not start until it does.

## Acceptance criteria

- [ ] AC1 (BC1): Varimax remains the default rotation for every engine; a fit with all defaults is byte-identical in its numerical outputs to the pre-change package, and the IP9 Forbes-reproduction settings and fidelity suite pass unchanged.
- [ ] AC2 (BC2): Under an oblique rotation, each level's stored `factor_cor` is the real within-level factor correlation — extracted from the engine, permuted by any column reordering (the `engine_esem.R` `ord` guard), and sign-flipped in lockstep with `align_signs` — and is surfaced in at least `summary()` and one tidier. No code path carries `factor_cor = I` for an oblique fit.
- [ ] AC3 (BC3): Oblique scoring weights are the correlation-preserving oblique variant (ten Berge et al. 1999) on the default path, with the oblique regression rule (`R^{-1}ΛΦ`) as the fallback; the orthogonal tenBerge formula is never applied to an oblique pattern matrix. The IP2 algebra-vs-scores oracle gains at least one oblique test case.
- [ ] AC4 (BC4): `variance` (and everything downstream of `.variance_explained()`) is computed by an oblique-valid formula when Φ ≠ I; the ESEM variance sort and the Φ permutation use the same ordering.
- [ ] AC5 (BC5): Output distinguishes marginal edges (`E`) from Φ-partialled lineage quantities (`B` or semipartials); no output surface labels an oblique marginal edge as a unique/lineage contribution. Which quantity drives primary-parent matching, sign anchoring, and the diagram is fixed by a design-session D-entry before implementation, not defaulted in code.
- [ ] AC6 (BC6): An oblique fit announces via cli (IP6) that edges are total correlations and that `redundancy_r`/`cut_show` conventions were calibrated under varimax; `prune()` under an oblique object either implements the design-session redundancy stance or warns that the default criterion assumes orthogonal levels.
- [ ] AC7 (BC7): The oblique path is oracle-backed (IP8) including a fidelity test against Forbes's `ExtendedBassAckwards` oblique branch with the standardization correspondence (`D ≠ I`) handled explicitly.
- [ ] AC8 (BC8): No new Imports: any rotation dependency (e.g. GPArotation for psych's oblique criteria) enters as a guarded Suggests only (D-011/GP5).
- [ ] AC9: `ackwards()` gains `rotation = "varimax"`. It also accepts `"oblimin"` and `"promax"` for `pca` and `efa` (psych, with GPArotation guarded by `rlang::check_installed()`), and `"oblimin"` and `"geomin"` (lavaan's oblique geomin) for `esem`. A test enumerates the full engine by rotation cross-product and asserts that each supported pair fits and each unsupported pair aborts with a cli error naming both. The chosen rotation is stored in the object's `rotation` field and shown by `print()` and `summary()`.
- [ ] AC10: `.variance_explained()` takes Φ and returns `diag(Φ Λ'Λ) / p`, which is psych's `Vaccounted` convention for oblique solutions and reduces to `colSums(Λ^2)/p` at Φ = I. A live-oracle test asserts equality with `psych::fa()`'s `Vaccounted` row on an oblimin fit within 1e-8. The ESEM sort key at `R/engine_esem.R:200` uses the same quantity, asserted by a test on a level where the oblique and orthogonal keys order differently.
- [ ] AC11: The oblique Forbes expected values ship as `tests/testthat/fixtures/forbes2023_oblique.rds` with a `provenance` attribute that names its md5-pinned `data-raw/` generator (OSF `7jfkw`, M53 pattern), registered in `cairn/ORACLES.md`.
- [ ] AC12: The `rotation` bullet in `?ackwards`, the `ackwards-engines` vignette (RR02 Q7 passage with the option-form ◆ sentences), NEWS.md, and DESIGN §9's rotation row (no longer "not a user argument") document the option. `Rscript vignettes/precompute.R` and `Rscript tools/check-vignette-freshness.R` pass.
- [ ] AC13: `Rscript tools/dod-gate.R` passes with GPArotation and lavaan installed. Both are required in dev and CI, with no coverage exemption on oblique branches.

Deviations from RR02 (narrowed readings, not softenings: the fresh-context criteria audit of 2026-09-16 found each universal unenumerable as written):

| BC | Reading this milestone verifies |
|---|---|
| BC1 "byte-identical" | Equal within 1e-12 to M89's master-generated baseline fixture (`baseline-m89.rds`, all three engines) plus the unchanged fidelity suite. Bitwise equality does not survive the refactor BC2 requires. |
| BC2 "No code path carries `factor_cor = I`" | Each engine's oblique fixture asserts that its stored `factor_cor` equals the engine's own Φ (psych `$Phi`, lavaan `cor.lv`) permuted and flipped, and is not the identity. |
| BC3 "never applied to an oblique pattern matrix" | `.tenBerge_weights()` gains a required `Phi` argument, so every caller passes Φ. A test asserts that the Φ = I call reproduces the current weights and that an oblique call reproduces `psych::factor.scores(method = "tenBerge")` weights within 1e-8. |
| BC5 "no output surface labels" | The enumerated surfaces `print()`, `summary()`, `tidy()` (`edges`, `variance`, `factor_cor`), and `autoplot()` edge labels are snapshotted on an oblique fit. They label `r` as a total correlation and `beta` as the partialled coefficient. |

## Coverage

- AC1 → T1, T9
- AC2 → T2, T3, T4
- AC3 → T3, T4
- AC4 → T5
- AC5 → T6
- AC6 → T7
- AC7 → T8
- AC8 → T2
- AC9 → T2
- AC10 → T5
- AC11 → T8
- AC12 → T9
- AC13 → T9

## Tasks

- [ ] T1: The pre-implementation gate reads the design-session D-entry (BC5 quantity, BC6 stance) and records both in this file's Decisions. Unblock only then. Re-run M89's baseline test before any change.
- [ ] T2: Add the `rotation` argument on `ackwards()` (`R/ackwards.R:311`), the per-engine validation table, and the object stamp (`R/ackwards.R:1054`). Add GPArotation to Suggests with `rlang::check_installed()` on the psych oblique paths. Write the cross-product test.
- [ ] T3: PCA and EFA oblique. Pass `rotate` through (`R/engine_pca.R:28`, `R/engine_efa.R:16`). Write `.tenBerge_weights(R, L, Phi)` as `A = Σ^{-1/2} C Φ^{1/2}`, `C = Σ^{-1/2} L (L'Σ^{-1}L)^{-1/2}`, `L = ΛΦ^{1/2}` (tenberge1999 Eq. 9, Thm 1), with the `R^{-1}ΛΦ` fallback. Live oracle against `psych::factor.scores(method = "tenBerge")`. IP2 algebra-vs-scores oblique case.
- [ ] T4: ESEM oblique. Pass lavaan `rotation` through (`R/engine_esem.R:74-79`). Weights through T3's helper with `cor.lv`. Regression fallback with Φ (`R/engine_esem.R:219-226`).
- [ ] T5: `.variance_explained(L, p, labels, Phi)` per AC10. ESEM sort key (`R/engine_esem.R:200`). psych `Vaccounted` oracle test.
- [ ] T6: Apply the D-entry to `match_parents()` (`R/utils.R:230`), `.align_signs()` (`R/utils.R:275`), and `ba_layout()` and `autoplot()` labels. Snapshot the enumerated surfaces on an oblique fit.
- [ ] T7: Fit-time cli advisory (IP6) and the `prune()` stance per the D-entry. Tests assert the message text and the warn-or-implement branch.
- [ ] T8: `data-raw/forbes2023-oblique.R` (md5-pinned OSF `7jfkw` functions, oblique branch). Reconcile her unstandardized `comp.corr` with standardized edges through `D`. Fixture plus provenance, `test-forbes-fidelity.R` oblique block, ORACLES rows.
- [ ] T9: Docs per AC12 (roxygen, `ackwards-engines.Rmd.orig` plus precompute, NEWS, DESIGN §9 row). Run `devtools::document()` and `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan and set `blocked` at creation. Blocker: no design-session D-entry fixes BC5's lineage quantity and BC6's redundancy stance (the D-034 gate). Collision sweep: absorbs the candidate row "Oblique rotation support (D-034; implementation gated)". The row stays until completion. RR02 BC1 to BC8 ingested verbatim as AC1 to AC8.
- 2026-09-16: criteria audit ran in full mode (fresh-context [O] reader, shared run with M89). M90 fixes: BC4's unnamed formula pinned to psych's `diag(Φ Λ'Λ)/p` convention with a live oracle (AC10). GPArotation Suggests made an explicit requirement (AC9, D-035). geomin pinned to its oblique variant. Error-branch test widened to the full cross-product. DESIGN §9 row added to the doc criterion. Fixture provenance made explicit (AC11). BC1, BC2, BC3, and BC5 universals narrowed in the Deviations table. Coverage stance fixed as packages required in dev and CI (AC13).
- 2026-09-16: plan gate chose psych's `diag(Φ Λ'Λ)/p` over the structure-matrix sum of squares for oblique variance because it is the convention the wrapped engine already reports (GP4, live oracle). Falsified by the design session naming a different published convention.
- 2026-09-16: cairn_validate's sizing tripwire fires (13 criteria). It stands: eight are RR02's binding criteria, which one milestone must carry verbatim, and the other five pin what the audit found unverifiable in them. The oblique semantics were split off instead: M89 carries the Φ plumbing and partialled reporting.
- 2026-09-16: plan gate chose GPArotation as a guarded Suggests over an ESEM-only oblique option because Forbes's reference implementation is oblique PCA and EFA, the oracle's own engine. Falsified by GPArotation leaving CRAN or breaking psych's oblique paths.

## Decisions

## Review
