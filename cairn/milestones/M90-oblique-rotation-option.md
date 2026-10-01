# M90: Oblique rotation as a documented non-default option

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M89
- **Driving RR:** RR02
- **Principles touched:** IP1, IP2, IP4, IP6, IP8, IP9, GP1, GP4, GP5
- **Resolves:** —
- **Surface tier:** user-facing, because it adds a `rotation` argument on `ackwards()` and changes its outputs
- **Branch/PR:** m090-oblique-rotation-option

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

- [x] T1: The pre-implementation gate reads the design-session D-entry (BC5 quantity, BC6 stance) and records both in this file's Decisions. Unblock only then. Re-run M89's baseline test before any change.
- [x] T2: Add the `rotation` argument on `ackwards()` (`R/ackwards.R:311`), the per-engine validation table, and the object stamp (`R/ackwards.R:1054`). Add GPArotation to Suggests with `rlang::check_installed()` on the psych oblique paths. Write the cross-product test. Pass the object's rotation through `.fit_levels_muffled()` so that `boot_edges()` refits with it.
- [x] T3: PCA and EFA oblique. Pass `rotate` through (`R/engine_pca.R:28`, `R/engine_efa.R:16`). Write `.tenBerge_weights(R, L, Phi)` as `A = Σ^{-1/2} C Φ^{1/2}`, `C = Σ^{-1/2} L (L'Σ^{-1}L)^{-1/2}`, `L = ΛΦ^{1/2}` (tenberge1999 Eq. 9, Thm 1), with the `R^{-1}ΛΦ` fallback. Live oracle against `psych::factor.scores(method = "tenBerge")`. IP2 algebra-vs-scores oblique case.
- [x] T4: ESEM oblique. Pass lavaan `rotation` through (`R/engine_esem.R:74-79`). Weights through T3's helper with `cor.lv`. Regression fallback with Φ (`R/engine_esem.R:219-226`).
- [x] T5: `.variance_explained(L, p, labels, Phi)` per AC10. ESEM sort key (`R/engine_esem.R:200`). psych `Vaccounted` oracle test.
- [x] T6: Apply the D-entry to `match_parents()` (`R/utils.R:230`), `.align_signs()` (`R/utils.R:275`), and `ba_layout()` and `autoplot()` labels. Snapshot the enumerated surfaces on an oblique fit.
- [ ] T7: Fit-time cli advisory (IP6) and the `prune()` stance per the D-entry. Tests assert the message text and the warn-or-implement branch.
- [ ] T8: `data-raw/forbes2023-oblique.R` (md5-pinned OSF `7jfkw` functions, oblique branch). Reconcile her unstandardized `comp.corr` with standardized edges through `D`. Fixture plus provenance, `test-forbes-fidelity.R` oblique block, ORACLES rows.
- [ ] T9: Docs per AC12 (roxygen, `ackwards-engines.Rmd.orig` plus precompute, NEWS, DESIGN §9 row). Run `devtools::document()` and `Rscript tools/dod-gate.R`.

## Work log

- 2026-09-16: created by /milestone-plan and set `blocked` at creation. Blocker: no design-session D-entry fixes BC5's lineage quantity and BC6's redundancy stance (the D-034 gate). Collision sweep: absorbs the candidate row "Oblique rotation support (D-034; implementation gated)". The row stays until completion. RR02 BC1 to BC8 ingested verbatim as AC1 to AC8.
- 2026-09-16: criteria audit ran in full mode (fresh-context [O] reader, shared run with M89). M90 fixes: BC4's unnamed formula pinned to psych's `diag(Φ Λ'Λ)/p` convention with a live oracle (AC10). GPArotation Suggests made an explicit requirement (AC9, D-035). geomin pinned to its oblique variant. Error-branch test widened to the full cross-product. DESIGN §9 row added to the doc criterion. Fixture provenance made explicit (AC11). BC1, BC2, BC3, and BC5 universals narrowed in the Deviations table. Coverage stance fixed as packages required in dev and CI (AC13).
- 2026-09-16: plan gate chose psych's `diag(Φ Λ'Λ)/p` over the structure-matrix sum of squares for oblique variance because it is the convention the wrapped engine already reports (GP4, live oracle). Falsified by the design session naming a different published convention.
- 2026-09-16: cairn_validate's sizing tripwire fires (13 criteria). It stands: eight are RR02's binding criteria, which one milestone must carry verbatim, and the other five pin what the audit found unverifiable in them. The oblique semantics were split off instead: M89 carries the Φ plumbing and partialled reporting.
- 2026-09-16: plan gate chose GPArotation as a guarded Suggests over an ESEM-only oblique option because Forbes's reference implementation is oblique PCA and EFA, the oracle's own engine. Falsified by GPArotation leaving CRAN or breaking psych's oblique paths.
- 2026-09-16: carried from the M89 review (four [O] findings routed here because they only bite once a real non-identity Phi exists): (1) the PCA/EFA carry assumes psych sorts `$Phi` with the loadings and no M89 test can reach a non-identity psych `$Phi`; M90's oblique fixtures must assert `factor_cor` matches the stored loadings' order and signs on all three engines. (2) The routing test cannot detect a wrong `ord`; the same fixtures cover it. (3) Docs do not distinguish the rotation correlation (`factor_cor`, shown by `tidy(what = "factor_cor")`) from the score correlation Phi_s behind `beta` and `r2`; name both once oblique makes them differ. (4) PCA/EFA store `factor_cor` unnamed while ESEM labels it; label all three here, where O13's attribute-for-attribute baseline is superseded by the oblique fixtures.
- 2026-09-30: /milestone-implement started. Owner override: the owner settled BC5 and BC6 without the Forbes design session that D-034 gated on, recorded as D-036 (`r` drives lineage, `prune()` warns). Status blocked to in-progress. Branch m090-oblique-rotation-option cut from master e0d227c.
- 2026-09-30: implement gate. The owner allowed the md5-pinned download of Forbes's OSF `7jfkw` script for T8. The GPArotation question rested on a wrong premise (that promax runs without it). The owner picked "oblimin only", then a check showed `psych::kaiser()` stops without GPArotation, so both rotations are guarded per D-035 as written.
- 2026-09-30: T1 done. D-036 recorded, both stances copied to Decisions, M89 baseline test re-run before any code change (5 fits, 95 expectations, 0 failures).
- 2026-09-30: minor amendment. T2 gains the `boot_edges()` pass-through, which Scope's "inheriting the fixed engines" needs. Without it, the bootstrap of an oblique object refits varimax.
- 2026-09-30: T2 done. `rotation` is the last named argument of `ackwards()`, after `correct`, so no positional call shifts. Added `.check_rotation()` and `.supported_rotations` (`R/utils.R`), the object stamp, GPArotation in Suggests, and the `boot_edges()` pass-through. `test-rotation.R` covers the cross-product (9 fits, 3 aborts) and the arg_match errors. It also covers the GPArotation routing, the print and summary header, and the boot routing, each routing test with a control. Suite: 2915 expectations, the 3 failures were one boot test whose `pca_levels` mock lacked `rotation`, fixed.
- 2026-09-30: discovered sub-task. The installed roxygen 8.1.0 writes grouped `importFrom()` blocks. Base R's `parseNamespaceFile()` then reads the backticked `%||%` with literal backticks, so the import is lost. The three `@importFrom rlang` tags now quote it. `Config/roxygen2/version` moved to 8.1.0.
- 2026-09-30: T3 done. PCA and EFA pass `rotation` to psych. `.tenBerge_weights(R, L, Phi)` takes a required Φ, and a Φ within 1e-12 of I runs the old orthogonal formula unchanged. EFA reads one Φ for its weights and `factor_cor`. PCA keeps psych's weights, which are R⁻¹ΛΦ under oblique. PCA and EFA `factor_cor` now carry level labels (M89 carried finding 4), so the baseline test compares its values only. ESEM passes the identity until T4.
- 2026-09-30: T3 tests in `test-oblique.R` cover the Φ = I formula and psych tenBerge oracles (1e-12, 1e-8) and the not-positive-definite error. They also check stored Φ against psych `$Phi` with order and signs, scores reproducing Φ, and oblique algebra-vs-scores on all pairs. Planting Φ-blind weights turned 4 tests red. Related files: 650 expectations, 0 failures.
- 2026-09-30: T4 done. ESEM passes `rotation` to lavaan. `cor.lv` is read by lavaan factor name before the sort and permuted with the loadings. It feeds the tenBerge weights, the R⁻¹ΛΦ fallback, and `factor_cor`. An oblique level whose Φ cannot be read now truncates with an error instead of storing the identity.
- 2026-09-30: T4 tests check stored Φ against `cor.lv` with order and signs (oblimin, geomin), scores reproducing Φ, and algebra-vs-scores. Dropping the permutation turned both Φ tests red, so M89's carried findings 1 and 2 are covered on all three engines. Baseline gap on this machine: 0 for PCA and EFA, 1.6e-16 for ESEM.
- 2026-09-30: T5 done. `.variance_key()` returns diag(ΦΛ'Λ). For a Φ within 1e-12 of I it returns colSums(Λ²), through `.near_identity()`, which the tenBerge helper shares. `.variance_explained()` takes a required Φ and divides that key by p. The ESEM sort orders by the same key in lavaan's order.
- 2026-09-30: T5 tests: PCA and EFA variance against psych's `Vaccounted` proportion row (1e-8, with a guard that the old formula misses by more than 1e-5). The ESEM sort test uses bfi25 with geomin at k = 5 and seed 1. There the oblique key orders the factors 2 4 1 5 3 and the squared loadings order them 2 4 1 3 5. Planting the old key turned all three red. Full suite: 3032 expectations, 0 failures.
- 2026-09-30: T6 done. Per D-036, `match_parents()` and `.align_signs()` keep r, with comments that cite D-036. Under oblique, print() and summary() add a note that r is a total correlation and beta the partialled coefficient, and autoplot() adds a caption. `?tidy.ackwards` and `?summary.ackwards` describe the oblique case and name `factor_cor` against the score correlation behind beta (M89 carried finding 3). `ba_layout()` is unchanged.
- 2026-09-30: T6 tests: snapshots of print, summary, and the three tidy tables on a PCA oblimin fit. Assertions check that `is_primary` marks the largest |r| with r > 0, that beta differs from r, and that the notes and caption appear. A varimax control shows none of them.

## Decisions

- 2026-09-30 (T1): BC5 and BC6 read from D-036. The marginal `r` drives `match_parents()`, `.align_signs()`, and the diagram under every rotation. `prune()` keeps its chase on `r` and warns on an oblique object.
- 2026-09-30: `rlang::check_installed("GPArotation")` guards both psych oblique rotations, `oblimin` and `promax`. psych's promax runs through `psych::kaiser()`, which stops without GPArotation. This follows D-035 as written.
- 2026-09-30: PCA and EFA keep the column order that psych returns. Under oblique, `psych::fa()` sorts by `diag(ΦΛ'Λ)`, the AC10 key, and `psych::principal()` sorts by `diag(Λ'Λ)`. A re-sort breaks column correspondence with psych and with Forbes's reference output. So a PCA level's `variance` vector can fall out of descending order in rare cases. ESEM sorts by the AC10 key.

## Review
