<!-- Entombed by M64 (2026-07-16): the embedded decision log formerly at DESIGN.md §14
     (lines 461–876 at extraction; source commit 0ce1095ae067a5cc39111df8a5dd00a623f0c23b), preserved VERBATIM below.
     Historical "§14.x" citations anywhere in this repo (DECISIONS.md D-001–D-015 source
     lines, milestone archives, legacy files) resolve against this file. Still-governing
     cross-cutting decisions are re-recorded in cairn/DECISIONS.md (D-001 onward); live
     known limitations moved to DESIGN.md "Known limitations". Do not edit. -->

## 14. Decisions resolved & remaining

**Resolved (design round):**
1. Rotation → **varimax (orthogonal) across all engines, hardcoded internal constant since M13**;
   oblique rotation is **out of scope** (resolved 2026-06 — it confounds the cross-level signal that
   is the method's whole point). `rotation` removed as a user argument in M13.
2. Correlation basis → **`pearson` default, `polychoric` opt-in**, with an ordinal-detection cli
   **warning** (no silent switching).
3. Within-level order → **primary-parent (recursive), variance tiebreak**; `f{j}` IDs follow it.
4. Ordinal ESEM scoring → **tenBerge on polychoric (linear → algebra)** as default; **EAP opt-in**.
   Default Likert path does **not** need `data` at edge time.
5. Package name → **`ackwards`** (verify via `available::available("ackwards")`).

**Resolved during build (M1–M3):**
6. Ordinal-detection heuristic: **≤ 7 distinct integer values** triggers the warning.
7. CF `kappa`: **κ = 1/p** (varimax-equivalent; Crawford & Ferguson 1970). The `kappa` argument
   was accepted and stored but never wired to any engine — all engines hardcoded `cfT → "varimax"`.
   Removed entirely in M13: CF(κ=1/p) ≡ varimax; the literature never varies kappa; exposing it
   implied quartimax/equamax were reasonable alternatives (they are not for this method).
8. `cor = "spearman"` added alongside `"pearson"` as a non-polychoric rank-based option.
9. `fm` argument exposed for EFA engine (`"minres"` default, `"ml"`, `"pa"`).
10. `cut_show` argument exposed (default 0.3) to control which edges are flagged `above_cut`.

**Resolved for M4:**
11. ESEM engine → **`lavaan::efa()`** (EFA-in-SEM, not full ESEM with structural paths). Gives
    per-level fit indices + rotation-aware SEs + WLSMV ordinal estimation. Does **not** require
    full ESEM block syntax; the bass-ackwards levels are independent EFA solutions.
12. ESEM scoring → self-compute tenBerge weights from lavaan's estimated loadings `Λ` and latent
    correlation matrix `R` via the existing `.tenBerge_weights(R, Λ)`. `lavPredict()` is not used
    for the default path (it lacks tenBerge); it is the eventual hook for EAP opt-in.
13. EAP scoring → **out of scope (declined M28).** EAP's shrinkage attenuates cross-level
    correlations — the primary signal bass-ackwards measures — making it theoretically inferior to
    tenBerge for this method. An EAP request still returns `cli_abort("not yet implemented")`.
    The scores-route seam in `compute_edges()` is preserved, but implementing EAP is not planned.
14. `cor = "polychoric"` → **general basis** for all engines (not ESEM-only). PCA/EFA compute `R`
    via `psych::polychoric()` in `ackwards()` then feed it to the engine as usual. ESEM uses
    lavaan's own latent correlation matrix as `R` for edges (do not mix psych-polychoric `R` with
    lavaan-WLSMV loadings).
15. ESEM ordinal estimator → **WLSMV** default; ULSMV a documented option.
16. `loadings_se` → added to the §4 level contract (p×k matrix, NULL for PCA/EFA).

**Completed in M4 (in addition to items 11–16 above):**
- NPD polychoric matrices: warn + smooth via `psych::cor.smooth()` (implemented).
- Convergence truncation tested for the ESEM engine (k=4+ on p=6 triggers lavaan error → warn + truncate → object builds to k_eff=3).
- Algebra-vs-scores cross-check documented as Pearson/continuous paths only; polychoric ESEM edges (algebra uses lavaan polychoric R; scores route uses Pearson standardization) diverge by design and are excluded from the oracle.

**Resolved for M5 (Forbes extension):**
17. API → ~~two orthogonal args: `pairs = c("adjacent","all")` and `prune = c("none","redundant","artefact")` (char vector; default `"none"`). `prune != "none"` auto-upgrades `pairs` to `"all"` with a loud `cli_inform()` (chains need all-levels edges to assess).~~ **Superseded by M34** (item 27 below): pruning is now a standalone `prune()` verb, not an `ackwards()` argument, and `pairs` no longer auto-upgrades.
18. Prune action → **flag-only, never remove.** Adds `pruned`/`prune_reason` annotation columns to the edge/node tidy structures; the object retains all levels (preserves Invariant 5 and the algebra-vs-scores oracle). Pruning is interpretive relabeling, **not** re-estimation — say so in `print`/docs (extends the §2 honesty caveat).
19. Redundancy → faithful to Forbes (2023): score-correlation chains `|r| ≥ .9` (default, tunable), retention rule = keep the bottom node if the chain reaches level k (most-specific, best-defined), else keep the topmost node (broadest manifestation). Tucker's φ (`> .95`, Lorenzo-Seva & ten Berge 2006) computed on aligned loadings as an **optional conjunctive** criterion. φ formula `Σaᵢbᵢ / sqrt(Σaᵢ² · Σbᵢ²)`; base R, no new dependency.
    - **M53 amendment — `redundancy_criterion` (default `"direct"`).** The "chain" is now traced by Forbes's actual `ChaseCorrPaths` rule: chase upward via the **direct (skip-level)** correlation between a factor and each ancestor level, contiguously while `|r| ≥ redundancy_r`. The pre-M53 implementation used **adjacent** primary-parent links only; the M41 "Σr² ≤ 1 ⟹ adjacent ≡ direct" argument holds on shallow/transitive hierarchies (all three simulations agree) but **fails under non-transitivity** in deep ones — on Forbes's 10-level AMH applied example the two diverge on 7 of 54 components. Direct is the honest operationalization of "same construct" (two scores share `≥ redundancy_r²` variance *directly*) and reproduces her published AMH chase exactly (M53, test-backed). `"adjacent"` is retained as an opt-in. This supersedes the implicit adjacent-chain assumption in items 19–20.
20. **Additive enrichments over the paper**: always *report* both `r` and `φ` for every redundancy candidate (report-first, flag-second — borderline cases like the paper's own `.89`/`.93` alcohol component stay visible); report endpoint `r` (direct root-to-leaf, from all-levels edges) alongside the chain as an at-a-glance cross-check. **(M53:** with the direct criterion now the default, the earlier "chain (perpetuates at every level) vs endpoint `r` (same construct)" framing describes the *adjacent* opt-in specifically; under `"direct"` the chain *is* the same-construct criterion.**)**
21. Artifact → **never auto-flagged.** `prune(x, "artifact")` surfaces φ for inspection; removal is a documented researcher judgment (Forbes is explicit this introduces researcher DoF / confirmation bias; cf. Wicherts et al. 2016).

**Resolved for M32 (API-shape & naming; owner-reviewed, no equivalent guidance elsewhere in this
document):**
22. `tidy(x, what = "fit")` long-format key column → **`statistic`**, replacing `index` (which read
    as a row position; it actually held fit-index names for EFA/ESEM and eigenvalue positions for
    PCA). Wide format (`format = "wide"`) unaffected — column names still come from the values of
    this key.
23. `k_max` naming collision (`ackwards()` = extraction depth vs. `suggest_k()` = max
    factors/components evaluated) → **keep the shared name in both functions.** They're genuinely
    the same dial applied at different workflow stages (`suggest_k(k_max = ...)` →
    `ackwards(k_max = ...)`); renaming one would create vocabulary drift without removing any real
    ambiguity. Resolved instead via roxygen: each function's `@param k_max` states its own meaning
    and cross-references the other's.
24. `tidy(what = "fit")` cutoff flags → **removed.** The `cutoffs = TRUE` argument and the
    `meets`/`{statistic}_meets` columns it produced (`.flag_fit()`) are dropped; a pass/fail boolean
    quietly endorsed Hu & Bentler (1999) thresholds that the package elsewhere (§9, this section)
    treats as conventional and contested — continuing the M28/M31 output-honesty trajectory.
    `.fit_cutoffs()` is retained as an internal helper: `autoplot(what = "fit")` still draws
    threshold reference lines and `summary()` still annotates inline with a check/cross mark: both
    are visual guides, not a returned pass/fail column a user could mistake for a computed
    judgment.
25. Variance scale → **proportion, 0-1** (`tidy(what = "variance")` columns renamed
    `variance_pct`/`cumulative_pct` → `proportion`/`cumulative`, values divided by 100). This
    aligns `tidy()` with the engine's internal `variance` slot (already 0-1; see `print.ackwards`,
    which reads it directly) and with broom/psych convention (`psych` reports "Proportion Var" as
    0-1). Percent formatting moves to the display layer (`print()`, `summary()`, vignette `gt`
    tables) rather than living in the tidy data.
26. **M31-deferred: effective ESEM estimator recorded in `$meta`.** M31 explicitly deferred this
    ("better bundled with M32's meta/column decisions than bolted on here"): `x$meta$estimator`
    now stores the effective estimator after auto-selection (`"ML"`/`"MLR"`/`"WLSMV"`/`"ULSMV"`;
    `NA` for PCA/EFA). `summary()` gains a one-line footnote naming lavaan's scaled-variant
    reporting (§14.M31 point 2/Post-review) whenever the effective estimator is
    `"WLSMV"`/`"ULSMV"`/`"MLR"` — not shown for `"ML"`, which has no scaled variant.

All five M32 changes are breaking with no deprecation path (pre-CRAN, no users; consistent with the
M34 pruning-verb precedent of clean moves over compatibility shims).

**Resolved for M34 (pruning verb; owner-reviewed, breaking, no deprecation path — pre-CRAN, no
users):**
27. Pruning extracted into a standalone, pipeable S3 generic `prune()` (`prune.ackwards`), not an
    `ackwards()` argument. The five prune-related args (`prune` → renamed `rules`, `redundancy_r`,
    `redundancy_phi`, `min_items`, `orphan_r`) leave `ackwards()` entirely:
    `ackwards(...) |> prune(...)`. Rationale: extraction is the expensive, deterministic step;
    pruning is cheap and interpretive. Separating them lets a researcher re-prune with new
    thresholds without re-extracting (re-running ESEM per level is the expensive path M26
    optimized). `prune()` is a generic (`UseMethod`), not a plain function, so it coexists with the
    `prune` generics already defined by recursive-partitioning packages (e.g. `rpart::prune`)
    regardless of package load order. Returns the same `ackwards` object with `$prune` populated
    (replacing any prior pruning) — no new class, so `print`/`summary`/`tidy`/`glance`/`augment`/
    `autoplot` all work unchanged.
28. Edges for redundancy chains and artifact φ are recomputed **fresh inside `prune()`** via
    `compute_edges(pairs = "all")` from the object's stored `levels`/`r` (Invariant 3), and never
    written back to `x$edges` (Invariant 1: one edge path). `ackwards()`'s `pairs` auto-upgrade to
    `"all"` when pruning was requested (M5 item 17) is removed along with it — `pairs` is now a
    pure display/storage setting on `ackwards()`, decoupled from pruning; `prune()` works correctly
    regardless of the fit-time `pairs` value.
29. Manual pruning: `prune(x, rules = "none", manual = c("m4f3", "m4f4"))` flags user-named nodes
    directly (standalone, no auto rule needed), or unions them onto an auto rule's flags
    (`prune(x, "redundant", manual = c(...))`). Unknown labels error. On overlap between an auto
    rule and `manual`, the auto rule's `prune_reason` wins (more informative); only
    otherwise-unflagged manual nodes get `prune_reason = "manual"`.
30. Naming: canonical rule name is **`"artifact"`** (US spelling, the owner's preference), with
    `"artefact"` accepted as an alias (nod to Commonwealth spelling and to Forbes' own usage) and
    normalized internally to `"artifact"`. Existing code passing `"artefact"` keeps working.
    `"tucker"` was considered and rejected as an alias/rename: the mode surfaces more than Tucker's
    φ (it also computes the `few_items`/`orphan`/`split_merge` structural signals), so naming it
    after the statistic would mislabel the umbrella — `"artifact"` names what the mode is *for*.
31. `ackwards()` explicitly rejects the five removed args if passed via `...` (rather than silently
    absorbing them) with a pointer to `prune()` — silent absorption would be a masked-argument
    footgun, not the clean break intended (Invariant 6: loud, not silent).

**Resolved for M38 (`missing = "fiml"` for PCA/EFA; reverses a resolved default —
owner-signed-off):**
32. **FIML promoted to a first-class PCA/EFA route.** The M16 contract "`missing = "fiml"`
    errors for PCA, EFA" is **reversed for the Pearson basis**: under `engine = "pca"/"efa"` with
    `cor = "pearson"`, `missing = "fiml"` now routes R through `psych::corFiml()` (full-information
    ML, MVN) and feeds it to the existing `W'RW` algebra — Invariant-1-clean (one edge path, one
    corFiml call per run since R is computed once and reused across all levels; no new dependency,
    `psych` already Imports). It **still errors** for a non-Pearson PCA/EFA basis (Spearman,
    polychoric — corFiml cannot honor them) and for WLSMV/ULSMV. `.resolve_missing()` gained a `cor`
    argument to enforce this guard matrix. Grew out of the M37 doc-planning observation that a user
    could already smuggle FIML in via the M22 correlation-matrix seam (`ackwards(corFiml(x), …)`);
    M38 makes the capability discoverable and owns the `n_obs` tradeoff.
33. **`n_obs` string option + default.** On the raw-data FIML PCA/EFA path, `n_obs` accepts
    `"total"` (**default**) or `"complete"`. `"total"` = every row contributing to the FIML
    likelihood, matching the convention a genuine FIML analysis reports (Enders 2010) and giving
    cross-engine parallelism with ESEM-ML FIML. Point estimates (loadings/edges) do not depend on
    the choice; only the EFA fit indices do, and those are *approximate* under this two-step
    (FIML matrix → normal-theory EFA) procedure regardless of N (Zhang & Savalei 2020) — so the
    route announces itself and the caveat via cli (Invariant 6). `"effective"` was considered and
    **dropped**: no canonical formula exists for a corFiml→EFA route, so it would be a
    package-invented convention masquerading as a standard. A string `n_obs` is rejected off this
    path; cor-matrix input still requires a numeric N.

**Resolved for M45 (out-of-sample scoring; owner-approved 2026-07-01):**
34. **Scoring standardizes by fit-time moments by default.** `augment()`'s `scaling` argument
    defaults to `"fit"`: supplied data are standardized by the training means/SDs (stored in
    `meta$item_means`/`item_sds` since M45) before the weight matrices apply. Rationale: this is
    what "applying the trained model" means — a test observation's score must not depend on which
    other observations share its split, and train/test scores must share one metric (the
    cross-validation use case that motivated the milestone). `"sample"` (the pre-M45 behaviour)
    remains as an explicit opt-in for deliberately re-standardizing a sample from a different
    population in its own metric, and is the only option for cor-matrix objects. Breaking for
    subset/new-data scoring values (pre-CRAN, no deprecation path; M32/M34 precedent). The
    companion `predict.ackwards()` export was chosen over an augment-only surface for
    discoverability by replicators (`psych`/`lavaan` users reach for `predict()`).

**Resolved for M46 (Girard extension — replicability gate; owner-approved 2026-07-01):**
35. **Split-half factor comparability as a standalone verb.** `comparability(data, k_max, engine,
    cor, fm, n_splits = 10, seed)` resurrects Everett (1983) factor comparability coefficients —
    the split-half replication gate of Goldberg's research program (Saucier 1997; Saucier et al.
    2005; *not* Goldberg 1990, which contains no split-halves — citation corrected 2026-07-16) —
    as the hierarchy-depth gate, extended to every bass-ackwards level.
    Mechanics: a full-sample `ackwards()` fit anchors labels (results report the same `m{k}f{j}`
    the user will see); per split, levels 1..k are fit independently in each random half via the
    engine internals; each half-solution is matched to the anchor by **greedy-with-removal**
    max-|r| bijection (square, same-k — well-posed here, unlike §7 parent matching; no `clue`
    dependency, and near-ties are themselves the instability being measured); the coefficient is
    the correlation between the two matched half-solution scores on the **pooled** correlation
    matrix, computed by `compute_edges()` on a two-element levels list (Invariant 1 — one edge
    path), signed after orienting both matched factors positively toward the anchor (a negative
    coefficient stays visible as a diagnostic). Tucker's φ on the matched loading columns is
    reported alongside (score agreement vs. pattern agreement). **Report-first, judge-never**:
    nothing auto-flagged; the conventional .90/.95 benchmarks appear only as reference lines /
    footer prose (M32 cutoff philosophy). Scope: **PCA/EFA only** (ESEM deferred — `2 * n_splits`
    lavaan hierarchies need their own performance treatment; logged in `ROADMAP.md`) and
    **pearson/spearman only** (the `suggest_k()` precedent: polychoric estimation in every
    half-sample is slow and NPD-prone; screen on Pearson, fit the final model polychoric).
    Repeated splits default `n_splits = 10` (Goldberg's practice; a single split misleads by luck
    of the draw). Per-half engine warnings are muffled and any convergence shortfall is
    summarised once (`summary$n_splits_ok`). The recommended workflow — suggest_k ceiling →
    comparability floor → fit → prune → interpret → out-of-sample validation — is documented in
    `vignette("ackwards-girard")` and named the *replicability-gated (or Girard) workflow*.

36. **Bootstrap edge CIs as a standalone verb.** `boot_edges(x, data, n_boot = 1000, conf = 0.95,
    seed)` attaches nonparametric bootstrap SEs + **percentile** CIs to every edge — the
    inferential-honesty companion to the point-estimate `prune()` rules and the Forbes strongest-
    edge practice (resurrects the §14 e4 deferral). A **standalone verb** (not an `ackwards()`
    argument), consistent with the M34 direction of extracting expensive/optional work into pipeable
    verbs (`prune()`, `comparability()`); the raw `data` is re-supplied because the light core does
    not store it (Invariant 3). Mechanics per replicate: resample rows with replacement → recompute
    `R` with the fit's `cor`/`missing` routine → refit levels 1..k_max (engine internals, muffled)
    → **anchor** each replicate level to the full-sample solution (the item-35 greedy-max-|r|
    matching + sign orientation, evaluated on the full-sample `R`) → edges via `compute_edges()` on
    the replicate `R` (Invariant 1). Anchoring is load-bearing: without it, factor label-switching
    and sign-flipping across replicates corrupt the pooled distributions. **All resample indices are
    drawn upfront from `seed`**, so each replicate is deterministic given its indices and the serial
    and `future.apply`-parallel paths (M26 pattern) agree bit-for-bit — the design choice that makes
    "serial ≡ parallel" test-assertable. Failed replicates drop to `NA` and are counted per edge
    (`n_boot_ok`), never aborting (Invariant 7). Output attaches to `x$boot`; `tidy(what = "edges")`
    gains `se`/`lo`/`hi`/`n_boot_ok`, and `print`/`summary` note coverage. Scope: **PCA/EFA +
    pearson/spearman only**, mirroring item 35 (ESEM = `n_boot` lavaan hierarchies; polychoric =
    per-resample polychoric estimation, slow + NPD-prone; both logged in `ROADMAP.md`). Cor-matrix
    input, ESEM, and polychoric objects error with pointers. **DESIGN §14's original "reuse the
    `loadings_se` infrastructure" phrasing is amended**: `loadings_se` is lavaan's analytic
    delta-method SE, which has no analogue for edges of varimax-rotated hierarchies — bootstrap *is*
    the SE method; what is reused is the storage/tidy display pattern. **Percentile (not Fisher-z or
    normal-approx) CIs**: percentile respects the [-1, 1] bound and the skew of `r` near the 0.9
    prune threshold. *Oracle note:* the full-pipeline bootstrap SE legitimately **exceeds** the
    Fisher-z analytic SE of the materialized-score correlation, because each replicate re-extracts
    the factors (loading uncertainty) while Fisher-z treats scores as fixed observed variables; the
    exact Fisher-z match holds only for a *fixed-weights* percentile bootstrap (both halves are
    test-asserted). The intervals are **per-edge error bars, not a familywise correction** — they
    make each edge's precision visible but do not undo the selection bias of scanning many edges for
    the strongest (stated in the object docs + girard vignette).

**Resolved for M50 (release polish; owner-approved 2026-07-02):**
37. **`bfi25` carries public-domain IPIP variable labels.** Each item column gets its IPIP marker
    stem (Goldberg, 1999; ipip.ori.org) as a plain `label` attribute, populating the M36 capture
    path so `top_items()` prints `code: label` with zero setup. Sourced directly as public-domain
    IPIP items (text is necessarily identical to `psych::bfi.dictionary` — same public pool — but
    not framed as a copy). Plain attributes are dropped by base row-subsetting (`na.omit()`, `[`),
    unlike the class-based `labelled`/`haven` vectors M36 targets, so the documented pattern is to
    fit the dataset directly (`missing = "listwise"`) rather than pre-filtering. A `label_items()`
    setter and a third teaching dataset were **declined** (logged in M49 Phase A to avoid
    double-logging); a persistent *factor*-label pipeline (distinct from these *item* labels) is
    deferred to 0.2.0.
38. **`suggest_k()` ordinal-detection warning (Invariant 6 symmetry).** `suggest_k()` now runs
    `detect_ordinal()` on raw-data input and warns once per session, matching `ackwards()` and
    `comparability()`. Because `suggest_k()` screens on the Pearson/Spearman basis *by design*, the
    wording points at the final `ackwards()` fit (`cor = "polychoric"`), not at `suggest_k()` — it
    does not (and cannot) switch basis itself. Guarded to the raw-data path (a correlation matrix
    has no items to inspect).
39. **Console/example polish (no contract change).** Engine name renders lowercase across all
    print methods; `summary()` per-level figures use fixed decimals with trailing zeros; pruned
    footers consolidate to one note. Roxygen examples split by intent: mechanics demos use the
    continuous `sim16` (no ordinal warning, faster checks), content demos keep `bfi25` on the
    polychoric basis. These are display/documentation choices, not behavioural defaults.

**Resolved for M49 Phase A (roadmap cleanup; owner-approved 2026-07-02):**
40. **Dual EFA chi-squares (review enhancement e2) — declined.** A proposed `chi_empirical` column
    (psych's residual-based `$chi`) reported alongside the likelihood-ratio `chi` the EFA fit row
    now carries (the M42/C1 fix). **Rejected as net-negative:** it re-opens the exact mispairing C1
    fixed (a second chi-square inviting users to pair it with the wrong p-value); it is an
    NA-heavy, EFA-only column in a fit schema shared with ESEM (which has no empirical chi-square),
    adding asymmetry + doc burden; and it has **zero downstream consumer** in the package (anyone
    who wants it can read it off a `keep_fits = TRUE` object). Same shape as the M40 `categorical`
    decline: a confusion surface for no capability. Removed from `ROADMAP.md`.
41. **Item-label ecosystem — one shipped, two declined, one deferred.** M50 shipped bfi25's
    public-domain IPIP variable labels (item 37). Three adjacent asks were resolved here to avoid
    double-logging: (a) a **`label_items()` setter** is **declined** — it would duplicate
    `labelled::var_label()` (the established idiom `.capture_item_labels()` already reads) and, worse,
    deepen the `label` overload between *item* labels and `label_template()`'s *factor* labels; the
    `attr(col, "label") <- ...` / `labelled::var_label()` route is documented instead. (b) A **third
    teaching dataset** is **declined** — `sim16` (idealized continuous) and `bfi25` (realistic
    ordinal) are a deliberate two-foil pair; a third dilutes that contrast without adding a distinct
    teaching case. (M54 clarification: this declines a third *teaching* dataset only. The exported
    `forbes2023` matrix is a *fidelity/reproduction* dataset — the published case the Forbes-fidelity
    contract reproduces, not a pedagogical foil — and is therefore outside this decision's scope.) (c) A persistent **factor-label pipeline** (a `set_factor_labels()`-style setter
    honored by `autoplot`/`print`/`summary`/`tidy`, with `label_template()` as scaffold) is
    **deferred to 0.2.0** (not declined): it is purely additive, so 0.1.0 loses nothing, and it
    deserves its own naming/storage/precedence design. The **item-vs-factor "label" vocabulary
    split** is load-bearing for all of the above — the two concepts must stay lexically distinct in
    docs and any future API. (Banked in `ROADMAP.md`.)

**Resolved for M49 (polychoric robustness fix; owner-approved 2026-07-02):**
42. **`correct` argument exposes the polychoric continuity correction.** Real-world ordinal data
    surfaced a hard failure: `psych::polychoric()` errors under its default continuity correction
    (`correct = 0.5`) when an item has a near-empty (singleton) response category or items with
    unequal category counts produce a sparse cross-cell, and — because psych runs the item pairs in
    parallel — one failure collapses the whole call into an opaque `supply both 'x' and 'y'` message.
    psych's own advice is `correct = 0`, but `ackwards()` gave no way to pass it. Resolution:
    `ackwards()` gains `correct = 0.5` (matching psych's default and name for discoverability),
    forwarded to `psych::polychoric()` on the **PCA/EFA polychoric path only** (ESEM computes its own
    polychoric inside lavaan; Pearson/Spearman ignore it). Additive, non-breaking — the default
    reproduces prior behaviour. The failure message was rewritten to name the `correct = 0` remedy
    instead of passing psych's opaque error through, and the NPD guard was hardened to catch a
    polychoric matrix with `NA`/`NaN` entries (naming the offending items) before `eigen()`, rather
    than crashing with base R's "missing value where TRUE/FALSE needed". Not a reimplementation of
    polychoric — a wrap-and-diagnose fix (§ "wrap established engines").
43. **`check_items()` — pre-analysis item screen (new export).** The `correct` fix (item 42) is
    reactive; it does not stop a *near-constant* item from silently producing a plausible-looking but
    meaningless factor, nor does it name the offender, nor does it catch a truly constant item (which
    `psych::polychoric()` silently deletes, crashing downstream with `subscript out of bounds`).
    `check_items(data, cor)` reports, one row per item, the stats that predict these failures
    (`n_valid`, `pct_missing`, `n_distinct`, `min_count`, `top_prop`) and a worst-case `flag`
    (`constant` / `near-constant` / `sparse category` / `high missing` / `ok`). Thresholds are
    deliberately conservative — `constant` = one distinct value; `near-constant` = a category holding
    ≥ 95% of responses (the case that actually breaks polychoric); `sparse category` = smallest
    observed category < 5, flagged **only** under `cor = "polychoric"` (rare-but-present Likert
    categories usually fit fine, so this is report-only). `ackwards()` runs the **same** internal
    screen (`.screen_items()`, shared — DRY): it **errors** on a constant item (naming it, before
    psych can delete it) and **warns once** on a near-constant item, but deliberately does **not**
    warn on a merely sparse category (avoids nagging ordinary ordinal data). Report-first, like
    `suggest_k()`/`comparability()`; never modifies data. New export, no new dependency.
44. **Trust guidance + durable near-singularity signal.** Fit-time warnings are transient (they
    scroll off; a reloaded/shared object shows none), so a user can trust a rank-deficient fit
    whose `CFI` is silently `NA`. Two additions: (a) a **"When to trust the result"** `@section` in
    `?ackwards` tiering every diagnostic ackwards raises as *fatal* (constant item, polychoric
    failure, non-convergence, near-singular matrix) / *caution* (Heywood, near-constant item,
    ordinal-on-Pearson) / *informational* (pairwise-missing, sparse category, ordinal-detection when
    Pearson was intended); (b) the near-singularity check (`.near_singular_check()`, shared by the
    raw-data and cor-matrix paths — **basis-agnostic**: it also catches redundant items on the
    Pearson/Spearman/FIML bases and user-supplied matrices, not just polychoric) records
    `meta$min_eigenvalue` and `meta$near_singular` (smallest eigenvalue `< 1e-4`) so
    `print()`/`summary()` re-surface a **durable** "near-singular -- fit indices and scores may be
    unreliable" caution on the object itself, not just once at fit time. The messaging avoids
    singling out CFI (an ESEM-only index — for EFA the tell is TLI/RMSEA from psych's residual
    fallback). Report-only; changes no estimate. (M49 review follow-up generalized this from the
    initial polychoric-only implementation.)

**Resolved for M51 (factor-label pipeline; first milestone of the 0.2.0 cycle; owner-approved
2026-07-02):**
45. **Persistent factor labels — new `set_factor_labels()` verb + `factor_labels()` getter.**
    Deferred-not-declined at §14.41(c); this resolves the naming/storage/precedence design. It is
    the *factor*-label counterpart to M50's *item* labels (`meta$item_labels`, item 37) and the two
    must stay lexically distinct — every symbol here says **factor** label; nothing touches
    `item_labels`.
    - **Storage.** `x$meta$factor_labels`: a named character vector, names = factor IDs
      (`m{k}f{j}`), values = user strings. Partial labeling is allowed (a factor with no entry
      falls back to its ID everywhere). Living in `meta` (not a top-level slot) means it rides
      through `prune()`, `boot_edges()`, `augment()`, and `predict()` unmodified with no per-verb
      code — those verbs already copy `meta`. `NULL`/absent when no label is set (the default),
      so an unlabeled object is byte-identical in behaviour to a pre-M51 one.
    - **API — a verb + a getter, not a replacement function.** `set_factor_labels(x, labels)`
      returns the modified `ackwards` object (pipeable, matching the package's verb grammar:
      `prune()`, `boot_edges()`). Repeated calls **merge/update** (so you can label incrementally);
      a `labels` value of `NULL` clears **all** labels; an `NA` or `""` value for a given ID
      **removes** just that one. The setter **errors** (not warns) on any name matching no factor ID
      — a typo'd ID is a mistake, and the object's IDs are knowable up front (consistent with the
      package's loud-validation posture, Invariant 6). It does **not** warn when two levels share a
      value (`m2f1` and `m3f1` both "Internalizing" is the ordinary hierarchical case). The
      companion getter `factor_labels(x)` returns the stored named vector (or `NULL`). Chose the
      verb+getter pair over a `factor_labels(x) <- value` replacement method to keep one pipeline
      idiom across the package and avoid the `names<-`-style in-place mutation that reads oddly next
      to the object's otherwise-functional API.
    - **Display form `label (id)`.** Every text surface that lists a factor shows
      `Neuroticism (m5f1)` when a label is set, bare `m5f1` otherwise. Keeping the ID visible is
      load-bearing: it is the handle a user needs to index into `tidy()`/`edges`, and mirrors how
      `autoplot()` keeps IDs discoverable. The rule is simply: **wherever a text surface names a
      factor, it shows `label (id)`.** Concretely: `summary()` (both the per-level variance listing
      and the lineage tree), `print()` (one status line noting how many of the factors are labeled,
      only when any are), and `top_items()` (the factor dimension either way — group headers under
      `by = "factor"`, body entries under `by = "item"`). **`autoplot()` is the exception** — a
      diagram node
      shows the substantive label *only* (no parenthetical ID), because a stored label is exactly a
      persistent default for the existing `node_labels` argument, which has always replaced the node
      text wholesale. Precedence: stored labels form the node-text baseline; a call-time
      `node_labels` entry overrides per node (call-time beats stored; every existing script is
      unaffected).
    - **`tidy()` — label columns appear only when set (Invariant 5).** The ID columns (`factor` in
      loadings/variance/scores; `from`/`to` in edges) are **never** mutated — lineage lives in IDs.
      When labels are set, `tidy()` gains a `factor_label` column (loadings/variance/scores) and
      `from_label`/`to_label` (edges), carrying the label for labeled factors and `NA` otherwise. An
      **unlabeled** object's `tidy()` output keeps its pre-M51 column set byte-for-byte — zero churn
      for existing consumers, which is why the columns are conditional rather than always-present.
    - **Out of scope for this verb.** `comparability()` fits fresh hierarchies from data (not from a
      labeled object), so it takes no labels. `augment()`/`predict()` score **column names** stay
      `m{k}f{j}` — they are syntactic identifiers a user binds to, not display text. `label_template()`
      is unchanged; its printed `c(...)` literal is now framed in the docs as the scaffold you fill
      in and feed to `set_factor_labels()`.
46. **Factorability & data-adequacy diagnostics — new `factorability()` export + an internal
    `ackwards()` screen (M52).** Adds the pre-fit "is this matrix worth factoring, and is N large
    enough" layer that item 43 (`check_items()`, per-item) and `suggest_k()` (depth) did not cover.
    No new dependency (`psych` already supplies `KMO()`/`cortest.bartlett()`; the Ledermann bound is
    base arithmetic), no signature or default change to `ackwards()`.
    - **What it reports.** KMO sampling adequacy (overall + per-item MSA), Bartlett's sphericity
      test, N / p / N:p, and the Ledermann bound (largest number of common factors identifiable from
      p variables, `.ledermann_bound(p)` = max k with `((p-k)^2 - (p+k))/2 >= 0`). Computed on the
      chosen `cor` basis, from raw data **or** a correlation matrix (N-based rows are `NA` when a
      matrix arrives without `n_obs`).
    - **Honesty stance (extends item 24).** Item 24 removed the Hu–Bentler fit pass/fail flags
      because the package refuses to endorse contested cutoffs as returned booleans. The same
      principle governs here: `factorability()` reports values and their **conventional bands**
      (Kaiser KMO labels, the 5:1/10:1 N:p rules, Bartlett at .05), each explicitly framed in the
      printout and roxygen as a *rule of thumb, not a settled threshold* (required N depends on
      communalities and overdetermination — MacCallum et al. 1999). No pass/fail column is returned.
    - **Internal `ackwards()` screen — engine-split, advisory-only (Invariant 6/7).** `ackwards()`
      runs `.factorability_screen()` on every fit and emits at most two `.frequency = "once"`
      warnings, never aborting and never altering the fit: (a) a **Ledermann** under-identification
      warning naming the first over-bound level, fired **only for EFA/ESEM** (`k_max` compared to
      the bound) — PCA is exempt because components are exact linear combinations with no
      latent-model df constraint (its existing `k_max <= p` ceiling suffices); and (b) a single
      **consolidated sampling-adequacy** warning when KMO `< 0.5`, N:p `< 5`, or N `< 100`,
      pointing to `factorability()` for the full report. The Ledermann warning compares against the
      *requested* `k_max` (the design error to flag) and warns-and-builds — a truncated level still
      warns separately, so the two are complementary.

**Known limitations / deferred to future milestones:**
- `factor_cor` in the ESEM engine is not permuted by the variance-sort `ord` vector. Safe permanently: only orthogonal rotation is supported (`factor_cor = I`; permutation of I is I), and oblique rotation is out of scope (§9, §14.1). The guard comment in `engine_esem.R` documents what *would* be required if that decision were ever reversed.
- Algebra-vs-scores cross-check does not cover `cor = "polychoric"` paths (see above), nor the
  `missing = "fiml"` PCA/EFA path (M38): there the algebra uses the `psych::corFiml()` matrix while
  the scores route standardizes the raw, NA-bearing data (pairwise Pearson SDs), so the two bases
  diverge under missingness by design — the same reason polychoric is excluded. The oracle tests
  therefore run only on complete-data linear engines; no FIML/polychoric object is fed to them, so
  there is no false-failure risk, but the cross-check does not *certify* those paths.
- `cor = "spearman"` + `engine = "esem"` is semantically inconsistent (lavaan fits Pearson ML on raw data while edges use Spearman R); a warning is now emitted (M10).
- ~~ESEM engine does not detect or warn on improper/Heywood solutions~~ — **resolved M10**: engine now inspects `lavaan::lavInspect(fit, "theta")` and warns when `any(diag(theta) <= 0)` (parity with EFA engine, Invariant 7).
- **ESEM ML/MLR with `missing = "pairwise"`**: lavaan uses listwise deletion for the model fit while edges are computed from a separately-computed pairwise `stats::cor()` — a minor inconsistency (fit statistics at complete-case N, edges at full pairwise N). Documented in `$meta$missing`; a per-call advisory warning fires whenever NAs are detected. Use `missing = "listwise"` or `"fiml"` to resolve. *(Added M16; §9 `missing` row cross-references this.)*
- **Forbes-extension improvements deferred past M5** (the published method has weaker spots worth strengthening later; M5 ships the faithful method + the §14.20 reporting enrichments):
  - ~~*Structural artefact signals.*~~ **Done M25 (Wave 2).** `prune = "artefact"` now populates `x$prune$structural` with per-factor `few_items` / `orphan` / `split_merge` signals (flag/report only; no auto-pruning). Args: `min_items = 3L`, `orphan_r = 0.5`. `split_merge` is `TRUE` when a factor's primary items came from multiple different primary factors at the preceding level. CLI and `print()`/`summary()` report the flagged count.
  - ~~*Factor-score-indeterminacy caveat for EFA/ESEM redundancy.*~~ **Done M25 (Wave 3).** `redundancy_phi = NULL` now auto-resolves: PCA → no φ filter; EFA/ESEM → `0.95`. `NA` is the explicit opt-out. Announces via cli when auto-applied (Invariant 6). See §9 defaults table.
  - ~~*Selection bias in the "strongest" edge.*~~ **Done M47 (item 36).** `boot_edges()` adds
    nonparametric bootstrap SEs + percentile CIs on every edge (PCA/EFA + pearson/spearman; upfront
    seeded indices → serial ≡ parallel; full-sample anchoring per replicate; `future.apply`
    dispatch; `tidy`/`print`/`summary` integration; girard-vignette honesty section). The intervals
    quantify *per-edge* precision but by design do **not** correct the multiplicity of scanning many
    edges for the strongest — that caveat is documented, not silently "fixed".
- **M40 spin-off (code/viz deferred out of the doc-only M39; shipped M40 except the declined flag):**
  - ~~*Ordinal `categorical` convenience flag.*~~ **Declined M40 (owner sign-off, 2026-07-01).** A proposed `categorical = TRUE/FALSE` argument on `ackwards()` flipping `cor` (pearson→polychoric) and the ESEM estimator (ML/MLR→WLSMV) together. **Rejected as redundant:** `cor = "polychoric"` *already* auto-selects WLSMV (the `estimator = NULL` auto-rule, §9), so `categorical = TRUE` would be a pure synonym for `cor = "polychoric"` — not a two-settings-in-one shortcut. It would add a conflict surface (`categorical = TRUE` + `cor = "pearson"`?) and a §9 defaults change for zero new capability. Discoverability is handled at the docs layer instead (the ordinal vignette already states that `cor = "polychoric"` alone flips the estimator; the ordinal-detection cli warning names the option at runtime). No `R/`/`DESCRIPTION`/§9 change.
  - ~~*Ordinal correlation-comparison visualization.*~~ **Done M40.** The two raw `round(x$r[1:5,1:5], 2)` matrix chunks in `vignettes/ackwards-ordinal.Rmd` were replaced with one dodged bar chart of the ten `N1`–`N5` lower-triangle item-pair correlations, `fill = basis` (Pearson vs polychoric), reshape code hidden. Viz-only; `ggplot2` already in Suggests; no package-code or dependency change. (The gt long-format Δ-table alternative was set aside — the vignette already carries two gt Δ-tables, so a chart adds variety and directly addresses "matrices are hard to compare".)
  - ~~*Forbes pruned-level axis-label styling.*~~ **Done M40.** A **fully-pruned** level (every factor flagged) now gets an *italic* `autoplot()` axis label in the normal (non-`drop_pruned`) render path, denoting its status alongside the existing grey node fill; partially-pruned levels keep a plain label. Automatic (no new argument). See §11.

