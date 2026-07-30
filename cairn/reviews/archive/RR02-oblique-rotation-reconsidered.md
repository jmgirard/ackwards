# RR02: Does oblique rotation actually confound the between-level signal? — Review Report

- **Date:** 2026-07-30
- **Brief:** `cairn/reviews/RB02-oblique-rotation-reconsidered.md` (no milestone)
- **Reviewer:** independent Fable-level review, per `/milestone-brief`
- **Materials read:** `cairn/DESIGN.md` §4, §5.1, §9, Design principles, Known limitations;
  `cairn/DECISIONS.md` D-002, D-004, D-007, D-009, D-010, D-017, D-031, D-032;
  `cairn/reviews/archive/RR01-oblique-algebra-claim.md` (in full);
  `R/compute_edges.R`; `R/utils.R` (`.score_var`, `.variance_explained`,
  `match_parents`, `.align_signs`, `flip_weights`, `.compute_scores`);
  `R/engine_pca.R`; `R/engine_efa.R`; `R/engine_esem.R`; `R/prune.R`;
  `R/ackwards.R` (roxygen rotation bullet; lineage/alignment assembly);
  `cairn/references/goldberg2006.md`, `forbes2023.md`, `waller2007.md`;
  shelf PDFs rendered as page images (not pdftotext): goldberg2006 pp. 348, 349,
  350, 356; kim2015 pp. 1066, 1067; forbes2023 p. 2 (main text + footnotes 1–2).

## Verdict in one line

The confounding claim is half-true: the oblique between-level correlation is a
different, fully interpretable quantity (a structure/total correlation) rather than
an uninterpretable one — what oblique genuinely breaks is not the *edge* but the
*lineage overlay* read onto the edges (primary parent, split narrative, additive
variance partitioning), and a Φ-partialled decomposition that restores the lineage
reading is exactly identified from what the package already stores; the DESIGN §9
attribution to Goldberg (2006) does not survive the page image, and Kim & Eaton
(2015) assert the claim bare. D-002's recorded rationale should be corrected
regardless of the decision; the decision itself is worth superseding to
oblique-as-documented-non-default, with the user-facing semantics settled at the
Forbes design session before implementation.

## Q1. Is the confounding claim true as stated?

**Formal statement.** Per the settled RR01 algebra (not re-derived here), for any
fixed linear scoring `S_x = Z W_x` the package's edge is exactly

```
E_ab = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2},   D_x = diag(W_x' R W_x)
```

(`R/compute_edges.R:95-101`; `.score_var()` `R/utils.R:104-106`). `E_ab[i, j]` **is**
the correlation between the standardized score composites `s_ai` and `s_bj` — under
oblique rotation exactly as under varimax. For PCA the composites are determinate,
so `E_ab[i, j]` is the correlation *between the oblique components themselves*, not
an estimate of it; for EFA/ESEM with correlation-preserving scoring it estimates the
model-implied cross-level factor correlation, subject to RR01 Q3's precision note
(score-level Φ_s equals the latent Φ only for correlation-preserving scoring rules).

**What changes under oblique is not the quantity but its relation to the lineage
question.** Let `Φ_sa = D_a^{-1/2} (W_a' R W_a) D_a^{-1/2}` be the within-level score
correlation at the shallower level (Φ_sa = I under the shipped varimax defaults).
Regress each level-b score on the full set of level-a scores. The standardized
regression coefficients are

```
B_ab = Φ_sa^{-1} E_ab        so that        E_ab = Φ_sa B_ab .
```

This is the classic structure-vs-pattern distinction transplanted to the
between-level edges:

- **`E_ab[i, j]` (structure / marginal / total).** Answers: *"how strongly does
  construct i, as defined at granularity a, correlate with construct j at
  granularity b?"* — a total association, mixing a_i's direct overlap with b_j and
  the indirect overlap routed through a_i's correlated within-level siblings
  (the off-diagonal of Φ_sa).
- **`B_ab[i, j]` (pattern / partial / unique).** Answers: *"how does b_j decompose
  over the level-a construct set — how much of ancestor i lives in descendant j,
  net of the other ancestors?"* — the lineage question. Likewise
  `R²(b_j | level a) = E_ab[, j]' Φ_sa^{-1} E_ab[, j]`, which reduces to
  `Σ_i E_ab[i, j]²` only when Φ_sa = I.

Under varimax the two readings **coincide** (Φ_sa = I ⇒ E = B), which is exactly
what licenses the bass-ackwards practice of drawing the marginal correlations as
quasi-path coefficients (Goldberg's own framing, with his own part-whole caveat:
goldberg2006 p. 350 "used as path coefficients", pp. 356–357 caveats). Under
oblique they come apart: `E_ab` remains a perfectly interpretable correlation
between well-defined constructs, but reading the *set* of marginal edges as a
parent–child diagram over-reads — a child can show a strong edge to an ancestor
purely because that ancestor is correlated with the child's actual parent.

**Conclusion Q1.** The claim as stated in DESIGN §9 (`cairn/DESIGN.md:418`) —
that only orthogonal rotation produces *interpretable* between-level correlations —
is **too strong**. The oblique edge is a different, still-interpretable quantity
answering the total-association question. What is genuinely confounded under
oblique is the **lineage reading** of those edges (primary-parent assignment,
split narratives, additive variance partitioning, path-coefficient display), i.e.
the identification of marginal correlation with unique contribution. RR01 Q3's
"mixing" framing was correct as far as it went; the further step — "therefore
uninterpretable" — does not follow, and RR01 was never asked it (brief, item 3).

## Q2. Is there a decomposition?

**Yes, and it is identified from what the package already stores.**

- The **within-level score correlation** `Φ_sx = D_x^{-1/2}(W_x' R W_x) D_x^{-1/2}`
  is computable from the stored weights and `R` alone (both always retained: IP3;
  `cairn/DESIGN.md:305-327`). The model-level Φ is additionally already reserved in
  the level contract as `factor_cor` (`cairn/DESIGN.md:203`), populated as `diag(k)`
  by PCA/EFA (`R/engine_pca.R:58`, `R/engine_efa.R:143`) and extracted for real from
  lavaan by ESEM (`R/engine_esem.R:294-307`).
- The **lineage component** is `B_ab = Φ_sa^{-1} E_ab` (equivalently, per-edge
  semipartial correlations if a correlation-metric version is preferred); the
  **within-level-overlap component** is the difference `E_ab − B_ab =
  (Φ_sa − I) B_ab`, the part of each marginal edge routed through correlated
  siblings. `Φ_sa` is invertible whenever the level is non-degenerate (the same
  condition the rank-deficiency guard in `.tenBerge_weights()` polices,
  `R/engine_efa.R:177-194`).
- An equivalent change-of-basis form exists: symmetric orthogonalization
  `S_a* = S_a Φ_sa^{-1/2}` gives an orthonormal frame *within* each level, and
  edges between starred frames are `Φ_sa^{-1/2} E_ab Φ_sb^{-1/2}`. It is identified
  from the same objects, but it renames the constructs (each starred factor is a
  blend), so the regression/semipartial form is the one that keeps the level-b
  constructs as the user defined them.

**Two honest caveats.** (i) The decomposition is *statistical*, not causal: "unique
given siblings" is one legitimate operationalization of lineage, the one that
reduces to the current edges under Φ = I; it is not the only conceivable one, and
which direction to condition (descendant-on-ancestors vs ancestor-on-descendants,
where Φ_sb becomes the confounder) is a design choice. (ii) For EFA/ESEM the
decomposition holds exactly at the level of *score* correlations; transporting it
to latent-factor claims inherits the rotation-independent determinacy caveat
(RR01 Q1, Waller §4) unchanged.

**Conclusion Q2.** A lineage/overlap decomposition exists, is closed-form, and
requires nothing beyond `W`, `R`, and (for the latent version) `factor_cor` — all
already stored. Oblique support would not need new stored state, only new derived
output.

## Q3. Is Forbes's position coherent?

**First half — "each level in isolation can be useful, so oblique better reflects
the intercorrelated constructs" — coherent and well-grounded.** It is standard EFA
doctrine (conceded even by Kim & Eaton in the sentence immediately before their
orthogonality claim: "While oblique rotations are typically preferable in EFA…",
kim2015 p. 1067). Psychological constructs at a given granularity are correlated;
forcing Φ = I distorts within-level simple structure (inflated cross-loadings,
variance spread). An oblique level (pattern + Φ) is a complete, standalone
measurement model, and her own paper frames the method as rotation-agnostic:
"with orthogonal or oblique rotations" (forbes2023 p. 2, main text) and footnote 1
says oblique and EFA "give closely similar results when examining large numbers of
variables (e.g., Goldberg, 1990)" with her `ExtendedBassAckwards` function "equipped
to do orthogonal or oblique PCA or EFA" (forbes2023 p. 2, fn. 1 — verified from the
page image).

**Second half — between-level correlations as a rotation-consistency/robustness
diagnostic — partially holds, but is in tension with the structural reading.**
The kernel of truth: in the typical bass-ackwards pattern (level k+1 ≈ level k plus
one new factor), near-unity between-level correlations for the persisting
constructs do indicate that the oblique criterion is landing on the same constructs
at both depths — a real stability signal, and one that forced orthogonality can
mask. But as a *diagnostic* it is confounded in exactly the Q1 way: a high oblique
edge can mean rotational consistency **or** within-level-overlap leakage (false
positive), and a low edge can mean rotational instability **or** genuine
differentiation (indistinguishable without extra information). The deeper problem
is role conflict: the same number cannot simultaneously serve as an *estimate of
structure* and a *diagnostic of the estimation procedure* — when an edge is
intermediate, the two readings license opposite responses (interpret the hierarchy
vs distrust the rotation). The tension is resolvable, and the package already owns
the resolution: **Tucker's φ is the rotation-consistency metric** (loading-space
congruence, rotation-sensitive, sample-association-free) and **r is the structural
association** — precisely the dual reporting `prune()` already does
(`R/prune.R:40-79`, D-017's report-first stance). Forbes's workflow is coherent
once the diagnostic duty is assigned to φ (or to explicit congruence checks) and
the edges keep the structural duty, marginal or partialled as Q1 requires.

Note also the revealed-preference datum cuts both ways: she endorses oblique in
correspondence, equipped her reference implementation for it, but published the
paper's examples under varimax (forbes2023 p. 2, fn. 1; the package reproduces the
AMH varimax results to 1.3e-14). Her "conceptually it can still work out" is,
on this review's Q1/Q2 analysis, correct — with the lineage overlays re-derived,
not inherited.

## Q4. Citation audit

**The attribution to Goldberg (2006) is not defensible.** Verified against page
images:

- **goldberg2006 (p. 356)**: "The author's **preference** for examining orthogonal
  factors at each hierarchical level, instead of using an oblique rotational
  algorithm such as promax, is **not without its critics, and many readers may opt
  for the oblique alternative**. However, orthogonal factor scores have the
  advantage of **parsimony when used in multiple-regression analyses** to predict
  any important criteria (Goldberg, 1993), and they **encourage the development of
  factor markers that are maximally unrelated** to each other (Saucier, 2002)."
  A stated preference with two grounds — regression parsimony and marker
  separation — **neither of which is the between-level-confounding argument**.
- **goldberg2006 (p. 348 — the brief cites p. 349, but the passage sits in the
  Introduction on p. 348)**: his one rejection of an oblique-derived route targets
  the oblique → Schmid–Leiman orthogonalized hierarchy: "what exactly can one make
  of factors that are independent of all factors at other levels?" — a different
  target, already out of scope here (D-007/D-002 territory; DESIGN §2). Nothing on
  pp. 348–350 argues that oblique bass-ackwards edges are uninterpretable; p. 350
  even accommodates readers "whose theoretical propensities favor principal
  factors."
- **kim2015 (p. 1067, left column, first paragraph)**: "While oblique rotations are
  typically preferable in EFA, **only orthogonal rotations produce interpretable
  between-level factor score correlations in the bass-ackwards method, and
  orthogonal rotations are therefore suggested by Goldberg (2006)**." This is the
  DESIGN §9 sentence's source, near-verbatim. It is **stated bare**: no derivation,
  no citation other than the (misplaced) attribution to Goldberg, and nothing
  elsewhere in the paper supports it (p. 1066 contains only the
  EFA-vs-PCA-divergence setup and the statement that they follow Goldberg's
  methodology). The citation chain is circular: DESIGN cites Kim & Eaton and
  Goldberg; Kim & Eaton cite Goldberg; Goldberg made a preference claim on other
  grounds.
- **forbes2023 (p. 2 + fn. 1)**: the reference literature is not even unanimous —
  the paper this package treats as its fidelity contract frames the method as
  supporting "orthogonal or oblique rotations" in the main text.

**Should the rationale be re-sourced to Kim & Eaton alone?** No — re-sourcing a
bare assertion does not repair it. The right fix is to state the rationale in the
package's own voice (the Q1 Φ = I ⇒ marginal-equals-unique argument, which RR01 Q3
already identified as the load-bearing reason and which is *correct as an argument
for the default*), cite Kim & Eaton as the origin of the strong claim, and cite
Goldberg only for what he actually said.

**Recommended replacement wording for the §9 rationale (valid whatever Q6's
outcome, per the brief):**

> **Why varimax is the default.** Varimax keeps the within-level factors
> uncorrelated (Φ = I), so each between-level edge equals the unique contribution
> of that ancestor to that descendant — the marginal correlation and the
> Φ-partialled regression coefficient coincide, which is what licenses reading the
> edges as a lineage diagram (Goldberg 2006 reads them as path coefficients, with a
> part-whole caveat, p. 350). Under an oblique rotation the two quantities come
> apart: an edge is then a *total* correlation that also carries within-level
> factor overlap, and the lineage overlays (primary parent, split narrative,
> additive variance partitioning) no longer read off the raw edges. Varimax also
> matches the published analyses this package reproduces (Goldberg 2006; Kim &
> Eaton 2015; Forbes 2023's examples; Forbush et al. 2024). Kim & Eaton (2015,
> p. 1067) state the stronger claim that *only* orthogonal rotations produce
> interpretable between-level correlations, attributing it to Goldberg (2006);
> Goldberg in fact offered orthogonality as a preference on other grounds
> (regression parsimony, marker separation; p. 356) and Forbes (2023, p. 2)
> frames the method as supporting oblique rotation. The varimax criterion itself
> originates with Kaiser (1958); CF(κ = 1/p) ≡ varimax (Crawford & Ferguson 1970;
> Browne 2001).

(If D-002 stands unchanged, append the existing "Oblique rotation is out of scope
(D-002)" sentence; if superseded per Q6, append the loud-advisory framing instead.
The M76/RR01 correction parenthetical about algebra-exactness should be retained
either way.)

## Q5. Blast radius if oblique were supported

Surveyed against the current source. Legend: **unaffected** / **needs
re-derivation** / **breaks silently** (wrong numbers without erroring).

1. **`W'RW` path + standardization — unaffected.** `R/compute_edges.R:95-101`
   divides by `sqrt(.score_var(W, R))` unconditionally (`R/utils.R:104-106`); exact
   for any fixed linear `W` (RR01, settled). The scores route
   (`R/compute_edges.R:102-112`, `.compute_scores()` `R/utils.R:324-365`) divides by
   the same real SDs. IP1/IP2 survive untouched.
2. **Engine rotation sites — mechanical plumbing.** Varimax is hardcoded at
   `R/engine_pca.R:28`, `R/engine_efa.R:16`, `R/engine_esem.R:74`, and stamped into
   the object at `R/ackwards.R:1023`. lavaan::efa has built-in oblique rotations
   (geomin etc.); psych's oblique rotations require **GPArotation**, removed in M21
   — it would need to return as a *guarded Suggests* (permissible under D-011/GP5;
   it must not enter Imports).
3. **`factor_cor = I` assumptions — breaks silently, three ways.**
   (a) PCA/EFA hardcode `factor_cor = diag(k)` (`R/engine_pca.R:58`,
   `R/engine_efa.R:143`) — must extract the real Φ (psych's `$Phi`).
   (b) ESEM extracts real Φ but does **not** permute it by the variance-sort `ord`
   (`R/engine_esem.R:193-197` — the documented Known-limitations guard,
   `cairn/DESIGN.md:580-584`): under oblique the stored Φ would be silently
   misaligned with the loading columns.
   (c) Sign alignment flips loadings (`R/ackwards.R:922-923`) and weights
   (`R/ackwards.R:925-929` via `flip_weights`, `R/utils.R:307-309`) but **never
   `factor_cor`** — a flip of factor j must also flip row/column j of Φ. Invariant
   for Φ = I; silently sign-wrong under oblique.
4. **tenBerge weights — needs re-derivation; would break silently if reused.**
   `.tenBerge_weights()` implements the orthogonal formula
   `W = R^{-1} L (L' R^{-1} L)^{-1/2}` (`R/engine_efa.R:159-204`). Fed an oblique
   *pattern* matrix, it still yields a valid linear `W` — and `compute_edges()`
   would still return genuine correlations of *those* composites — but they would
   not be the correlation-preserving oblique tenBerge scores (ten Berge, Krijnen,
   Wansbeek & Shapiro 1999 give the oblique variant, which carries Φ into the
   scores; cf. DESIGN §9 scores row). Edges would be right-for-the-wrong-scores:
   the worst failure mode. Same for the ESEM regression fallback
   `W = R^{-1} L` (`R/engine_esem.R:219-226`), which omits the Φ factor of the
   oblique regression rule `W = R^{-1} Λ Φ` (cf. `cairn/DESIGN.md:230-231`).
5. **Variance explained — breaks silently.** `.variance_explained()` computes
   `colSums(L^2)/p` (`R/utils.R:97-100`) — valid only for orthogonal factors. Under
   oblique the SS-pattern-loadings are not variance contributions and do not sum
   meaningfully; the per-level `variance` vector, `tidy(what = "variance")`, and
   the ESEM variance-sort key (`R/engine_esem.R:200`, which also fixes the `ord`
   permutation) would all be silently wrong.
6. **Sign alignment to primary parent (D-010, IP4) — mechanics survive.** One sign
   DoF per factor is rotation-independent (RR01 Q5c); `.align_signs()`
   (`R/utils.R:275-304`) needs no change beyond item 3(c)'s Φ flip.
7. **Primary-parent greedy argmax (D-009) — needs re-derivation of the
   criterion, not the code.** `match_parents()` (`R/utils.R:230-239`) takes the
   column argmax of |marginal r|. Under oblique the marginal argmax can assign a
   child to whichever ancestor best *correlates with the true parent* (Q1); the
   lineage question is answered by the argmax of `|B| = |Φ_sa^{-1} E|` instead.
   Whether primary-parent should switch to the partialled criterion under oblique
   is a design decision (it changes the diagram, ordering, and sign anchoring) —
   flagged for the Forbes session.
8. **Redundancy chase (D-017) — threshold semantics shift silently.** The chase
   code (`.strong_links_direct()` `R/prune.R:232-277`; near band
   `R/prune.R:105-159`) runs unchanged, but `|r| ≥ .9` was calibrated on varimax
   edges (Forbes's published practice); oblique marginal edges are generically
   inflated by within-level overlap, so the chase would over-flag. Needs either a
   partialled criterion, a documented threshold reinterpretation, or an oblique
   advisory. Same applies to `cut_show = 0.3` (`R/compute_edges.R:50`), which
   inherits Goldberg's ≥ .30 display convention (goldberg2006 p. 351) — a
   display-only, lower-stakes version of the same shift.
9. **Tucker's φ filter — unaffected mechanically, one doc nuance.** `.tucker_phi()`
   (`R/utils.R:113-119`) on pattern columns remains the standard congruence usage
   (Lorenzo-Seva & ten Berge 2006); under oblique it should be documented as
   pattern-matrix congruence. Its role arguably *grows* under oblique (Q3: φ is the
   clean rotation-consistency metric).
10. **`redundancy_phi` auto-resolve (D-008 rule; `R/prune.R:869-889`) —
    unaffected.** Its rationale is determinacy (PCA scores determinate, EFA/ESEM
    not), which is rotation-independent — re-attaching it to rotation would
    recommit the M43 conflation (RR01 Q5c). Oblique PCA components remain
    determinate; the PCA → no-φ rule stands.
11. **ESEM engine rotation handling — plumbing plus items 3(b)/4.** lavaan side is
    the easy one (native oblique rotations, real Φ already extracted); the `ord`
    permutation and the tenBerge/regression weight formulas are the substantive
    changes.
12. **Downstream verbs.** `boot_edges()` and `comparability()` refit through the
    same engines and `compute_edges()` (D-022/D-023) — they inherit correctness
    once the engines are fixed, but comparability's greedy max-|r| matching and the
    bootstrap's full-sample anchoring carry the same marginal-vs-partial caveat as
    item 7. `ba_layout()`/`autoplot()` are edge-driven and mechanically unaffected;
    the *meaning* of what they draw is item 7/8's question. `augment()`/`predict()`
    are weight-driven and follow item 4.
13. **Forbes fidelity correspondence — constrain, don't break.** The M44
    correspondence relies on `W'RW = I` for varimax PCA making her unstandardized
    `comp.corr` equal our standardized edges (`cairn/references/forbes2023.md:62-65`;
    RR01 Q5a). Under oblique `D ≠ I`: any future fidelity test against her
    `ExtendedBassAckwards` oblique branch must reconcile her unstandardized
    products with our standardized edges explicitly. IP9 is untouched — varimax
    reproduction settings stay available.

**Silent-wrong-numbers summary (the flag the brief asked for):** items 3(a–c), 4,
5, and 8 would all produce plausible-looking wrong output rather than erroring if
oblique were enabled by plumbing alone. None of them errors today because the
`factor_cor = I` invariant (D-002's consequence line, `cairn/DECISIONS.md:25`)
makes them all vacuously correct — which is precisely why D-002 must not be
superseded by a rotation argument alone, without the accompanying kit.

## Q6. Recommendation on the decision itself

**Supersede D-002: keep varimax as the sole default, allow oblique as a
documented, loudly-advisory, non-default option — with the user-facing lineage
semantics settled at the Forbes design session before implementation.** Argued
against D-002's recorded rationale, as the brief requires:

- D-002's rationale had two clauses. The algebra clause was already struck by
  RR01/M76. The surviving interpretive clause — "oblique confounds the cross-level
  signal that is the method's whole point" — is, per Q1, an overstatement: the
  oblique edge is an interpretable total correlation; what it confounds is the
  *lineage overlay*, and Q2 shows the overlay is recoverable by a closed-form,
  already-identified partialling. The rationale's citation basis (Q4) is a bare
  assertion in Kim & Eaton misattributed to Goldberg, against a reference
  literature whose founding author calls orthogonality a preference (goldberg2006
  p. 356) and whose fidelity-contract paper declares the method
  orthogonal-or-oblique (forbes2023 p. 2).
- GP1's published-method bar is met: oblique bass-ackwards is not a
  package-invented convention — Forbes's reference implementation is "equipped to
  do orthogonal or oblique PCA or EFA" (forbes2023 p. 2, fn. 1), giving both the
  published lineage and a fidelity oracle for the new path (IP8).
- The method's author dissents from the hard exclusion in practice, and the
  package is her cited reference implementation — a standing wall against the
  rotation family she herself uses needs a better foundation than a broken
  citation chain. That foundation does not exist (Q4).

**What must ship alongside oblique for the result to remain honest** (Q5 items in
parentheses):

1. Real Φ end-to-end: extracted (3a), `ord`-permuted (3b), sign-flipped (3c),
   and *reported* — in `summary()` and a tidier — never silently carried as I.
2. Correct oblique scoring: the ten Berge et al. (1999) correlation-preserving
   oblique weights (and `R^{-1}ΛΦ` for the regression fallback) (4); the
   algebra-vs-scores oracle (IP2) extended to an oblique case.
3. Both edge readings available: marginal `E` (default display, labeled as total
   correlations) and Φ-partialled `B`/semipartials (the lineage quantities), never
   conflated in output labels. Which one drives primary-parent matching and the
   diagram is the headline design-session question (7).
4. A loud advisory (IP6) at fit time under oblique stating what an edge means and
   that redundancy/display thresholds were calibrated under varimax (8).
5. Variance-explained corrected for oblique (5).
6. Redundancy-chase stance under oblique decided with Forbes (partialled criterion
   vs documented reinterpretation) (8); `redundancy_phi` auto-resolve unchanged (10).
7. Fidelity test against Forbes's oblique `ExtendedBassAckwards` branch, with the
   standardization reconciliation made explicit (13); varimax defaults untouched,
   IP9 settings untouched.
8. GPArotation only ever as a guarded Suggests (2); no Imports change (D-011).

**What would change my mind** (if the recommendation were instead "stand"): a
demonstration that the Φ-partialled edges are not stably estimable in realistic
bass-ackwards regimes (e.g., ill-conditioned Φ_sa at deep levels making `B`
noise-dominated), or a design-session outcome where Forbes herself concludes the
oblique hierarchy overlays cannot be given publishable semantics. Either would
justify keeping oblique out on *engineering* or *design* grounds — honest grounds,
unlike the current recorded ones.

Independent of all of the above: **correct D-002's recorded rationale and the
DESIGN §9/goldberg2006.md citations now** (Q4), whatever happens to the decision —
the append-a-parenthetical pattern M76 used for RR01 applies directly.

## Q7. What should users be told?

Proposed vignette substance (~230 words; sentences flagged ◆ depend on Q6's
landing):

> **Orthogonal or oblique?** The bass-ackwards edges are correlations between
> factor scores at different depths, and the rotation you choose decides what a
> single edge means. Under varimax, the factors within each level are uncorrelated,
> so an edge from an ancestor to a descendant carries only their direct
> relationship — the correlation *is* the ancestor's unique contribution, which is
> what lets the edges be drawn as a lineage diagram and summed, squared, into
> variance accounted for. Under an oblique rotation, each level's factors are
> allowed to correlate — usually a more realistic measurement model for that level
> taken on its own — but a between-level edge then becomes a *total* correlation:
> it mixes the direct ancestor–descendant relationship with overlap routed through
> the ancestor's correlated neighbours. That quantity is meaningful (it answers
> "how strongly do these two constructs correlate?"), but it is no longer, by
> itself, a lineage statement, and thresholds tuned to varimax edges — redundancy
> chases, display cutoffs — read differently against it. ◆ `ackwards` therefore
> fits varimax by default and treats it as the only rotation. ◆ How to decide:
> if your question is the hierarchy — what splits into what, and how strongly —
> keep varimax; if your question is a single level's structure on its own terms,
> an oblique EFA of that level *outside* the hierarchy answers it without
> re-purposing the bass-ackwards edges.

If Q6 lands on oblique-as-option, the two ◆ sentences become: "`ackwards` fits
varimax by default; `rotation = "<oblique>"` is available, reports the
within-level factor correlations alongside the edges, and labels marginal and
partialled edges separately" and "…an oblique fit answers it — read its edges as
total correlations, and use the partialled edges for lineage claims." Every other
sentence is true under either outcome.

## Beyond the brief

1. **The brief's own page cite is off by one:** Goldberg's Schmid–Leiman rejection
   is on **p. 348** (Introduction), not p. 349 (which is Fig. 1 and the unrotated
   discussion). Verified from page images of both pages.
2. **`cairn/references/goldberg2006.md:37-39` carries the same over-attribution**
   the brief indicts in DESIGN §9: "Orthogonal (varimax) rotation at every level —
   matches ackwards' default (DESIGN §9; only orthogonal rotations make
   between-level score correlations cleanly interpretable)" — presented under "Key
   methodological positions taken in the paper," which the paper does not take.
   Should be corrected in the same pass as the §9 rewrite.
3. **`cairn/references/forbes2023.md` does not record footnote 1** (the
   orthogonal-or-oblique framing and the oblique-equipped reference
   implementation) — the single most decision-relevant fact in her paper for this
   question. Worth a "Rotation stance" bullet with the fn. 1 quote and page.
4. **`factor_cor` is produced by every engine but consumed by nothing** (grep of
   `R/`: no reader outside the engines). Under D-002 it is a vacuous slot; under
   oblique it becomes the single most load-bearing new output. The §4 contract
   already reserving it (`cairn/DESIGN.md:203`) means oblique support changes no
   object shape — a genuinely small structural footprint hiding the real work
   (Q5 items 3–5).
5. **Kim & Eaton's concessive clause is a usable doc asset:** "While oblique
   rotations are typically preferable in EFA…" (p. 1067) — quoting it alongside
   their strong claim inoculates the vignette against the same over-reading the
   package inherited.
6. **D-033 is a precedent in miniature** for Q6's shape: a Forbes-correspondence-
   driven change where the published convention became the default and the more
   permissive behaviour stayed available behind an argument — the same
   default-vs-capability split recommended here (and the same IP6 loud-NEWS
   handling, pre-1.0).

## Recommendations

1. **Apply now (no milestone needed beyond a tracking/doc pass):** rewrite the
   DESIGN §9 rotation-row rationale per Q4's replacement text; append an RR02
   correction parenthetical to D-002's Context (M76/RR01 house pattern — do not
   rewrite the entry); fix `goldberg2006.md:37-39`; add the fn. 1 rotation-stance
   bullet to `forbes2023.md`.
2. **Apply (decision):** supersede D-002 per Q6 — varimax-only *default*, oblique
   as documented non-default, implementation gated on the Forbes design session
   resolving the marginal-vs-partialled lineage semantics (Q5 item 7, Q6 item 3)
   — recorded as a new D-entry; register the implementation as a ROADMAP candidate
   carrying this report's Binding criteria.
3. **Apply (docs):** the Q7 vignette passage, with the ◆ sentences set by
   recommendation 2's outcome.
4. **Consider:** shipping the Φ-partialled edge decomposition (`B = Φ_s^{-1} E`,
   `R²` per node) as a *reporting* enhancement even before/without oblique
   rotation — under varimax it is an identity check (B ≡ E), and its presence
   makes the eventual oblique semantics legible rather than novel.
5. **Consider:** having the design session also decide whether Tucker's φ is
   formally promoted to the rotation-consistency diagnostic in docs (Q3's
   division of labor), so Forbes's robustness reading gets a sanctioned home that
   is not the score edges.
6. **Reject — restating the exclusion with better citations** (i.e., keeping
   D-002 but re-sourcing the same "uninterpretable" claim to Kim & Eaton alone):
   Q4 shows the claim is asserted bare there; re-sourcing would launder, not
   repair, it. If the decision is retained unchanged, retain it on the honest
   grounds (default fidelity to the published workflows + the unbuilt Q5 kit),
   not on the interpretability claim.
7. **Reject — enabling oblique by plumbing alone** (rotation argument passed
   through to engines without Q6's accompaniments): Q5 items 3–5 and 8 would
   produce silently wrong Φ, scores, variance, and redundancy flags. This is the
   one implementation shape this review rules out entirely.

## Binding criteria

For any future milestone implementing oblique rotation support (recommendation 2;
these do not bind the doc-only recommendations 1 and 3):

- **BC1.** Varimax remains the default rotation for every engine; a fit with all
  defaults is byte-identical in its numerical outputs to the pre-change package,
  and the IP9 Forbes-reproduction settings and fidelity suite pass unchanged.
- **BC2.** Under an oblique rotation, each level's stored `factor_cor` is the real
  within-level factor correlation — extracted from the engine, permuted by any
  column reordering (the `engine_esem.R` `ord` guard), and sign-flipped in lockstep
  with `align_signs` — and is surfaced in at least `summary()` and one tidier.
  No code path carries `factor_cor = I` for an oblique fit.
- **BC3.** Oblique scoring weights are the correlation-preserving oblique variant
  (ten Berge et al. 1999) on the default path, with the oblique regression rule
  (`R^{-1}ΛΦ`) as the fallback; the orthogonal tenBerge formula is never applied to
  an oblique pattern matrix. The IP2 algebra-vs-scores oracle gains at least one
  oblique test case.
- **BC4.** `variance` (and everything downstream of `.variance_explained()`) is
  computed by an oblique-valid formula when Φ ≠ I; the ESEM variance sort and the
  Φ permutation use the same ordering.
- **BC5.** Output distinguishes marginal edges (`E`) from Φ-partialled lineage
  quantities (`B` or semipartials); no output surface labels an oblique marginal
  edge as a unique/lineage contribution. Which quantity drives primary-parent
  matching, sign anchoring, and the diagram is fixed by a design-session D-entry
  before implementation, not defaulted in code.
- **BC6.** An oblique fit announces via cli (IP6) that edges are total
  correlations and that `redundancy_r`/`cut_show` conventions were calibrated
  under varimax; `prune()` under an oblique object either implements the
  design-session redundancy stance or warns that the default criterion assumes
  orthogonal levels.
- **BC7.** The oblique path is oracle-backed (IP8) including a fidelity test
  against Forbes's `ExtendedBassAckwards` oblique branch with the
  standardization correspondence (`D ≠ I`) handled explicitly.
- **BC8.** No new Imports: any rotation dependency (e.g. GPArotation for psych's
  oblique criteria) enters as a guarded Suggests only (D-011/GP5).
