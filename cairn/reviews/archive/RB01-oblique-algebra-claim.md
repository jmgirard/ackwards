# RB01: Does orthogonal rotation make the W'RW edge algebra "exact", or is that a conflation? (M76)

- **Date:** 2026-07-23
- **Output required:** write findings to `cairn/reviews/RR01-oblique-algebra-claim.md`

You are performing an independent expert review. This brief is fully
self-contained — do not assume any conversation context. Read only what this
brief directs you to read, answer the numbered questions, and write your
findings to the output path above using the same numbering.

## Background

`ackwards` is an R package implementing Goldberg's (2006) bass-ackwards method
and Forbes's (2023) extension: extract factor/component solutions from 1..k
factors, then characterize the hierarchy via **between-level factor-score
correlations**. Those between-level correlations are computed by a closed-form
identity (the "`W'RW` edge algebra") rather than by materializing scores.

Within each level, factors are rotated with **varimax** (orthogonal). Varimax
is the only supported rotation; oblique is out of scope by design decision
**D-002**. Milestone **M76** is a documentation-clarity pass. In the course of
it, two existing pieces of user-facing prose were found to justify the
varimax-only choice partly on the grounds that **orthogonality is what makes
the edge algebra exact**:

- `vignettes/ackwards-intro.Rmd.orig:123-132` (a teaching aside): "Orthogonality
  is what gives the between-level score correlations a closed form: an
  orthonormal rotation satisfies T′ = T⁻¹, which is exactly what makes the
  `W′RW` edge algebra **exact rather than an approximation**."
- `R/ackwards.R:11-17` (roxygen): "the `T'=T^-1` property of orthogonal rotation
  **enables the closed-form `W'RW` edge algebra** and keeps within-level factors
  uncorrelated so cross-level edges reflect only the hierarchical signal."

The concern: the repo's own design doc (`cairn/DESIGN.md` §5.1) states the
identity "holds for **any linear W**", and a prior correction (M43, recorded in
the §9 defaults table, `redundancy_phi` row) explicitly flagged an earlier
"`W′RW` algebra is exact" phrasing as a **conflation** of algebra-exactness
(true of any linear scoring) with a different, scoring-specific property. If the
algebra is exact for oblique linear scores too, then orthogonality is an
**interpretive** choice (it keeps within-level factors uncorrelated), not a
**numerical** prerequisite for the closed form — and the two prose passages
above overstate orthogonality's role.

This review must independently settle whether that concern is correct before the
reviewed text is rewritten, because it touches the mathematical rationale of a
standing design decision (D-002).

## Materials

Read these:

- `cairn/DESIGN.md` §5.1 "Why the algebra generalizes" (lines ~217-235) — the
  identity `cor(S_a,S_b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}`, its stated
  generality ("holds for any linear W … regression (Thurstone), Bartlett, and
  tenBerge"), and the "standardization is the trap" note.
- `cairn/DESIGN.md` §9 defaults table, the `redundancy_phi` row (search
  "redundancy_phi" / "conflated") — the M43 correction wording.
- `cairn/DECISIONS.md` D-002 (varimax only; oblique out of scope) — its stated
  context/reason.
- `vignettes/ackwards-intro.Rmd.orig:120-135` — the "Why varimax?" aside.
- `R/ackwards.R:8-17` — the `@section Defaults and why` bullets for `engine`
  and `rotation`.
- `cairn/references/waller2007.md` and `cairn/references/forbes2023.md` — the
  source summaries (Waller derived the identity for orthogonal components;
  Forbes uses `comp.corr = t(W_a) R W_b`).

Optional cross-check: `R/*.R` for `compute_edges()` to see how `W`, `R`, and
the score-variance standardization `D = diag(W'RW)` are actually used.

## Questions

1. **Exactness under oblique.** For scores that are a linear map `S = Z W` of
   the standardized observed variables, is the identity
   `cor(S_a,S_b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}` (with
   `D_x = diag(W_x' R W_x)`) **algebraically exact regardless of whether the
   within-level rotation is orthogonal or oblique**? Specifically: if an oblique
   rotation is used and scores are formed by a linear scoring rule (regression /
   Bartlett / tenBerge), does the closed form remain exact, or does obliqueness
   introduce an approximation? Show the algebra.

2. **Accuracy of the existing prose.** Given (1), are the two quoted passages
   ("orthogonality … is exactly what makes the `W′RW` edge algebra exact rather
   than an approximation"; "the `T'=T^-1` property … enables the closed-form
   `W'RW` edge algebra") **accurate, an overstatement, or wrong**? Distinguish
   the two claims if they differ in accuracy.

3. **The actual reason for varimax-only.** Is the correct justification for
   D-002 (varimax only) the **within-level orthogonality of the factors**
   (Φ = I under varimax vs Φ ≠ I under oblique), so that a between-level edge
   isolates the cross-level relationship rather than mixing in within-level
   factor intercorrelation — an **interpretive** reason — rather than any
   numerical/exactness property of the algebra? Or is that framing itself
   imprecise?

4. **Proposed replacement text.** Are the two replacements below accurate and
   complete? Flag any error, missing nuance, or overcorrection.

   *Intro vignette replacement:* "The reason is not numerical: the `W′RW`
   between-level correlation identity is exact for any linear scoring, oblique
   included. The reason is interpretive. Varimax leaves the factors within a
   level uncorrelated (Φ = I), so a between-level score correlation reflects
   only the cross-level relationship. Under an oblique rotation the within-level
   factors are themselves correlated (Φ ≠ I), and that within-level correlation
   leaks into every between-level edge — confounding the within-versus-between
   question the method exists to answer."

   *Roxygen replacement:* "varimax keeps the within-level factors uncorrelated
   (Φ = I), so each between-level edge reflects only the cross-level
   relationship; an oblique rotation's correlated within-level factors (Φ ≠ I)
   would leak into every edge and confound the between-level signal. The
   closed-form `W'RW` edge algebra is itself exact for any linear scoring,
   orthogonal or not, so the choice is interpretive, not a numerical necessity."

5. **Any legitimate numerical role for orthogonality** the "interpretive, not
   numerical" framing would wrongly erase? Consider: (a) the score-variance
   standardization `D_x` (does it silently rely on orthogonality?); (b) Waller's
   (2007) original derivation being stated for orthogonal components; (c) sign
   alignment / determinacy. If any such role exists, state precisely what the
   corrected prose must preserve.

## Constraints

- **D-002 is fixed** (varimax only; oblique out of scope). This review does
  **not** relitigate that decision — only the accuracy of the *stated
  mathematical rationale* for it. If you find the decision's real justification
  differs from both the old and proposed prose, say so; do not propose
  supporting oblique rotation.
- **IP1** (one edge path via `compute_edges()`) and the §5.1 identity as
  written are the repo's ground truth; if you believe §5.1 itself is wrong,
  flag that explicitly rather than working around it.
- The deliverable is prose accuracy for user-facing docs, not a code change to
  the algebra.

## Output format

In `RR01-oblique-algebra-claim.md`: answer each question by number with your
reasoning and evidence (show the algebra for Q1). List any additional findings
separately under "Beyond the brief"; end with concrete recommendations, each
marked apply / consider / reject-with-reason. If your findings bind the M76
rewrite (e.g., specific wording that must or must not appear), add a
`## Binding criteria` section: numbered `BC1…`, each a measurable assertion
checkable against the rewritten text, with any numeric projection stating its
tolerance. Binding criteria are ingested verbatim into M76's acceptance
criteria.
