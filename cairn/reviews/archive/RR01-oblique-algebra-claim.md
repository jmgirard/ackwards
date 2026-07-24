# RR01: Does orthogonal rotation make the W'RW edge algebra "exact"? — Review Report

- **Date:** 2026-07-23
- **Brief:** `cairn/reviews/RB01-oblique-algebra-claim.md` (M76)
- **Reviewer:** independent Fable-level review, per `/milestone-brief`
- **Materials read:** `cairn/DESIGN.md` §5.1 + §9 (rotation, scores, `redundancy_phi` rows);
  `cairn/DECISIONS.md` D-002; `vignettes/ackwards-intro.Rmd.orig:110-144`;
  `R/ackwards.R:1-40`; `cairn/references/waller2007.md`; `cairn/references/forbes2023.md`;
  cross-checked `R/compute_edges.R` (lines 95-101) and `.score_var()` (`R/utils.R:104-106`).

## Verdict in one line

The concern is correct: the `W'RW` identity is algebraically exact for **any fixed linear
scoring map**, oblique included — orthogonality is an interpretive choice (Φ = I), not a
numerical prerequisite for the closed form. Both quoted passages overstate orthogonality's
role; the vignette passage is flatly wrong, the roxygen passage half-wrong. The proposed
replacements are accurate with two small nuances flagged under Q4/Q5.

## Q1. Exactness under oblique — the algebra

Let `Z` be the n×p matrix of standardized observed variables (columns mean 0, sd 1), so
that the sample covariance of `Z` **is** the correlation matrix: `cov(Z) = R`. Let scores
at level x be a fixed linear map `S_x = Z W_x` (p×k_x weight matrix `W_x`). Then by
bilinearity of covariance alone:

```
cov(S_a, S_b) = cov(Z W_a, Z W_b) = W_a' cov(Z) W_b = W_a' R W_b
var(S_x[, j]) = [W_x' R W_x]_{jj}       =>  D_x = diag(W_x' R W_x)
cor(S_a, S_b) = D_a^{-1/2} (W_a' R W_b) D_b^{-1/2}
```

Every step uses only (i) `S = Z W` is linear with fixed `W`, and (ii) `cov(Z) = R`. **No
property of `W` is used anywhere** — not orthonormality of a rotation, not `T' = T^{-1}`,
not `Φ = I`, not unit score variance. Rotation enters only through the *construction* of
`W`, and any rotation (orthogonal or oblique) composed with any linear scoring rule
(regression `W = R^{-1}ΛΦ`, Bartlett, tenBerge and its oblique correlation-preserving
variant) still yields a fixed linear `W`. Therefore the closed form remains **exact** under
oblique rotation; obliqueness introduces **no approximation whatsoever**.

What obliqueness *does* change is `D_x`. For varimax-rotated standardized PC scores,
`W = Q_k Λ_k^{-1/2} T` with `T'T = I` gives

```
W' R W = T' Λ_k^{-1/2} (Q_k' R Q_k) Λ_k^{-1/2} T = T' Λ_k^{-1/2} Λ_k Λ_k^{-1/2} T = T'T = I,
```

so `D = I` and the raw product `W_a' R W_b` is already a correlation matrix. Under oblique
rotation `D ≠ I` in general and the standardization step is mandatory — which is exactly
the "standardization is the trap" note in §5.1, and exactly what `compute_edges()` does
unconditionally (`R/compute_edges.R:98-101` divides by `sqrt(.score_var(W, R))` on every
algebra-path call; `.score_var()` is `diag(crossprod(W, R %*% W))` with no orthogonality
assumption). So orthogonality buys a *convenience* (D = I for PCA), never *exactness*.

Independent confirmation from the primary source: Waller (2007) himself extends the
identity to oblique rotations in his §3 (`T_i^{-1} S T_j'^{-1}`; recorded in
`cairn/references/waller2007.md:37-38`). The primary source for the orthogonal derivation
thus *also* supplies the oblique closed form — the claim that orthogonality is what makes
the closed form exist is contradicted by the very paper cited for the algebra.

The only genuine "exact vs. approximate" distinction in this neighborhood is Waller's §4
caveat: for the **factor** (not component) model, correlations among *estimated* factor
scores are not the correlations among the *latent factors* (indeterminacy; Grice 2001).
That is an **engine/determinacy** issue, identical under orthogonal and oblique rotation —
precisely the conflation M43 already corrected once in the `redundancy_phi` row.

**Conclusion Q1:** Exact for any linear `W`, oblique included. §5.1 as written is correct;
no flag against it.

## Q2. Accuracy of the two existing passages

**Vignette passage** (`ackwards-intro.Rmd.orig:124-127`): "Orthogonality is what gives the
between-level score correlations a closed form … T′ = T⁻¹ … is exactly what makes the
`W′RW` edge algebra exact rather than an approximation." — **Wrong**, not merely
overstated. It asserts two false propositions: (a) that the closed form exists *because of*
orthogonality (Waller §3 gives the oblique closed form; the Q1 derivation needs no rotation
property at all), and (b) that without orthogonality the algebra would be "an
approximation" (it would not; nothing in the oblique case is approximate). It also
implies the exact/approximate boundary runs along the orthogonal/oblique axis, when it
actually runs along the component/estimated-factor-score axis (Waller §4) — the same
conflation M43 fixed in DESIGN §9. The passage's *last* sentence ("correlated factors
would confound the cross-level signal") is correct and is the real rationale.

**Roxygen passage** (`R/ackwards.R:11-17`): "the `T'=T^-1` property … **enables** the
closed-form `W'RW` edge algebra **and** keeps within-level factors uncorrelated…" — a
conjunction of a wrong clause and a right clause. "Enables" is milder than the vignette's
"exact rather than an approximation" but still asserts orthogonality as a prerequisite of
the closed form, which is false. The second clause (Φ = I so cross-level edges reflect
only the hierarchical signal) is correct and sufficient on its own.

**Conclusion Q2:** Vignette — wrong. Roxygen — half-wrong (the "enables" clause);
its interpretive half is accurate and should survive the rewrite.

## Q3. The actual reason for varimax-only

Yes, the interpretive framing is the correct justification, with one precision note.

The load-bearing reason for D-002 is: under varimax, the within-level factor correlation
is Φ = I, so a between-level edge is not contaminated by within-level factor overlap;
under oblique rotation Φ ≠ I, and the within-level intercorrelation propagates into the
cross-level score correlations, mixing "how much of ancestor i lives in descendant j"
with "how much descendants overlap each other" — confounding the exact question
bass-ackwards exists to answer. A secondary, also legitimate reason is **literature
fidelity**: Goldberg (2006), Kim & Eaton (2015), Forbush et al. (2024), and Forbes (2023)
all use varimax (CF(κ=1/p) ≡ varimax), so varimax-only keeps results commensurable with
the entire reference literature and the fidelity suite. Neither reason is numerical.

**Precision note:** "varimax keeps the within-level factors uncorrelated" is exact at the
*factor/component* level (Φ = I). Whether the within-level **scores** are uncorrelated
additionally depends on the scoring rule: PCA component scores and tenBerge scores
(correlation-preserving by construction — the stated reason tenBerge is the default,
DESIGN §9 scores row) carry Φ = I through to the scores; regression or Bartlett scores of
orthogonal factors generally have a non-diagonal score correlation (`Λ'R^{-1}Λ` etc. need
not be diagonal). The shipped defaults make the interpretive story airtight; prose should
attribute uncorrelatedness to the **factors (Φ)** — as both proposed replacements already
do — rather than unconditionally to estimated scores.

## Q4. Proposed replacement text

Both replacements are accurate in substance and correctly relocate the rationale from
numerical to interpretive. Findings:

- "exact for any linear scoring, oblique included" — correct; matches §5.1 and the Q1
  derivation. Consider "any *fixed* linear scoring" for precision; not required.
- "Varimax leaves the factors within a level uncorrelated (Φ = I)" — correct as stated
  (factor-level claim; see Q3 precision note). Do not strengthen it to a claim about
  arbitrary estimated scores.
- "leaks into every between-level edge" (vignette) / "would leak into every edge"
  (roxygen) — "every" is a strong universal; in degenerate configurations a particular
  edge can be numerically unaffected. Generically true and rhetorically fine for a
  teaching aside; softening to "into the between-level edges" would be strictly safer.
  Minor — not a blocker either way.
- Neither replacement claims the algebra recovers *latent factor* correlations for EFA —
  good; do not add such a claim (Waller §4 caveat is rotation-independent and engine-level).
- No overcorrection detected: both texts immediately follow the "exact for oblique too"
  concession with the interpretive reason, so neither reads as an invitation to support
  oblique (D-002 untouched, per the brief's constraint).
- The current vignette aside's surrounding content — the Goldberg lineage sentence, the
  CF-VARIMAX ≡ varimax equivalence, and the closing "correlated factors would confound
  the cross-level signal" — is correct and should be preserved around the replacement.

**Conclusion Q4:** Apply both, optionally with the two micro-tweaks above.

## Q5. Legitimate numerical roles for orthogonality

Checked the three candidates; none contradicts the "interpretive, not numerical" framing,
but one narrow fact must not be erased:

(a) **Score-variance standardization `D_x`.** Does *not* rely on orthogonality — the code
computes `diag(W'RW)` unconditionally and §5.1 says "NOT assumed 1". However,
orthogonality is exactly what makes `D = I` for varimax PCA (algebra in Q1), and that fact
is load-bearing **for the Forbes correspondence**: Forbes's `comp.corr = t(W_a) R W_b` is
*unstandardized*, and it equals our standardized edges entrywise only because varimax PCA
gives `W'RW = I` (recorded in `cairn/references/forbes2023.md`, "Correspondence
conventions"). This is a numerical role for orthogonality in the *fidelity
correspondence*, not in the package's own algebra. The corrected prose need not mention
it, but the rewrite must not delete or contradict the `forbes2023.md` correspondence note
or the §5.1 standardization-trap note.

(b) **Waller's derivation.** His main derivation is orthogonal, but his own §3 supplies
the oblique closed form — so his paper is evidence *for* the correction, not a numerical
role to preserve. Citing "Waller (2007) derived this for orthogonal components, but the
identity holds for any linear W" (§5.1's existing sentence) remains the right framing.

(c) **Sign alignment / determinacy.** Column sign flips of `W` flip the corresponding
edge rows/columns identically under orthogonal and oblique linear scoring; `align_signs`
does not depend on orthogonality for correctness. Score determinacy (Grice 2001) is
engine-level (PCA determinate, EFA not), rotation-independent — re-attaching it to
rotation would recommit the M43 conflation.

**Conclusion Q5:** Nothing the "interpretive, not numerical" framing wrongly erases,
provided the rewrite (i) keeps the standardize-by-real-SDs trap note intact and (ii)
leaves the Forbes `W'RW = I` correspondence fact in `forbes2023.md` untouched.

## Beyond the brief

1. **DESIGN §9 rotation row carries the same conflation.** The `rotation` row's rationale
   reads "T'T = I so T' = T^-1, **enabling closed-form W'RW algebra** (Waller 2007)
   without materialising scores" — directly contradicting §5.1 ("holds for any linear W")
   two hundred lines earlier in the same document. Since §9 is a living source-of-truth
   table, this should be corrected in the same pass, M43-parenthetical style (the
   `redundancy_phi` row shows the house pattern for annotating a corrected rationale).
2. **DECISIONS.md D-002 Context line has it too**: "Only orthogonal rotation (T' = T^-1)
   yields interpretable between-level score correlations **and the closed-form W'RW
   algebra**." The *decision* is fine; the recorded *context* conflates. Recommend an
   appended correction parenthetical (as M43 did) rather than rewriting the entry — the
   decision record's history should stay legible.
3. `cairn/references/waller2007.md` already records the §3 oblique extension and needs no
   change; it is the internal evidence base for this correction.
4. The two flagged passages disagree with each other in strength (vignette: "exact rather
   than an approximation"; roxygen: "enables") — a sign the claim drifted during
   paraphrase rather than being traced to a source. Post-rewrite, both should trace to
   §5.1's wording.

## Recommendations

1. **Apply** — replace the vignette aside's exactness sentences with the proposed intro
   replacement (optionally softening "every"), preserving the surrounding
   Goldberg/CF-VARIMAX/confound sentences.
2. **Apply** — replace the roxygen `rotation` bullet with the proposed roxygen replacement;
   keep the existing literature-matching sentence ("Matches Goldberg (2006)…") and the
   "varimax is the only supported rotation" sentence.
3. **Apply** — correct the DESIGN §9 `rotation` row's "enabling closed-form W'RW algebra"
   rationale in the same pass, with an M43-style correction parenthetical.
4. **Consider** — append a one-line correction parenthetical to D-002's Context in
   `cairn/DECISIONS.md` (do not rewrite the entry).
5. **Consider** — the micro-tweaks: "any *fixed* linear scoring"; "leaks into the
   between-level edges" instead of "every … edge".
6. **Reject** — any rewrite that moves the exact/approximate distinction onto the
   estimated-factor-scores axis *inside these two passages*; that caveat (Waller §4) is
   real but engine-level, belongs where determinacy is discussed (e.g. `redundancy_phi`
   docs), and importing it here would clutter a rotation rationale with a non-rotation
   issue.

## Binding criteria

- **BC1.** The rewritten vignette aside contains no claim that orthogonality (or
  `T' = T^{-1}`) is necessary for, enables, or is what makes the `W'RW` closed form exact;
  it states affirmatively that the identity is exact for any (fixed) linear scoring,
  oblique included.
- **BC2.** The rewritten roxygen `rotation` bullet contains no "enables the closed-form"
  (or equivalent necessity) claim about orthogonality, and states the interpretive Φ = I
  rationale (within-level factors uncorrelated, so between-level edges are not confounded
  by within-level factor intercorrelation).
- **BC3.** Both rewritten texts retain the substantive point that an oblique rotation's
  within-level factor correlation would contaminate/confound the between-level edges —
  the interpretive rationale must not be weakened while removing the numerical one.
- **BC4.** Neither rewritten text suggests oblique rotation is supported, supportable, or
  desirable (D-002 unchanged).
- **BC5.** Uncorrelatedness claims in the rewritten texts are attributed to the factors
  (Φ = I), not asserted unconditionally of estimated factor scores.
- **BC6.** The rewrite does not remove or contradict (i) DESIGN §5.1's "holds for any
  linear W" and "standardization is the trap" notes, or (ii) the `forbes2023.md`
  correspondence note that `W'RW = I` for varimax PCA is what makes Forbes's
  unstandardized `comp.corr` equal the standardized edges.
- **BC7.** The DESIGN §9 `rotation` row no longer asserts that `T' = T^{-1}` enables the
  closed-form algebra; the corrected row carries an M43-style parenthetical noting the
  wording correction and this review (RR01) as its source.
