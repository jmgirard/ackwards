# RB02: Does oblique rotation actually confound the between-level signal? (no milestone)

- **Date:** 2026-07-30
- **Output required:** write findings to `cairn/reviews/RR02-oblique-rotation-reconsidered.md`

You are performing an independent expert review. This brief is fully self-contained — do
not assume any conversation context. Read only what this brief directs you to read, answer
the numbered questions, and write your findings to the output path above using the same
numbering.

No milestone is gated on this brief. It feeds a `candidate` ROADMAP row and a design
session with Forbes. (M84 is `blocked` for an unrelated reason.)

## Background

`ackwards` is an R package implementing Goldberg's (2006) bass-ackwards method and Forbes's
(2023) extension. A hierarchy is built by extracting solutions from 1..k factors and
characterizing the structure via **between-level factor-score correlations**: for levels
`a` (shallower) and `b` (deeper), the edge is `cor(S_a, S_b)` where `S_x = Z W_x` are
standardized factor/component scores.

**Rotation is fixed at varimax and is not a user argument** (D-002, implemented M13:
hardcoded internal constant; the `kappa` argument was removed since CF(κ = 1/p) ≡ varimax).
Oblique rotation is declared out of scope. The stated rationale in `cairn/DESIGN.md` §9 is:

> **Only orthogonal rotation produces interpretable between-level factor score
> correlations** in the bass-ackwards method (Goldberg 2006; Kim & Eaton 2015): varimax
> makes the within-level factors orthogonal, so a between-level edge reflects only the
> cross-level relationship instead of being confounded by within-level factor
> intercorrelation.

Three developments have put that rationale under question, and this brief exists to test
it. None of them is a numerical claim — **RR01 (archived, 2026-07-23) already settled the
numerics** and must not be re-derived (see Constraints).

**(1) The method's author dissents in practice.** M. K. Forbes — author of Forbes (2023),
whose footnote 3 names this package as the reference implementation of the extended method
— wrote to the maintainer on 2026-07-30:

> one of the things I like about this framework is that each level in isolation can be
> useful, so I like to use oblique rotations to better reflect the intercorrelated
> constructs we're interested in. Then the correlations between levels also tell us about
> whether those rotations are consistent/robust. I see what you mean about the
> shared/overlapping variance within and between levels, but I think conceptually it can
> still work out to have the hierarchy derived with oblique rotations.

Her published analyses nonetheless used varimax: Forbes (2023) footnote 1 states she
focuses on orthogonal (varimax) PCA in the paper's examples, notes oblique rotations and
EFA are also possible, and says her supplemental `ExtendedBassAckward` function is equipped
to do orthogonal or oblique PCA or EFA. The package reproduces her published AMH component
correlations to 1.3e-14 across all 45 level-pairs, confirming varimax for the paper.

**(2) The citation chain behind DESIGN §9 is broken at one link.** The sentence above is
Kim & Eaton (2015, p. 1067) almost verbatim — they write that while oblique rotations are
typically preferable in EFA, only orthogonal rotations produce interpretable between-level
factor score correlations in the bass-ackwards method, and that orthogonal rotations are
therefore suggested by Goldberg (2006). That attribution to Goldberg does not hold up.
Goldberg (2006, p. 356) calls examining orthogonal factors at each level *the author's
preference*, states it is not without its critics and that many readers may opt for the
oblique alternative, and justifies it on two grounds unrelated to between-level
confounding: parsimony of orthogonal factor scores in multiple-regression prediction
(Goldberg 1993) and encouraging factor markers maximally unrelated to each other
(Saucier 2002). His one rejection of an oblique route (p. 349) targets the oblique →
Schmid–Leiman orthogonalized hierarchy, on the interpretive ground that factors independent
of every other level are hard to make sense of — a different target, and Schmid–Leiman is
already out of scope here (D-007). Goldberg's abstract does frame the method in terms of
correlations among *orthogonal* factor scores, so the practice is his; the specific
argument attributed to him is not.

**(3) RR01 asserted the interpretive rationale but was not asked to test it.** RR01's Q3
stated that under oblique rotation the within-level intercorrelation propagates into the
cross-level score correlations, mixing "how much of ancestor i lives in descendant j" with
"how much descendants overlap each other". That was offered as the correct *framing* for a
wording fix, not as a tested proposition, and RR01 was never asked whether the resulting
quantity is genuinely uninterpretable or merely a different interpretable quantity.

## Materials

Read these. Paths are repo-relative.

- `cairn/DESIGN.md` — §5.1 (edge algebra, standardization), §9 (the `rotation`, `scores`,
  and `redundancy_phi` rows), and the "Design principles" block (IP1–IP9, GP1–GP6).
- `cairn/DECISIONS.md` — D-002 (varimax only), D-004 (one shared edge path), D-007
  (tenBerge default, EAP out of scope), D-009 (primary-parent greedy argmax), D-010 (sign
  alignment anchors to the primary parent), D-017 (redundancy chased by direct skip-level
  correlation), D-031 (IP/GP adoption; IP9 Forbes capability).
- `cairn/reviews/archive/RR01-oblique-algebra-claim.md` — **read in full before answering.**
  Its Q1 (exactness under oblique), Q3 (the actual reason for varimax-only, plus its
  precision note on Φ vs. score correlation), and Q5 are directly upstream of this brief.
- `R/compute_edges.R` — the single edge path; note the unconditional standardization by
  `sqrt(.score_var(W, R))` (`.score_var()` in `R/utils.R`).
- `R/engine_pca.R`, `R/engine_efa.R`, `R/engine_esem.R` — where varimax is applied per
  engine, and any place `factor_cor = I` is assumed.
- `R/prune.R` — the redundancy chase and Tucker's φ congruence filter.
- `cairn/references/goldberg2006.md`, `cairn/references/forbes2023.md`,
  `cairn/references/waller2007.md` — committed source notes.
- Shelf PDFs (gitignored, local): `cairn/references/sources/goldberg2006.pdf` (p. 356 and
  p. 349), `cairn/references/sources/kim2015.pdf` (p. 1067, left column, first paragraph),
  `cairn/references/sources/forbes2023.pdf` (footnote 1). **Verify quotations by rendering
  the page image, not by `pdftotext` alone** — extraction drops inter-word spaces and
  de-hyphenates line-wrapped terms in this corpus.

## Questions

1. **Is the confounding claim true as stated?** Under oblique rotation at levels `a` and
   `b` (within-level factor correlations `Φ_a ≠ I`, `Φ_b ≠ I`), state formally what
   `cor(S_a, S_b)` estimates. Is the between-level correlation genuinely *confounded* — so
   that no lineage claim can be read off it — or is it a different, still-interpretable
   quantity that answers a different question? If the latter, say precisely which question
   each version answers.

2. **Is there a decomposition?** Can the oblique between-level correlation be separated
   into a lineage component and a within-level-overlap component (for instance by
   partialling `Φ`, or by a change of basis to an orthogonal frame within each level)? If
   such a decomposition exists, is it identified from what the package already stores
   (loadings, `Φ`, weights, `R`)? If it does not exist, say why not.

3. **Is Forbes's position coherent?** She argues that each level is useful in isolation
   (so oblique better reflects the constructs at that level), and that the between-level
   correlations then additionally indicate whether the rotations are consistent or robust.
   Assess this on its merits. In particular: does the second half — between-level
   correlations as a robustness diagnostic for the rotation — hold up, and is it
   compatible with reading the same correlations as hierarchical structure, or are the two
   readings in tension?

4. **Citation audit.** Is DESIGN §9's attribution of the orthogonality requirement to
   Goldberg (2006) defensible on the text at p. 356 and p. 349? Should the rationale be
   re-sourced to Kim & Eaton (2015) alone? Is Kim & Eaton's own assertion supported by
   argument or evidence anywhere in their paper, or is it stated bare? Recommend exact
   replacement wording for the §9 rationale that is accurate to the sources, **whatever
   the answer to Q1** — this question is about attribution, not about the decision.

5. **Blast radius if oblique were supported.** Enumerate the code-level consequences,
   citing files and lines. At minimum consider: the `W'RW` path and its standardization
   (`R/compute_edges.R`); sign alignment to the primary parent (D-010, IP4); primary-parent
   assignment by greedy argmax (D-009); the redundancy chase and its `|r| ≥ .9` threshold
   (D-017); Tucker's φ as a congruence filter; the `redundancy_phi` auto-resolve rule and
   its determinacy rationale; anywhere `factor_cor = I` is assumed; and the ESEM engine's
   rotation handling. Which of these break, which need re-derivation, and which are
   unaffected? Flag any that would silently produce wrong numbers rather than erroring.

6. **Recommendation on the decision itself.** Given Q1–Q5, should D-002 stand, be
   superseded to allow oblique as a documented non-default option, or be superseded to
   something else? If you recommend allowing oblique, state what must ship alongside it for
   the result to remain honest (for example: reporting `Φ` in the output, a loud advisory,
   a different default redundancy criterion). If you recommend it stand, state what would
   change your mind.

7. **What should users be told?** Independently of Q6, the package intends to add a
   vignette passage teaching users how to think about the orthogonal/oblique choice in this
   method. Sketch its substance in 150–300 words: what the choice actually trades off, what
   a between-level edge means under each, and how a user should decide. Write it to be true
   under whichever answer Q6 lands on, flagging any sentence that would need to change.

## Constraints

Flag disagreement with any of these explicitly rather than silently working around them.

- **RR01's algebra finding is settled and must not be re-derived.** The `W'RW` identity is
  exact for any fixed linear scoring map, oblique included; orthogonality is an interpretive
  choice (Φ = I), not a numerical prerequisite. Waller (2007) §3 gives the oblique closed
  form. Do not spend the review re-establishing this; assume it.
- **D-002 is a standing decision and stands unless superseded.** A recommendation to change
  it is in scope and welcome, but it must be argued against D-002's recorded rationale,
  which this brief quotes above — not around it.
- **Out of scope, do not propose:** Schmid–Leiman / orthogonalized hierarchical solutions
  and higher-order SEM (D-007, DESIGN §2); EAP scoring (D-007); adding any package to
  `Imports` (D-011 — `psych` is the only heavy Imports dependency; everything else is a
  guarded Suggests); Rcpp.
- **IP9 / D-031:** exact reproduction of Forbes (2023) is a permanent test-backed
  *capability*, not a default lock-in. Defaults may adopt a better method with loud
  documentation (IP6) and a D-entry. So "it would change the published AMH results" is not
  by itself an objection to a default change — but it is a constraint that the reproducing
  settings stay available.
- The package is pre-1.0 and deliberately fluid; a breaking change is permissible with a
  NEWS announcement. Do not weight backward compatibility heavily.

## Output format

In `RR02-oblique-rotation-reconsidered.md`: answer each question by number with your
reasoning and evidence; cite file:line and `citekey (p. N)` for every claim that rests on a
source. List additional findings separately under "Beyond the brief". End with concrete
recommendations, each marked apply / consider / reject-with-reason.

Emit a `## Binding criteria` section **only if** you recommend a code change concrete
enough to constrain a future milestone; if your recommendation is "D-002 stands" or "take
it to the design session", omit that section rather than inventing criteria.
