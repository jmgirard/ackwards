# everett1983

**Full citation.** Everett, J. E. (1983). Factor comparability as a means of
determining the number of factors and their rotation. *Multivariate
Behavioral Research, 18*(2), 197–218.
https://doi.org/10.1207/s15327906mbr1802_5

**PDF.** `pdf/everett1983.pdf` (local only; gitignored).

## Why this is a primary source

The basis of `comparability()` (M46; DESIGN §14.35): split-half factor
comparability as the criterion for how many factors are *reliably
recoverable* — the quality gate [[goldberg1990]] applied and the modern
HiTOP-era literature largely dropped.

## The procedure (as specified)

1. Split the `t` respondents into halves `1` and `2`; run the **identical**
   factor analysis in each half, obtaining score-coefficient matrices `S₁`,
   `S₂`.
2. Apply **both** halves' coefficients to the **total** sample:
   `F₁ₜ = S₁·Vₜ`, `F₂ₜ = S₂·Vₜ` — duplicate factor scores for every
   respondent.
3. Cross-correlate the duplicate scores → **comparability coefficients** (one
   per matched factor).

Key positions:

- **Comparability ≠ congruence.** Tucker's φ measures stability of
  *loadings*; comparability measures stability of the *scores* — the right
  target when factors are used as summary/classificatory measures (the
  "taxonomic" view, which Everett argues suits PCA). The two agree closely
  but not exactly in practice. (Complement: [[lorenzoseva2006]] calibrates
  the loading-side index.)
- **Threshold ≥ .90**: split-half factors diverging < 26°, ≥ 81% shared
  variance.
- Choose the number of factors **and the rotation** by maximizing split-half
  comparability — his synthetic demos show comparability collapsing exactly
  where extraction passes the true k (e.g., his Table 2: 3 true factors,
  comparabilities .98/.97/.96 at k = 3 vs .99/.94/.82/**.04** at k = 4).
- For non-homogeneous populations, split by subpopulation as well as at
  random — factors should be stable across both.

## Relation to our implementation

`comparability(cor, fm, n_splits = 10, seed)` (DESIGN §14 / M46) resurrects
this procedure with repeated random splits (Everett used single splits);
scores-side comparability is the headline diagnostic, with Tucker's φ on
matched loadings kept visible as the secondary index. ESEM/basis extension is
demand-gated (ROADMAP candidate). His "comparability collapses past true k"
result is the interpretive story for the vignette: a level whose factors
aren't split-half-stable is below the reliable depth of the hierarchy.
