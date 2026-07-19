# tenberge1999

**Full citation.** ten Berge, J. M. F., Krijnen, W. P., Wansbeek, T., &
Shapiro, A. (1999). Some new results on correlation-preserving factor scores
prediction methods. *Linear Algebra and its Applications, 289*(1–3), 311–318.
https://doi.org/10.1016/S0024-3795(97)10007-6

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `229fdd6`) from `cairn/references/sources/tenberge1999.pdf` (local only;
gitignored). Pagination: journal pages (311–318).
Extraction: verified 2026-07-19 (M67) — every standing fact read directly against the source (pp. 311–318): Eq. 3's LCP form, Eqs. 5/7/9, the three loss functions, Thm 1 (Eq. 9 is uniquely minimizing when Ψ is nonsingular, despite being introduced in Lemma 1's singular case), Thm 2, Thm 3, Cor. 1, and the p. 318 conclusion that all LCP methods coincide under orthogonality — all confirmed exactly, no corrections. The issue number `289(1–3)` is **not** printed in the source (header reads "289 (1999) 311–318"); the printed PII `S0024-3795(97)10007-6` matches the cited DOI — observed 2026-07-19.

## Why this is a primary source

The paper behind ackwards' **default scoring method** (`method = "tenBerge"`,
DESIGN §9): it defines the class of **linear correlation-preserving (LCP)**
factor-score predictors and gives the closed-form solution that
`psych::factor.scores(method = "tenBerge")` (and our own
`.tenBerge_weights(R, Λ)` for the ESEM engine) implements.

## The result

LCP predictors are linear predictors `f̂ = A'x` constrained to reproduce the
factor correlation matrix: `E(f̂f̂') = E(ff') = Φ` — equivalently
`A = Σ^{-1/2} C Φ^{1/2}` for some columnwise-orthonormal `C` (their Eq. 3).
Three historical LCP methods differ only in the loss they minimize:

| Method | Loss | Solution |
|---|---|---|
| Anderson–Rubin (1956) / McDonald (1981) | trace of variable-residual covariance | SVD of `Σ^{1/2}Ψ^{-1}ΛΦ^{1/2}` (Eq. 5); undefined when `Ψ` singular |
| Green (1969) | trace of factor MSE matrix | SVD of `Σ^{-1/2}ΛΦ^{3/2}` (Eq. 7) |
| Krijnen, Wansbeek, & ten Berge (1996) | determinant of factor MSE | closed form `C_d = Σ^{-1/2}L(L'Σ^{-1}L)^{-1/2}`, `L = ΛΦ^{1/2}` (Eq. 9, Thm 1 — the paper's contribution; previously iterative) |

Equivalences proved: determinant solution = McDonald's whenever `Ψ` is
nonsingular (Thm 2), and = Green's when `Φ` commutes with `Λ'Σ^{-1}Λ` —
in particular whenever the factors are **orthogonal** (Thm 3 + Cor. 1). So
**with orthogonal factors, all LCP methods coincide**; there is effectively
one "tenBerge" score.

## Relation to our implementation

- **Why it's the default** (DESIGN §9; DECISIONS D-007): correlation preservation is *the*
  property bass-ackwards cares about — edges are correlations among scores,
  and LCP scores reproduce the within-level factor correlation structure
  exactly instead of shrinking it (regression scores) or distorting it
  (Bartlett under some `Ψ`). Linear → algebra-eligible (Invariant 1).
- **Varimax default ⇒ orthogonal factors ⇒ the coincidence case (Thm 3)** —
  the method's behavior is unambiguous in our default configuration.
- Even correlation-preserving scores are still *predictions*, not the factors
  (indeterminacy) — the reason `redundancy_phi` auto-resolves to a φ guard for
  EFA/ESEM but not PCA (DESIGN §9 `redundancy_phi` row, M25/M43 wording), and a companion caveat
  to [[waller2007]] §4.
- ESEM engine computes tenBerge weights itself from lavaan's `Λ` and latent
  `R` (`.tenBerge_weights`), since `lavPredict()` lacks the method (DECISIONS
  D-016).
