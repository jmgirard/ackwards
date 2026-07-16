# lorenzoseva2006

**Full citation.** Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's
congruence coefficient as a meaningful index of factor similarity.
*Methodology, 2*(2), 57–64. https://doi.org/10.1027/1614-2241.2.2.57
*(NB: the PDF's printed DOI uses the print ISSN, `10.1027/1614-1881.2.2.57`,
which does **not** resolve; the registered DOI above — online ISSN — is the
one `prune()`'s roxygen correctly cites. Verified 2026-07-16.)*

**PDF.** `pdf/lorenzoseva2006.pdf` (local only; gitignored).

## Why this is a primary source

The empirical basis for every **Tucker's-φ threshold** in the package:

- `prune("redundant")`'s φ guard (`redundancy_phi` auto-resolves to **.95**
  for EFA/ESEM — DESIGN §14, M25);
- the optional conjunctive φ > **.95** criterion in the Forbes-faithful
  redundancy rule (DESIGN §14.19);
- φ interpretation in `comparability()` diagnostics (DESIGN §14 / M46).

φ formula (their Eq. 1): `φ(x,y) = Σxᵢyᵢ / sqrt(Σxᵢ² · Σyᵢ²)` — the cosine
between loading vectors; scale-invariant, sign-of-variable-invariant,
sensitive to additive constants.

## The study and the thresholds

Prior thresholds were unsourced rules of thumb (.80 Horn; .90 Mulaik/Bentler;
Tucker's own tiered bands via MacCallum et al. 1999). This study calibrates φ
against **expert judgment**: 56 factor-analysis practitioners rated the
similarity of 448 loading-column pairs with constructed φ ∈ {.62 … .97};
subjective ratings tracked φ nearly linearly (r = .974). Calibration:

- **φ .85–.94 → "fair" similarity** — and, importantly, **φ < .85 should not
  be read as indicating any factor similarity at all**;
- **φ > .95 → factors/components can be considered equal** — the source of
  our .95 default.

Secondary findings worth remembering: individual judges varied widely (single-
judge rules of thumb are hazardous — the whole point of the empirical
calibration); fewer variables and having variable labels made judges slightly
*more* conservative, not less.

## Relation to our implementation

The .95 default is the "factors are interchangeable" bound — exactly the
semantics `prune("redundant")` needs (a child that is *the same factor* as its
ancestor is redundant). The .85 floor is the right talking point for why we do
**not** expose a permissive φ default: below .85 congruence is meaningless as
similarity evidence. Both thresholds are quotable in docs/vignettes and the
BRM manuscript when a reviewer asks "why .95?".
