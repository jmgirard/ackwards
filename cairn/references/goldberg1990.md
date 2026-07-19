# goldberg1990

**Full citation.** Goldberg, L. R. (1990). An alternative "description of
personality": The Big-Five factor structure. *Journal of Personality and
Social Psychology, 59*(6), 1216–1229. https://doi.org/10.1037/0022-3514.59.6.1216

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `b85bee0`) from `cairn/references/sources/goldberg1990.pdf` (local only;
gitignored). Pagination: journal pages (1216–1229).
Extraction: verified 2026-07-19 (M67) — every standing fact read directly against the source (pp. 1216–1224): Study 1's 1,431 terms / 75 Norman categories / N = 187, the 10 method combinations (5 extractions × varimax/oblimin), the .950–.996 inter-method range, the 30-of-3,750 (< 1%) modal-factor result, the verbatim "nearly invariant across rotations of up to 13 factors" (p. 1221), Study 2's 479 terms / 133 clusters / four samples with congruence .86–.94 averaging .91 (p. 1222), and Study 3's 100 clusters from 339 terms (p. 1223) — all confirmed exactly. The ⚠ block's negative claim was re-confirmed against the full text. Two corrections: the stale "docs currently over-attribute" observation (the repointing has landed) and the "established the label-set" overclaim. Issue number `59(6)` **is** printed in the source — observed 2026-07-19.

## Why this is a primary source

Cited alongside [[everett1983]] as the precedent for `comparability()`: the
bass-ackwards inventor's own replication-based quality gate for deciding
which factors are real, applied at scale years before the 2006 method paper.

## What it shows (three studies, English trait adjectives)

- **Study 1** (1,431 terms in Norman's 75 categories; N = 187): the same five
  varimax factors emerge from **10 method combinations** (5 extraction
  procedures × orthogonal/oblique rotation) — mean factor-score correlations
  between methods .950–.996; < 1% of 3,750 loadings had a different modal
  factor. The first 7 factors stay "nearly invariant across rotations of up
  to 13 factors"; everything past 7 is defined by only one or two variables.
- **Study 2** (479 terms, 133 synonym clusters; 4 samples, self + peer
  ratings): Big Five replicate in every sample (cross-sample Tucker
  congruence .86–.94, mean .91); **no factor beyond the Big Five replicated
  across the four samples** — the retention decision is made by replication,
  not by an eigenvalue rule.
- **Study 3**: 100 clusters from 339 terms proposed as Big-Five markers.

## ⚠ Citation caution (found at ingest, 2026-07-16)

**The paper contains no Everett-style comparability coefficients, no
split-halves, and no "repeated random splits" anywhere** — re-confirmed
against the full text at M67 (2026-07-19). Its quality gate is replication
across *methods* (Study 1: 10 extraction × rotation combinations) and across
*whole, independent samples* (Study 2's four samples A–D; Study 3's reuse of
two of them), measured by inter-method score correlations and cross-sample
Tucker congruence — never by partitioning one sample.

**Resolved 2026-07-16 after checking five lab papers** (this one, Saucier &
Goldberg 1996, Goldberg & Somer 2000 Turkish, Ashton/Lee/Goldberg 2004
hierarchy — all negative): the real lineage is **[[saucier1997]]** (split-half
stability gate; tried Everett's method and de-emphasized it; drop-off
criterion ~.70) and **[[saucier2005]]** (random participant halves,
congruence, "conventional threshold of .90", used to choose the hierarchical
level). The .90 floor is *also* independently traceable to [[everett1983]]
(his 26°/81% passage).

**The repointing is done** — observed 2026-07-19. The four over-attributing
sites this block once listed (`comparability()`'s "Goldberg's own lab"
roxygen, the `(Everett, 1983; Goldberg, 1990)` benchmark lines in
`print.comparability()` + `autoplot()`, the `n_splits` "following Goldberg's
practice" doc, and the `ackwards-girard` vignette's "repeated random
split-halves" sentence) were all corrected by M61/M63 and none survives in
`R/`, `man/`, or `vignettes/`. *(Corrected M67: this block previously said
the docs "currently" over-attribute — a dated observation left standing after
the fix landed. Goldberg 1990 is now cited only for the generic
replication-gate ethos, as intended.)*

## Relation to our implementation

- The operative principle — *retain a level's factors only if they reproduce
  under resampling/measurement perturbation* — is what `comparability()`
  packages (per-split score comparability; φ secondary).
- The Study 1 method-invariance result (components vs. factors, orthogonal
  vs. oblique: near-identical scores) is the empirical backing for
  Goldberg's "virtually identical either way" claim in [[goldberg2006]] and
  useful in the manuscript when defending engine-agnostic defaults.
- Historical note: the source of the Big Five label-set used throughout
  ackwards' bfi25 examples — Table 2's factor headings are exactly Surgency,
  Agreeableness, Conscientiousness, Emotional Stability, Intellect. Note the
  paper presents I–IV as *traditional* labels it inherits, not ones it coins
  ("these 'Big-Five' factors have **traditionally** been numbered and labeled
  as follows", p. 1217); its own move is adopting **Intellect** for Factor V,
  where the traditional label was *Culture*. *(Corrected M67: previously said
  this paper "established" the label-set, which it explicitly disclaims.)*
