# goldberg1990

**Full citation.** Goldberg, L. R. (1990). An alternative "description of
personality": The Big-Five factor structure. *Journal of Personality and
Social Psychology, 59*(6), 1216–1229. https://doi.org/10.1037/0022-3514.59.6.1216

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `b85bee0`) from `cairn/references/sources/goldberg1990.pdf` (local only;
gitignored). Pagination: journal pages (1216–1229).
Extraction: unverified — first pass, values not re-read against the source — observed 2026-07-19.

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
split-halves, and no "repeated random splits" anywhere.** Its quality gate is
replication across *methods* (Study 1) and across *whole samples* (Studies
2–3), measured by inter-method score correlations and cross-sample Tucker
congruence. The package docs currently attribute more to it than it contains:
`comparability()`'s roxygen ("the depth gate Goldberg's own lab applied to
its hierarchies"), the `.95 comfortable / .90 floor (Everett, 1983; Goldberg,
1990)` benchmark lines in `print.comparability()` + `autoplot()`, the
`n_splits` doc ("following Goldberg's practice"), and
`vignettes/ackwards-girard.Rmd.orig` line ~51 ("Goldberg (1990) retained
factor solutions only when they held up across repeated random
split-halves" — **false as stated**). **Resolved 2026-07-16 after checking
five lab papers** (this one, Saucier & Goldberg 1996, Goldberg & Somer 2000
Turkish, Ashton/Lee/Goldberg 2004 hierarchy — all negative): the real
lineage is **[[saucier1997]]** (split-half stability gate; tried Everett's
method and de-emphasized it; drop-off criterion ~.70) and **[[saucier2005]]**
(random participant halves, congruence, "conventional threshold of .90",
used to choose the hierarchical level). The .90 floor is *also* independently
traceable to [[everett1983]] (his 26°/81% passage). Repoint the package
citations to Everett 1983 + Saucier 1997 + Saucier et al. 2005; keep
Goldberg 1990 only for the generic replication-gate ethos.

## Relation to our implementation

- The operative principle — *retain a level's factors only if they reproduce
  under resampling/measurement perturbation* — is what `comparability()`
  packages (per-split score comparability; φ secondary).
- The Study 1 method-invariance result (components vs. factors, orthogonal
  vs. oblique: near-identical scores) is the empirical backing for
  Goldberg's "virtually identical either way" claim in [[goldberg2006]] and
  useful in the manuscript when defending engine-agnostic defaults.
- Historical note: the paper that established the Big Five label-set
  (Surgency, Agreeableness, Conscientiousness, Emotional Stability,
  Intellect) used throughout ackwards' bfi25 examples.
