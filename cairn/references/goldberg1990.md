# goldberg1990

**Full citation.** Goldberg, L. R. (1990). An alternative "description of
personality": The Big-Five factor structure. *Journal of Personality and
Social Psychology, 59*(6), 1216–1229. https://doi.org/10.1037/0022-3514.59.6.1216

**PDF.** `pdf/goldberg1990.pdf` (local only; gitignored).

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
