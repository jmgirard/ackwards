# saucier1997

**Full citation.** Saucier, G. (1997). Effects of variable selection on the
factor structure of person descriptors. *Journal of Personality and Social
Psychology, 73*(6), 1296–1312. https://doi.org/10.1037/0022-3514.73.6.1296

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `254e023`) from `cairn/references/sources/saucier1997.pdf` (local only;
gitignored; fetched 2026-07-16 from Saucier's public page,
`pages.uoregon.edu/gsaucier/`). Pagination: journal pages (1296–1312).
Extraction: verified 2026-07-19 (M67) — standing facts read directly against the source (pp. 1296–1305): the 500 descriptors, N = 700 / N = 201 samples, PCA + varimax at 2–10 factors, the "single estimate of factor reliability for each factor solution" averaging, the maximize-magnitude matching, the variables-split-not-participants design, the Everett-method quotation (p. 1304, ellipsis correct), the "substantially decrease, well below .70" drop-off after seven factors (five for dispositions), and footnote 14's parallel-analysis quotation — all confirmed exactly; one correction ("nested" variable selections). M61's DOI and fn-14 evidence absorbed by citation. Issue number `73(6)` **is** printed in the source — observed 2026-07-19.

## Why this is a primary source

**The Goldberg-lab paper that actually used a split-half stability depth
gate** — the fact-checked replacement for the mis-attributed "Goldberg
(1990)" in `comparability()`'s docs (see the ⚠ block in [[goldberg1990]];
G1990 contains no split-halves). 500 familiar English person descriptors,
four variable selections of differing breadth (all 500 terms, 455 nonphysical,
broad dispositions, dispositions-and-states); self (N = 700) and acquaintance
(N = 201) samples; PCA + varimax at 2–10 factors. *(Corrected M67: "four
nested variable selections" — the paper describes four selections varying in
breadth and never calls them nested, and strict containment does not hold
between the two disposition selections.)*

## The stability procedure (Method, "Analyses")

Two indices, averaged into "a single estimate of factor reliability for each
factor solution":

1. **Between-sample**: Tucker congruence between self- and acquaintance-
   sample solutions (matched to maximize coefficient magnitude — the
   greedy-match ancestor of our `.match_square()`).
2. **Within-sample split-half**: solutions fit on random halves, factor
   scores generated, and the two sets of scores correlated. NB: the primary
   split was of the **variables**; splitting **participants** — Everett's
   (1983) method, i.e., what `comparability()` implements — was computed
   "for comparison purposes" and **de-emphasized**: "Everett's method
   generated coefficients with a similar pattern … but they were
   systematically higher, impeding precise comparisons; this method is not
   emphasized in this article."

**Decision rule**: a drop-off pattern, not a fixed gate — coefficients
"substantially decrease, well below .70" beyond the stable depth (after 7
factors for wide selections, after 5 for disposition-only selections). The
lab's explicit **.90 threshold** for subsample replication appears later, in
[[saucier2005]].

## Two useful side findings

- **Footnote 14**: parallel analysis "suggests a large number of factors when
  the number of variables is large, estimating as many as 30 factors in the
  present analyses" — in-the-wild support for `suggest_k()`'s "PA-PC tends to
  overextract; treat as an upper bound."
- Abstract's framing — "stable" factors defined by replication, wider
  selections buying extra stable factors (Attractiveness, Negative Valence) —
  is the variable-selection caution to cite when docs discuss what enters the
  battery.

## Relation to our implementation

`comparability()` implements **Everett's participant-split** version (scores
correlated on the pooled R), not Saucier's variable-split. Saucier's
observation that Everett coefficients run *systematically higher* than
variable-split/between-sample coefficients is worth echoing in the docs:
comparability coefficients are a *within-population, same-items* bound on
replicability — an easier test than cross-sample congruence, so a factor that
fails even this gate is decisively unstable (and a passing one is not thereby
cross-sample-proof).
