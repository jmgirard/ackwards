# saucier2005

**Full citation.** Saucier, G., Georgiades, S., Tsaousis, I., & Goldberg,
L. R. (2005). The factor structure of Greek personality adjectives. *Journal
of Personality and Social Psychology, 88*(5), 856–875.
https://doi.org/10.1037/0022-3514.88.5.856

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `254e023`) from `cairn/references/sources/saucier2005.pdf` (local only;
gitignored; fetched 2026-07-16 from Goldberg's ORI page,
`projects.ori.org/lrg/`). Pagination: journal pages (856–875).
Extraction: verified 2026-07-19 (M67) — all three quoted passages confirmed verbatim against the source (the hierarchical-level plan, the 495/496 subsample split, and the "conventional threshold of .90" gate, whose ellipsis is correctly placed); 3,302 adjectives and the emic six-factor stability claim also exact. **Two errors found and corrected in the surrounding summary line**, which the page's earlier "verified against the PDF" note had never covered: a sample size of "201" that appears nowhere in the paper (Samples are 991 / 429 / 538), and an asserted varimax rotation the paper never names — observed 2026-07-19.

## Why this is a primary source

The Goldberg-coauthored paper where the lab's split-half **.90 replication
gate for choosing the hierarchical level** is explicit — the second half of
the fact-checked replacement for the mis-attributed "Goldberg (1990)" (see
[[goldberg1990]] ⚠ block and [[saucier1997]]).

Key passages (verification status is recorded once, in the Provenance block above):

- Plan: "to choose the optimal hierarchical level (number of factors) **based
  primarily on replication across subsamples**."
- Procedure: "we randomly divided Sample 1 in half (i.e., 495 in one
  subsample, 496 in the other) and examined whether the principal components
  generated in one half, for each number of factors, had high matching
  **congruence coefficients** with those from the other half."
- Gate: for the solutions retained, "all congruence coefficients were above
  the **conventional threshold of .90**. These solutions can be considered
  highly robust across subsamples … as was true in an earlier study of
  English adjectives (Saucier, 1997)."

Greek lexical study: 3,302 adjectives extracted from a dictionary of modern
Greek → 400-term high-frequency and high-clarity subsets (and a reduced 248
common-term set). Samples: **1 = 991, 2 = 429, 3 = 538** participants; the
replication analyses focus on Samples 1 and 2. Solutions are
**principal-components, "1 unrotated and 2 to 10 rotated factors"** — the paper
never names the rotation method. One- and two-factor structures the most
replicable; "an emic 6-factor structure showed relative stability" (abstract).

*(Corrected M67, two errors in this summary — the three quotations above had
been verified at ingest but this line had not: it read "N ≈ 991 + 201", and no
"201" appears anywhere in the paper (it is saucier1997's acquaintance-sample N,
cross-contaminated from the page ingested in the same commit); and it read
"PCA + varimax", though "varimax" — like "oblique" and "promax" — occurs zero
times in the source.)*

## Relation to our implementation

- This is the citable Goldberg-lab precedent for the **.90 benchmark** and
  for using split-half replication as the **hierarchy-depth floor** — exactly
  `comparability()`'s workflow role. Note the index difference: they gated on
  **loading congruence** (our reported φ column); Everett's score
  comparability (our `r` column) is the stricter-lineage index the verb
  headlines. Citing this paper supports the .90 line for either column.
- Single split here; `comparability()`'s `n_splits = 10` repeated-splits
  design is our own robustness choice and should be owned as such in the docs
  (not attributed to lab practice).
