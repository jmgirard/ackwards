# asparouhov2009

**Full citation.** Asparouhov, T., & Muthén, B. (2009). Exploratory
structural equation modeling. *Structural Equation Modeling, 16*(3),
397–438. https://doi.org/10.1080/10705510903008204

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `b85bee0`) from `cairn/references/sources/asparouhov2009.pdf` (local only;
gitignored). Pagination: journal pages (397–438).
Extraction: verified 2026-07-19 (M67) — standing facts read directly against the source (pp. 397–401): the m² identification restrictions on Λ and Ψ (pp. 398, 401), the Jennrich-extending rotated-parameter standard errors and overall fit tests (p. 399), and the categorical-indicator formulation — Eq. 3's τ thresholds, the underlying normal variable, the probit equivalence, and limited-information weighted least squares following Muthén (1984) (p. 400) — all confirmed exactly; one correction (the Browne quotation, previously compressed with an unmarked elision). The issue number `16(3)` is not printed in the source (header reads "16:397–438, 2009") — observed 2026-07-19.

## Why this is a primary source

The founding paper of **ESEM** — the third ackwards engine
(`engine = "esem"`, implemented on lavaan). Everything the engine needs from
the framework is defined here.

## The framework

- An EFA measurement block inside SEM: for `m` factors, only `m²`
  identification restrictions on `Λ` and `Ψ` (vs. the forest of fixed zero
  cross-loadings in CFA). Motivated by Browne's (2001) argument, block-quoted
  at p. 398, that "the discovery of misspecified loadings … is more direct
  through rotation of the factor matrix than through the examination of model
  modification indices" (Browne 2001, p. 113; [[rotation-and-k]] → browne2001a).
  *(Quotation corrected M67: the page previously compressed this to "…than
  through modification indices" with the elision unmarked.)*
- Rotation applied *within* the SEM: standard errors for all rotated
  parameters (extending Jennrich), overall model-fit tests — the properties
  that let bass-ackwards levels carry per-level fit and loading SEs.
- **Categorical indicators** via the underlying-normal-variable / probit
  formulation (their Eq. 3, following Muthén 1984) with limited-information
  weighted least squares — the basis for ackwards' **WLSMV ordinal default**
  (DESIGN §9).
- Discusses geomin and target rotation; Mplus implementation.

## Relation to our implementation

- The ESEM engine wraps **lavaan's** ESEM (`efa()` blocks), not Mplus; the
  model class is this paper's.
- Fit-bearing levels + non-convergence handling: an ESEM level that fails to
  converge is data, not an error (Invariant 7) — the in-the-wild precedent
  is Kim & Eaton's NESARC 8-factor level (`applications.md` → kim2015, which
  cites this paper for its ESEM machinery).
- Scoring: ackwards does **not** use ESEM-native score estimation for edges —
  tenBerge weights are self-computed from lavaan's `Λ` and latent `R`
  (DECISIONS D-016; [[tenberge1999]]) so the one-edge-path algebra (Invariant
  1) applies unchanged.
- Rotation default stays varimax across engines for between-level score-
  correlation interpretability ([[goldberg2006]]); geomin/oblique is out of
  scope (DESIGN §2).
