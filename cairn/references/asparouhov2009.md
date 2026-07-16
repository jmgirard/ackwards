# asparouhov2009

**Full citation.** Asparouhov, T., & Muthén, B. (2009). Exploratory
structural equation modeling. *Structural Equation Modeling, 16*(3),
397–438. https://doi.org/10.1080/10705510903008204

**PDF.** `pdf/asparouhov2009.pdf` (local only; gitignored).

## Why this is a primary source

The founding paper of **ESEM** — the third ackwards engine
(`engine = "esem"`, implemented on lavaan). Everything the engine needs from
the framework is defined here.

## The framework

- An EFA measurement block inside SEM: for `m` factors, only `m²`
  identification restrictions on `Λ` and `Ψ` (vs. the forest of fixed zero
  cross-loadings in CFA). Motivated by Browne's (2001) argument that
  discovering misspecified loadings "is more direct through rotation … than
  through modification indices" ([[rotation-and-k]] → browne2001a).
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
  (DESIGN §14.12; [[tenberge1999]]) so the one-edge-path algebra (Invariant
  1) applies unchanged.
- Rotation default stays varimax across engines for between-level score-
  correlation interpretability ([[goldberg2006]]); geomin/oblique is out of
  scope (DESIGN §2).
