# goldberg2006

**Full citation.** Goldberg, L. R. (2006). Doing it all Bass-Ackwards: The
development of hierarchical factor structures from the top down. *Journal of
Research in Personality, 40*(4), 347–358.
https://doi.org/10.1016/j.jrp.2006.01.001

**PDF.** `pdf/goldberg2006.pdf` (local only; gitignored).

## Why this is a primary source

The founding paper of the method this package implements. Everything in
`ackwards()` descends from the procedure described here; Forbes (2023)
([[forbes2023]]) is the modern extension we treat as the fidelity contract.

## The original procedure (as specified)

1. Extract and rotate solutions with 1, 2, 3, … components until "no variable
   has its highest factor loading" on a new component (his stopping criterion;
   we require an explicit `k_max` instead — DESIGN §9).
2. Compute and save factor scores at every level; correlate the full set of
   scores in one large correlational analysis.
3. Read the correlations between **adjoining** levels as "path coefficients"
   in a top-down diagram, with the first unrotated principal component (FUPC)
   at the apex. Only paths ≥ .30 are displayed in his figures.

Key methodological positions taken in the paper:

- **Components over factors.** Goldberg recommends PCA because component
  scores are directly computable while factor scores must be estimated; he
  notes the resulting structures are "virtually identical" either way. (Waller
  2007 ([[waller2007]]) later showed the score-correlation step itself is
  unnecessary; Kim & Eaton 2015 argue the EFA side — see `applications.md`.)
- **Orthogonal (varimax) rotation at every level** — matches ackwards' default
  (DESIGN §9; only orthogonal rotations make between-level score correlations
  cleanly interpretable).
- **Size-ordered labels within level** ("4/1" = largest component of the
  four-component solution) — the ancestor of our `m{k}f{j}` labels; like ours,
  the labels carry no lineage (CLAUDE.md Invariant 5).
- **Honest caveats** (§5): the representation is "sequential" rather than
  truly hierarchical in the Yung/Thissen/McLeod sense, and the "path
  coefficients" are strictly speaking akin to part-whole correlations. Useful
  framing to echo in the manuscript/vignettes.

## Illustrations in the paper

1,710 English trait adjectives (Fig. 1–2, from Ashton, Lee, & Goldberg 2004);
435 familiar English adjectives (Fig. 3); 440 Turkish adjectives (Fig. 4);
31 dissociation symptoms (Fig. 5); 48 food-preference items (Fig. 6);
33 musical-behavior items (Fig. 7). None ship with usable matrices — these are
narrative illustrations, not oracle material.

## Relation to our implementation

- No published correlation matrix or reproducible numeric output → **not an
  oracle source**; fidelity is anchored to Forbes (2023) instead
  (`cairn/ORACLES.md`).
- Footnote credits Waller with "a set of procedures, and a computer program,
  for reproducing these structures without having to compute factor scores" —
  the algebraic route that became [[waller2007]] and, ultimately, our
  `compute_edges()` (Invariant 1).
