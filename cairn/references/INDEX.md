# References index

_One line per source summary: `citekey — full title — which milestones/oracles use it`._

## Method sources (one note per source)

- goldberg2006.md — Goldberg, L. R. (2006). Doing it all Bass-Ackwards: The development of hierarchical factor structures from the top down. *Journal of Research in Personality.* — founding paper of the method; no oracle material.
- waller2007.md — Waller, N. (2007). A general method for computing hierarchical component structures by Goldberg's Bass-Ackwards method. *Journal of Research in Personality.* — mathematical basis of the algebraic edge route (Invariants 1–2).
- forbes2023.md — Forbes, M. K. (2023). Improving hierarchical models of individual differences: An extension of Goldberg's bass-ackward method. *Psychological Methods.* — oracles O1 (AMH) + O2 (simulations).
- tenberge1999.md — ten Berge, Krijnen, Wansbeek, & Shapiro (1999). Some new results on correlation-preserving factor scores prediction methods. *Linear Algebra and its Applications.* — basis of the `tenBerge` scoring default.
- lorenzoseva2006.md — Lorenzo-Seva & ten Berge (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. *Methodology.* — source of the φ > .95 / .85 thresholds (`prune()`, `comparability()`).
- everett1983.md — Everett, J. E. (1983). Factor comparability as a means of determining the number of factors and their rotation. *Multivariate Behavioral Research.* — basis of `comparability()` (M46).
- goldberg1990.md — Goldberg, L. R. (1990). An alternative "description of personality": The Big-Five factor structure. *JPSP.* — the inventor's replication quality gate (`comparability()` precedent); method-invariance evidence.
- asparouhov2009.md — Asparouhov & Muthén (2009). Exploratory structural equation modeling. *Structural Equation Modeling.* — foundation of the ESEM engine + WLSMV ordinal default.
- saucier1997.md — Saucier, G. (1997). Effects of variable selection on the factor structure of person descriptors. *JPSP.* — the Goldberg-lab split-half stability depth gate (verified `comparability()` lineage; Everett variant de-emphasized).
- saucier2005.md — Saucier, Georgiades, Tsaousis, & Goldberg (2005). The factor structure of Greek personality adjectives. *JPSP.* — the lab's explicit split-half .90 gate for choosing the hierarchical level.

## Rotation & number-of-factors (collapsed in rotation-and-k.md)

- crawford1970 — Crawford & Ferguson (1970), general rotation criterion — CF family; Kim & Eaton's CF-varimax κ = 1/p.
- browne2001a — Browne (2001), overview of analytic rotation — rotation-choice rationale; motivates ESEM.
- horn1965 — Horn (1965), parallel analysis — PA rows of the k-advice table.
- velicer1976 — Velicer (1976), MAP criterion — MAP row of the k-advice table.
- lim2019 — Lim & Jahng (2019), PA and its variants — "PA estimate as ±1 range" = our advisory stance.
- achim2021 — Achim (2021), comment on Lim & Jahng — counterpoint (NEST); cite the pair together.

## Applications (collapsed in applications.md)

- markon2005 — Markon, Krueger, & Watson (2005), normal + abnormal personality structure — the influential early application; Big Two→Five unfolding.
- kim2015 — Kim & Eaton (2015), hierarchical structure of common mental disorders — ESEM/WLSMV precedent; in-the-wild Invariant-7 case.
- wright2014a — Wright & Simms (2014), structure of PD traits (CAT-PD/PID-5/NEO-PI-3) — canonical personality-pathology application.
- forbush2018 — Forbush et al. (2018), Hi-TIDE eating-disorder classification — self-report application.
- forbush2024 — Forbush et al. (2024), lumpers vs. splitters ED taxonomy — clinician-rating ESEM application.
- cowan2024 — Cowan et al. (2024), psychosis risk onto HiTOP — pure-Goldberg-recipe application; public code.

## Framing & contrasts (collapsed in background.md)

- kotov2017 — Kotov et al. (2017), HiTOP founding paper — Intro framing; the taxonomy the applications target.
- saucier1996 — Saucier & Goldberg (1996), Big Five in familiar adjectives — robustness/anti-prestructuring evidence; ruled out as comparability source.
- schmid1957 — Schmid & Leiman (1957), hierarchical factor solutions — the bottom-up contrast; out-of-scope citation (DESIGN §2).
- yung1999 — Yung, Thissen, & McLeod (1999), higher-order vs. hierarchical factor models — the "sequential, not truly hierarchical" caveat.
- widiger2019 — Widiger et al. (2019), Criterion A of the AMPD in HiTOP — background review only; no bass-ackwards analysis.
