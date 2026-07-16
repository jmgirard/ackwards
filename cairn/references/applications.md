# applications — published uses of the bass-ackwards method

_Collapsed entries: these are **examples of the method in use** (plus one
background review), not method or oracle sources. Kept for manuscript/vignette
citation and as precedents for engine/design choices. PDFs in `pdf/`
(gitignored)._

## kim2015 — Kim & Eaton (2015)

Kim, H., & Eaton, N. R. (2015). The hierarchical structure of common mental
disorders: Connecting multiple levels of comorbidity, bifactor models, and
predictive validity. *Journal of Abnormal Psychology, 124*(4), 1064–1078.
https://doi.org/10.1037/abn0000113

**The precedent for ackwards' ESEM/ordinal engine.** NESARC (N = 43,093),
12 DSM-IV diagnoses as **categorical indicators**: ESEM in Mplus 7.11 with
**WLSMV** and orthogonal CF-varimax (κ = 1/p), bass-ackwards across levels
1–7. Argues EFA over Goldberg's PCA (comorbidity = covariance, not variance).
Their **8-factor solution failed to converge and they stopped at 7** — the
in-the-wild case for Invariant 7 (convergence is data, not an error). Also
links hierarchy levels to a bifactor model (bass-ackwards p ↔ bifactor p at
r = .99) and runs within/between-level predictive-validity comparisons.

## wright2014a — Wright & Simms (2014)

Wright, A. G. C., & Simms, L. J. (2014). On the structure of personality
disorder traits: Conjoint analyses of the CAT-PD, PID-5, and NEO-PI-3 trait
models. *Personality Disorders: Theory, Research, and Treatment, 5*(1),
43–54. https://doi.org/10.1037/per0000037

Conjoint EFA + bass-ackwards "unfolding" (citing Goldberg 2006) of CAT-PD-SF,
PID-5, and NEO-PI-3FH facet scales in n = 628 psychiatric outpatients (Mplus,
FIML for planned missingness). Five-factor level optimal; higher levels
recover Internalizing/Externalizing-like spectra. Canonical personality-
pathology application; the "a" suffix distinguishes it from other Wright 2014
papers.

## forbush2018 — Forbush et al. (2018)

Forbush, K. T., Chen, P.-Y., Hagan, K. E., Chapa, D. A. N., Gould, S. R.,
Eaton, N. R., & Krueger, R. F. (2018). A new approach to eating-disorder
classification: Using empirical methods to delineate diagnostic dimensions
and inform care. *International Journal of Eating Disorders, 51*(7), 710–721.
https://doi.org/10.1002/eat.22891

Hi-TIDE model: hierarchical factor analysis of self-reported eating/mood/
anxiety symptoms in N = 243 community adults with DSM-5 EDs; Internalizing →
distress / fear–avoidance / body dissatisfaction → 15 lowest-level dimensions;
ESEM shows the hierarchy out-predicts DSM categories (60.1% vs 35.8% of
6-month outcome variance).

## forbush2024 — Forbush et al. (2024)

Forbush, K. T., Chen, Y., Chen, P.-Y., Bohrer, B. K., Hagan, K. E., Chapa,
D. A. N., … Wildes, J. E. (2024). Integrating "lumpers" versus "splitters"
perspectives: Toward a hierarchical dimensional taxonomy of eating disorders
from clinician ratings. *Clinical Psychological Science, 12*(4), 625–643.
https://doi.org/10.1177/21677026231186803

Clinician-rating follow-up to forbush2018 (N = 252 community adults with
EDs): "modified version of Goldberg's method" — sequential extraction via
**ESEM** — yielding a 10-factor hierarchical model; hierarchical dimensions
out-predict DSM diagnoses/severity indicators at 6 months and 1 year.

## cowan2024 — Cowan et al. (2024)

Cowan, H. R., Williams, T. F., Schiffman, J., Ellman, L. M., & Mittal, V. A.
(2024). Mapping psychosis risk states onto the Hierarchical Taxonomy of
Psychopathology using hierarchical symptom dimensions. *Clinical
Psychological Science, 12*(1), 3–21.
https://doi.org/10.1177/21677026221146178

Self-report symptom items, N = 3,460 community young adults: factor analysis
at 1–10 factors with **varimax**, cross-level **Pearson correlations of
estimated factor scores**, paths ≤ .30 omitted (pure Goldberg recipe). Level
count chosen by parallel analysis + "≥ 2 items loading ≥ .30" — a concrete
published stopping rule to cite alongside `k_max`. Psychosis-risk (CHR)
mapping onto HiTOP; analysis code public at
`github.com/hrcowan/HiTOP-MAP`.

## widiger2019 — Widiger et al. (2019) _(background review; no bass-ackwards analysis)_

Widiger, T. A., Bach, B., Chmielewski, M., Clark, L. A., DeYoung, C.,
Hopwood, C. J., … Thomas, K. M. (2019). Criterion A of the AMPD in HiTOP.
*Journal of Personality Assessment, 101*(4), 345–355.
https://doi.org/10.1080/00223891.2018.1465431

Consortium review of whether DSM-5 AMPD Criterion A (self/interpersonal
deficits) belongs in HiTOP, largely via the g-PD/p-factor literature. No
hierarchical structural analysis of its own — context for HiTOP framing in
the manuscript, omit from method citations.
