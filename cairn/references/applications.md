# applications — published uses of the bass-ackwards method

_Collapsed entries: these are **examples of the method in use** (plus one
background review), not method or oracle sources. Kept for manuscript/vignette
citation and as precedents for engine/design choices. PDFs in `sources/`
(gitignored)._

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `351a916`) from the shelf PDFs named per entry below, under
`cairn/references/sources/` (local only; gitignored). Pagination: journal pages,
per entry.
Extraction: verified 2026-07-19 (M68) — all six entries re-read
against their shelf PDFs (kim2015 pp. 1064–1078; markon2005 pp. 139–157;
wright2014a pp. 43–54; forbush2018 pp. 710–721; forbush2024 pp. 625–643;
cowan2024 pp. 3–21), plus [[goldberg2006]] p. 356 for the markon2005
attribution. Three corrections, marked inline below — observed 2026-07-19.

## kim2015 — Kim & Eaton (2015)

Kim, H., & Eaton, N. R. (2015). The hierarchical structure of common mental
disorders: Connecting multiple levels of comorbidity, bifactor models, and
predictive validity. *Journal of Abnormal Psychology, 124*(4), 1064–1078.
https://doi.org/10.1037/abn0000113

**The precedent for ackwards' ESEM/ordinal engine.** NESARC Wave 1
(N = 43,093; Wave 2 N = 34,653), 12 DSM-IV diagnoses as **categorical
indicators**: ESEM in Mplus 7.11 with **WLSMV** and orthogonal CF-varimax,
bass-ackwards across levels 1–7 (pp. 1066–1067). Argues EFA over Goldberg's
PCA — EFA "produces interpretable factors that account for multivariate
*covariance* (comorbidity), while PCA produces components … that account for
*variance*" (p. 1066). Their **8-factor solution failed to converge and they
stopped at 7** ("the eight-factor solution failed to converge, and thus we
used it as a reasonable stopping point," p. 1067) — the in-the-wild case for
Invariant 7 (convergence is data, not an error). Also links hierarchy levels
to a bifactor model (bass-ackwards p ↔ bifactor p at r = .99, p. 1067) and
runs within/between-level predictive-validity comparisons (p. 1066).

**p. 1067 is the verbatim source for two DESIGN §9 positions**: "only
orthogonal rotations produce interpretable between-level factor score
correlations in the bass-ackwards method, and orthogonal rotations are
therefore suggested by Goldberg (2006)"; and "we used the orthogonal
Crawford-Ferguson family varimax rotation (where κ = 1/p; for more
information, see Browne, 2001; Crawford & Ferguson, 1970)" — the κ = 1/p
claim in [[rotation-and-k]] traces here.

## markon2005 — Markon, Krueger, & Watson (2005)

Markon, K. E., Krueger, R. F., & Watson, D. (2005). Delineating the
structure of normal and abnormal personality: An integrative hierarchical
approach. *Journal of Personality and Social Psychology, 88*(1), 139–157.
https://doi.org/10.1037/0022-3514.88.1.139

An early top-down application: sequential **varimax** factor solutions at
**2–5 factors** ("only factor solutions including two to five factors were
examined further" — note there is no one-factor level) integrating normal and
abnormal personality measures, replicated across a meta-analytic correlation
matrix and an empirical sample (abstract, p. 139). The two-factor solution
replicated less well under varimax, so they rotated the meta-analytic loadings
to Study 2's via a partially specified target criterion, giving congruences of
.97 and .98. The 2→5 unfolding is what wright2014a and the psychopathology
applications take as the expected top-of-hierarchy pattern; "Big Two → Big
Three → Big Four → Big Five" is *our* shorthand for it, not the paper's
language (those terms do not appear in it), and the paper never uses the word
"bass-ackwards" — it predates [[goldberg2006]]'s publication.

*(Corrected M68: the page previously said [[goldberg2006]] "§5 singles it
out" as the most influential early application. It does not. At p. 356
Goldberg names **Saucier (1997)** — "one of the most influential uses of the
method may have been Saucier's (1997) comprehensive analyses" — and lists
Markon, Krueger, and Watson (2005) as one of four "recent reports" providing
"additional examples of this top down procedure," alongside Mlačić &
Ostendorf (2005), Roberts et al. (2005), and Saucier et al. (2005). See
[[saucier1997]].)*

## wright2014a — Wright & Simms (2014)

Wright, A. G. C., & Simms, L. J. (2014). On the structure of personality
disorder traits: Conjoint analyses of the CAT-PD, PID-5, and NEO-PI-3 trait
models. *Personality Disorders: Theory, Research, and Treatment, 5*(1),
43–54. https://doi.org/10.1037/per0000037

Conjoint EFA + hierarchical "unfolding" (citing Markon et al. 2005 and
Goldberg 2006) of the **primary scales** of CAT-PD-SF, PID-5, and NEO-PI-3FH
in n = 628 psychiatric outpatients (64% female; 695 before exclusions), in
Mplus 7, with FIML for the planned missingness of a balanced incomplete block
design. Five-factor solution gave "conceptually coherent alignment" across the
three measures (abstract, p. 43); higher levels recover
Internalizing/Externalizing-like spectra. Canonical personality-pathology
application; the "a" suffix distinguishes it from other Wright 2014 papers.

⚠ **Note the rotation: oblique geomin, not varimax** — "due to its desirable
weighting of factor complexity and interpretability (Browne, 2001; Sass &
Schmitt, 2010)." This is a published bass-ackwards-style application that
departs from the orthogonal-only stance DESIGN §9 rests on (Goldberg 2006;
kim2015 above), so do not cite it as support for that stance.

## forbush2018 — Forbush et al. (2018)

Forbush, K. T., Chen, P.-Y., Hagan, K. E., Chapa, D. A. N., Gould, S. R.,
Eaton, N. R., & Krueger, R. F. (2018). A new approach to eating-disorder
classification: Using empirical methods to delineate diagnostic dimensions
and inform care. *International Journal of Eating Disorders, 51*(7), 710–721.
https://doi.org/10.1002/eat.22891

Hi-TIDE model: hierarchical factor analysis of self-reported eating/mood/
anxiety symptoms in N = 243 community-recruited adults with a DSM-5 ED
(82.2% women), assessed at baseline, 6 months, and 1 year. A broad
Internalizing factor branches into three subfactors — distress, fear–avoidance,
and body dissatisfaction (the last **nested within distress**) — over a lowest
level of 15 factors. ESEM shows the hierarchy out-predicts DSM categories:
60.1% of outcome variance at 6-month follow-up vs 35.8% for all DSM eating,
mood, and anxiety disorders combined (abstract, p. 710). Verified exact.

## forbush2024 — Forbush et al. (2024)

Forbush, K. T., Chen, Y., Chen, P.-Y., Bohrer, B. K., Hagan, K. E., Chapa,
D. A. N., … Wildes, J. E. (2024). Integrating "lumpers" versus "splitters"
perspectives: Toward a hierarchical dimensional taxonomy of eating disorders
from clinician ratings. *Clinical Psychological Science, 12*(4), 625–643.
https://doi.org/10.1177/21677026231186803

Clinician-rating follow-up to forbush2018 (N = 252 community-recruited adults
with an ED, 81.9% female): a "modified version of Goldberg's method, which
involved sequentially extracting latent factors using exploratory structural
equation modeling" — yielding a 10-factor hierarchical-dimensional model
(abstract, p. 625).

**Method detail (p. 633).** 73 indicators; parallel analysis set the maximum,
so they ran ESEMs from 1 to 11 factors and correlated factor scores between
*consecutive* levels. They state the lineage explicitly: "Whereas Goldberg
extracted the factors using the principal component analysis, we followed the
modification from Kim and Eaton (2015) and used exploratory structural
equation models (ESEMs) for factor extraction." ⚠ **Their estimator is
ULSMV, not WLSMV** — "a mean-and-variance-adjusted *unweighted* least squares
estimator with the orthogonal Crawford-Ferguson rotation ('crawfer' in
Mplus)." kim2015 (above) is the WLSMV precedent; this paper is not. (κ is not
stated here, so this supports *orthogonal CF* generally, not κ = 1/p
specifically.)

Two *distinct* results, easily conflated: (a) the dimensions predicted
**92.4%** and **58.7%** of the variance in recovery outcomes at 6 months and
1 year respectively — that is their own predictive power, not a comparison;
(b) against other illness indicators (DSM diagnoses, dimensional ED impairment
scores, weight/shape overvaluation, DSM ED-severity specifiers), hierarchical
dimensions predicted **0.88 to 334 times** more variance in ED behaviors at
baseline and **1.95 to 80.8 times** more in psychiatric impairment at 1-year
follow-up. *(Corrected M68: the page previously read the two as one claim that
the dimensions "out-predict DSM diagnoses/severity indicators at 6 months and
1 year." The comparison's low end is 0.88× — below parity — so the advantage
is not uniform, and the two comparisons are at baseline and 1 year, not
6 months and 1 year.)*

## cowan2024 — Cowan et al. (2024)

Cowan, H. R., Williams, T. F., Schiffman, J., Ellman, L. M., & Mittal, V. A.
(2024). Mapping psychosis risk states onto the Hierarchical Taxonomy of
Psychopathology using hierarchical symptom dimensions. *Clinical
Psychological Science, 12*(1), 3–21.
https://doi.org/10.1177/21677026221146178

Self-report symptom items, N = 3,460 community young adults (mean age 20.3):
factor analysis at 1–10 factors with **varimax**, cross-level **Pearson
correlations of estimated factor scores** between *adjacent* levels, and
"Correlations *r* ≤ .30 are omitted for clarity" (Fig. 2 caption, p. 9 —
note this omits *at* .30, where Goldberg displays paths ≥ .30). Level count
chosen by parallel analysis + all factors having "at least two items loaded
≥ .30" (11- and 12-factor solutions failed the second test) — a concrete
published stopping rule to cite alongside `k_max`. Psychosis-risk (CHR)
mapping onto HiTOP; analysis code public at
`github.com/hrcowan/HiTOP-MAP`. p. 8 also states the rationale we rely on:
"Orthogonal rotation is typical in unfolding analyses (e.g., Kim & Eaton,
2015) because it produces interpretable cross-level correlations between
latent variables (Goldberg, 2006)."

*(Corrected M68: the page called this a "pure Goldberg recipe." The
**procedure** is Goldberg's, but the **engine is not** — Cowan et al. used
minimum residual factor analysis with regression-based factor scores, where
Goldberg prescribes PCA with directly computed component scores. Relevant to
us twice over: it is an EFA-engine precedent, and its regression scores are
the method ackwards deliberately does not default to (`tenBerge`, DESIGN §9).)*

_(widiger2019 moved to `background.md` 2026-07-16.)_
