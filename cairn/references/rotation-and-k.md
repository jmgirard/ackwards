# rotation-and-k — rotation criteria & number-of-factors sources

_Collapsed entries: methodological supports for the varimax default and the
k-selection advisory table (DESIGN §9 + the PA/MAP table). None is an oracle
source. PDFs in `sources/` (gitignored)._

**Provenance.** Ingested 2026-07-16 by a cairn hygiene pass (no milestone;
commit `b85bee0`) from the shelf PDFs named per entry below, under
`cairn/references/sources/` (local only; gitignored). Pagination: journal pages,
per entry — **except `lim2019`, whose shelf PDF is the APA online-first version
paginated from 1**; its 452–467 range is taken from the published record (as
cited in `achim2021`'s reference list, p. 73), not from the shelf copy.
Extraction: verified 2026-07-19 (M68) — all six entries re-read against their
shelf PDFs (crawford1970 pp. 321–332; browne2001a pp. 111–150 incl. Table 1
p. 119; horn1965 pp. 179–185; velicer1976 pp. 321–327; lim2019 abstract;
achim2021 pp. 69–73). One correction, marked inline below.

## Rotation

### crawford1970 — Crawford & Ferguson (1970)

Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and
its use in orthogonal rotation. *Psychometrika, 35*(3), 321–332.
https://doi.org/10.1007/BF02310792

The Crawford–Ferguson family: minimizing the weighted sum of *test parsimony*
and *factor parsimony* (Eq. 7, `G(φ) = K₁T + K₂F`, p. 323) gives a general
criterion for either oblique or orthogonal rotation, with quartimax, varimax,
and equamax as special cases (abstract, p. 321). In the **orthogonal** case the
CF family is equivalent to the orthomax family (pp. 324–326). The CF
parameterization is what modern software (incl. Mplus/lavaan rotation
machinery) implements; **CF-varimax with κ = 1/p is the rotation Kim & Eaton
(2015) used** for the ESEM bass-ackwards precedent (`applications.md`) — the
κ indexing comes from browne2001a's Table 1 (p. 119, below), not Crawford &
Ferguson's own `K₁`/`K₂` notation.

Their closing argument (Discussion, p. 331) is the direct support for our
varimax default: "The most important factor bearing on which of these two
criteria should be used for a particular problem is the amount of information
that is available on the number of factors that should be rotated," closing
with — "The varimax criterion should usually be used whenever an accurate
estimate of the number of factors is not available." Bass-ackwards extracts a
*range* of levels precisely because k is not known, so this is the paper's own
recommendation rather than an analogy (DESIGN §9 `rotation`).

### browne2001a — Browne (2001)

Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
analysis. *Multivariate Behavioral Research, 36*(1), 111–150.
https://doi.org/10.1207/S15327906MBR3601_05

The standard modern survey of analytic rotation (CF family, Yates' geomin,
oblique target rotation; behavior under complex patterns — abstract, p. 111).
**Table 1, p. 119** ("The Orthomax Family of Rotation Criteria") is the source
of the κ indexing DESIGN §9 uses: κ = 0 quartimax, **κ = 1/p varimax**,
κ = m/(2p) equamax, κ = (m−1)/(p+m−2) parsimax, κ = 1 factor parsimony — p
variables, m factors. Cited in DESIGN's key refs for the rotation-choice
rationale and quoted by [[asparouhov2009]] as the motivation for ESEM: the
"discovery of misspecified loadings … more direct through rotation" of the
factor matrix than through modification indices (p. 113; the elision is
`, however, is` — see [[asparouhov2009]] for the full quotation, corrected
M67). The "a" suffix keeps room for Browne's other 2001 papers.

## Number of factors (advisory only — `k_max` stays user-set)

### horn1965 — Horn (1965)

Horn, J. L. (1965). A rationale and test for the number of factors in factor
analysis. *Psychometrika, 30*(2), 179–185. https://doi.org/10.1007/BF02289447

Parallel analysis: retain factors whose observed eigenvalues exceed the mean
eigenvalues of random data of the same order and N — Kaiser/Guttman's
latent-root-one corrected for sampling error and least-squares
"capitalization" on it (abstract, p. 179). Operationally: generate K random
`m × N` matrices, average their kth roots to get the curve `R_a`, and read k
off where `R_a` intersects the real-data curve (p. 182). The PA-PC / PA-FA
rows of DESIGN's k-advice table (`psych::fa.parallel(fa = "both")`).

### velicer1976 — Velicer (1976)

Velicer, W. F. (1976). Determining the number of components from the matrix
of partial correlations. *Psychometrika, 41*(3), 321–327.
https://doi.org/10.1007/BF02293557

MAP: `f_m` is the average squared partial correlation after the first m
components are partialled out (Eq. 9, p. 323); extract while `f_m` falls and
stop at its minimum — "This approach gives an exact stopping point"
(abstract, p. 321). Components-oriented, and conservative by the author's own
report: "the proposed stopping rule would suggest extracting fewer factors
than previous analyses have done" (p. 324). The MAP row of the k-advice table
(`psych::vss()`).

### lim2019 — Lim & Jahng (2019)

Lim, S., & Jahng, S. (2019). Determining the number of factors using
parallel analysis and its recent variants. *Psychological Methods, 24*(4),
452–467. https://doi.org/10.1037/met0000230

Monte Carlo of 13 PA variants: traditional PA on the full correlation matrix
performs best in most conditions, **but** accuracy degrades with weak
loadings/high factor correlations; they recommend treating the PA estimate
as a **±1 range** resolved by interpretability — the closest published
statement of ackwards' own stance (k selection is advisory; extract a range
of levels and inspect). All of this is abstract-verbatim (first page; the
shelf PDF is online-first, so no journal page number is available here).
⚠ **"Performs best" is scoped to the 13 PA variants compared here, not to
number-of-factors criteria in general** — M61 found and removed a vignette
clause that made the wider claim. Do not restate it unscoped.

### achim2021 — Achim (2021)

Achim, A. (2021). Determining the number of factors using parallel analysis
and its recent variants: Comment on Lim and Jahng (2019). *Psychological
Methods, 26*(1), 69–73. https://doi.org/10.1037/met0000269

Counterpoint: across lim2019's simulations **involving noise factors**, in
"nearly 17%" (abstract, p. 69; body rate 16.8%, p. 70) the first noise-factor
eigenvalue exceeded the kth main-factor eigenvalue — so a nominally correct k
could mean retaining a noise dimension at the expense of a signal one, and
those cases "did not deserve qualifying as successes." Hence "traditional PA
is best" overreaches ("their conclusion about traditional PA being the method
to prefer is not valid," p. 70); argues for formal sequential tests of
k-factor sufficiency (his NEST — introduced in Achim, 2017, promoted here).
*(Corrected M68: the page previously gave the 17% as a share of "correct"
retentions. The denominator is the noise-factor simulations, not the correct
retentions.)* Cite the pair together whenever the docs
touch PA — the disagreement itself is the argument for `k_max` staying a
user decision with loud advisories rather than an auto-resolved default.
