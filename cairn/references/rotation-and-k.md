# rotation-and-k — rotation criteria & number-of-factors sources

_Collapsed entries: methodological supports for the varimax default and the
k-selection advisory table (DESIGN §9 + the PA/MAP table). None is an oracle
source. PDFs in `pdf/` (gitignored)._

## Rotation

### crawford1970 — Crawford & Ferguson (1970)

Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and
its use in orthogonal rotation. *Psychometrika, 35*(3), 321–332.
https://doi.org/10.1007/BF02310792

The Crawford–Ferguson family: a weighted sum of *test parsimony* and *factor
parsimony* whose special cases include quartimax, varimax, and equamax. The
CF parameterization is what modern software (incl. Mplus/lavaan rotation
machinery) implements; **CF-varimax with κ = 1/p is the rotation Kim & Eaton
(2015) used** for the ESEM bass-ackwards precedent (`applications.md`).
Their closing argument — the right criterion depends on how much is known
about the number of factors — reads well next to bass-ackwards' extract-many-
levels stance.

### browne2001a — Browne (2001)

Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
analysis. *Multivariate Behavioral Research, 36*(1), 111–150.
https://doi.org/10.1207/S15327906MBR3601_05

The standard modern survey of analytic rotation (CF family, geomin, target
rotation; behavior under complex patterns). Cited in DESIGN's key refs for
the rotation-choice rationale and quoted by [[asparouhov2009]] as the
motivation for ESEM ("discovery of misspecified loadings … more direct
through rotation"). The "a" suffix keeps room for Browne's other 2001 papers.

## Number of factors (advisory only — `k_max` stays user-set)

### horn1965 — Horn (1965)

Horn, J. L. (1965). A rationale and test for the number of factors in factor
analysis. *Psychometrika, 30*(2), 179–185. https://doi.org/10.1007/BF02289447

Parallel analysis: retain factors whose observed eigenvalues exceed the mean
eigenvalues of random data of the same order and N — Kaiser/Guttman's
latent-root-one corrected for sampling "capitalization." The PA-PC / PA-FA
rows of DESIGN's k-advice table (`psych::fa.parallel(fa = "both")`).

### velicer1976 — Velicer (1976)

Velicer, W. F. (1976). Determining the number of components from the matrix
of partial correlations. *Psychometrika, 41*(3), 321–327.
https://doi.org/10.1007/BF02293557

MAP: extract components while the average squared partial correlation of the
residual matrix keeps falling; the minimum is an exact stopping point.
Components-oriented, usually conservative — the MAP row of the k-advice table
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
of levels and inspect).

### achim2021 — Achim (2021)

Achim, A. (2021). Determining the number of factors using parallel analysis
and its recent variants: Comment on Lim and Jahng (2019). *Psychological
Methods, 26*(1), 69–73. https://doi.org/10.1037/met0000269

Counterpoint: in lim2019's noise-factor conditions, ~17% of "correct"
retentions actually kept a noise dimension while dropping a signal one, so
"traditional PA is best" overreaches; argues for formal sequential tests of
k-factor sufficiency (his NEST). Cite the pair together whenever the docs
touch PA — the disagreement itself is the argument for `k_max` staying a
user decision with loud advisories rather than an auto-resolved default.
