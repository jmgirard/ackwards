# Screen a dataset for factorability and sampling adequacy

Before extracting a hierarchy, it is worth asking two prior questions:
*is this correlation matrix even worth factoring*, and *is the sample
large enough to trust the answer*. `factorability()` reports the
standard diagnostics for both, one object you can print or pull values
from:

- **KMO** – the Kaiser-Meyer-Olkin measure of sampling adequacy, overall
  and per variable (via
  [`psych::KMO()`](https://rdrr.io/pkg/psych/man/KMO.html)). It
  contrasts correlations with partial correlations: a low value means
  the variables share little common variance once other variables are
  partialled out, so a factor model has little to recover.

- **Bartlett's test of sphericity** – tests whether the correlation
  matrix differs from the identity (via
  [`psych::cortest.bartlett()`](https://rdrr.io/pkg/psych/man/cortest.bartlett.html)).
  Needs `N`; a non-significant result says the correlations are too weak
  to factor.

- **Sample size** – `N`, the number of variables `p`, and the `N:p`
  ratio.

- **Ledermann bound** – the largest number of *common factors*
  identifiable from `p` variables (a hard limit for EFA/ESEM, though not
  for PCA).

**Read these as conventions, not verdicts.** Every cutoff here –
Kaiser's KMO bands, the 5:1 / 10:1 `N:p` rules of thumb, Bartlett
significance at `.05` – is a widely repeated *rule of thumb*, not a
settled threshold, and each is contested in the methodological
literature (the required `N` in particular depends on communalities and
factor overdetermination far more than on any fixed ratio; MacCallum et
al. 1999). `factorability()` deliberately reports the numbers and their
conventional bands rather than returning a pass/fail flag.
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
runs the same screen internally and warns only at the genuinely
consequential extreme (KMO `< 0.5`, `N:p < 5`).

## Usage

``` r
factorability(data, cor = c("pearson", "spearman", "polychoric"), n_obs = NULL)
```

## Arguments

- data:

  A data frame or numeric matrix of item responses, **or** a correlation
  matrix (detected automatically). KMO and Bartlett are computed on the
  correlation matrix in the basis you name via `cor`.

- cor:

  The correlation basis to screen on: `"pearson"` (default),
  `"spearman"`, or `"polychoric"`. Ignored when `data` is already a
  correlation matrix. Use the same basis you plan to pass to
  [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md).

- n_obs:

  Number of observations. Taken from `nrow(data)` for raw data (a
  supplied value is ignored with a warning). **Required** for the
  N-based diagnostics (Bartlett, `N:p`) when `data` is a correlation
  matrix; those are reported as `NA` when it is absent.

## Value

An object of class `factorability` (a named list): `cor`, `n_obs`,
`n_vars`, `np_ratio`, `kmo_overall`, `kmo_items` (a data frame of
per-item MSA), `bartlett` (a list of `chisq`/`df`/`p_value`, or `NULL`
when `N` is unknown), and `ledermann` (the bound). Print it for a banded
summary; index the list for the raw values.

## References

Kaiser, H. F. (1974). An index of factorial simplicity. *Psychometrika*,
39(1), 31–36.

MacCallum, R. C., Widaman, K. F., Zhang, S., & Hong, S. (1999). Sample
size in factor analysis. *Psychological Methods*, 4(1), 84–99.

## See also

[`check_items()`](https://jmgirard.github.io/ackwards/reference/check_items.md)
(per-item screening),
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
(how many factors), and
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
(split-half replicability) – the other pre-analysis diagnostics;
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md),
which runs this screen internally.

## Examples

``` r
# Continuous data, Pearson basis
factorability(sim16)
#> 
#> ── Factorability screen (ackwards) ─────────────────────────────────────────────
#> Basis: pearson
#> Observations: 1,000
#> Variables: 16
#> N:p ratio: 62.5:1
#> 
#> ── Sampling adequacy ──
#> 
#> Overall KMO: 0.86 (meritorious)
#> 
#> Bartlett's test of sphericity: chi-square(120) = 6476.3, p < .001
#> 
#> ── Identifiability ──
#> 
#> Ledermann bound: at most 10 common factors are identifiable from 16 variables
#> (EFA/ESEM; PCA is unbounded).
#> ────────────────────────────────────────────────────────────────────────────────
#> KMO bands (Kaiser 1974), the N:p >= 5/10 rules, and Bartlett at .05 are widely
#> used *rules of thumb*, not settled thresholds -- required N depends on
#> communalities and factor overdetermination (MacCallum et al. 1999). Read the
#> numbers, not a pass/fail.

# From a correlation matrix: supply n_obs for the N-based diagnostics
R <- cor(sim16)
factorability(R, n_obs = nrow(sim16))
#> 
#> ── Factorability screen (ackwards) ─────────────────────────────────────────────
#> Basis: pearson
#> Observations: 1,000
#> Variables: 16
#> N:p ratio: 62.5:1
#> 
#> ── Sampling adequacy ──
#> 
#> Overall KMO: 0.86 (meritorious)
#> 
#> Bartlett's test of sphericity: chi-square(120) = 6476.3, p < .001
#> 
#> ── Identifiability ──
#> 
#> Ledermann bound: at most 10 common factors are identifiable from 16 variables
#> (EFA/ESEM; PCA is unbounded).
#> ────────────────────────────────────────────────────────────────────────────────
#> KMO bands (Kaiser 1974), the N:p >= 5/10 rules, and Bartlett at .05 are widely
#> used *rules of thumb*, not settled thresholds -- required N depends on
#> communalities and factor overdetermination (MacCallum et al. 1999). Read the
#> numbers, not a pass/fail.
```
