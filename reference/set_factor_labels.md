# Attach persistent factor labels to an ackwards object

Store substantive names for factors (e.g. `"Neuroticism"` for `"m5f1"`)
on the object itself, so that
[print()](https://jmgirard.github.io/ackwards/reference/print.ackwards.md),
[summary()](https://jmgirard.github.io/ackwards/reference/summary.ackwards.md),
[tidy()](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md),
[autoplot()](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md),
and
[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
display them without re-supplying the labels each time. This is the
*factor*-label counterpart to item / variable labels (see
[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
and
[`?ackwards`](https://jmgirard.github.io/ackwards/reference/ackwards.md));
the two are distinct and never interchanged.

## Usage

``` r
set_factor_labels(x, labels)
```

## Arguments

- x:

  An `ackwards` object.

- labels:

  A named character vector mapping factor IDs (`"m{k}f{j}"`) to label
  strings, or `NULL` to clear all stored labels. Every name must match a
  factor ID in `x` (an unknown ID is an error, not a warning – the
  object's IDs are knowable up front; see
  [`factor_labels()`](https://jmgirard.github.io/ackwards/reference/factor_labels.md)
  to read the current set).

## Value

The `ackwards` object, with `meta$factor_labels` updated. Pipeable.

## Details

Labels are display only – they never change a factor's stable ID
(`m{k}f{j}`), and every lineage / edge / score column continues to key
on the ID. A factor with no label falls back to its ID everywhere. The
stored labels ride along through
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md),
[`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md),
[augment()](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md),
and
[predict()](https://jmgirard.github.io/ackwards/reference/predict.ackwards.md)
unchanged.

## Updating and clearing

`set_factor_labels()` **merges** into any labels already stored, so you
can build them up incrementally. Within a single call:

- a normal string sets (or overwrites) that factor's label;

- an `NA` or `""` value **removes** just that factor's label;

- passing `labels = NULL` clears **all** labels at once.

The scaffold printed by
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
is a convenient starting point: copy its `c(...)` literal, fill in the
substantive names, and pass it here.

## See also

[`factor_labels()`](https://jmgirard.github.io/ackwards/reference/factor_labels.md)
to read them back,
[`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
for a ready-to-edit scaffold,
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
(`node_labels` overrides a stored label per node).

## Examples

``` r
x <- ackwards(sim16, k_max = 4, engine = "pca")

# Start from a scaffold, fill in names, attach them:
x <- set_factor_labels(x, c(m4f1 = "Alpha", m4f2 = "Beta"))
factor_labels(x)
#>    m4f1    m4f2 
#> "Alpha"  "Beta" 

# Merge in more; remove one; everything downstream now shows the labels:
x <- set_factor_labels(x, c(m2f1 = "Broad", m4f2 = NA))
summary(x)
#> 
#> ── Summary: Bass-Ackwards Analysis (ackwards) ──────────────────────────────────
#> Engine: pca
#> Rotation: varimax
#> Basis: pearson
#> n: 1,000
#> k (max): 4
#> 
#> ── Levels ──
#> 
#> k = 1: 1 factor (28.2% cumulative variance)
#> m1f1 28.2% eigenvalue 4.51
#> 
#> k = 2: 2 factors (46.5% cumulative variance)
#> Broad (m2f1) 23.3% eigenvalue 4.51
#> m2f2 23.2% eigenvalue 2.93
#> 
#> k = 3: 3 factors (57.5% cumulative variance)
#> m3f1 23.0% eigenvalue 4.51
#> m3f2 17.5% eigenvalue 2.93
#> m3f3 16.9% eigenvalue 1.76
#> 
#> k = 4: 4 factors (67.7% cumulative variance)
#> Alpha (m4f1) 17.2% eigenvalue 4.51
#> m4f2 16.9% eigenvalue 2.93
#> m4f3 16.8% eigenvalue 1.76
#> m4f4 16.8% eigenvalue 1.63
#> 
#> ── Lineage (primary parents) ──
#> 
#> m1f1 → Broad (m2f1), m2f2
#> Broad (m2f1) → m3f2, m3f3
#> m2f2 → m3f1
#> m3f1 → m4f3, m4f4
#> m3f2 → Alpha (m4f1)
#> m3f3 → m4f2
#> ────────────────────────────────────────────────────────────────────────────────
#> Note: This is a series of linked solutions, not a fitted hierarchical model.
#> Cross-level edges are descriptive score correlations. Per-level fit indices
#> (EFA/ESEM) describe how well a k-factor model fits the items at that level --
#> they do not validate the edges or the hierarchy itself.
```
