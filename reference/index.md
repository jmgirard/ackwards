# Package index

## Fit the model

The main modelling function and a helper for choosing k.

- [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  : Bass-ackwards hierarchical structural analysis
- [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  : Suggest a maximum number of factors for bass-ackwards analysis
- [`print(`*`<suggest_k>`*`)`](https://jmgirard.github.io/ackwards/reference/print.suggest_k.md)
  : Print a suggest_k object
- [`autoplot(`*`<suggest_k>`*`)`](https://jmgirard.github.io/ackwards/reference/autoplot.suggest_k.md)
  : Plot a suggest_k diagnostic

## Interpret factors

Helpers for understanding what each factor means and labelling it.

- [`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md)
  : Display the salient items for each factor
- [`print(`*`<top_items>`*`)`](https://jmgirard.github.io/ackwards/reference/print.top_items.md)
  : Print a top_items object
- [`label_template()`](https://jmgirard.github.io/ackwards/reference/label_template.md)
  : Generate a node-label scaffold for autoplot

## Extract results

Tidy, glance, augment, summary, and print methods for ackwards objects.

- [`tidy(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
  : Tidy an ackwards object into a long data frame
- [`glance(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)
  : Glance at an ackwards object
- [`augment(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
  : Augment data with factor scores from an ackwards object
- [`summary(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/summary.ackwards.md)
  : Summarise an ackwards object
- [`print(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
  : Print an ackwards object
- [`print(`*`<summary_ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/print.summary_ackwards.md)
  : Print a summary_ackwards object

## Low-level utilities

Functions for advanced use or building on top of ackwards results.

- [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md)
  : Compute between-level factor-score correlations

## Visualize

Diagram and layout helpers for the bass-ackwards hierarchy.

- [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  : autoplot generic
- [`autoplot(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  [`plot(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  : Plot a bass-ackwards diagram or per-level fit index chart
- [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)
  : Compute a layered layout for a bass-ackwards diagram

## Data

Example dataset bundled with the package.

- [`bfi25`](https://jmgirard.github.io/ackwards/reference/bfi25.md) :
  Big Five Inventory – 25-item IPIP example dataset

## Re-exported generics

Generics re-exported from the generics package so they work after
[`library(ackwards)`](https://jmgirard.github.io/ackwards/) without also
loading broom.

- [`reexports`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`augment`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`tidy`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`glance`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  : Objects exported from other packages
