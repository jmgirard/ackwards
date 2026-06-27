# Package index

## Fit the model

The main modelling function and a helper for choosing k.

- [`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
  : Bass-ackwards hierarchical structural analysis
- [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
  : Suggest a maximum number of factors for bass-ackwards analysis
- [`print(`*`<suggest_k>`*`)`](https://jmgirard.github.io/ackwards/reference/print.suggest_k.md)
  : Print a suggest_k object

## Extract results

Tidy, glance, augment, and print methods for ackwards objects.

- [`tidy(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/tidy.ackwards.md)
  : Tidy an ackwards object into a long data frame
- [`glance(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/glance.ackwards.md)
  : Glance at an ackwards object
- [`augment(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/augment.ackwards.md)
  : Augment data with factor scores from an ackwards object
- [`print(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/print.ackwards.md)
  : Print an ackwards object

## Visualize

Diagram and layout helpers for the bass-ackwards hierarchy.

- [`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md)
  : autoplot generic
- [`autoplot(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  [`plot(`*`<ackwards>`*`)`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md)
  : Plot a bass-ackwards diagram
- [`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)
  : Compute a layered layout for a bass-ackwards diagram

## Re-exported generics

Generics re-exported from the generics package so they work after
[`library(ackwards)`](https://jmgirard.github.io/ackwards/) without also
loading broom.

- [`reexports`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`augment`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`tidy`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  [`glance`](https://jmgirard.github.io/ackwards/reference/reexports.md)
  : Objects exported from other packages
