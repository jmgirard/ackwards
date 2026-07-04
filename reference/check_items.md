# Screen items for problems before factor analysis

Real-world item sets often contain columns that break or *silently*
degrade a factor analysis: an item everyone answered the same way (no
variance), an item dominated by one response with a couple of stray
answers, or an item with heavy missingness. On the polychoric basis
these are especially costly – a near-empty response category can make
[`psych::polychoric()`](https://rdrr.io/pkg/psych/man/tetrachor.html)
fail outright (see the `correct` argument of
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)),
and a near-constant item can produce a plausible-looking but meaningless
factor with no warning.

`check_items()` reports these problems *before* you fit, one row per
item, so you can collapse rare categories, drop degenerate items, or set
`cor`/`correct` deliberately rather than debugging a cryptic error. It
only reports; it never changes your data.
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
runs the same screen internally: it **errors** on a constant item and
**warns** on a near-degenerate one, naming the offenders.

## Usage

``` r
check_items(data, cor = c("polychoric", "pearson", "spearman"))
```

## Arguments

- data:

  A data frame or numeric matrix (items in columns).

- cor:

  The correlation basis you plan to use: `"polychoric"` (default),
  `"pearson"`, or `"spearman"`. Only affects which problems are flagged
  – sparse response categories destabilise the polychoric basis but not
  the others.

## Value

A data frame (class `check_items`) with one row per item and columns
`item`, `n_valid`, `pct_missing`, `n_distinct`, `min_count` (smallest
observed category count), `top_prop` (proportion of valid responses in
the most common value), and `flag` – one of `"ok"`, `"constant"`,
`"near-constant"`, `"sparse category"`, or `"high missing"`. Print it
for a grouped summary with guidance; treat it as a plain data frame for
the full per-item table.

## See also

[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
and its `correct` argument (the polychoric failure this screens for);
[`factorability()`](https://jmgirard.github.io/ackwards/reference/factorability.md),
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md),
and
[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md),
the other pre-analysis diagnostics.

## Examples

``` r
# sim16 is clean continuous data -> nothing flagged
check_items(sim16, cor = "pearson")
#> 
#> ── Item quality check (ackwards) ───────────────────────────────────────────────
#> Basis: pearson
#> Items: 16
#> Flagged: 0
#> ✔ No item problems detected.
#> ────────────────────────────────────────────────────────────────────────────────
#> Constant items must be dropped (no variance). A near-constant item (one
#> response dominates) can yield a meaningless factor; a sparse category can make
#> `cor = "polychoric"` fail -- collapse rare categories, try `correct = 0`, or
#> drop the item. Full per-item table: treat this object as a data frame.

# An item dominated by one response with a couple of stray answers is flagged
d <- sim16
d$bad <- c(rep(1L, nrow(d) - 2L), 2L, 3L)
check_items(d)
#> 
#> ── Item quality check (ackwards) ───────────────────────────────────────────────
#> Basis: polychoric
#> Items: 17
#> Flagged: 17
#> 
#> ── Flagged items ──
#> 
#> ! near-constant: "bad"
#> ! sparse category: "i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10",
#> "i11", "i12", "i13", "i14", "i15", and "i16"
#> ────────────────────────────────────────────────────────────────────────────────
#> Constant items must be dropped (no variance). A near-constant item (one
#> response dominates) can yield a meaningless factor; a sparse category can make
#> `cor = "polychoric"` fail -- collapse rare categories, try `correct = 0`, or
#> drop the item. Full per-item table: treat this object as a data frame.
```
