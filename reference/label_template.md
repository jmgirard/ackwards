# Generate a node-label scaffold for autoplot

Returns a named character vector covering all factor IDs in the object,
ready to pass to `autoplot(x, node_labels = ...)`. Printing the result
shows an editable `c(...)` literal so you can copy it into a script,
fill in substantive labels, and use it directly; assigning the result
(`labs <- label_template(x)`) produces no console output.

## Usage

``` r
label_template(x, style = c("id", "forbes", "blank"))

# S3 method for class 'ackwards_labels'
print(x, ...)
```

## Arguments

- x:

  An `ackwards` object.

- style:

  One of `"id"` (default), `"forbes"`, or `"blank"`. See **Style
  options** above.

- ...:

  Ignored; included for S3 method consistency.

## Value

A named character vector of class `"ackwards_labels"`: names are factor
IDs (`"m{k}f{j}"`), values are the label strings for the chosen `style`.
The vector is suitable for direct use as the `node_labels` argument to
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md).
Its [`print()`](https://rdrr.io/r/base/print.html) method renders the
vector as an editable `c(...)` literal for copy-paste.

## Details

The factor IDs are returned in the same left-to-right, top-to-bottom
order used by
[`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)
and
[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md),
so the printed literal maps directly onto the diagram.

## Style options

- `"id"` *(default)* – every value equals the factor ID (`"m1f1"`,
  `"m2f1"`, ...). This is a round-trip no-op: passing the result to
  `node_labels` without editing reproduces the default labels exactly.
  Useful as the starting point for adding substantive labels.

- `"forbes"` – values follow the Forbes (2023) convention:
  level-letter + within-level index (`"A1"`, `"B1"`, `"B2"`, ...). Level
  1 -\> `A`, level 2 -\> `B`, level 3 -\> `C`, and so on. Within-level
  indices are assigned in canonical layout order (left to right).
  Requires `k_max <= 26` (LETTERS has 26 entries); an error is raised
  for deeper objects.

- `"blank"` – all values are empty strings. Useful as a starting
  scaffold when you want to supply every label from scratch with no
  defaults showing through.

## See also

[`autoplot.ackwards()`](https://jmgirard.github.io/ackwards/reference/autoplot.ackwards.md),
[`top_items()`](https://jmgirard.github.io/ackwards/reference/top_items.md),
[`ba_layout()`](https://jmgirard.github.io/ackwards/reference/ba_layout.md)

## Examples

``` r
x <- ackwards(sim16, k_max = 5)

# Start from ID defaults, then fill in your own labels:
labs <- label_template(x)
labs["m5f1"] <- "My factor name"

# Forbes letter convention:
label_template(x, style = "forbes")
#> `label_template()` scaffold (forbes style):
#> c(
#>   "m1f1" = "A1",
#>   "m2f1" = "B1",
#>   "m2f2" = "B2",
#>   "m3f1" = "C1",
#>   "m3f2" = "C2",
#>   "m3f3" = "C3",
#>   "m4f1" = "D1",
#>   "m4f2" = "D2",
#>   "m4f3" = "D3",
#>   "m4f4" = "D4",
#>   "m5f1" = "E1",
#>   "m5f2" = "E2",
#>   "m5f3" = "E3",
#>   "m5f4" = "E4",
#>   "m5f5" = "E5"
#> )
```
