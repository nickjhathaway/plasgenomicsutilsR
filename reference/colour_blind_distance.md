# Smallest colour difference within (or between) palettes, under colour blindness

The worst-case pair: for one palette, the two colours a viewer would
find hardest to tell apart; for two, the closest a colour in one gets to
a colour in the other. Reported as CIEDE2000 distance – roughly, below
~5 is a confusion risk, above ~10 is comfortable.

## Usage

``` r
colour_blind_distance(
  x,
  y = NULL,
  vision = c("normal", "deutan", "protan", "tritan")
)

color_blind_distance(
  x,
  y = NULL,
  vision = c("normal", "deutan", "protan", "tritan")
)
```

## Arguments

- x:

  A vector of colours.

- y:

  Optional second vector. When given, distances are measured *between*
  `x` and `y` (e.g. region colours against ancestry-component colours)
  rather than within `x`.

- vision:

  Which vision types to report, from `"normal"`, `"deutan"`, `"protan"`,
  `"tritan"`.

## Value

A named numeric vector, one worst-case distance per vision type.

## Examples

``` r
colour_blind_distance(color_palette(6))
#>    normal    deutan    protan    tritan 
#> 22.238274 12.521241  4.903496  8.349361 
# ColorBrewer "Paired", 6 colours: its light/dark pairs collapse for deutan vision
paired6 <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C")
colour_blind_distance(paired6)
#>    normal    deutan    protan    tritan 
#> 21.345256  4.731983 18.735687  8.176823 
```
