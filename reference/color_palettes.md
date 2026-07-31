# Colour-blind-friendly categorical palettes

Three categorical palettes (8, 12, 15 colours) for colouring metadata
levels.

## Usage

``` r
colorPalette_08

colorPalette_12

colorPalette_15
```

## Format

Character vectors of hex colours.

## Examples

``` r
if (requireNamespace("scales", quietly = TRUE)) scales::show_col(colorPalette_12)
```
