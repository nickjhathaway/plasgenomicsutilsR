# Assign colours to the levels of metadata columns

Builds a named list (one entry per column) of `level -> hex colour`
vectors, so a single mapping can colour points in one plot and bars in
another. Levels are taken in factor order (or sorted); a
colour-blind-friendly palette is chosen by level count (interpolated
past 15). `overrides` replaces specific colours.

## Usage

``` r
meta_colors(meta, cols = NULL, overrides = NULL)
```

## Arguments

- meta:

  A data frame of metadata.

- cols:

  Columns to build colours for (default: all but a `sample` column).

- overrides:

  Optional named list `column -> (level -> colour)` to override
  individual assignments (e.g.
  `list(country = c(Ethiopia = "#FFDC3D"))`).

## Value

A named list, one `level -> colour` named vector per column.

## Examples

``` r
meta_colors(data.frame(sample = 1:3, region = c("A", "B", "A")))
#> $region
#>         A         B 
#> "#2271B2" "#3DB7E9" 
#> 
```
