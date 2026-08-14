# PCA scatter plot

PCA scatter plot

## Usage

``` r
plot_pca(
  x,
  pcs = c(1, 2),
  colour = NULL,
  colors = NULL,
  point_size = 1.6,
  point_alpha = 0.8,
  legend_point_size = 3,
  color = NULL,
  colours = NULL
)
```

## Arguments

- x:

  A `pop_structure` object or a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  R6 object.

- pcs:

  Which two PCs to plot (default `c(1, 2)`).

- colour, color:

  Metadata column to colour points by (needs `meta`).

- colors, colours:

  Optional named `level -> colour` vector for the colour scale (e.g.
  from
  [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md));
  a `PopStructure` supplies its stored map.

- point_size, point_alpha:

  Point aesthetics.

- legend_point_size:

  Size of the coloured dots in the legend (default `3`, larger than the
  plotted points so the key is easy to read); `NULL` leaves it as-is.

## Value

A ggplot object.
