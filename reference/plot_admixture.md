# Admixture (STRUCTURE) bar plot from a Q matrix

Stacked ancestry-proportion bars, one per sample, optionally faceted by
a grouping column, with an optional group colour strip whose colours can
match another plot (e.g. UMAP points).

## Usage

``` r
plot_admixture(
  q,
  samples = NULL,
  meta = NULL,
  group = NULL,
  order_within = TRUE,
  sample_order = NULL,
  colours = NULL,
  group_bar = FALSE,
  group_colours = NULL,
  border = TRUE,
  border_colour = "black",
  border_linewidth = 0.15
)
```

## Arguments

- q:

  A samples-by-K ancestry matrix (e.g. from
  [`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)).

- samples:

  Sample ids (defaults to `rownames(q)`).

- meta:

  Optional metadata data frame with a `sample` column.

- group:

  Optional metadata column to facet by (e.g. `"region"`).

- order_within:

  Order samples within each group by clustering (ignored when
  `sample_order` is supplied).

- sample_order:

  Optional explicit sample order (see
  [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md));
  keeps bars in the same place across different K.

- colours:

  Optional fill colours for the K clusters.

- group_bar:

  Draw a coloured strip above the bars keyed by `group`.

- group_colours:

  Named `level -> colour` vector for the group strip (see
  [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md));
  pass the same mapping you use for the UMAP to match colours.

- border:

  Outline each sample's bar (default `TRUE`) so neighbours with nearly
  identical ancestry stay distinct; set `FALSE` for borderless bars.

- border_colour, border_linewidth:

  Colour and width of the per-sample outline.

## Value

A ggplot object.
