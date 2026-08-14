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
  border_linewidth = 0.15,
  legend_position = c("right", "bottom", "top", "left", "none"),
  legend_rows = NULL,
  colors = NULL,
  group_colors = NULL,
  border_color = NULL
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

- colours, colors:

  Optional fill colours for the K clusters.

- group_bar:

  Draw a coloured strip above the bars keyed by `group`.

- group_colours, group_colors:

  Named `level -> colour` vector for the group strip (see
  [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md));
  pass the same mapping you use for the UMAP to match colours.

- border:

  Outline each sample's bar (default `TRUE`, matching
  [`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md))
  so neighbours with nearly identical ancestry stay distinct; `FALSE`
  for borderless bars. With very many samples in a narrow render the
  outlines can swamp the fills – either render wider
  ([`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
  uses the attached width) or set `border = FALSE` / a thinner
  `border_linewidth`.

- border_colour, border_color, border_linewidth:

  Colour and width of the per-sample outline.

- legend_position:

  Where the legends go: `"right"` (default), `"bottom"`, `"top"`,
  `"left"`, or `"none"`. A large K makes for a tall legend stack, so
  `"bottom"` often fits better on a wide, short admixture panel.

- legend_rows:

  Keys per legend column (or per row when the legend is horizontal)
  before wrapping into another column/row. `NULL` (default) wraps a side
  legend every 10 keys and splits a horizontal one over two rows, which
  keeps `K` = 15 plus a group strip on the page. The suggested output
  height accounts for whatever this works out to.

## Value

A ggplot object.
