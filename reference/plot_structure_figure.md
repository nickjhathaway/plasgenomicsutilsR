# Combined UMAP + admixture figure

Draws a UMAP scatter and an sNMF admixture bar plot as **one** figure
that shares a single theme (matching font sizes), a single region colour
map (so UMAP points and the admixture colour strips match), and
collected legends. The admixture is faceted by a metadata column with a
colour strip – and no text label – over each region, samples ordered
once and reused, laid out over one or more rows.

## Usage

``` r
plot_structure_figure(
  x,
  group = NULL,
  colour = group,
  K = NULL,
  rows = NULL,
  orientation = c("vertical", "horizontal"),
  sample_order = NULL,
  umap_colour = colour,
  region_colours = NULL,
  cluster_colours = NULL,
  region_label = NULL,
  base_size = 11,
  border = TRUE,
  border_colour = "black",
  border_linewidth = 0.15,
  legend = "right",
  legend_point_size = 3.5,
  point_size = 1.6,
  point_alpha = 0.8,
  umap_ratio = 1,
  file = NULL,
  width = NULL,
  height = NULL,
  color = NULL,
  border_color = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  with a UMAP (`run_umap()`) and an sNMF fit
  ([`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)).

- group:

  Metadata column to facet/colour the admixture by (default the first
  non-`sample` metadata column).

- colour, color:

  Metadata column supplying the shared region colours (default `group`).

- K:

  Number of ancestral populations: an integer for a specific K, or
  `NULL` / `"best_k"` to use the cross-entropy best K
  ([PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)'s
  `best_k()`).

- rows:

  A list of character vectors giving which `group` levels sit on each
  admixture row, e.g.
  `list(c("DRC", "Kenya"), c("Tanzania", "Uganda"))`. Default: all
  levels on one row. Levels not listed are dropped.

- orientation:

  `"vertical"` puts the UMAP above the admixture; `"horizontal"` puts it
  to the left.

- sample_order:

  Optional explicit sample order (see
  [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md));
  by default computed once (within group) and reused across the figure.

- umap_colour:

  Metadata column colouring the UMAP points (default `colour`).

- region_colours, cluster_colours:

  Optional named colour vectors overriding the region strip / K-cluster
  fills.

- region_label:

  Legend title for the region colours (default `colour`).

- base_size:

  Base font size shared by every panel.

- border:

  Outline each sample's admixture bar (default `TRUE`) so neighbours
  with nearly identical ancestry stay distinct.

- border_colour, border_color, border_linewidth:

  Colour and width of the per-sample outline.

- legend:

  Where to collect the shared legends (`"right"`, `"left"`, `"bottom"`,
  `"top"`).

- legend_point_size:

  Size of the region dots in the UMAP legend (default `3.5`, so the key
  reads clearly next to the small plotted points).

- point_size, point_alpha:

  Size and opacity of the UMAP scatter points.

- umap_ratio:

  Relative size of the UMAP vs the admixture block (a single number;
  default 1 means roughly equal).

- file:

  Optional path to save to (via
  [`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md))
  at the auto-computed size.

- width, height:

  Optional output size in inches (default auto from sample counts, row
  count, and orientation); also attached to the result as attributes.

## Value

A patchwork object (invisibly if `file` is given), with `width`/`height`
attributes carrying the suggested output size.

## Examples

``` r
if (FALSE) { # \dontrun{
ps <- example_pop_structure("africa")
ps$run_snmf(K = 1:9)
plot_structure_figure(ps, group = "site",
  rows = list(c("DRC", "Ethiopia", "Sudan"),
              c("Kenya_East", "Kenya_West", "Tanzania_East", "Tanzania_West")))
} # }
```
