# Triangle heatmap of group x group differentiation

A lower-triangle heatmap of a group summary (see
[`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)),
styled like the drug-gene triangles: the fill legend sits in the empty
upper triangle. Optionally adds a clustering dendrogram across the top
and a metadata annotation strip (e.g. colour each site by its country).

## Usage

``` r
plot_diff_heatmap(
  pd,
  stat = c("mean", "median", "top_mean", "max"),
  top = 0.05,
  cluster = TRUE,
  dendrogram = TRUE,
  triangle = TRUE,
  legend_inside = triangle,
  annotate = NULL,
  meta = NULL,
  annotate_colours = NULL,
  label = FALSE,
  digits = 2,
  colors = c("white", "#fde0dd", "#fa9fb5", "#c51b8a", "#7a0177"),
  trans = "identity",
  base_size = 11
)

plot_jost_d_heatmap(pd, stat = "mean", ...)
```

## Arguments

- pd:

  A
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  result.

- stat, top:

  Summary passed to
  [`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  (`"mean"`, `"median"`, `"top_mean"`, `"max"`; `top` for the
  top-percentile mean).

- cluster:

  Order groups by hierarchical clustering of the matrix (default
  `TRUE`), which is what puts similar groups side by side. With
  `cluster = FALSE` the axes follow the grouping column's levels instead
  – its factor order if it has one (e.g. from
  `PopStructure$set_levels()`), otherwise a natural sort. Annotation
  strips and the dendrogram tips take their colours in level order
  either way, so a group keeps the same colour whichever ordering is
  drawn.

- dendrogram:

  Draw the clustering dendrogram across the top (needs `cluster`).

- triangle:

  Show only the lower triangle (with the legend in the empty corner).

- legend_inside:

  Put the legends in the empty half of the matrix rather than in a
  margin column (default: on whenever `triangle` is `TRUE`, since that
  is what leaves the space). The annotation strips' legends move onto
  the matrix too, so there is one legend area. Turn it off when there
  are few groups: the empty corner is then small and the legends would
  sit over the cells.

- annotate:

  One or more annotations drawn as colour strips above the columns: a
  metadata column name (looked up in `meta`), a vector of several column
  names, a named `group -> value` vector, or a named list mixing these
  (the names label the strips/legends). Each strip aligns to the
  clustered column order.

- meta:

  Metadata data frame for resolving `annotate` column names.

- annotate_colours:

  Named list `annotation -> (value -> colour)` giving custom colours per
  annotation (unlisted annotations get an automatic palette).

- label:

  Print the value in each cell (default `FALSE`; the text can distract).

- digits:

  Cell-label digits.

- colors:

  Fill ramp (low -\> high differentiation).

- trans:

  Fill transform, e.g. `"identity"` (default) or `"sqrt"` to lift a
  scale dominated by near-zero values.

- base_size:

  Base font size.

- ...:

  For `plot_jost_d_heatmap()`, arguments forwarded to
  `plot_diff_heatmap()`.

## Value

A ggplot (or patchwork) object.
