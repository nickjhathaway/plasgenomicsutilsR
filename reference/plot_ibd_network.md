# IBD network at a gene or locus

A sample-level network for one gene or locus: each node is a sample and
an edge joins two samples whose pair shares an **IBD block overlapping**
the interval. Optionally colours and/or shapes nodes by metadata
columns, and optionally drops isolated nodes for a cleaner picture – or
keeps them so the total sample count is visible.

## Usage

``` r
plot_ibd_network(
  x,
  gene = NULL,
  locus = NULL,
  color_group = NULL,
  colors = NULL,
  shape_group = NULL,
  shapes = NULL,
  na_shape = .NA_SHAPE,
  na_colour = "grey70",
  within = 0,
  sharing = c("overlap", "complete"),
  include_isolated = FALSE,
  layout = "fr",
  spread = 1.5,
  node_size = 3,
  node_alpha = 0.9,
  edge_colour = "grey65",
  edge_alpha = 0.5,
  edge_width = 1,
  title = NULL,
  subtitle = TRUE,
  seed = 42,
  colour_group = NULL,
  colours = NULL,
  na_color = NULL,
  edge_color = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  with IBD `blocks`.

- gene:

  A single gene name from the object's track (its interval is used).

- locus:

  Alternatively, a locus: `"chr:pos"`, `"chr:start-end"`, or a one-row
  data frame with `chr`, `start`, `end`. Give exactly one of `gene` /
  `locus`. Coordinates are 0-based half-open (see
  [plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)),
  so `"chr:1000-2000"` is the 1000 bases from 0-based 1000 up to but not
  including 2000, and a bare `"chr:1000"` is the single base at 0-based
  1000.

- color_group, colour_group:

  Optional metadata column to colour nodes by (needs `meta`). To colour
  by the graph's own connected components, add them to the metadata
  first with
  [`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)
  and name the column it creates.

- colors, colours:

  Colours for `color_group`. A **named** `level -> colour` vector maps
  by name and may cover only some levels (the rest keep their automatic
  colour); an unnamed vector is taken positionally, in the column's
  level order. `NULL` (default) uses
  [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md),
  so the same mapping can be shared with other plots.

- shape_group:

  Optional metadata column to set node **shape** by, so a second
  variable can be read off the same plot (needs `meta`). Independent of
  `color_group`: use either, both, or neither.

- shapes:

  Shapes for `shape_group`, named or positional exactly like `colors`
  (values are ggplot2 shape codes). `NULL` picks distinguishable
  defaults and errors if the column has more levels than there are
  distinct shapes.

- na_shape, na_colour, na_color:

  What a sample with **no value** in `shape_group` / `color_group` gets:
  by default a hollow circle (shape `1`) and grey. The default shape is
  deliberately outside the automatic palette so it cannot be mistaken
  for a real level; if you set it to one that a level also uses, the
  plot says so rather than letting two things look alike.

  The legends stack in a fixed order – colour, then shape – so two plots
  of the same cohort are comparable. ggplot2 otherwise orders guides by
  an internal hash that changes with the labels. Override per plot with
  `+ ggplot2::guides(shape = ggplot2::guide_legend(order = 1))`.

- within:

  Pad the interval by this many bp on both sides (default `0`).

- sharing:

  What an edge requires of a pair's IBD segment:

  `"overlap"`

  :   (default) the segment touches the interval anywhere – the pair
      shares *some* of the gene/locus.

  `"complete"`

  :   the segment spans the whole interval – the pair shares the
      *entire* gene/locus. A stricter, usually much sparser graph.

  `within` applies either way, so with padding `"complete"` asks the
  segment to cover the padded interval. Use
  [`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
  to see per-pair which of the two each segment satisfies (its
  `coverage` column) and how much of the gene is covered.

- include_isolated:

  Keep samples with no IBD edge (default `FALSE`). `TRUE` also shows
  every analyzed sample; the isolated samples are drawn in a grid below
  the connected component (sorted by `group` when colouring by one).
  `FALSE` shows only the connected nodes.

- layout:

  A ggraph/igraph layout name (default `"fr"`, Fruchterman- Reingold).
  Other options include `"kk"` (Kamada-Kawai), `"stress"`, `"drl"`,
  `"circle"`, `"nicely"`, `"graphopt"`, and `"lgl"`. The connected
  component is laid out on its own, with its aspect preserved.

- spread:

  How strongly to open up densely inter-connected groups (default
  `1.5`). Edge attraction is weighted by `(1 - J)^spread`, where `J` is
  the Jaccard overlap of the two samples' IBD neighbourhoods: edges
  inside a near-clique pull only weakly, so the group expands enough to
  resolve individual samples, while edges linking otherwise separate
  groups keep their full pull and cluster separation is preserved. `0`
  leaves edges unweighted; raise it above the default to loosen the
  densest groups further, at the cost of them taking up more of the
  canvas. Applies to weight-aware layouts (`"fr"`, `"kk"`, `"drl"`,
  `"stress"`, ...); others are unaffected.

- node_size, node_alpha:

  Node point aesthetics.

- edge_colour, edge_color, edge_width:

  Edge aesthetics (default width `1`).

- edge_alpha:

  Edge opacity (default `0.5`). `NULL` instead scales opacity down with
  edge count (`~120 / n_edges`, clamped to \[0.06, 0.6\]).

- title:

  Plot title: `NULL` (default) uses `"IBD network: <label>"`, a string
  sets a custom title, and `NA`/`FALSE` draws no title.

- subtitle:

  Show the sample / IBD-pair count line under the title (default
  `TRUE`).

- seed:

  Random seed for stochastic layouts (e.g. `"fr"`).

## Value

A ggplot (ggraph) object.

## Details

Needs an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
built with `blocks =` (and `meta =` for grouping / the full sample set).
See
[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
for the pairwise-fraction summary of the same data.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", meta = meta, genes = PF_EXAMPLE_DRUG_GENES)
plot_ibd_network(ibd, gene = "pfcrt", color_group = "region")
# colour by region, shape by year: two variables on one plot
plot_ibd_network(ibd, gene = "pfcrt", color_group = "region",
                 shape_group = "collection_year")
# a named vector pins chosen levels; the rest keep their automatic colour
plot_ibd_network(ibd, gene = "pfcrt", color_group = "region",
                 colors = c(North = "#1f78b4", East = "#33a02c"))
plot_ibd_network(ibd, locus = "Pf3D7_07_v3:403500", include_isolated = TRUE)
} # }
```
