# Genome-wide IBD relatedness network

Nodes are samples; an edge joins a pair sharing more than `min_ibd` of
the genome IBD, and its width is how much.
[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
asks who shares one locus, where every edge means the same thing; this
asks who is related overall, so the amount is the point.

## Usage

``` r
plot_ibd_pair_network(
  pairs,
  meta = NULL,
  weight = NULL,
  min_ibd = 0.01,
  samples = NULL,
  color_group = NULL,
  colors = NULL,
  shape_group = NULL,
  shapes = NULL,
  na_shape = .NA_SHAPE,
  na_colour = "grey70",
  include_isolated = TRUE,
  layout = "fr",
  spread = 1.5,
  node_size = 3,
  node_alpha = 0.9,
  edge_colour = "grey65",
  edge_alpha = 0.6,
  weight_range = c(0.15, 2.6),
  weight_breaks = NULL,
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

- pairs:

  The per-pair table from
  `plasgenomicsutils ibd_fraction_and_snp_density`
  (`*.pair_ibd_fraction.tsv.gz`): a path, a data frame, or an
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  carrying one. Needs `sample1`/`sample2`, which that command writes.

- meta:

  Sample metadata for `color_group` / `shape_group`; taken from the
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  when `pairs` is one.

- weight:

  Column holding the fraction. Defaults to `ibd_fraction_accessible`,
  the callable-genome denominator; `ibd_fraction_full_genome` divides by
  the whole genome instead, so it reads lower for the same pair.

- min_ibd:

  Draw an edge only above this fraction (default `0.01`). Every pair is
  in the table, most of them sharing essentially nothing, so without a
  cutoff the graph is complete.

- samples:

  Optional subset of samples to keep.

- color_group, colour_group, shape_group:

  Metadata columns for node colour and shape.

- colors, colours, shapes:

  Values for those scales, named or positional (see
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)).

- na_shape, na_colour, na_color:

  What a sample with no value in those columns gets.

  Legends stack in a fixed order – colour, shape, then the IBD width –
  so plots stay comparable; override with
  `+ ggplot2::guides(linewidth = ggplot2::guide_legend(order = 1))`.

- include_isolated:

  Show samples with no edge (default `TRUE`).

- layout, spread, seed:

  Layout algorithm, clique spreading, and the seed that makes it
  reproducible.

- node_size, node_alpha, edge_colour, edge_color, edge_alpha:

  Node and edge aesthetics.

- weight_range:

  Narrowest and widest edge, in `linewidth` units.

- weight_breaks:

  Legend breaks; defaults to powers of two spanning the data, since
  sharing runs over orders of magnitude.

- title:

  Plot title; `NULL` (default) writes one, `FALSE` or `NA` drops it.

- subtitle:

  `TRUE` (default) counts samples, edges and unconnected samples; a
  string replaces it, `FALSE` drops it.

## Value

A ggplot object.

## Details

Samples with no edge above the cutoff are drawn as a grid underneath,
separated by a dashed rule and counted – an unrelated sample is a
result, and dropping it silently would overstate how connected the
cohort is.

## See also

[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
for one gene or locus.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_ibd_pair_network("ibd_fraction.pair_ibd_fraction.tsv.gz", meta = meta,
                      color_group = "region", min_ibd = 0.03)
} # }
```
